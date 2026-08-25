"""
cudaq_core.py - Native NVIDIA CUDA-Q implementation of the SSH-Hubbard
adiabatic-evolution circuit.

Unlike the old backend (which built a qiskit circuit, transpiled it, and
replayed the gates through CUDA-Q's deprecated `make_kernel` builder), this
module expresses the whole algorithm directly in CUDA-Q:

  * the circuit is executed by a single `@cudaq.kernel` (JIT-compiled to MLIR
    and run on the GPU statevector simulator via `cudaq.get_state`);
  * all gates are CUDA-Q native operations: x, h, rx, rz, s / s.adj,
    x.ctrl (CNOT), ry.ctrl (controlled-RY for Givens rotations) and
    r1.ctrl (controlled-phase for the Hubbard interaction);
  * qiskit is NOT needed anywhere on this path.

Because CUDA-Q kernels are statically compiled, the circuit is encoded
classically as a flat instruction program (parallel lists of opcode / qubit
indices / angle) and interpreted inside one kernel. This keeps a single
compiled kernel for every (N, T_A, steps, U_A, U_B, ...) parameter set, so
sweeps do not pay a re-JIT cost per point.

Gate conventions match `core.build_annealing_circuit` exactly (same order,
same angles), so a qiskit run and a cudaq run of the same experiment produce
identical statevectors -- verified by `aqs selftest --backend cudaq`.

Decompositions used (all exact, standard identities):

  RZZ(phi) = CX . (I x RZ(phi)) . CX                    = exp(-i phi/2 ZZ)
  RXX(phi) = (H x H)       RZZ(phi) (H x H)             = exp(-i phi/2 XX)
  RYY(phi) = (RX(pi/2)^x2) RZZ(phi) (RX(-pi/2)^x2)      = exp(-i phi/2 YY)

  R(theta) = RXX(theta/2) RYY(theta/2)          (real hopping, cf. core.py)
  G(theta) = Sdg_0 RXX(theta/2) S_0 Sdg_1 RXX(-theta/2) S_1   (imag hopping)
  CP(theta) = r1.ctrl(theta)                     (on-site Hubbard n_up n_dn)
  Givens(theta, phi) = RZ(phi)_k RZ(-phi)_j CX_{jk} CRY(-2 theta)_{k->j} CX_{jk}
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------- opcodes
OP_X, OP_H, OP_RX, OP_RZ, OP_CX, OP_CRY, OP_CR1, OP_S, OP_SDG, OP_RY = range(10)


@dataclass
class InstructionProgram:
    """Flat gate program consumed by the CUDA-Q interpreter kernel."""

    n_qubits: int
    op: list = field(default_factory=list)     # opcode per instruction
    q0: list = field(default_factory=list)     # first qubit (or control)
    q1: list = field(default_factory=list)     # second qubit (or target; -1 if unused)
    ang: list = field(default_factory=list)    # rotation angle (0.0 if unused)

    # ---- primitive emitters ------------------------------------------------
    def _emit(self, op, q0, q1=-1, ang=0.0):
        self.op.append(int(op)); self.q0.append(int(q0))
        self.q1.append(int(q1)); self.ang.append(float(ang))

    def x(self, q):            self._emit(OP_X, q)
    def h(self, q):            self._emit(OP_H, q)
    def s(self, q):            self._emit(OP_S, q)
    def sdg(self, q):          self._emit(OP_SDG, q)
    def rx(self, ang, q):      self._emit(OP_RX, q, ang=ang)
    def ry(self, ang, q):      self._emit(OP_RY, q, ang=ang)
    def rz(self, ang, q):      self._emit(OP_RZ, q, ang=ang)
    def cx(self, c, t):        self._emit(OP_CX, c, t)
    def cry(self, ang, c, t):  self._emit(OP_CRY, c, t, ang)
    def cr1(self, ang, c, t):  self._emit(OP_CR1, c, t, ang)

    # ---- composite gates (exact decompositions) ----------------------------
    def rzz(self, phi, a, b):
        self.cx(a, b); self.rz(phi, b); self.cx(a, b)

    def rxx(self, phi, a, b):
        self.h(a); self.h(b)
        self.rzz(phi, a, b)
        self.h(a); self.h(b)

    def ryy(self, phi, a, b):
        self.rx(math.pi / 2, a); self.rx(math.pi / 2, b)
        self.rzz(phi, a, b)
        self.rx(-math.pi / 2, a); self.rx(-math.pi / 2, b)

    def R(self, theta, a, b):
        """Real-hopping gate: RXX(theta/2) RYY(theta/2) (== core.create_R_gate)."""
        self.rxx(theta / 2.0, a, b)
        self.ryy(theta / 2.0, a, b)

    def G(self, theta, a, b):
        """Imaginary-hopping gate (== core.create_G_gate)."""
        self.sdg(a); self.rxx(theta / 2.0, a, b); self.s(a)
        self.sdg(b); self.rxx(-theta / 2.0, a, b); self.s(b)

    def CP(self, theta, a, b):
        """On-site Hubbard interaction (== core.create_CP_gate)."""
        self.cr1(theta, a, b)

    def givens(self, theta, phi, j, k):
        """Givens rotation (== core._givens_instruction on qubits [j, k])."""
        self.rz(phi, k)
        self.rz(-phi, j)
        self.cx(j, k)
        self.cry(-2.0 * theta, k, j)
        self.cx(j, k)


# =====================================================================
# Circuit construction (mirrors core.build_annealing_circuit 1:1)
# =====================================================================
def build_annealing_program(model, Q_up, Q_dn,
                            U_A=0.0, U_B=0.0, T_A=1.0, steps=0,
                            ramp_U=True) -> InstructionProgram:
    """Native CUDA-Q version of `core.build_annealing_circuit`.

    Emits the same gate sequence with the same angles; the returned
    InstructionProgram is executed by `statevector(program)` below.
    Slater-determinant Givens parameters come from OpenFermion (classical
    preprocessing, identical to the qiskit path).
    """
    import openfermion  # classical pre-processing only

    N_cells = model.N_cells
    v, w = model.v, model.w
    PBC = model.PBC
    num_up, num_dn = model.num_up, model.num_dn
    L = 2 * N_cells
    prog = InstructionProgram(n_qubits=2 * L)

    # ---- 1. occupation pattern ------------------------------------------
    for i in range(num_up):
        prog.x(i)
    for i in range(num_dn):
        prog.x(L + i)

    # ---- 2. Slater-determinant preparation (Givens network) --------------
    for parallel_ops in openfermion.circuits.slater_determinant_preparation_circuit(Q_up):
        for j, k, theta, phi in parallel_ops:
            prog.givens(theta, phi, j, k)
    for parallel_ops in openfermion.circuits.slater_determinant_preparation_circuit(Q_dn):
        for j, k, theta, phi in parallel_ops:
            prog.givens(theta, phi, j + L, k + L)

    # ---- 3. Trotterized adiabatic evolution -------------------------------
    if steps > 0:
        tau = T_A / steps
        parity_up = (-1) ** num_up
        parity_dn = (-1) ** num_dn
        for step in range(1, steps + 1):
            s = (step - 0.5) / steps

            def apply_hopping_layer(start_i, t_val):
                for i in range(start_i, L, 2):
                    is_pbc_bond = i == L - 1
                    if not PBC and is_pbc_bond:
                        continue
                    j = (i + 1) % L
                    w_R, w_I = float(np.real(t_val)), float(np.imag(t_val))
                    theta_R = -2.0 * tau * w_R
                    theta_I = -2.0 * tau * w_I
                    f_R_up = -parity_up * theta_R if is_pbc_bond else theta_R
                    f_I_up = -parity_up * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        prog.R(f_R_up, i, j)
                    if abs(w_I) > 1e-8:
                        prog.G(f_I_up, i, j)
                    f_R_dn = -parity_dn * theta_R if is_pbc_bond else theta_R
                    f_I_dn = -parity_dn * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        prog.R(f_R_dn, L + i, L + j)
                    if abs(w_I) > 1e-8:
                        prog.G(f_I_dn, L + i, L + j)

            apply_hopping_layer(0, v)
            apply_hopping_layer(1, w)

            if U_A != 0 or U_B != 0:
                ramp = s if ramp_U else 1.0
                for i in range(L):
                    current_U = U_A if i % 2 == 0 else U_B
                    if current_U != 0:
                        theta_U = -tau * ramp * current_U
                        prog.CP(theta_U, i, L + i)

    return prog


# =====================================================================
# Kernel execution
# =====================================================================
_KERNEL = None


def _get_kernel():
    """Compile (once) the CUDA-Q interpreter kernel for InstructionPrograms."""
    global _KERNEL
    if _KERNEL is not None:
        return _KERNEL

    import cudaq

    @cudaq.kernel
    def aqs_program_kernel(n: int, op: list[int], q0: list[int],
                           q1: list[int], ang: list[float]):
        q = cudaq.qvector(n)
        for k in range(len(op)):
            o = op[k]
            if o == 0:                       # X
                x(q[q0[k]])
            elif o == 1:                     # H
                h(q[q0[k]])
            elif o == 2:                     # RX
                rx(ang[k], q[q0[k]])
            elif o == 3:                     # RZ
                rz(ang[k], q[q0[k]])
            elif o == 4:                     # CNOT
                x.ctrl(q[q0[k]], q[q1[k]])
            elif o == 5:                     # controlled-RY (Givens)
                ry.ctrl(ang[k], q[q0[k]], q[q1[k]])
            elif o == 6:                     # controlled-phase (Hubbard)
                r1.ctrl(ang[k], q[q0[k]], q[q1[k]])
            elif o == 7:                     # S
                s(q[q0[k]])
            elif o == 8:                     # S-dagger
                s.adj(q[q0[k]])
            elif o == 9:                     # RY
                ry(ang[k], q[q0[k]])

    _KERNEL = aqs_program_kernel
    return _KERNEL


def statevector(program: InstructionProgram) -> np.ndarray:
    """Run an InstructionProgram with `cudaq.get_state` and return the final
    statevector in qiskit's little-endian amplitude ordering (qubit k = bit k),
    which is what aqs.observables expects.

    cudaq orders qubit 0 as the most-significant bit, so logical qubit i is
    mapped onto CUDA-Q qubit (n-1-i); this is a pure wire relabelling and
    makes the returned amplitude indices match bit-for-bit.
    """
    import cudaq  # noqa: F401  (ensures cudaq is importable before JIT)

    n = program.n_qubits
    q0 = [n - 1 - i for i in program.q0]
    q1 = [(n - 1 - i) if i >= 0 else 0 for i in program.q1]

    kernel = _get_kernel()
    state = cudaq.get_state(kernel, n, program.op, q0, q1, program.ang)
    return np.array(state, copy=True).astype(complex).ravel()
