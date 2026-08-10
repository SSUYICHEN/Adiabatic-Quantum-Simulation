"""
core.py - Built-in SSH-Hubbard model: gates and the annealing circuit.

Qubit layout (spin-block mapping)
---------------------------------
L = 2 * N_cells  spatial orbitals per spin.
Qubits [0, L)        -> spin-up   block
Qubits [L, 2L)       -> spin-down block
Within a block, even index = A sublattice, odd index = B sublattice.

Hopping amplitudes follow the standard SSH naming:
    v = intra-cell hopping (even bonds)
    w = inter-cell hopping (odd bonds)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import openfermion
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate


# =====================================================================
# 1. Physical model
# =====================================================================
@dataclass
class SSHHModel:
    """One-body SSH-Hubbard formulation on N_cells unit cells.

    Parameters
    ----------
    N_cells : number of unit cells (each cell = one A + one B site).
    v, w    : intra-cell and inter-cell hopping amplitudes.
    num_up, num_dn : electron counts per spin sector.
    PBC     : periodic (True) or open (False) boundary conditions.
    """

    N_cells: int
    v: complex
    w: complex
    num_up: int
    num_dn: int
    PBC: bool = True

    def __post_init__(self):
        self.L = 2 * self.N_cells
        self.number_of_electrons = self.num_up + self.num_dn
        self.H = self._build_hamiltonian()

    @classmethod
    def from_total_electrons(cls, N_cells, v, w, number_of_electrons, PBC=True):
        num_up = (number_of_electrons // 2) + (number_of_electrons % 2)
        num_dn = number_of_electrons // 2
        return cls(N_cells, v, w, num_up, num_dn, PBC)

    def _build_hamiltonian(self) -> np.ndarray:
        dim = 2 * self.L
        H = np.zeros((dim, dim), dtype=complex)
        for i in range(self.L):
            if not self.PBC and i == self.L - 1:
                continue
            j = (i + 1) % self.L
            t = self.v if i % 2 == 0 else self.w
            H[i, j] = -t
            H[j, i] = -np.conj(t)
            H[self.L + i, self.L + j] = -t
            H[self.L + j, self.L + i] = -np.conj(t)
        return H

    def slater_Q_matrices(self):
        """Occupied single-particle orbitals for Slater-determinant prep."""
        H_single = self.H[: self.L, : self.L]
        eigenvalues, eigenvectors = np.linalg.eigh(H_single)
        idx = np.argsort(eigenvalues)
        wf = eigenvectors[:, idx]
        Q_up = wf[:, : self.num_up].T
        Q_dn = wf[:, : self.num_dn].T
        return Q_up, Q_dn, self.num_up, self.num_dn


# =====================================================================
# 2. Modular two-qubit gates
# =====================================================================
def create_R_gate(theta):
    """R(theta) = exp(-i theta/2 (XX + YY)) : real (Hermitian) hopping term.

    Since rxx(p) = exp(-i p/2 XX) and XX commutes with YY, this is exactly
    rxx(theta) . ryy(theta) -- no factor of two anywhere.
    """
    qc = QuantumCircuit(2, name=f"R({theta:.3f})")
    qc.rxx(theta, 0, 1)
    qc.ryy(theta, 0, 1)
    return qc.to_instruction()


def create_G_gate(theta):
    """G(theta) = exp(-i theta/2 (X0 Y1 - Y0 X1)) : imaginary hopping term.

    Note the operator ordering: X on the FIRST qubit passed, Y on the second.
    G is antisymmetric under swapping the two qubits, unlike R.
    """
    qc = QuantumCircuit(2, name=f"G({theta:.3f})")
    qc.sdg(0); qc.rxx(-theta, 0, 1); qc.s(0)
    qc.sdg(1); qc.rxx(theta, 0, 1); qc.s(1)
    return qc.to_instruction()


def create_CP_gate(theta):
    """CP(theta) = exp(-i theta n_i n_j) : density-density / on-site Hubbard term.

    qiskit's cp(l) = diag(1, 1, 1, e^{+i l}), so the argument is negated to
    realise exp(-i theta n_i n_j) = diag(1, 1, 1, e^{-i theta}).
    """
    qc = QuantumCircuit(2, name=f"CP({theta:.3f})")
    qc.cp(-theta, 0, 1)
    return qc.to_instruction()


def _givens_instruction(theta, phi):
    """Two-qubit Givens rotation in OpenFermion's convention.

    slater_determinant_preparation_circuit emits (i, j, theta, phi) tuples whose
    single-particle rotation is

        G = [[cos t,  -e^{i p} sin t],
             [sin t,   e^{i p} cos t]],        det G = e^{i p}

    Lifted to the two-qubit occupation basis this leaves |00> alone, applies G
    on the single-occupancy subspace, and multiplies |11> by det G.

    Previously decomposed as rz(phi,1), rz(-phi,0), cx, cry(-2*theta), cx, which
    had two defects:

      * the cry angle was negated. For real hoppings that is exactly
        G(-theta, 0), so the circuit prepared the determinant built from the
        HIGHEST single-particle orbitals. On a chiral-symmetric chain those
        energies are the exact negatives of the lowest ones, so preparation
        landed on +E instead of the ground energy -E.

      * for phi != 0 it was not a mis-parameterised Givens rotation at all: it
        left |11> with phase 1 instead of det G = e^{i phi}, and no choice of
        (theta, phi) reproduces the correct gate (global fit residual 0.64).

    Verified at fidelity 1.0 against Slater determinants constructed directly
    with OpenFermion, for real and complex hoppings, OBC and PBC, N = 2..4.
    """
    c, s = np.cos(theta), np.sin(theta)
    e = np.exp(1j * phi)
    U = np.zeros((4, 4), dtype=complex)
    U[0, 0] = 1.0
    U[1, 1] = c
    U[2, 2] = c * e
    U[1, 2] = s
    U[2, 1] = -s * e
    U[3, 3] = e
    return UnitaryGate(U, label="Givens")


# =====================================================================
# 3. Annealing circuit (parity-corrected PBC, staggered Hubbard U)
# =====================================================================
def build_annealing_circuit(model: SSHHModel, Q_up, Q_dn,
                            U_A=0.0, U_B=0.0, T_A=1.0, steps=0, ramp_U=True):
    """Build the Trotterized adiabatic-evolution circuit (see README)."""
    N_cells = model.N_cells
    v, w = model.v, model.w
    PBC = model.PBC
    num_up, num_dn = model.num_up, model.num_dn
    L = 2 * N_cells
    qc = QuantumCircuit(2 * L)

    for i in range(num_up):
        qc.x(i)
    for i in range(num_dn):
        qc.x(L + i)

    for parallel_ops in openfermion.circuits.slater_determinant_preparation_circuit(Q_up):
        for j, k, theta, phi in parallel_ops:
            qc.append(_givens_instruction(theta, phi), [j, k])
    for parallel_ops in openfermion.circuits.slater_determinant_preparation_circuit(Q_dn):
        for j, k, theta, phi in parallel_ops:
            qc.append(_givens_instruction(theta, phi), [j + L, k + L])
    qc.barrier()

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
                    # exp(-i tau H_bond) with H_bond = -Re(t)/2 (XX+YY)
                    #                              + Im(t)/2 (X0Y1 - Y0X1)
                    # in qiskit's little-endian basis, so R takes -tau*Re and
                    # G takes +tau*Im. (Beware: openfermion's big-endian matrix
                    # conjugates the imaginary part -- always compare against a
                    # little-endian H when re-deriving the G sign.)
                    theta_R = -tau * w_R
                    theta_I = tau * w_I
                    f_R_up = -parity_up * theta_R if is_pbc_bond else theta_R
                    f_I_up = -parity_up * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        qc.append(create_R_gate(f_R_up), [i, j])
                    if abs(w_I) > 1e-8:
                        qc.append(create_G_gate(f_I_up), [i, j])
                    f_R_dn = -parity_dn * theta_R if is_pbc_bond else theta_R
                    f_I_dn = -parity_dn * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        qc.append(create_R_gate(f_R_dn), [L + i, L + j])
                    if abs(w_I) > 1e-8:
                        qc.append(create_G_gate(f_I_dn), [L + i, L + j])

            apply_hopping_layer(0, v)
            apply_hopping_layer(1, w)

            if U_A != 0 or U_B != 0:
                ramp = s if ramp_U else 1.0
                for i in range(L):
                    current_U = U_A if i % 2 == 0 else U_B
                    if current_U != 0:
                        theta_U = tau * ramp * current_U
                        qc.append(create_CP_gate(theta_U), [i, L + i])

    return qc
