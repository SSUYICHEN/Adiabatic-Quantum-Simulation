"""spinless.py - Extended spinless SSH model with nearest-neighbor interaction.

Implements the adiabatic quantum simulation of
``.claude/specs/reconstruct_goals.md`` (the "spec"):

    H = H_SSH + H_NN,   ramped as H(t) = H_SSH + (t/T) H_NN,

on 2N qubits with the layout q_{2j} = A_j, q_{2j+1} = B_j (spec sec. 1).

Normative conventions
---------------------
The spec's *Pauli-operator representation* (sec. 2) and Trotter angle table
(sec. 3.2) are the normative sources (see .claude/rules/gates-and-circuits.md);
they are mutually consistent and are implemented literally here.

The spec's fermionic display (`+v b^dag a` etc.) is NOT consistent with its own
Pauli representation (the same cosmetic tension exists between Eq. 12 and
Eq. 30 of the paper). The fermionic Hamiltonian actually equivalent to the
normative Pauli form is

    intracell:            -v       b_j^dag a_j       + h.c.
    intercell and wrap:   -conj(w) b_j^dag a_{j+1}   + h.c.
    interaction:          +V_v n_Aj n_Bj  +V_w n_Bj n_Aj+1   (diagonal)

For real v, w (all planned numerics) the conjugation is invisible; for complex
w it is load-bearing and is pinned by the H_eff round-trip tests in
tests/test_spinless_model.py. Do not "simplify" the conjugation away.

PBC decisions (spec-silent, recorded here):
* The wrap hopping bond carries the Jordan-Wigner string; within the fixed
  particle-number sector it reduces to the parity factor -(-1)^{n_particles}
  on the R angle (and the opposite sign on the G angle, because the wrap
  bond's creation operator sits on the *higher* site index, unlike the bulk
  intercell bonds). Pinned by the fixed-sector wrap test.
* The wrap interaction V_w n_{B,N-1} n_{A,0} is included under PBC (diagonal,
  string-free, NO parity factor) so the interacting ring stays translation
  invariant.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable, Sequence

import numpy as np
import openfermion
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from scipy.sparse.linalg import expm_multiply

from aqs.core import (
    create_CP_gate,
    create_G_gate,
    create_R_gate,
    _givens_instruction,
)
from aqs.observables import _bit, _nonzero_probs, density_profile_qubits


# =====================================================================
# 1. Model definition
# =====================================================================
@dataclass
class SpinlessSSHModel:
    """Extended spinless SSH + NN interaction on N unit cells (spec sec. 1-2).

    Parameters
    ----------
    N_cells     : number of unit cells; 2*N_cells qubits/sites.
    v, w        : intracell / intercell hopping amplitudes.
    V_v, V_w    : intracell / intercell nearest-neighbor interaction strengths.
    PBC         : periodic (True) or open (False) boundary conditions.
    n_particles : fermion count; defaults to N (PBC, half filling) or
                  N+1 (OBC, edge modes occupied) per spec sec. 4.
    """

    N_cells: int
    v: complex
    w: complex
    V_v: float = 0.0
    V_w: float = 0.0
    PBC: bool = True
    n_particles: int | None = None
    n_sites: int = field(init=False)

    def __post_init__(self) -> None:
        if self.N_cells < 2:
            raise ValueError("N_cells must be >= 2")
        self.n_sites = 2 * self.N_cells
        if self.n_particles is None:
            self.n_particles = self.N_cells if self.PBC else self.N_cells + 1
        if not 0 < self.n_particles <= self.n_sites:
            raise ValueError(
                f"n_particles must be in (0, {self.n_sites}], got {self.n_particles}"
            )

    def hopping_bonds(self) -> list[tuple[int, int, complex, bool]]:
        """Bonds as (low_qubit, high_qubit, coupling, is_wrap).

        Bulk bonds are (2j, 2j+1, v) and (2j+1, 2j+2, w); under PBC the wrap
        bond is (0, 2N-1, w) with is_wrap=True.
        """
        bonds: list[tuple[int, int, complex, bool]] = []
        for j in range(self.N_cells):
            bonds.append((2 * j, 2 * j + 1, self.v, False))
        for j in range(self.N_cells - 1):
            bonds.append((2 * j + 1, 2 * j + 2, self.w, False))
        if self.PBC:
            bonds.append((0, self.n_sites - 1, self.w, True))
        return bonds

    def interaction_bonds(self) -> list[tuple[int, int, float]]:
        """Diagonal NN bonds as (qubit_a, qubit_b, strength). Never carry parity."""
        bonds = [(2 * j, 2 * j + 1, self.V_v) for j in range(self.N_cells)]
        bonds += [(2 * j + 1, 2 * j + 2, self.V_w) for j in range(self.N_cells - 1)]
        if self.PBC:
            bonds.append((self.n_sites - 1, 0, self.V_w))
        return bonds

    def single_particle_hamiltonian(self) -> np.ndarray:
        """2N x 2N matrix M with H_hopping = sum_mn M[m, n] c_m^dag c_n.

        Derived from the normative Pauli representation (module docstring):
        M[B_j, A_j] = -v, M[B_j, A_{j+1}] = -conj(w) (wrap included under PBC).
        Guarded by the one-particle-sector spectrum test.
        """
        M = np.zeros((self.n_sites, self.n_sites), dtype=complex)
        for j in range(self.N_cells):
            a, b = 2 * j, 2 * j + 1
            M[b, a] = -self.v
            M[a, b] = -np.conj(self.v)
        last = self.N_cells - 1 if self.PBC else self.N_cells - 2
        for j in range(last + 1):
            b, a_next = 2 * j + 1, (2 * j + 2) % self.n_sites
            M[b, a_next] = -np.conj(self.w)
            M[a_next, b] = -self.w
        return M

    def slater_Q_matrix(self) -> np.ndarray:
        """Rows = the n_particles lowest single-particle orbitals (site basis)."""
        eigenvalues, eigenvectors = np.linalg.eigh(self.single_particle_hamiltonian())
        order = np.argsort(eigenvalues)
        return eigenvectors[:, order][:, : self.n_particles].T


# =====================================================================
# 2. Reference Hamiltonians (little-endian, qiskit convention)
# =====================================================================
def _label(n_qubits: int, placements: dict[int, str]) -> str:
    """Full-length qiskit Pauli label with explicit little-endian placement.

    qiskit label strings put qubit n-1 leftmost; building the label from an
    explicit {qubit: char} map keeps the bit order impossible to get wrong.
    """
    chars = ["I"] * n_qubits
    for q, c in placements.items():
        chars[q] = c
    return "".join(reversed(chars))


def spec_pauli_hamiltonian(model: SpinlessSSHModel) -> SparsePauliOp:
    """The normative Pauli Hamiltonian of spec sec. 2 (plus PBC wrap terms).

    Bulk bonds (lower qubit lo, higher qubit hi, coupling t):
        -Re(t)/2 (X_lo X_hi + Y_lo Y_hi) - Im(t)/2 (X_lo Y_hi - Y_lo X_hi)
    PBC wrap bond (fermionic -conj(w) c_{2N-1}^dag c_0 + h.c., creation on the
    higher site) with its Jordan-Wigner Z-string on qubits 1..2N-2:
        [-Re(w)/2 (X_0 X_m + Y_0 Y_m) + Im(w)/2 (X_0 Y_m - Y_0 X_m)] * string
    Interaction: V n_a n_b = V/4 (I - Z_a)(I - Z_b), string-free.
    """
    n = model.n_sites
    terms: list[tuple[str, complex]] = []
    for lo, hi, t, is_wrap in model.hopping_bonds():
        string: dict[int, str] = (
            {k: "Z" for k in range(lo + 1, hi)} if is_wrap else {}
        )
        im_sign = +1.0 if is_wrap else -1.0
        terms += [
            (_label(n, {lo: "X", hi: "X", **string}), -np.real(t) / 2),
            (_label(n, {lo: "Y", hi: "Y", **string}), -np.real(t) / 2),
            (_label(n, {lo: "X", hi: "Y", **string}), im_sign * np.imag(t) / 2),
            (_label(n, {lo: "Y", hi: "X", **string}), -im_sign * np.imag(t) / 2),
        ]
    for qa, qb, V in model.interaction_bonds():
        if V == 0.0:
            continue
        terms += [
            (_label(n, {}), V / 4),
            (_label(n, {qa: "Z"}), -V / 4),
            (_label(n, {qb: "Z"}), -V / 4),
            (_label(n, {qa: "Z", qb: "Z"}), V / 4),
        ]
    return SparsePauliOp.from_list(terms).simplify()


def exact_hamiltonian(model: SpinlessSSHModel) -> np.ndarray:
    """Dense little-endian many-body Hamiltonian (independent test reference)."""
    return spec_pauli_hamiltonian(model).to_matrix()


def _sector_indices(n_qubits: int, n_particles: int) -> np.ndarray:
    return np.array(
        [i for i in range(2**n_qubits) if bin(i).count("1") == n_particles]
    )


def exact_ground_state_fixed_n(model: SpinlessSSHModel) -> np.ndarray:
    """Ground state of H restricted to the model's particle-number sector.

    The restriction is essential: the circuits conserve particle number, and
    the *global* interacting ground state may live at a different filling
    (see .claude/rules/models-and-observables.md).
    """
    H = exact_hamiltonian(model)
    idx = _sector_indices(model.n_sites, model.n_particles)
    _, vecs = np.linalg.eigh(H[np.ix_(idx, idx)])
    psi = np.zeros(2**model.n_sites, dtype=complex)
    psi[idx] = vecs[:, 0]
    return psi


# =====================================================================
# 3. Simulation driver (spec sec. 6)
# =====================================================================
class SpinlessSSHHSim:
    """Adiabatic quantum simulation driver for the spinless model (spec sec. 6).

    Builds the Givens-rotation initial state (spec sec. 4), the Trotterized
    adiabatic evolution (spec sec. 3.2), and the observable estimators
    (spec sec. 5), all in qiskit's little-endian convention.
    """

    def __init__(self, model: SpinlessSSHModel, T: float = 1.0, L: int = 40):
        if L < 1:
            raise ValueError("L must be >= 1")
        self.model = model
        self.T = float(T)
        self.L = int(L)

    # -- circuits ------------------------------------------------------
    def preparation_circuit(self) -> QuantumCircuit:
        """Slater-determinant ground state of H_SSH via OpenFermion Givens.

        X gates set the occupation pattern, then the Givens network rotates to
        the lowest-orbital determinant. Guarded by the lowest-fill-energy test
        (the historical failure mode prepares the highest band at exactly -E).
        """
        model = self.model
        qc = QuantumCircuit(model.n_sites)
        for q in range(model.n_particles):
            qc.x(q)
        circuit_description = openfermion.circuits.slater_determinant_preparation_circuit(
            model.slater_Q_matrix()
        )
        for parallel_ops in circuit_description:
            for op in parallel_ops:
                if not (isinstance(op, tuple) and len(op) == 4):
                    raise NotImplementedError(
                        f"unsupported Givens operation from OpenFermion: {op!r}"
                    )
                j, k, theta, phi = op
                qc.append(_givens_instruction(theta, phi), [j, k])
        return qc

    def trotter_circuit(self) -> QuantumCircuit:
        """L Trotter layers of the lambda-ramped evolution (spec sec. 3.2).

        Hopping angles are constant, theta_R = -dt Re(t), theta_G = -dt Im(t),
        gates on (lower, higher) qubits. The PBC wrap bond multiplies theta_R
        by p = -(-1)^{n_particles} and theta_G by -p (creation on the higher
        site flips the G sign relative to bulk intercell bonds; see module
        docstring). Interaction CP angles are dt * V * (2l-1)/(2L), parity-free.
        """
        model = self.model
        dt = self.T / self.L
        parity_factor = -((-1) ** model.n_particles)
        qc = QuantumCircuit(model.n_sites)
        for step in range(1, self.L + 1):
            lam = (2 * step - 1) / (2 * self.L)
            for lo, hi, t, is_wrap in model.hopping_bonds():
                theta_R = -dt * float(np.real(t))
                theta_G = -dt * float(np.imag(t))
                if is_wrap:
                    theta_R *= parity_factor
                    theta_G *= -parity_factor
                if abs(np.real(t)) > 1e-12:
                    qc.append(create_R_gate(theta_R), [lo, hi])
                if abs(np.imag(t)) > 1e-12:
                    qc.append(create_G_gate(theta_G), [lo, hi])
            for qa, qb, V in model.interaction_bonds():
                if abs(V) > 1e-12:
                    qc.append(create_CP_gate(dt * V * lam), [qa, qb])
        return qc

    def full_circuit(self) -> QuantumCircuit:
        """Initial-state preparation followed by the adiabatic evolution."""
        qc = self.preparation_circuit()
        qc.barrier()
        return qc.compose(self.trotter_circuit())

    def evolve(self) -> np.ndarray:
        """Statevector after preparation + adiabatic evolution."""
        return Statevector(self.full_circuit()).data

    # -- observables (spec sec. 5) ------------------------------------
    def berry_phase(self, statevector: np.ndarray) -> float:
        """Many-body Berry phase gamma = Im ln <exp(i 2 pi/N X)> (spec 5.2).

        X = sum_q floor(q/2) n_q. PBC only; under OBC the phase is not a
        topological indicator. Absolute values carry a constant filling
        offset -- only the pi jump across the transition is physical.
        """
        if not self.model.PBC:
            raise ValueError("the Berry phase is defined for PBC models only")
        z = 0.0 + 0.0j
        for idx, p in zip(*_nonzero_probs(np.asarray(statevector))):
            X_b = sum(
                (q // 2) * _bit(idx, q) for q in range(self.model.n_sites)
            )
            z += p * np.exp(2j * np.pi * X_b / self.model.N_cells)
        return float(np.angle(z))

    def polarization_profile(self, statevector: np.ndarray) -> np.ndarray:
        """Sublattice polarization P_j = <n_{2j} - n_{2j+1}> (spec 5.1).

        OBC only: under PBC the polarization is identically zero and a
        "result" computed there is meaningless.
        """
        if self.model.PBC:
            raise ValueError("the sublattice polarization requires an OBC model")
        density = density_profile_qubits(
            np.asarray(statevector), self.model.n_sites
        )
        return np.array(
            [density[2 * j] - density[2 * j + 1] for j in range(self.model.N_cells)]
        )

    # -- spec sec. 6 drivers ------------------------------------------
    def trotter_benchmark(
        self, L_values: Sequence[int]
    ) -> list[dict[str, float]]:
        """Fidelity vs L at fixed T (spec sec. 6.3, paper Fig. 2 protocol).

        Noninteracting models compare against the exact final state; with
        interaction the reference is the previous L's final state, which
        diagnoses stability under refinement only, not absolute accuracy.
        """
        model = self.model
        records: list[dict[str, float]] = []
        psi0 = Statevector(self.preparation_circuit()).data
        interacting = model.V_v != 0.0 or model.V_w != 0.0
        if not interacting:
            H_sparse = spec_pauli_hamiltonian(model).to_matrix(sparse=True)
            reference = expm_multiply(-1j * self.T * H_sparse, psi0)
        else:
            reference = psi0
        for L in L_values:
            sim = SpinlessSSHHSim(model, T=self.T, L=L)
            evolved = Statevector(sim.full_circuit()).data
            records.append(
                {"L": L, "fidelity": float(abs(np.vdot(reference, evolved)) ** 2)}
            )
            if interacting:
                reference = evolved
        return records

    def interaction_sweep(
        self, V_points: Iterable[tuple[float, float]]
    ) -> list[dict[str, float]]:
        """Sweep (V_v, V_w) and record the topological indicator (spec 6.4).

        PBC models record the Berry phase; OBC models record the polarization
        profile (edge cell first).
        """
        records: list[dict[str, float]] = []
        for V_v, V_w in V_points:
            model = replace(self.model, V_v=V_v, V_w=V_w)
            sim = SpinlessSSHHSim(model, T=self.T, L=self.L)
            psi = sim.evolve()
            record: dict[str, float] = {"V_v": V_v, "V_w": V_w}
            if model.PBC:
                record["berry_phase"] = sim.berry_phase(psi)
            else:
                profile = sim.polarization_profile(psi)
                record["edge_polarization"] = float(profile[0])
                record["polarization"] = profile.tolist()
            records.append(record)
        return records
