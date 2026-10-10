"""spin_ssh.py - SSH-Hubbard chain with the spin as the sublattice.

Implements the adiabatic quantum simulation of
``.claude/specs/reconstruct_goals.md`` (revision 2026-10-06, the "spec"):

    H = H_spin + H_Hubbard,   ramped as H(t) = H_spin + (t/T) H_Hubbard,

    H_spin    = sum_j ( v c_{j,up}^dag c_{j,dn} + h.c. )
              + sum_j ( w c_{j,dn}^dag c_{j+1,up} + h.c. )   (+ PBC term)
    H_Hubbard = U sum_j n_{j,up} n_{j,dn}

on 2N qubits with q_{2j} = (j, up), q_{2j+1} = (j, dn) (spec sec. 1). Each
site's two spin modes play the role of the SSH sublattices A (up) and B (dn):
v is the SSH intracell hopping (here an on-site spin flip), w the intercell
hopping, U the intracell interaction. There is no spin-conserving hopping,
so only the total particle number is conserved.

Normative conventions
---------------------
The fermionic Hamiltonian of spec sec. 2.1 is normative; its Jordan-Wigner
image (spec sec. 2.2) is what `pauli_hamiltonian` builds and what the Trotter
angles follow:

    t c_a^dag c_b + h.c.  (a < b neighbours)
        = Re(t)/2 (X_a X_b + Y_a Y_b) - Im(t)/2 (X_a Y_b - Y_a X_b)
    exp(-i dt H_t) = R(dt Re t) G(-dt Im t)       (exact for real t)

The previous revision of the spec used the gauge image v -> -v, w -> -w
(fermionic form -v b^dag a + h.c.); the two differ by c_{j,dn} -> -c_{j,dn}
and give identical spectra and observables. Angle signs here follow the
current spec and are pinned by the H_eff round-trip tests in
tests/test_spin_ssh_model.py with complex v, w. Do not "simplify" them.

PBC decisions (recorded here):
* The wrap term w c_{N-1,dn}^dag c_{0,up} + h.c. creates on the HIGHER mode,
  so its imaginary part has the opposite sign to the bulk form, and it
  carries the Jordan-Wigner string Z_1 ... Z_{2N-2}. Inside the sector with
  n particles the string is the number p = -(-1)^n, applied to the R angle;
  the G angle gets -p (sign flip from the creation side, times p). Pinned by
  the fixed-sector wrap test.
* The Hubbard term is on-site, so the ring has no interaction wrap term.
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
class SpinSSHModel:
    """SSH-Hubbard chain with spin as the sublattice, N sites (spec sec. 1-2).

    Parameters
    ----------
    N           : number of sites; 2*N modes (qubits), two per site.
    v, w        : on-site spin-flip and intersite hopping amplitudes
                  (the SSH intracell / intercell hoppings).
    U           : on-site Hubbard repulsion between the two spin modes.
    PBC         : periodic (True) or open (False) boundary conditions.
    n_particles : fermion count; defaults to N (PBC, half filling) or
                  N+1 (OBC, edge modes occupied) per spec sec. 4.
    """

    N: int
    v: complex
    w: complex
    U: float = 0.0
    PBC: bool = True
    n_particles: int | None = None
    n_modes: int = field(init=False)

    def __post_init__(self) -> None:
        if self.N < 2:
            raise ValueError("N must be >= 2")
        self.n_modes = 2 * self.N
        if self.n_particles is None:
            self.n_particles = self.N if self.PBC else self.N + 1
        if not 0 < self.n_particles <= self.n_modes:
            raise ValueError(
                f"n_particles must be in (0, {self.n_modes}], got {self.n_particles}"
            )

    def hopping_bonds(self) -> list[tuple[int, int, complex, bool]]:
        """Bonds as (low_qubit, high_qubit, coupling, is_wrap).

        Spin flips are (2j, 2j+1, v); intersite hops (2j+1, 2j+2, w); under
        PBC the wrap bond is (0, 2N-1, w) with is_wrap=True.
        """
        bonds: list[tuple[int, int, complex, bool]] = []
        for j in range(self.N):
            bonds.append((2 * j, 2 * j + 1, self.v, False))
        for j in range(self.N - 1):
            bonds.append((2 * j + 1, 2 * j + 2, self.w, False))
        if self.PBC:
            bonds.append((0, self.n_modes - 1, self.w, True))
        return bonds

    def hubbard_pairs(self) -> list[tuple[int, int]]:
        """The (up, down) qubit pair of every site; diagonal, never carries parity."""
        return [(2 * j, 2 * j + 1) for j in range(self.N)]

    def single_particle_hamiltonian(self) -> np.ndarray:
        """2N x 2N matrix M with H_spin = sum_mn M[m, n] c_m^dag c_n (spec 2.1).

        M[2j, 2j+1] = v, M[2j+1, 2j+2] = w (mod 2N for the PBC term), plus
        Hermitian conjugates. Guarded by the one-particle-sector spectrum test.
        """
        M = np.zeros((self.n_modes, self.n_modes), dtype=complex)
        for j in range(self.N):
            up, dn = 2 * j, 2 * j + 1
            M[up, dn] = self.v
            M[dn, up] = np.conj(self.v)
        last = self.N - 1 if self.PBC else self.N - 2
        for j in range(last + 1):
            dn, up_next = 2 * j + 1, (2 * j + 2) % self.n_modes
            M[dn, up_next] = self.w
            M[up_next, dn] = np.conj(self.w)
        return M

    def slater_Q_matrix(self) -> np.ndarray:
        """Rows = the n_particles lowest single-particle orbitals (mode basis)."""
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


def pauli_hamiltonian(model: SpinSSHModel) -> SparsePauliOp:
    """The Jordan-Wigner Pauli form of H (spec sec. 2.2).

    Bulk bonds (lower qubit lo, higher qubit hi, coupling t):
        Re(t)/2 (X_lo X_hi + Y_lo Y_hi) - Im(t)/2 (X_lo Y_hi - Y_lo X_hi)
    PBC wrap bond (creation on the higher mode) with its Jordan-Wigner
    Z-string on qubits 1..2N-2:
        [Re(w)/2 (X_0 X_m + Y_0 Y_m) + Im(w)/2 (X_0 Y_m - Y_0 X_m)] * string
    Hubbard: U n_up n_dn = U/4 (I - Z_up)(I - Z_dn), string-free.
    """
    n = model.n_modes
    terms: list[tuple[str, complex]] = []
    for lo, hi, t, is_wrap in model.hopping_bonds():
        string: dict[int, str] = (
            {k: "Z" for k in range(lo + 1, hi)} if is_wrap else {}
        )
        im_sign = +1.0 if is_wrap else -1.0
        terms += [
            (_label(n, {lo: "X", hi: "X", **string}), np.real(t) / 2),
            (_label(n, {lo: "Y", hi: "Y", **string}), np.real(t) / 2),
            (_label(n, {lo: "X", hi: "Y", **string}), im_sign * np.imag(t) / 2),
            (_label(n, {lo: "Y", hi: "X", **string}), -im_sign * np.imag(t) / 2),
        ]
    if model.U != 0.0:
        for qa, qb in model.hubbard_pairs():
            terms += [
                (_label(n, {}), model.U / 4),
                (_label(n, {qa: "Z"}), -model.U / 4),
                (_label(n, {qb: "Z"}), -model.U / 4),
                (_label(n, {qa: "Z", qb: "Z"}), model.U / 4),
            ]
    return SparsePauliOp.from_list(terms).simplify()


def exact_hamiltonian(model: SpinSSHModel) -> np.ndarray:
    """Dense little-endian many-body Hamiltonian (independent test reference)."""
    return pauli_hamiltonian(model).to_matrix()


def _sector_indices(n_qubits: int, n_particles: int) -> np.ndarray:
    return np.array(
        [i for i in range(2**n_qubits) if bin(i).count("1") == n_particles]
    )


def exact_ground_state_fixed_n(model: SpinSSHModel) -> np.ndarray:
    """Ground state of H restricted to the model's particle-number sector.

    The restriction is essential: the circuits conserve particle number, and
    the *global* interacting ground state may live at a different filling
    (see .claude/rules/models-and-observables.md).
    """
    H = exact_hamiltonian(model)
    idx = _sector_indices(model.n_modes, model.n_particles)
    _, vecs = np.linalg.eigh(H[np.ix_(idx, idx)])
    psi = np.zeros(2**model.n_modes, dtype=complex)
    psi[idx] = vecs[:, 0]
    return psi


# =====================================================================
# 3. Simulation driver (spec sec. 6)
# =====================================================================
class SpinSSHHSim:
    """Adiabatic quantum simulation driver for the spin-sublattice chain.

    Builds the Givens-rotation initial state (spec sec. 4), the Trotterized
    adiabatic evolution (spec sec. 3.2), and the observable estimators
    (spec sec. 5), all in qiskit's little-endian convention.
    """

    def __init__(self, model: SpinSSHModel, T: float = 1.0, L: int = 40):
        if L < 1:
            raise ValueError("L must be >= 1")
        self.model = model
        self.T = float(T)
        self.L = int(L)

    # -- circuits ------------------------------------------------------
    def preparation_circuit(self) -> QuantumCircuit:
        """Slater-determinant ground state of H_spin via OpenFermion Givens.

        X gates set the occupation pattern, then the Givens network rotates to
        the lowest-orbital determinant. Guarded by the lowest-fill-energy test
        (the historical failure mode prepares the highest band at exactly -E).
        """
        model = self.model
        qc = QuantumCircuit(model.n_modes)
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

        Hopping angles are constant, theta_R = +dt Re(t), theta_G = -dt Im(t),
        gates on (lower, higher) qubits. The PBC wrap bond multiplies theta_R
        by p = -(-1)^{n_particles} and theta_G by -p (creation on the higher
        mode flips the G sign relative to bulk bonds; see module docstring).
        Hubbard CP angles are dt * U * (2l-1)/(2L), parity-free.
        """
        model = self.model
        dt = self.T / self.L
        parity_factor = -((-1) ** model.n_particles)
        qc = QuantumCircuit(model.n_modes)
        for step in range(1, self.L + 1):
            lam = (2 * step - 1) / (2 * self.L)
            for lo, hi, t, is_wrap in model.hopping_bonds():
                theta_R = dt * float(np.real(t))
                theta_G = -dt * float(np.imag(t))
                if is_wrap:
                    theta_R *= parity_factor
                    theta_G *= -parity_factor
                if abs(np.real(t)) > 1e-12:
                    qc.append(create_R_gate(theta_R), [lo, hi])
                if abs(np.imag(t)) > 1e-12:
                    qc.append(create_G_gate(theta_G), [lo, hi])
            if abs(model.U) > 1e-12:
                for qa, qb in model.hubbard_pairs():
                    qc.append(create_CP_gate(dt * model.U * lam), [qa, qb])
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
        """Many-body Berry phase gamma = arg <exp(i 2 pi/N X)> (spec 5.2).

        X = sum_q floor(q/2) n_q = sum_j j (n_{j,up} + n_{j,dn}). PBC only;
        under OBC the phase is not a topological indicator. Absolute values
        carry a constant filling offset -- only the pi jump is physical.
        """
        if not self.model.PBC:
            raise ValueError("the Berry phase is defined for PBC models only")
        z = 0.0 + 0.0j
        for idx, p in zip(*_nonzero_probs(np.asarray(statevector))):
            X_b = sum(
                (q // 2) * _bit(idx, q) for q in range(self.model.n_modes)
            )
            z += p * np.exp(2j * np.pi * X_b / self.model.N)
        return float(np.angle(z))

    def magnetisation_profile(self, statevector: np.ndarray) -> np.ndarray:
        """Site magnetisation m_j = <n_{j,up} - n_{j,dn}> (spec 5.1).

        This is the SSH sublattice polarization. OBC only: under PBC it is
        identically zero and a "result" computed there is meaningless.
        """
        if self.model.PBC:
            raise ValueError("the magnetisation profile requires an OBC model")
        density = density_profile_qubits(
            np.asarray(statevector), self.model.n_modes
        )
        return np.array(
            [density[2 * j] - density[2 * j + 1] for j in range(self.model.N)]
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
        interacting = model.U != 0.0
        if not interacting:
            H_sparse = pauli_hamiltonian(model).to_matrix(sparse=True)
            reference = expm_multiply(-1j * self.T * H_sparse, psi0)
        else:
            reference = psi0
        for L in L_values:
            sim = SpinSSHHSim(model, T=self.T, L=L)
            evolved = Statevector(sim.full_circuit()).data
            records.append(
                {"L": L, "fidelity": float(abs(np.vdot(reference, evolved)) ** 2)}
            )
            if interacting:
                reference = evolved
        return records

    def interaction_sweep(
        self, U_values: Iterable[float]
    ) -> list[dict[str, float]]:
        """Sweep U and record the topological indicator (spec 6.4).

        PBC models record the Berry phase; OBC models record the magnetisation
        profile (edge site first).
        """
        records: list[dict[str, float]] = []
        for U in U_values:
            model = replace(self.model, U=U)
            sim = SpinSSHHSim(model, T=self.T, L=self.L)
            psi = sim.evolve()
            record: dict[str, float] = {"U": U}
            if model.PBC:
                record["berry_phase"] = sim.berry_phase(psi)
            else:
                profile = sim.magnetisation_profile(psi)
                record["edge_magnetisation"] = float(profile[0])
                record["magnetisation"] = profile.tolist()
            records.append(record)
        return records
