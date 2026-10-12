"""hubbard.py - 1D Hubbard chain: the spinful SSHH model without sublattices.

    H = sum_{i, s} ( -t c_{i,s}^dag c_{i+1,s} + h.c. ) + U sum_i n_{i,up} n_{i,dn}

This is the paper's SSH-Hubbard model (``aqs.core.SSHHModel``) with the A/B
distinction switched off: v = w = t and U_A = U_B = U. The adiabatic circuit is
therefore ``aqs.core.build_annealing_circuit`` unchanged; this module adds the
model definition, an exact reference, and the Hubbard observables.

Layout (same as the spinful SSHH model): n_sites sites on 2 * n_sites qubits,
site i spin-up on qubit i, spin-down on qubit n_sites + i, qiskit little-endian.

Conventions
-----------
* Single-particle matrix M[i, i+1] = -t, i.e. the term -t c_i^dag c_{i+1},
  identical to ``SSHHModel.H``. For complex t this fixes which direction
  carries the phase; pinned against OpenFermion in tests/test_hubbard_model.py.
* Closed shells only. The circuit starts from the U = 0 Slater determinant,
  which is unique only when the highest filled and lowest empty single-particle
  levels differ. A ring of 4m sites at half filling is an open shell (use a
  ring of 4m + 2 sites, or open boundaries).
* There is no sublattice, hence no SSH topological phase: the Berry-phase jump
  and edge polarization of the SSH models do not apply here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector

from aqs.models.sshh import SSHHModel
from aqs.circuits.trotter import build_annealing_circuit
from aqs.observables import _bit, _nonzero_probs


# =====================================================================
# 1. Model definition
# =====================================================================
@dataclass
class HubbardModel:
    """Uniform Hubbard chain on n_sites sites (2 * n_sites qubits).

    Parameters
    ----------
    n_sites        : number of lattice sites; must be even (two per SSHH cell).
    t              : hopping amplitude.
    U              : on-site interaction between the two spin species.
    PBC            : periodic (True) or open (False) boundary conditions.
    num_up, num_dn : electrons per spin; default n_sites // 2 each (half filling).
    """

    n_sites: int
    t: complex = 1.0
    U: float = 0.0
    PBC: bool = True
    num_up: int | None = None
    num_dn: int | None = None

    def __post_init__(self) -> None:
        if self.n_sites < 2 or self.n_sites % 2:
            raise ValueError(
                f"n_sites must be an even number >= 2, got {self.n_sites}")
        if self.num_up is None:
            self.num_up = self.n_sites // 2
        if self.num_dn is None:
            self.num_dn = self.n_sites // 2
        for name, n in (("num_up", self.num_up), ("num_dn", self.num_dn)):
            if not 0 < n <= self.n_sites:
                raise ValueError(f"{name} must be in (0, {self.n_sites}], got {n}")

    def sshh(self) -> SSHHModel:
        """The equivalent spinful SSHH model: v = w = t."""
        return SSHHModel(self.n_sites // 2, self.t, self.t,
                         self.num_up, self.num_dn, self.PBC)

    def bonds(self) -> list[tuple[int, int, bool]]:
        """Hopping bonds as (lower_site, higher_site, is_wrap)."""
        bonds = [(i, i + 1, False) for i in range(self.n_sites - 1)]
        if self.PBC:
            bonds.append((0, self.n_sites - 1, True))
        return bonds

    def single_particle_hamiltonian(self) -> np.ndarray:
        """n_sites x n_sites matrix M with H_hop = sum_ij M[i, j] c_i^dag c_j (one spin)."""
        n = self.n_sites
        return self.sshh().H[:n, :n].copy()

    def shell_gap(self) -> float:
        """Smallest gap between the highest filled and lowest empty level.

        Zero means the U = 0 ground state is degenerate (open shell).
        """
        levels = np.sort(np.linalg.eigvalsh(self.single_particle_hamiltonian()))
        gaps = [levels[n] - levels[n - 1]
                for n in (self.num_up, self.num_dn) if n < self.n_sites]
        return float(min(gaps)) if gaps else float("inf")


# =====================================================================
# 2. Reference Hamiltonian (little-endian, qiskit convention)
# =====================================================================
def _label(n_qubits: int, placements: dict[int, str]) -> str:
    """Full-length qiskit Pauli label from an explicit {qubit: char} map."""
    chars = ["I"] * n_qubits
    for q, c in placements.items():
        chars[q] = c
    return "".join(reversed(chars))


def hubbard_pauli_hamiltonian(model: HubbardModel) -> SparsePauliOp:
    """Jordan-Wigner Pauli form of the Hubbard Hamiltonian.

    Bulk bond (lo, hi = lo + 1), term -t c_lo^dag c_hi + h.c.:
        -Re(t)/2 (X_lo X_hi + Y_lo Y_hi) + Im(t)/2 (X_lo Y_hi - Y_lo X_hi)
    PBC wrap bond (lo = 0, hi = n - 1), term -t c_hi^dag c_lo + h.c. (creation
    on the higher site, so the Im sign flips), times the Jordan-Wigner Z-string
    on the sites strictly between, within the same spin block:
        [-Re(t)/2 (X_lo X_hi + Y_lo Y_hi) - Im(t)/2 (X_lo Y_hi - Y_lo X_hi)] * string
    Interaction: U n_up n_dn = U/4 (I - Z_i)(I - Z_{n+i}), string-free.
    """
    n = model.n_sites
    nq = 2 * n
    re, im = float(np.real(model.t)), float(np.imag(model.t))
    terms: list[tuple[str, complex]] = []
    for offset in (0, n):
        for lo, hi, is_wrap in model.bonds():
            a, b = offset + lo, offset + hi
            string = {k: "Z" for k in range(a + 1, b)} if is_wrap else {}
            im_sign = -1.0 if is_wrap else +1.0
            terms += [
                (_label(nq, {a: "X", b: "X", **string}), -re / 2),
                (_label(nq, {a: "Y", b: "Y", **string}), -re / 2),
                (_label(nq, {a: "X", b: "Y", **string}), im_sign * im / 2),
                (_label(nq, {a: "Y", b: "X", **string}), -im_sign * im / 2),
            ]
    if model.U != 0.0:
        for i in range(n):
            terms += [
                (_label(nq, {}), model.U / 4),
                (_label(nq, {i: "Z"}), -model.U / 4),
                (_label(nq, {n + i: "Z"}), -model.U / 4),
                (_label(nq, {i: "Z", n + i: "Z"}), model.U / 4),
            ]
    return SparsePauliOp.from_list(terms).simplify()


def hubbard_hamiltonian(model: HubbardModel):
    """Sparse little-endian many-body Hamiltonian (2**(2 n_sites) square)."""
    return hubbard_pauli_hamiltonian(model).to_matrix(sparse=True).tocsr()


def sector_indices(model: HubbardModel) -> np.ndarray:
    """Basis indices with exactly num_up up-spins and num_dn down-spins."""
    n = model.n_sites
    idx = np.arange(1 << (2 * n), dtype=np.int64)
    up = np.zeros_like(idx)
    dn = np.zeros_like(idx)
    for q in range(n):
        up += (idx >> q) & 1
        dn += (idx >> (n + q)) & 1
    return idx[(up == model.num_up) & (dn == model.num_dn)]


def hubbard_ground_state_fixed_n(model: HubbardModel) -> np.ndarray:
    """Ground state of H restricted to the model's (num_up, num_dn) sector.

    The circuit conserves both spin counts, so this is the state it can reach;
    the global ground state may sit at a different filling.
    """
    idx = sector_indices(model)
    H = hubbard_hamiltonian(model)[idx][:, idx].toarray()
    _, vecs = np.linalg.eigh(H)
    psi = np.zeros(1 << (2 * model.n_sites), dtype=complex)
    psi[idx] = vecs[:, 0]
    return psi


# =====================================================================
# 3. Simulation driver
# =====================================================================
class HubbardSim:
    """Adiabatic preparation of the Hubbard ground state.

    Stage 1 prepares the U = 0 Slater determinant; stage 2 applies L Trotter
    layers over total time T while U ramps linearly from 0 to its target.
    """

    def __init__(self, model: HubbardModel, T: float = 1.0, L: int = 40):
        if L < 1:
            raise ValueError("L must be >= 1")
        if model.shell_gap() < 1e-9:
            raise ValueError(
                "open shell: the U = 0 ground state is degenerate for "
                f"n_sites={model.n_sites}, num_up={model.num_up}, "
                f"num_dn={model.num_dn}, PBC={model.PBC}, so the adiabatic "
                "starting state is not defined. Use a ring of 4m + 2 sites or "
                "open boundaries.")
        self.model = model
        self.T = float(T)
        self.L = int(L)
        self._H = None

    # -- circuits ------------------------------------------------------
    def _circuit(self, steps: int) -> QuantumCircuit:
        model = self.model
        sshh = model.sshh()
        Q_up, Q_dn, _, _ = sshh.slater_Q_matrices()
        return build_annealing_circuit(
            sshh, Q_up, Q_dn, U_A=model.U, U_B=model.U, T_A=self.T, steps=steps)

    def preparation_circuit(self) -> QuantumCircuit:
        """Stage 1 only: the U = 0 Slater determinant."""
        return self._circuit(0)

    def full_circuit(self) -> QuantumCircuit:
        """Stage 1 followed by the L-layer U ramp."""
        return self._circuit(self.L)

    def evolve(self, backend=None) -> np.ndarray:
        """Statevector after preparation + adiabatic evolution.

        backend: None for qiskit's Statevector, or a name / instance from
        aqs.backends (e.g. "mlx" for the Apple-silicon GPU).
        """
        if backend is None:
            return Statevector(self.full_circuit()).data
        from aqs.backends import get_backend
        bk = backend if hasattr(backend, "statevector") else get_backend(backend)
        return bk.statevector(self.full_circuit())

    # -- observables ---------------------------------------------------
    def energy(self, statevector: np.ndarray) -> float:
        """<H> at the model's target U."""
        if self._H is None:
            self._H = hubbard_hamiltonian(self.model)
        psi = np.asarray(statevector)
        return float(np.real(np.vdot(psi, self._H @ psi)))

    def _occupations(self, statevector: np.ndarray):
        """(probabilities, n_up[site, state], n_dn[site, state]) on the support."""
        n = self.model.n_sites
        idx, p = _nonzero_probs(np.asarray(statevector))
        n_up = np.array([_bit(idx, i) for i in range(n)], dtype=float)
        n_dn = np.array([_bit(idx, n + i) for i in range(n)], dtype=float)
        return p, n_up, n_dn

    def double_occupancy(self, statevector: np.ndarray) -> float:
        """D = (1 / n_sites) sum_i <n_{i,up} n_{i,dn}>."""
        p, n_up, n_dn = self._occupations(statevector)
        return float(np.sum(p * (n_up * n_dn)) / self.model.n_sites)

    def spin_correlation(self, statevector: np.ndarray) -> np.ndarray:
        """C[i, j] = <S^z_i S^z_j>, S^z_i = (n_{i,up} - n_{i,dn}) / 2."""
        p, n_up, n_dn = self._occupations(statevector)
        sz = 0.5 * (n_up - n_dn)
        return (sz * p) @ sz.T

    def spin_structure_factor(self, statevector: np.ndarray, q: float = np.pi) -> float:
        """S(q) = (1 / n_sites) sum_ij exp(i q (i - j)) <S^z_i S^z_j>."""
        C = self.spin_correlation(statevector)
        phase = np.exp(1j * q * np.arange(self.model.n_sites))
        return float(np.real(np.conj(phase) @ C @ phase) / self.model.n_sites)
