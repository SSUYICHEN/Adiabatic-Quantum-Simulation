"""sshh.py - The paper's spinful SSH-Hubbard model (paper Eq. 9-15).

Qubit layout (spin-block mapping, paper Eq. 36)
-----------------------------------------------
L = 2 * N_cells  spatial orbitals per spin.
Qubits [0, L)        -> spin-up   block
Qubits [L, 2L)       -> spin-down block
Within a block, even index = A sublattice, odd index = B sublattice.

Hopping amplitudes follow the standard SSH naming:
    v = intra-cell hopping (even bonds)
    w = inter-cell hopping (odd bonds)

Convention: the single-particle matrix carries H[i, j] = -t, the Hermitian
conjugate of the paper's `v b^dag a` form, so code v, w = conj(paper v, w)
(.claude/rules/models-and-observables.md). Moot for real hoppings.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


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
