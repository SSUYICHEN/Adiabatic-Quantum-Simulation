"""
core.py - Physics core for the SSH-Hubbard adiabatic quantum simulation.

This module is the single source of truth for:
  * the SSH-Hubbard one-body model and Slater-determinant ground state,
  * the modular two-qubit gates (R, G, CP, Givens),
  * the parity-corrected (PBC) / open (OBC) Trotterized annealing circuit.

The numerics here are kept byte-for-byte equivalent to the original
standalone scripts (TrotterforPBC.py / deltaU.py / halffillingedgewithdeltaU.py)
so that the packaged tool reproduces the published results exactly.

Qubit layout (spin-block mapping)
---------------------------------
L = 2 * N_cells  spatial orbitals per spin.
Qubits [0, L)        -> spin-up   block
Qubits [L, 2L)       -> spin-down block
Within a block, even index = A sublattice, odd index = B sublattice.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import openfermion
from qiskit import QuantumCircuit


# =====================================================================
# 1. Physical model (spin-block ordering)
# =====================================================================
@dataclass
class SSHHModel:
    """One-body SSH-Hubbard formulation on N_cells unit cells.

    Parameters
    ----------
    N_cells : number of unit cells (each cell = one A + one B site).
    t1, t2  : intra-cell (v) and inter-cell (w) hopping amplitudes.
    num_up, num_dn : electron counts per spin sector.
    PBC     : periodic (True) or open (False) boundary conditions.
    """

    N_cells: int
    t1: complex
    t2: complex
    num_up: int
    num_dn: int
    PBC: bool = True

    def __post_init__(self):
        self.L = 2 * self.N_cells
        self.number_of_electrons = self.num_up + self.num_dn
        self.H = self._build_hamiltonian()

    # -- constructors -------------------------------------------------
    @classmethod
    def from_total_electrons(cls, N_cells, t1, t2, number_of_electrons, PBC=True):
        """Split a total electron count into (up, dn) using the standard rule
        num_up = ceil(N/2), num_dn = floor(N/2)."""
        num_up = (number_of_electrons // 2) + (number_of_electrons % 2)
        num_dn = number_of_electrons // 2
        return cls(N_cells, t1, t2, num_up, num_dn, PBC)

    # -- hamiltonian --------------------------------------------------
    def _build_hamiltonian(self) -> np.ndarray:
        dim = 2 * self.L
        H = np.zeros((dim, dim), dtype=complex)
        for i in range(self.L):
            if not self.PBC and i == self.L - 1:
                continue
            j = (i + 1) % self.L
            t = self.t1 if i % 2 == 0 else self.t2
            # Up-spin block
            H[i, j] = -t
            H[j, i] = -np.conj(t)
            # Down-spin block
            H[self.L + i, self.L + j] = -t
            H[self.L + j, self.L + i] = -np.conj(t)
        return H

    # -- Slater determinant -------------------------------------------
    def slater_Q_matrices(self):
        """Return (Q_up, Q_dn, num_up, num_dn): the occupied single-particle
        orbitals (lowest eigenvectors) used to prepare the Slater determinant."""
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
    """exp(-i theta/2 (XX + YY)) : real (Hermitian) hopping term."""
    qc = QuantumCircuit(2, name=f"R({theta:.3f})")
    qc.rxx(theta / 2.0, 0, 1)
    qc.ryy(theta / 2.0, 0, 1)
    return qc.to_instruction()


def create_G_gate(theta):
    """exp(-i theta/2 (YX - XY)) : imaginary hopping term."""
    qc = QuantumCircuit(2, name=f"G({theta:.3f})")
    qc.sdg(0); qc.rxx(theta / 2.0, 0, 1); qc.s(0)
    qc.sdg(1); qc.rxx(-theta / 2.0, 0, 1); qc.s(1)
    return qc.to_instruction()


def create_CP_gate(theta):
    """Controlled-phase : the on-site Hubbard interaction n_up n_dn."""
    qc = QuantumCircuit(2, name=f"CP({theta:.3f})")
    qc.cp(theta, 0, 1)
    return qc.to_instruction()


def _givens_instruction(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1)
    gv.rz(-phi, 0)
    gv.cx(0, 1)
    gv.cry(-2.0 * theta, 1, 0)
    gv.cx(0, 1)
    return gv.to_instruction()


# =====================================================================
# 3. Annealing circuit (parity-corrected PBC, staggered Hubbard U)
# =====================================================================
def build_annealing_circuit(
    model: SSHHModel,
    Q_up,
    Q_dn,
    U_A=0.0,
    U_B=0.0,
    T_A=1.0,
    steps=0,
    ramp_U=True,
):
    """Build the Trotterized adiabatic-evolution circuit.

    Parameters
    ----------
    model        : SSHHModel providing N_cells, t1, t2, boundary, occupations.
    Q_up, Q_dn   : occupied orbitals for Slater-determinant state preparation.
    U_A, U_B     : Hubbard interaction on the A (even) and B (odd) sublattices.
                   Pass U_A == U_B for a uniform interaction.
    T_A          : total annealing time.
    steps        : number of Trotter steps (0 = state preparation only).
    ramp_U       : if True the interaction is ramped adiabatically with the
                   schedule s = (step-0.5)/steps (used for measurements);
                   if False the interaction is held constant (used for the
                   Trotter-convergence fidelity check).

    Returns
    -------
    qiskit.QuantumCircuit on 2L qubits.
    """
    N_cells = model.N_cells
    t1, t2 = model.t1, model.t2
    PBC = model.PBC
    num_up, num_dn = model.num_up, model.num_dn
    L = 2 * N_cells
    num_qubits = 2 * L

    qc = QuantumCircuit(num_qubits)

    # 1. Occupation: fill the lowest orbitals.
    for i in range(num_up):
        qc.x(i)
    for i in range(num_dn):
        qc.x(L + i)

    # 2. Slater-determinant preparation (Givens-rotation network).
    circ_up = openfermion.circuits.slater_determinant_preparation_circuit(Q_up)
    circ_dn = openfermion.circuits.slater_determinant_preparation_circuit(Q_dn)
    for parallel_ops in circ_up:
        for j, k, theta, phi in parallel_ops:
            qc.append(_givens_instruction(theta, phi), [j, k])
    for parallel_ops in circ_dn:
        for j, k, theta, phi in parallel_ops:
            qc.append(_givens_instruction(theta, phi), [j + L, k + L])
    qc.barrier()

    # 3. Trotterized evolution.
    if steps > 0:
        tau = T_A / steps
        parity_up = (-1) ** num_up
        parity_dn = (-1) ** num_dn

        for step in range(1, steps + 1):
            s = (step - 0.5) / steps  # adiabatic schedule for the interaction

            def apply_hopping_layer(start_i, t_val):
                for i in range(start_i, L, 2):
                    is_pbc_bond = i == L - 1
                    if not PBC and is_pbc_bond:
                        continue
                    j = (i + 1) % L

                    w_R, w_I = float(np.real(t_val)), float(np.imag(t_val))
                    theta_R = -2.0 * tau * w_R
                    theta_I = -2.0 * tau * w_I

                    # Up-spin block (parity correction on the wrap-around bond)
                    f_R_up = -parity_up * theta_R if is_pbc_bond else theta_R
                    f_I_up = -parity_up * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        qc.append(create_R_gate(f_R_up), [i, j])
                    if abs(w_I) > 1e-8:
                        qc.append(create_G_gate(f_I_up), [i, j])

                    # Down-spin block
                    f_R_dn = -parity_dn * theta_R if is_pbc_bond else theta_R
                    f_I_dn = -parity_dn * theta_I if is_pbc_bond else theta_I
                    if abs(w_R) > 1e-8:
                        qc.append(create_R_gate(f_R_dn), [L + i, L + j])
                    if abs(w_I) > 1e-8:
                        qc.append(create_G_gate(f_I_dn), [L + i, L + j])

            apply_hopping_layer(0, t1)
            apply_hopping_layer(1, t2)

            # On-site Hubbard interaction (staggered U_A / U_B).
            if U_A != 0 or U_B != 0:
                ramp = s if ramp_U else 1.0
                for i in range(L):
                    current_U = U_A if i % 2 == 0 else U_B
                    if current_U != 0:
                        theta_U = -tau * ramp * current_U
                        qc.append(create_CP_gate(theta_U), [i, L + i])

    return qc
