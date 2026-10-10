"""trotter.py - Trotterized adiabatic evolution of the spinful SSHH model
(paper Eq. 33-35, 40-41).

Per step: hopping layers (v bonds, then w bonds) on both spin blocks, then the
Hubbard CP layer ramped at the layer midpoint s_l = (2l-1)/(2L). The PBC wrap
bond carries the Jordan-Wigner parity factor -(-1)^{n_s}; the interaction does
not (.claude/rules/gates-and-circuits.md). The spin-sublattice model keeps its
own, structurally identical layer in aqs.models.spin_ssh.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit

from aqs.circuits.gates import create_R_gate, create_G_gate, create_CP_gate
from aqs.circuits.state_prep import prepare_slater_determinant
from aqs.models.sshh import SSHHModel


def build_annealing_circuit(model: SSHHModel, Q_up, Q_dn,
                            U_A=0.0, U_B=0.0, T_A=1.0, steps=0, ramp_U=True):
    """Build the Trotterized adiabatic-evolution circuit (see README)."""
    N_cells = model.N_cells
    v, w = model.v, model.w
    PBC = model.PBC
    num_up, num_dn = model.num_up, model.num_dn
    L = 2 * N_cells
    qc = QuantumCircuit(2 * L)

    # Slater determinant of each spin block (X gates + OpenFermion Givens network)
    prepare_slater_determinant(qc, Q_up, list(range(L)))
    prepare_slater_determinant(qc, Q_dn, list(range(L, 2 * L)))
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
