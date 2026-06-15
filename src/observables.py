"""
observables.py - Topological observables extracted from the final state.

  * twist_invariant   : the many-body Resta "twist" expectation z_N, from which
                        the Berry phase (bulk topological invariant) is read off.
  * density_profile   : site-resolved electron density, used for the real-space
                        electron polarization < n_A > - < n_B > (edge response).
"""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import Statevector


def _probabilities(sv_np, cutoff=1e-8):
    sv = Statevector(sv_np)
    return {bs: p for bs, p in sv.probabilities_dict().items() if p > cutoff}


def twist_invariant(sv_np, L_sites, mode="up_spin"):
    """Many-body twist expectation z_N = <exp(i 2pi/N_cells * X)>.

    X = sum_j (j+1) * n_cell_j, where n_cell_j is the particle number in cell j
    (spin-up only for mode='up_spin', total for mode='total').

    Returns
    -------
    z_N : complex
        |z_N| is the twist amplitude; arg(z_N) gives the Berry phase.
    """
    N_cells = L_sites // 2
    probs = _probabilities(sv_np)

    z_N = 0.0 + 0.0j
    for b, p in probs.items():
        rev_b = b[::-1]  # Qiskit bitstrings are little-endian
        X_b = 0.0
        for i in range(N_cells):
            n_2i, n_2i_1 = int(rev_b[2 * i]), int(rev_b[2 * i + 1])
            if mode == "up_spin":
                cell_particles = n_2i + n_2i_1
            else:  # total (up + down)
                n_dn_2i = int(rev_b[L_sites + 2 * i])
                n_dn_2i_1 = int(rev_b[L_sites + 2 * i + 1])
                cell_particles = n_2i + n_2i_1 + n_dn_2i + n_dn_2i_1
            X_b += (i + 1) * cell_particles
        z_N += p * np.exp(1j * (2.0 * np.pi / N_cells) * X_b)
    return z_N, probs


def density_profile(sv_np, L_sites):
    """Site-resolved total electron density <n_i> = <n_{i,up}> + <n_{i,dn}>."""
    probs = _probabilities(sv_np)
    density = np.zeros(L_sites)
    for b, p in probs.items():
        rev_b = b[::-1]
        for i in range(L_sites):
            n_up = int(rev_b[i])
            n_dn = int(rev_b[L_sites + i])
            density[i] += p * (n_up + n_dn)
    return density


def wrap_berry_phase(z_N):
    """Berry phase in units of pi, wrapped to (-1/2, 3/2] like the plots."""
    raw = np.imag(np.log(z_N))
    wrapped = np.mod(raw + np.pi / 2.0, 2.0 * np.pi) - np.pi / 2.0
    return wrapped / np.pi
