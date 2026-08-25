"""
observables.py - Topological observables + an extensible property registry.

Everything works directly on the raw statevector amplitudes with integer bit
masks (NumPy-vectorised). We deliberately AVOID Statevector.probabilities_dict(),
which materialises a label string for every one of the 2^n basis states
(~5 GiB at n=24) and runs out of memory; here we only touch the non-zero
amplitudes.

Adding a new property
---------------------
Write a function f(sv_np, layout) -> dict and register it in PROPERTIES. A
`layout` is a dict describing how qubits map to unit cells / sublattices:

    layout = {
        "n_qubits":   int,
        "n_cells":    int,
        "cell_qubits": [[q,...], ...],   # qubit indices counted in each cell
        "A_qubits":   [q, ...],          # (optional) sublattice-A site qubits
        "B_qubits":   [q, ...],          # (optional) sublattice-B site qubits
    }
"""

from __future__ import annotations

import numpy as np


# ------------------------------------------------------------------ helpers
def _nonzero_probs(sv_np, cutoff=1e-8):
    probs = np.abs(np.asarray(sv_np)) ** 2
    idx = np.nonzero(probs > cutoff)[0]
    return idx, probs[idx]


def _bit(idx, pos):
    """Occupation (0/1) of qubit `pos` for each integer basis index (little-endian)."""
    return (idx >> pos) & 1


# ------------------------------------------------------------------ layouts
def spin_block_layout(N_cells, mode="up_spin"):
    """Layout for the built-in SSH-Hubbard spin-block mapping."""
    L = 2 * N_cells
    cell_qubits = []
    A_qubits, B_qubits = [], []
    for j in range(N_cells):
        qs = [2 * j, 2 * j + 1]
        A_qubits.append(2 * j); B_qubits.append(2 * j + 1)
        if mode != "up_spin":  # total: include the down-spin block
            qs = qs + [L + 2 * j, L + 2 * j + 1]
            A_qubits.append(L + 2 * j); B_qubits.append(L + 2 * j + 1)
        cell_qubits.append(qs)
    return {
        "n_qubits": 2 * L, "n_cells": N_cells, "cell_qubits": cell_qubits,
        "A_qubits": A_qubits, "B_qubits": B_qubits, "L_sites": L,
    }


# ------------------------------------------------------------------ core observables
def twist_invariant_layout(sv_np, layout):
    """Many-body twist z_N = <exp(i 2pi/N_cells * X)>, X = sum_j (j+1) n_cell_j."""
    n_cells = layout["n_cells"]
    cell_qubits = layout["cell_qubits"]
    idx, p = _nonzero_probs(sv_np)
    if idx.size == 0:
        return 0.0 + 0.0j
    X = np.zeros(idx.shape, dtype=np.float64)
    for j, qubits in enumerate(cell_qubits):
        cell = np.zeros(idx.shape, dtype=np.int64)
        for q in qubits:
            cell = cell + _bit(idx, q)
        X += (j + 1) * cell
    return complex(np.sum(p * np.exp(1j * (2.0 * np.pi / n_cells) * X)))


def density_profile_qubits(sv_np, n_qubits):
    """Per-qubit occupation expectation <n_q>."""
    idx, p = _nonzero_probs(sv_np)
    dens = np.zeros(n_qubits)
    if idx.size == 0:
        return dens
    for q in range(n_qubits):
        dens[q] = np.sum(p * _bit(idx, q))
    return dens


# ------------------------------------------------------------------ legacy/simple wrappers
def twist_invariant(sv_np, L_sites, mode="up_spin"):
    """SSH convenience wrapper (spin-block layout)."""
    return twist_invariant_layout(sv_np, spin_block_layout(L_sites // 2, mode))


def density_profile(sv_np, L_sites):
    """Total electron density per spatial site, <n_i> = <n_{i,up}> + <n_{i,dn}>."""
    idx, p = _nonzero_probs(sv_np)
    density = np.zeros(L_sites)
    if idx.size == 0:
        return density
    for i in range(L_sites):
        density[i] = np.sum(p * (_bit(idx, i) + _bit(idx, L_sites + i)))
    return density


def wrap_berry_phase(z_N, convention="symmetric"):
    """Berry phase gamma in units of pi.

    convention="symmetric" (default): wrap to [-1, 1)  i.e. gamma in [-pi, pi).
        Standard branch; anything above +pi folds to negative.
    convention="zero_one": wrap to [-1/2, 3/2)  i.e. gamma in [-pi/2, 3*pi/2).
        The paper's phase-diagram convention where the topological/trivial
        plateaus sit near gamma=1 and gamma=0.
    """
    raw = np.imag(np.log(z_N))
    if convention == "zero_one":
        return (np.mod(raw + np.pi / 2.0, 2.0 * np.pi) - np.pi / 2.0) / np.pi
    # symmetric [-pi, pi)
    return (np.mod(raw + np.pi, 2.0 * np.pi) - np.pi) / np.pi


# ------------------------------------------------------------------ property registry
def _prop_berry(sv_np, layout):
    z = twist_invariant_layout(sv_np, layout)
    return {
        "property": "berry",
        "Z_N_real": float(np.real(z)), "Z_N_imag": float(np.imag(z)),
        "Twist_Amplitude": float(np.abs(z)),
        "Berry_Phase_pi": float(np.mod(np.imag(np.log(z)), 2.0 * np.pi) / np.pi),
        "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z)),
    }


def _prop_polarization(sv_np, layout):
    """Sublattice polarization per cell, <n_A,j> - <n_B,j>, from layout sublattices."""
    dens = density_profile_qubits(sv_np, layout["n_qubits"])
    n_cells = layout["n_cells"]
    # Map A/B site qubits back to cells (assumes A_qubits[j], B_qubits[j] per cell;
    # for the 'total' spin layout there are 2 A and 2 B qubits per cell -> summed).
    per_cell = max(1, len(layout.get("A_qubits", [])) // n_cells)
    A = layout.get("A_qubits", [])
    B = layout.get("B_qubits", [])
    nA, nB, diff = [], [], []
    for j in range(n_cells):
        a = sum(dens[A[j * per_cell + k]] for k in range(per_cell))
        b = sum(dens[B[j * per_cell + k]] for k in range(per_cell))
        nA.append(round(float(a), 10)); nB.append(round(float(b), 10))
        diff.append(round(float(a - b), 10))
    return {"property": "polarization", "unit_cell_j": list(range(n_cells)),
            "n_A": nA, "n_B": nB, "n_A_minus_n_B": diff}


# name -> measurement function. Extend this to add new properties.
PROPERTIES = {
    "berry": _prop_berry,
    "polarization": _prop_polarization,
}


# ======================================================================
# Counts-based estimators (real hardware / shot-based simulators)
# ======================================================================
# On real QPUs there is no statevector -- only measurement shots. The paper's
# protocol (Sec. IV-C) is built for exactly this: both observables are diagonal
# in the computational basis after Jordan-Wigner, so ONE set of Z-basis
# measurements suffices and no extra circuit depth is needed.
#
# BIT ORDERING -- the single most common source of silent errors here.
# Qiskit's get_counts() returns bitstrings written with qubit 0 on the RIGHT
# (most-significant character = highest qubit index). The statevector code above
# indexes qubits little-endian via `_bit(idx, pos) = (idx >> pos) & 1`.
# _bits_of_key() below reverses the string so that bits[q] is the occupation of
# logical qubit q, matching that convention exactly.
#
# measure_all() appends a fresh classical register; when a circuit has no other
# classical bits the key length equals num_qubits. Any whitespace (present when
# several registers are measured) is stripped.

def _bits_of_key(key, n_qubits):
    """Convert a get_counts() bitstring into a list `bits` with bits[q] = n_q.

    Qiskit prints qubit 0 last, so the string is reversed. Shorter keys are
    zero-padded on the left (i.e. high qubit indices) before reversing.
    """
    s = str(key).replace(" ", "")
    if len(s) < n_qubits:
        s = s.rjust(n_qubits, "0")
    elif len(s) > n_qubits:
        s = s[-n_qubits:]            # keep the low-index qubits
    return [int(c) for c in reversed(s)]


def counts_to_arrays(counts, n_qubits):
    """Vectorise a counts dict into (bit_matrix, weights).

    bit_matrix[m, q] = occupation of qubit q in the m-th distinct outcome,
    weights[m]       = normalised frequency of that outcome (sums to 1).
    """
    keys = list(counts.keys())
    shots = float(sum(counts.values()))
    if shots <= 0:
        raise ValueError("empty counts dictionary")
    bit_matrix = np.array([_bits_of_key(k, n_qubits) for k in keys], dtype=np.int8)
    weights = np.array([counts[k] for k in keys], dtype=float) / shots
    return bit_matrix, weights


def twist_invariant_from_counts(counts, layout):
    """Estimator of z_N from measurement shots -- paper Eqs. (47)-(49).

        X^(m)   = sum_j (j+1) * n_cell_j        (evaluated per bitstring)
        z_N^(m) = exp(i 2pi/N_cells * X^(m))
        z_bar_N = (1/M) sum_m z_N^(m)

    Averaging the phase factor per shot avoids ever materialising the full
    empirical distribution p(b), which is what makes the protocol scalable.
    """
    n_qubits = layout["n_qubits"]
    n_cells = layout["n_cells"]
    bits, wts = counts_to_arrays(counts, n_qubits)

    X = np.zeros(bits.shape[0], dtype=float)
    for j, qubits in enumerate(layout["cell_qubits"]):
        cell_occ = bits[:, qubits].sum(axis=1)
        X += (j + 1) * cell_occ
    return complex(np.sum(wts * np.exp(1j * (2.0 * np.pi / n_cells) * X)))


def density_profile_from_counts(counts, L_sites):
    """Total electron density per spatial site from shots, <n_i> = <n_i,up> + <n_i,dn>.

    Mirrors density_profile() but reads occupations out of bitstrings.
    """
    n_qubits = 2 * L_sites
    bits, wts = counts_to_arrays(counts, n_qubits)
    density = np.zeros(L_sites)
    for i in range(L_sites):
        density[i] = float(np.sum(wts * (bits[:, i] + bits[:, L_sites + i])))
    return density


def shot_noise_estimate(counts):
    """Rough 1/sqrt(M) statistical uncertainty scale for the estimators above."""
    shots = float(sum(counts.values()))
    return 1.0 / np.sqrt(shots) if shots > 0 else float("nan")