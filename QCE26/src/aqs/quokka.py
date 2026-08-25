"""
quokka.py - Compute SSHH observables with Quokka-Sharp (model counting).

How this differs from the other backends
----------------------------------------
Statevector/ddsim/hardware sample the FULL distribution then post-process.
Quokka-Sharp answers ONE joint probability per query via weighted model
counting (GPMC). It never materialises the state.

Verified facts (established interactively):
  * FULLY-specified joint probabilities are EXACT: sim(measurement={0:1,1:0,...})
    with every qubit fixed matches statevector to ~1e-3 (0.0130 vs 0.0130).
  * MARGINAL queries are BROKEN in this Quokka version: a partial dict like
    {i:1} returns 0 in the computational basis and segfaults in the pauli
    basis. So we CANNOT ask for <n_i> directly.
  * QASM export is faithful (reloading the QASM reproduces the exact density),
    so the discrepancy is purely Quokka's marginal handling, not our circuit.

Strategy: occupations by enumeration
------------------------------------
Because the circuit conserves particle number per spin sector, only
C(2N, num_up) * C(2N, num_dn) computational-basis configurations are non-zero.
We enumerate those, get each one's EXACT joint probability from Quokka, and sum
to obtain <n_i>. This sidesteps the broken marginal path entirely.

Scaling: configs = C(2N,num_up)*C(2N,num_dn).
  N=2 -> 36, N=4 -> 3136, N=6 -> ~630k.
So this is practical for N=2, borderline for N=4, and infeasible for N=6 --
the enumeration count IS the feasibility boundary, a result worth reporting
(analogous to ddsim's decision-diagram blow-up).

Polarization only. Berry phase needs a complex non-local phase expectation the
joint-probability interface cannot provide.
"""

from __future__ import annotations

import itertools
import json
import os

import numpy as np

from .core import SSHHModel, build_annealing_circuit

QUOKKA_BASIS_GATES = ["h", "s", "t", "rx", "ry", "rz", "cx", "cz"]


# =====================================================================
# Circuit -> QASM in Quokka's gate set
# =====================================================================
def circuit_to_qasm(qc, path):
    """Decompose to Quokka's gate set, strip barriers, zero global phase,
    write OpenQASM. (Verified faithful: reloading reproduces the density.)"""
    from qiskit import transpile
    from qiskit.transpiler.passes import RemoveBarriers

    qc2 = transpile(qc, basis_gates=QUOKKA_BASIS_GATES, optimization_level=1)
    qc2 = RemoveBarriers()(qc2)
    qc2.global_phase = 0

    try:
        from qiskit.qasm2 import dumps
        qasm_str = dumps(qc2)
    except Exception:
        qasm_str = qc2.qasm()
    with open(path, "w", encoding="utf-8") as f:
        f.write(qasm_str)
    return qc2.num_qubits


def build_circuit(N_cells, v, w, electrons, U_A, U_B,
                  T_A=1.0, steps=20, PBC=True):
    model = SSHHModel.from_total_electrons(N_cells, v, w, electrons, PBC=PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B,
                                 T_A=T_A, steps=steps, ramp_U=True)
    return qc, model


# =====================================================================
# Occupations by enumeration of particle-number-conserving configs
# =====================================================================
def occupations_by_enumeration(qasmfile, N_cells, num_up, num_dn,
                               progress=False):
    """<n_i> for all qubits via exact joint probabilities. See module docstring.
    Returns (occupations list | 'TIMEOUT', n_configs, total_prob)."""
    import quokka_sharp as qk

    L = 2 * N_cells
    n_qubits = 2 * L
    occ = np.zeros(n_qubits)
    total_p = 0.0
    n_cfg = 0

    combos = [(u, d)
              for u in itertools.combinations(range(L), num_up)
              for d in itertools.combinations(range(L), num_dn)]
    if progress:
        print(f"    enumerating {len(combos)} configs ...")

    for up_occ, dn_occ in combos:
        meas = {q: (1 if q in up_occ else 0) for q in range(L)}
        meas.update({L + q: (1 if q in dn_occ else 0) for q in range(L)})
        p = qk.functionalities.sim(qasmfile=qasmfile, basis="comp",
                                   measurement=meas)
        if isinstance(p, str):
            return p, n_cfg, total_p
        p = float(np.real(p))
        n_cfg += 1
        total_p += p
        for q in up_occ:
            occ[q] += p
        for q in dn_occ:
            occ[L + q] += p

    return occ.tolist(), n_cfg, total_p


# =====================================================================
# Polarization
# =====================================================================
def run_polarization(N_cells=2, electrons=None, v=0.5, w=1.5,
                     U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                     T_A=1.0, steps=20, out_dir="quokka_results/polarization",
                     qasm_dir="/tmp"):
    """Sublattice polarization via Quokka-Sharp enumeration. OBC, n=2N+2."""
    electrons = electrons if electrons is not None else 2 * N_cells + 2
    L_sites = 2 * N_cells
    data_dir = os.path.join(out_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    import math
    n_combos = math.comb(L_sites, (electrons + 0) // 2)  # rough; real uses model
    print(f"[quokka polarization] N={N_cells} Ne={electrons} v={v} w={w} "
          f"L={steps} OBC")

    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc, model = build_circuit(N_cells, v, w, electrons, U_A, U_B,
                                  T_A, steps, PBC=False)
        qasm = os.path.join(qasm_dir, f"sshh_pol_N{N_cells}_dU{delta_U:g}.qasm")
        circuit_to_qasm(qc, qasm)

        occ, ncfg, tot = occupations_by_enumeration(
            qasm, N_cells, model.num_up, model.num_dn, progress=True)
        if isinstance(occ, str):
            print(f"  dU={delta_U:<5g} -> {occ} (enumeration did not finish)")
            continue

        density = [occ[i] + occ[L_sites + i] for i in range(L_sites)]
        n_A = [round(density[2 * j], 10) for j in range(N_cells)]
        n_B = [round(density[2 * j + 1], 10) for j in range(N_cells)]
        diff = [round(n_A[j] - n_B[j], 10) for j in range(N_cells)]

        rec = {
            "metadata": {"N_cells": N_cells, "L_sites": L_sites,
                         "Total_Electrons": electrons, "v": v, "w": w,
                         "t1": v, "t2": w, "U_A": U_A, "U_B": U_B,
                         "delta_U": float(delta_U), "T_A": T_A, "steps": steps,
                         "PBC": False, "backend": "quokka-sharp.gpmc",
                         "n_configs": ncfg, "total_prob": round(tot, 6)},
            "data": {"unit_cell_j": list(range(N_cells)),
                     "n_A": n_A, "n_B": n_B, "n_A_minus_n_B": diff},
        }
        fname = f"Polarization_UA_{U_A:.1f}_UB_{U_B:.4f}_N{N_cells}.json"
        with open(os.path.join(data_dir, fname), "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=2)
        print(f"  dU={delta_U:<5g} edge(cell1)={diff[0]:+.4f}  "
              f"configs={ncfg} total_p={tot:.4f} (should ~1.0)")
        records.append(rec)

    print(f"  -> {data_dir}")
    return records, data_dir


# =====================================================================
# Verify against statevector (small N)
# =====================================================================
def verify_against_statevector(N_cells=2, electrons=4, v=0.5, w=1.5,
                               U_A=0.01, U_B=0.01, T_A=1.0, steps=1):
    from qiskit.quantum_info import Statevector

    qc, model = build_circuit(N_cells, v, w, electrons, U_A, U_B, T_A, steps, PBC=False)
    n_qubits = qc.num_qubits

    sv = Statevector(qc)
    probs = sv.probabilities()
    exact = [sum(probs[idx] for idx in range(len(probs)) if (idx >> q) & 1)
             for q in range(n_qubits)]

    qasm = "/tmp/sshh_verify.qasm"
    circuit_to_qasm(qc, qasm)
    occ, ncfg, tot = occupations_by_enumeration(qasm, N_cells,
                                                model.num_up, model.num_dn)

    print(f"configs={ncfg}  total_prob={tot:.4f} (should ~1.0)")
    print(f"{'qubit':>5} {'statevector':>12} {'quokka':>12} {'match':>6}")
    ok = True
    for q in range(n_qubits):
        m = abs(exact[q] - occ[q]) < 1e-3
        ok = ok and m
        print(f"{q:>5} {exact[q]:>12.6f} {occ[q]:>12.6f} {'OK' if m else 'XX':>6}")
    print("ALL MATCH" if ok else "MISMATCH")
    return ok