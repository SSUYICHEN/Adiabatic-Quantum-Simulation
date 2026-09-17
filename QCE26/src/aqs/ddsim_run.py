"""
ddsim_run.py - Run the SSH-Hubbard experiments on the MQT DDSIM decision-diagram
simulator (local, noiseless).

Why DDSIM (vs the statevector path you already have)
----------------------------------------------------
DDSIM represents the quantum state as a decision diagram instead of a dense
2^n vector. For circuits with structure and symmetry -- which the SSHH model
has (fixed particle number per spin sector, sublattice structure) -- the
diagram can stay far smaller than the full statevector, so DDSIM can reach
system sizes where a dense statevector runs out of memory. That makes it the
natural tool for pushing to N = 8, 10, ... beyond what Aer/CUDA-Q statevector
handles.

It is a LOCAL, synchronous, NOISELESS simulator: no cloud, no auth, no queue,
no cost. A run returns results immediately, so -- unlike ibm.py / azure.py --
there is no manifest and no separate fetch step. Sampling `shots` from the
final state gives the same counts-based estimators used everywhere else, so
the numbers should agree with your statevector results at N=6 (cross-check),
then extend to larger N.

Everything downstream (Eqs. 47-52 estimators, JSON schema, plotting) is reused
unchanged from aqs.observables / aqs.plotting.
"""

from __future__ import annotations

import json
import os

import numpy as np

from .core import SSHHModel, build_annealing_circuit
from .observables import (
    spin_block_layout,
    twist_invariant_from_counts,
    density_profile_from_counts,
    wrap_berry_phase,
    shot_noise_estimate,
)


# =====================================================================
# Backend
# =====================================================================
def get_backend(name="qasm_simulator"):
    """MQT DDSIM qiskit backend. 'qasm_simulator' = sampling (counts)."""
    try:
        from mqt import ddsim
    except ImportError as e:
        raise RuntimeError("MQT DDSIM not installed:  pip install mqt.ddsim") from e
    return ddsim.DDSIMProvider().get_backend(name)


def build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                           T_A=1.0, steps=20, PBC=True):
    model = SSHHModel.from_total_electrons(N_cells, v, w, electrons, PBC=PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B,
                                 T_A=T_A, steps=steps, ramp_U=True)
    qc.measure_all()
    return qc


def _run_counts(qc, shots, backend):
    job = backend.run(qc, shots=shots)
    result = job.result()
    try:
        return result.get_counts(qc)
    except Exception:
        return result.get_counts()


# =====================================================================
# Berry phase sweep (PBC)
# =====================================================================
def run_berry(N_cells=6, electrons=None, v=1.0,
              w_values=(0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0),
              U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
              T_A=1.0, steps=20, shots=8192,
              out_dir="ddsim_results/berry"):
    """Berry phase vs w at half filling (n = 2N). Writes one JSON per (dU, w)
    in the schema aqs.plotting.plot_berry_phase understands."""
    electrons = electrons if electrons is not None else 2 * N_cells
    backend = get_backend("qasm_simulator")
    data_dir = os.path.join(out_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    layout = spin_block_layout(N_cells, mode="up_spin")

    print(f"[ddsim berry] N={N_cells} Ne={electrons} v={v} L={steps} "
          f"shots={shots} PBC  ({len(delta_U_values)}x{len(w_values)} circuits)")
    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for w in w_values:
            qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                        T_A, steps, PBC=True)
            counts = _run_counts(qc, shots, backend)
            z_N = twist_invariant_from_counts(counts, layout)
            data = {
                "L_sites": 2 * N_cells, "N_cells": N_cells,
                "U_A": U_A, "U_B": U_B, "delta_U": float(delta_U),
                "v": float(v), "w": float(w), "t1": float(v), "t2": float(w),
                "Z_N_real": float(np.real(z_N)), "Z_N_imag": float(np.imag(z_N)),
                "Twist_Amplitude": float(np.abs(z_N)),
                "Berry_Phase_pi": float(np.mod(np.imag(np.log(z_N)), 2 * np.pi) / np.pi),
                "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z_N)),
                "shots": int(sum(counts.values())),
                "shot_noise": shot_noise_estimate(counts),
                "backend": "ddsim.qasm_simulator",
            }
            fname = f"UA_{U_A:.2f}_UB_{U_B:.4f}_w_{w:.2f}_N{N_cells}.json"
            with open(os.path.join(data_dir, fname), "w", encoding="utf-8") as f:
                json.dump({"data": data}, f, indent=2)
            print(f"  dU={delta_U:<5g} w={w:<5g} "
                  f"|z|={abs(z_N):.4f} gamma={data['Berry_Phase_pi_wrapped']:+.4f} pi")
            records.append(data)
    print(f"  -> {data_dir}")
    return records, data_dir


# =====================================================================
# Polarization sweep (OBC, n = 2N+2)
# =====================================================================
def run_polarization(N_cells=6, electrons=None, v=0.5, w=1.5,
                     U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                     T_A=1.0, steps=20, shots=8192,
                     out_dir="ddsim_results/polarization"):
    """Sublattice polarization vs unit cell for each dU (OBC). Writes one JSON
    per dU in the schema aqs.plotting.plot_polarization understands."""
    electrons = electrons if electrons is not None else 2 * N_cells + 2
    backend = get_backend("qasm_simulator")
    data_dir = os.path.join(out_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    L_sites = 2 * N_cells

    print(f"[ddsim polarization] N={N_cells} Ne={electrons} v={v} w={w} "
          f"L={steps} shots={shots} OBC  ({len(delta_U_values)} circuits)")
    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                    T_A, steps, PBC=False)
        counts = _run_counts(qc, shots, backend)
        density = density_profile_from_counts(counts, L_sites)
        n_A = [round(float(density[2 * j]), 10) for j in range(N_cells)]
        n_B = [round(float(density[2 * j + 1]), 10) for j in range(N_cells)]
        diff = [round(n_A[j] - n_B[j], 10) for j in range(N_cells)]
        rec = {
            "metadata": {"N_cells": N_cells, "L_sites": L_sites,
                         "Total_Electrons": electrons, "v": v, "w": w,
                         "t1": v, "t2": w, "U_A": U_A, "U_B": U_B,
                         "delta_U": float(delta_U), "T_A": T_A, "steps": steps,
                         "PBC": False, "backend": "ddsim.qasm_simulator",
                         "shots": int(sum(counts.values())),
                         "shot_noise": shot_noise_estimate(counts)},
            "data": {"unit_cell_j": list(range(N_cells)),
                     "n_A": n_A, "n_B": n_B, "n_A_minus_n_B": diff},
        }
        fname = f"Polarization_UA_{U_A:.1f}_UB_{U_B:.4f}_N{N_cells}.json"
        with open(os.path.join(data_dir, fname), "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=2)
        print(f"  dU={delta_U:<5g} edge(cell 1)={diff[0]:+.4f}")
        records.append(rec)
    print(f"  -> {data_dir}")
    return records, data_dir