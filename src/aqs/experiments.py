"""
experiments.py - High-level experiment runners for the built-in SSH-H model.

Each runner builds circuits via the shared core, simulates them through the
chosen backend (qiskit CPU or cudaq GPU), computes the observable, writes JSON,
and returns the in-memory result.

  * fidelity_scan        : Trotter convergence -> choose T_A and steps/L.
  * berry_sweep          : Berry phase vs inter-cell hopping w, PBC.
  * polarization_sweep   : electron polarization vs Delta U, OBC half-filling.

For arbitrary Hamiltonians (not the SSH model) see hamiltonian.measure_hamiltonian.
"""

from __future__ import annotations

import json
import os

import numpy as np

from .core import SSHHModel, build_annealing_circuit
from .observables import twist_invariant, density_profile, wrap_berry_phase
from .backends import get_backend


def _ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


def _resolve_backend(backend):
    return backend if hasattr(backend, "statevector") else get_backend(backend or "qiskit")


# =====================================================================
# Experiment 1 : Trotter-convergence fidelity (choose T_A and steps/L)
# =====================================================================
def fidelity_scan(N_cells, v, w, number_of_electrons, U=0.0,
                  T_A_values=(1, 5, 15, 80),
                  steps_list=(1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150),
                  PBC=True, reference="ground", out_json=None, backend="qiskit"):
    bk = _resolve_backend(backend)
    model = SSHHModel.from_total_electrons(N_cells, v, w, number_of_electrons, PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()

    qc0 = build_annealing_circuit(model, Q_up, Q_dn, U, U, T_A=0, steps=0, ramp_U=False)
    sv0 = bk.statevector(qc0)

    def fid(a, b):
        return float(np.abs(np.vdot(a, b)) ** 2)

    curves = {}
    for T_A in T_A_values:
        fidelities = []
        sv_prev = sv0
        for s in steps_list:
            qc = build_annealing_circuit(model, Q_up, Q_dn, U, U, T_A=T_A, steps=s, ramp_U=False)
            sv = bk.statevector(qc)
            fidelities.append(fid(sv0 if reference == "ground" else sv_prev, sv))
            sv_prev = sv
        curves[float(T_A)] = fidelities

    result = {
        "experiment": "fidelity", "backend": bk.name,
        "params": {"N_cells": N_cells, "v": float(v), "w": float(w),
                   "number_of_electrons": number_of_electrons, "U": float(U),
                   "PBC": PBC, "reference": reference,
                   "T_A_values": list(map(float, T_A_values)),
                   "steps_list": list(steps_list)},
        "curves": curves,
    }
    if out_json:
        _ensure_dir(os.path.dirname(out_json))
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=4)
    return result


# =====================================================================
# Experiment 2 : Berry phase vs inter-cell hopping (PBC, bulk invariant)
# =====================================================================
def berry_sweep(N_cells, number_of_electrons, v=1.0,
                w_values=(0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0),
                U_A=0.01, delta_U_values=(0.0,), T_A=1.0, steps=40,
                mode="up_spin", out_dir=None, backend="qiskit"):
    bk = _resolve_backend(backend)
    L_sites = 2 * N_cells
    total_orbitals = 4 * N_cells
    _ensure_dir(out_dir)

    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for w in w_values:
            model = SSHHModel.from_total_electrons(N_cells, v, w, number_of_electrons, PBC=True)
            Q_up, Q_dn, _, _ = model.slater_Q_matrices()
            qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B, T_A=T_A, steps=steps, ramp_U=True)
            z_N = twist_invariant(bk.statevector(qc), L_sites, mode=mode)
            data = {
                "L_sites": L_sites, "total_orbitals": total_orbitals,
                "N_cells": N_cells, "U_A": float(U_A), "U_B": float(U_B),
                "delta_U": float(U_B - U_A), "v": float(v), "w": float(w),
                # 'w' is the inter-cell hopping (was 't2'); kept as 't2' too for
                # backward-compatible plotting of older JSON files.
                "t1": float(v), "t2": float(w),
                "Z_N_real": float(np.real(z_N)), "Z_N_imag": float(np.imag(z_N)),
                "Twist_Amplitude": float(np.abs(z_N)),
                "Berry_Phase_pi": float(np.mod(np.imag(np.log(z_N)), 2.0 * np.pi) / np.pi),
                "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z_N)),
            }
            records.append(data)
            if out_dir:
                fname = f"UA_{U_A:.2f}_UB_{U_B:.4f}_w_{w:.2f}.json"
                with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                    json.dump({"data": data}, f, indent=4)
    return records


# =====================================================================
# Experiment 3 : Electron polarization vs Delta U (OBC half-filling)
# =====================================================================
def polarization_sweep(N_cells, total_electrons, v=0.5, w=1.5, U_A=1.0,
                       delta_U_values=(0, 0.0001, 0.0003, 0.001, 0.003, 0.01,
                                       0.03, 0.1, 0.3, 1, 3),
                       T_A=1.0, steps=40, out_dir=None, backend="qiskit"):
    bk = _resolve_backend(backend)
    L_sites = 2 * N_cells
    _ensure_dir(out_dir)

    model = SSHHModel.from_total_electrons(N_cells, v, w, total_electrons, PBC=False)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    j_list = list(range(N_cells))

    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B, T_A=T_A, steps=steps, ramp_U=True)
        density = density_profile(bk.statevector(qc), L_sites)
        n_A = [round(float(density[2 * j]), 10) for j in range(N_cells)]
        n_B = [round(float(density[2 * j + 1]), 10) for j in range(N_cells)]
        delta_n = [round(n_A[j] - n_B[j], 10) for j in range(N_cells)]
        record = {
            "metadata": {"N_cells": N_cells, "L_sites": L_sites,
                         "Total_Electrons": total_electrons, "v": float(v),
                         "w": float(w), "t1": float(v), "t2": float(w),
                         "U_A": float(U_A), "U_B": float(U_B),
                         "delta_U": float(delta_U), "T_A": float(T_A),
                         "steps": steps, "PBC": False, "backend": bk.name},
            "data": {"unit_cell_j": j_list, "n_A": n_A, "n_B": n_B,
                     "n_A_minus_n_B": delta_n},
        }
        records.append(record)
        if out_dir:
            fname = f"Polarization_UA_{U_A:.1f}_UB_{U_B:.4f}.json"
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump(record, f, indent=4)
    return records
