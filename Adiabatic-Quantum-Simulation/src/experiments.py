"""
experiments.py - High-level experiment runners.

Each runner builds circuits via the shared core, computes the relevant
observable, writes a JSON record, and returns an in-memory result dict so the
caller can plot immediately without re-reading from disk.

  * fidelity_scan        : Trotter convergence -> choose T_A and steps.
  * berry_sweep          : Berry phase vs inter-cell hopping w (= t2), PBC.
  * polarization_sweep   : electron polarization vs Delta U, OBC half-filling.
"""

from __future__ import annotations

import json
import os

import numpy as np
from qiskit.quantum_info import Statevector

from .core import SSHHModel, build_annealing_circuit
from .observables import twist_invariant, density_profile, wrap_berry_phase


def _statevector(qc):
    return np.array(Statevector(qc).data)


def _ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


# =====================================================================
# Experiment 1 : Trotter-convergence fidelity (choose T_A and steps/L)
# =====================================================================
def fidelity_scan(
    N_cells,
    t1,
    t2,
    number_of_electrons,
    U=0.0,
    T_A_values=(1, 5, 15, 80),
    steps_list=(1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150),
    PBC=True,
    reference="ground",
    out_json=None,
):
    """Scan Trotter `steps` for several total times `T_A` and report fidelity.

    reference='ground'   : fidelity against the prepared state psi(0)
                           (Trotter error of the *non-interacting* evolution;
                           the interaction is held constant, ramp_U=False).
    reference='adjacent' : fidelity between consecutive step counts
                           (self-convergence; use with U != 0).

    The interaction is held constant (ramp_U=False) for both modes, matching
    TrotterforPBC.py / TrotterforPBCvslaststep.py.
    """
    model = SSHHModel.from_total_electrons(N_cells, t1, t2, number_of_electrons, PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()

    # Reference state psi(0): pure state preparation, no evolution.
    qc0 = build_annealing_circuit(model, Q_up, Q_dn, U, U, T_A=0, steps=0, ramp_U=False)
    sv0 = Statevector(qc0)

    curves = {}
    for T_A in T_A_values:
        fidelities = []
        sv_prev = sv0
        for s in steps_list:
            qc = build_annealing_circuit(
                model, Q_up, Q_dn, U, U, T_A=T_A, steps=s, ramp_U=False
            )
            sv = Statevector(qc)
            ref = sv0 if reference == "ground" else sv_prev
            fidelities.append(float(np.abs(ref.inner(sv)) ** 2))
            sv_prev = sv
        curves[float(T_A)] = fidelities

    result = {
        "experiment": "fidelity",
        "params": {
            "N_cells": N_cells, "t1": float(t1), "t2": float(t2),
            "number_of_electrons": number_of_electrons, "U": float(U),
            "PBC": PBC, "reference": reference,
            "T_A_values": list(map(float, T_A_values)),
            "steps_list": list(steps_list),
        },
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
def berry_sweep(
    N_cells,
    number_of_electrons,
    t1=1.0,
    t2_values=(0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0),
    U_A=0.01,
    delta_U_values=(0.0,),
    T_A=1.0,
    steps=40,
    mode="up_spin",
    out_dir=None,
):
    """Adiabatically ramp the (staggered) Hubbard interaction, then measure the
    Berry phase / twist amplitude for every (delta_U, t2) point. PBC.

    One JSON file is written per (delta_U, t2) point (compatible with the
    original deltaU.py output and with plotting.plot_berry_phase).
    """
    L_sites = 2 * N_cells
    total_orbitals = 4 * N_cells
    _ensure_dir(out_dir)

    records = []  # flat list of per-point dicts (for immediate plotting)
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for t2 in t2_values:
            model = SSHHModel.from_total_electrons(
                N_cells, t1, t2, number_of_electrons, PBC=True
            )
            Q_up, Q_dn, _, _ = model.slater_Q_matrices()
            qc = build_annealing_circuit(
                model, Q_up, Q_dn, U_A, U_B, T_A=T_A, steps=steps, ramp_U=True
            )
            sv_np = _statevector(qc)
            z_N, probs = twist_invariant(sv_np, L_sites, mode=mode)

            data = {
                "L_sites": L_sites, "total_orbitals": total_orbitals,
                "N_cells": N_cells, "U_A": float(U_A), "U_B": float(U_B),
                "delta_U": float(U_B - U_A), "t1": float(t1), "t2": float(t2),
                "Z_N_real": float(np.real(z_N)), "Z_N_imag": float(np.imag(z_N)),
                "Twist_Amplitude": float(np.abs(z_N)),
                "Berry_Phase_pi": float(
                    np.mod(np.imag(np.log(z_N)), 2.0 * np.pi) / np.pi
                ),
                "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z_N)),
            }
            records.append(data)

            if out_dir:
                fname = f"UA_{U_A:.2f}_UB_{U_B:.4f}_t2_{t2:.2f}.json"
                with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                    json.dump({"data": data, "state_probabilities": probs},
                              f, indent=4)
    return records


# =====================================================================
# Experiment 3 : Electron polarization vs Delta U (OBC half-filling, edge)
# =====================================================================
def polarization_sweep(
    N_cells,
    total_electrons,
    t1=0.5,
    t2=1.5,
    U_A=1.0,
    delta_U_values=(0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3),
    T_A=1.0,
    steps=40,
    out_dir=None,
):
    """Adiabatically ramp the staggered interaction and measure the real-space
    sublattice polarization < n_A,j > - < n_B,j > per unit cell. OBC.

    One JSON file is written per delta_U value (compatible with
    halffillingedgewithdeltaU.py).
    """
    L_sites = 2 * N_cells
    _ensure_dir(out_dir)

    model = SSHHModel.from_total_electrons(N_cells, t1, t2, total_electrons, PBC=False)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    j_list = list(range(N_cells))

    records = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_annealing_circuit(
            model, Q_up, Q_dn, U_A, U_B, T_A=T_A, steps=steps, ramp_U=True
        )
        sv_np = _statevector(qc)
        density = density_profile(sv_np, L_sites)

        n_A = [round(float(density[2 * j]), 10) for j in range(N_cells)]
        n_B = [round(float(density[2 * j + 1]), 10) for j in range(N_cells)]
        delta_n = [round(n_A[j] - n_B[j], 10) for j in range(N_cells)]

        record = {
            "metadata": {
                "N_cells": N_cells, "L_sites": L_sites,
                "Total_Electrons": total_electrons, "t1": float(t1),
                "t2": float(t2), "U_A": float(U_A), "U_B": float(U_B),
                "delta_U": float(delta_U), "T_A": float(T_A), "steps": steps,
                "PBC": False,
            },
            "data": {
                "unit_cell_j": j_list, "n_A": n_A, "n_B": n_B,
                "n_A_minus_n_B": delta_n,
            },
        }
        records.append(record)

        if out_dir:
            fname = f"Polarization_UA_{U_A:.1f}_UB_{U_B:.4f}.json"
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump(record, f, indent=4)
    return records
