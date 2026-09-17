"""
ibm.py - Run the SSH-Hubbard experiments on IBM Quantum hardware.

Why this is a separate path from backends.py
--------------------------------------------
The qiskit/cudaq backends return a STATEVECTOR. Real QPUs cannot: you only get
measurement shots. The paper's protocol (Sec. IV-C) is designed for exactly
this -- both the many-body Berry phase and the sublattice polarization are
diagonal in the computational basis after Jordan-Wigner, so a single set of
Z-basis measurements suffices, with the estimators of Eqs. (47)-(52) applied
classically to the sampled bitstrings.

Workflow (mirrors the four Qiskit-patterns steps)
-------------------------------------------------
  1. map      : build_annealing_circuit(...) + measure_all()   [aqs.core]
  2. optimize : generate_preset_pass_manager(backend).run(qc)  -> ISA circuit
  3. execute  : SamplerV2.run([...]) -> job (async; keep the job ID!)
  4. analyze  : counts -> twist_invariant_from_counts / density_profile_from_counts

Jobs are asynchronous and hardware queues are long, so every submission is
recorded to a manifest JSON. `fetch_results()` can reconstruct all observables
later from job IDs alone, without re-running anything.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np

from .core import SSHHModel, build_annealing_circuit
from .observables import (
    spin_block_layout,
    twist_invariant_from_counts,
    density_profile_from_counts,
    wrap_berry_phase,
    shot_noise_estimate,
)

DEFAULT_BACKEND = "ibm_kingston"


# =====================================================================
# Service / backend helpers
# =====================================================================
def get_service():
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError as e:
        raise RuntimeError(
            "qiskit-ibm-runtime is required for the IBM path.\n"
            "  pip install 'qiskit[all]' qiskit-ibm-runtime"
        ) from e
    return QiskitRuntimeService()


def median_2q_error(backend):
    """Median two-qubit gate error of a backend (native cz/ecr/cx)."""
    errs = []
    for inst, qargs_map in backend.target.items():
        if inst in ("cz", "ecr", "cx"):
            for _, props in qargs_map.items():
                if props is not None and props.error is not None:
                    errs.append(props.error)
    return float(np.median(errs)) if errs else None


def pick_backend(service, name=None, min_qubits=24):
    if name:
        return service.backend(name)
    return service.least_busy(operational=True, simulator=False,
                              min_num_qubits=min_qubits)


# =====================================================================
# Circuit construction + transpilation
# =====================================================================
def build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                           T_A=1.0, steps=20, PBC=True):
    """SSHH annealing circuit with the paper's Z-basis measurement appended."""
    model = SSHHModel.from_total_electrons(N_cells, v, w, electrons, PBC=PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B,
                                 T_A=T_A, steps=steps, ramp_U=True)
    qc.measure_all()
    return qc


def transpile_for(backend, qc, optimization_level=3):
    from qiskit.transpiler import generate_preset_pass_manager
    pm = generate_preset_pass_manager(backend=backend,
                                      optimization_level=optimization_level)
    return pm.run(qc)


def circuit_report(qc, isa, backend, err=None):
    n2q = isa.num_nonlocal_gates()
    logical_2q = max(qc.num_nonlocal_gates(), 1)
    rep = {
        "logical_depth": qc.depth(),
        "logical_2q": qc.num_nonlocal_gates(),
        "isa_depth": isa.depth(),
        "isa_2q": n2q,
        "expansion": round(n2q / logical_2q, 2),
        "backend": backend.name,
    }
    if err:
        rep["median_2q_error"] = err
        rep["p_no_error"] = float((1 - err) ** n2q)
    return rep


# =====================================================================
# Submission
# =====================================================================
def _sampler(backend, shots, dynamical_decoupling=True, twirling=True):
    from qiskit_ibm_runtime import SamplerV2, SamplerOptions
    options = SamplerOptions()
    options.default_shots = shots
    if dynamical_decoupling:
        # Error SUPPRESSION: idle qubits are echoed so they dephase less.
        # Worth enabling for deep circuits like ours.
        options.dynamical_decoupling.enable = True
        options.dynamical_decoupling.sequence_type = "XY4"
    if twirling:
        # Pauli twirling turns coherent gate errors into stochastic ones,
        # which the shot average handles far more gracefully.
        options.twirling.enable_gates = True
        options.twirling.enable_measure = True
    return SamplerV2(mode=backend, options=options)


def submit_berry(N_cells=6, electrons=12, v=1.0,
                 w_values=(0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0),
                 U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                 T_A=1.0, steps=20, shots=8192,
                 backend_name=DEFAULT_BACKEND, optimization_level=3,
                 out_dir="ibm_results/berry", dry_run=False):
    """Submit the Berry-phase sweep (PBC) to IBM hardware.

    One job per (delta_U, w) pair. With the defaults that is 4 x 10 = 40 jobs,
    so consider trimming w_values first -- see the module docstring.
    Set dry_run=True to transpile and report sizes without submitting.
    """
    service = get_service()
    backend = pick_backend(service, backend_name, min_qubits=4 * N_cells)
    err = median_2q_error(backend)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[berry] backend={backend.name} qubits={backend.num_qubits} "
          f"median2Qerr={err:.2e} queue={backend.status().pending_jobs}")
    print(f"        N={N_cells} Ne={electrons} v={v} T_A={T_A} L={steps} "
          f"shots={shots} PBC")
    print(f"        dU={list(delta_U_values)}  w={list(w_values)}")
    print(f"        -> {len(delta_U_values) * len(w_values)} jobs"
          f"{' (DRY RUN)' if dry_run else ''}\n")

    manifest = {
        "experiment": "berry", "backend": backend.name,
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
        "params": {"N_cells": N_cells, "electrons": electrons, "v": v,
                   "U_A": U_A, "T_A": T_A, "steps": steps, "shots": shots,
                   "PBC": True, "optimization_level": optimization_level,
                   "median_2q_error": err},
        "jobs": [],
    }

    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for w in w_values:
            qc = build_measured_circuit(N_cells, v, w, electrons,
                                        U_A, U_B, T_A, steps, PBC=True)
            isa = transpile_for(backend, qc, optimization_level)
            rep = circuit_report(qc, isa, backend, err)

            entry = {"delta_U": float(delta_U), "U_A": float(U_A),
                     "U_B": float(U_B), "v": float(v), "w": float(w),
                     "report": rep}

            if dry_run:
                print(f"  dU={delta_U:<6g} w={w:<5g} isa2q={rep['isa_2q']:<6d} "
                      f"depth={rep['isa_depth']:<6d} P={rep['p_no_error']:.2e}")
            else:
                sampler = _sampler(backend, shots)
                job = sampler.run([isa])
                entry["job_id"] = job.job_id()
                print(f"  dU={delta_U:<6g} w={w:<5g} isa2q={rep['isa_2q']:<6d} "
                      f"P={rep['p_no_error']:.2e}  job={job.job_id()}")
            manifest["jobs"].append(entry)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    tag = f"L{steps}_{backend.name}_{stamp}"
    path = os.path.join(out_dir, f"manifest_berry_{tag}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  manifest -> {path}")
    if not dry_run:
        print("  KEEP THIS FILE: it maps job IDs back to physical parameters.")
        print(f"  Later:  aqs-ibm fetch --manifest {path}")
    return manifest


def submit_polarization(N_cells=6, electrons=14, v=0.5, w=1.5,
                        U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                        T_A=1.0, steps=20, shots=8192,
                        backend_name=DEFAULT_BACKEND, optimization_level=3,
                        out_dir="ibm_results/polarization", dry_run=False):
    """Submit the sublattice-polarization sweep (OBC, n = 2N+2) to IBM hardware.

    One job per delta_U -- only 4 jobs with the defaults.
    """
    service = get_service()
    backend = pick_backend(service, backend_name, min_qubits=4 * N_cells)
    err = median_2q_error(backend)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[polarization] backend={backend.name} median2Qerr={err:.2e} "
          f"queue={backend.status().pending_jobs}")
    print(f"        N={N_cells} Ne={electrons} v={v} w={w} T_A={T_A} "
          f"L={steps} shots={shots} OBC")
    print(f"        dU={list(delta_U_values)} -> {len(delta_U_values)} jobs"
          f"{' (DRY RUN)' if dry_run else ''}\n")

    manifest = {
        "experiment": "polarization", "backend": backend.name,
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
        "params": {"N_cells": N_cells, "electrons": electrons, "v": v, "w": w,
                   "U_A": U_A, "T_A": T_A, "steps": steps, "shots": shots,
                   "PBC": False, "optimization_level": optimization_level,
                   "median_2q_error": err},
        "jobs": [],
    }

    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_measured_circuit(N_cells, v, w, electrons,
                                    U_A, U_B, T_A, steps, PBC=False)
        isa = transpile_for(backend, qc, optimization_level)
        rep = circuit_report(qc, isa, backend, err)

        entry = {"delta_U": float(delta_U), "U_A": float(U_A),
                 "U_B": float(U_B), "v": float(v), "w": float(w),
                 "report": rep}

        if dry_run:
            print(f"  dU={delta_U:<6g} isa2q={rep['isa_2q']:<6d} "
                  f"depth={rep['isa_depth']:<6d} P={rep['p_no_error']:.2e}")
        else:
            sampler = _sampler(backend, shots)
            job = sampler.run([isa])
            entry["job_id"] = job.job_id()
            print(f"  dU={delta_U:<6g} isa2q={rep['isa_2q']:<6d} "
                  f"P={rep['p_no_error']:.2e}  job={job.job_id()}")
        manifest["jobs"].append(entry)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    tag = f"L{steps}_{backend.name}_{stamp}"
    path = os.path.join(out_dir, f"manifest_polarization_{tag}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  manifest -> {path}")
    if not dry_run:
        print("  KEEP THIS FILE: it maps job IDs back to physical parameters.")
        print(f"  Later:  aqs-ibm fetch --manifest {path}")
    return manifest


# =====================================================================
# Retrieval + post-processing (paper Eqs. 47-52)
# =====================================================================
def _counts_of(job):
    """Extract the counts dict from a SamplerV2 result (measure_all -> 'meas')."""
    res = job.result()[0]
    data = res.data
    for field in ("meas", "c", "c0"):
        if hasattr(data, field):
            return getattr(data, field).get_counts()
    # fall back to whatever single register exists
    name = list(data.keys())[0] if hasattr(data, "keys") else None
    if name:
        return data[name].get_counts()
    raise RuntimeError("could not locate a classical register in the result")


def fetch_results(manifest_path, out_dir=None, save_counts=False):
    """Retrieve finished jobs and apply the paper's estimators.

    Writes JSON files in exactly the schema aqs.plotting already understands,
    so the hardware data can be plotted with plot_berry_phase /
    plot_polarization / the cross-section functions with no changes.
    """
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    service = get_service()
    exp = manifest["experiment"]
    P = manifest["params"]
    N_cells = P["N_cells"]
    L_sites = 2 * N_cells
    out_dir = out_dir or os.path.join(os.path.dirname(manifest_path), "data")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[fetch] {exp}  {len(manifest['jobs'])} jobs  backend={manifest['backend']}")
    records, pending = [], 0

    for entry in manifest["jobs"]:
        jid = entry.get("job_id")
        if not jid:
            continue
        job = service.job(jid)
        status = str(job.status())
        if "DONE" not in status.upper():
            print(f"  {jid}  status={status}  (skipped)")
            pending += 1
            continue

        counts = _counts_of(job)
        shots = int(sum(counts.values()))
        sigma = shot_noise_estimate(counts)

        if exp == "berry":
            layout = spin_block_layout(N_cells, mode="up_spin")
            z_N = twist_invariant_from_counts(counts, layout)
            data = {
                "L_sites": L_sites, "N_cells": N_cells,
                "U_A": entry["U_A"], "U_B": entry["U_B"],
                "delta_U": entry["delta_U"],
                "v": entry["v"], "w": entry["w"],
                "t1": entry["v"], "t2": entry["w"],
                "Z_N_real": float(np.real(z_N)),
                "Z_N_imag": float(np.imag(z_N)),
                "Twist_Amplitude": float(np.abs(z_N)),
                "Berry_Phase_pi": float(np.mod(np.imag(np.log(z_N)), 2 * np.pi) / np.pi),
                "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z_N)),
                "shots": shots, "shot_noise": sigma,
                "job_id": jid, "backend": manifest["backend"],
                "isa_2q": entry["report"]["isa_2q"],
            }
            if save_counts:
                data["counts"] = counts
            fname = (f"UA_{entry['U_A']:.2f}_UB_{entry['U_B']:.4f}"f"_w_{entry['w']:.2f}_L{P['steps']}_{manifest['backend']}.json")
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump({"data": data}, f, indent=2)
            print(f"  dU={entry['delta_U']:<6g} w={entry['w']:<5g} "
                  f"|z|={abs(z_N):.4f} gamma={data['Berry_Phase_pi_wrapped']:+.4f} pi")
            records.append(data)

        else:  # polarization
            density = density_profile_from_counts(counts, L_sites)
            n_A = [round(float(density[2 * j]), 10) for j in range(N_cells)]
            n_B = [round(float(density[2 * j + 1]), 10) for j in range(N_cells)]
            diff = [round(n_A[j] - n_B[j], 10) for j in range(N_cells)]
            rec = {
                "metadata": {"N_cells": N_cells, "L_sites": L_sites,
                             "Total_Electrons": P["electrons"],
                             "v": entry["v"], "w": entry["w"],
                             "t1": entry["v"], "t2": entry["w"],
                             "U_A": entry["U_A"], "U_B": entry["U_B"],
                             "delta_U": entry["delta_U"],
                             "T_A": P["T_A"], "steps": P["steps"],
                             "PBC": False, "backend": manifest["backend"],
                             "shots": shots, "shot_noise": sigma,
                             "job_id": jid},
                "data": {"unit_cell_j": list(range(N_cells)),
                         "n_A": n_A, "n_B": n_B, "n_A_minus_n_B": diff},
            }
            if save_counts:
                rec["counts"] = counts
            fname = (f"Polarization_UA_{entry['U_A']:.1f}_UB_{entry['U_B']:.4f}"f"_L{P['steps']}_{manifest['backend']}.json")
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump(rec, f, indent=2)
            print(f"  dU={entry['delta_U']:<6g} edge(cell 1)={diff[0]:+.4f}")
            records.append(rec)

    print(f"\n  {len(records)} results -> {out_dir}")
    if pending:
        print(f"  {pending} job(s) still running; re-run fetch later.")
    return records, out_dir


def job_status(manifest_path):
    """Print the status of every job in a manifest (no post-processing)."""
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)
    service = get_service()
    print(f"[status] {manifest['experiment']} on {manifest['backend']}")
    for entry in manifest["jobs"]:
        jid = entry.get("job_id")
        if not jid:
            continue
        label = (f"dU={entry['delta_U']:<6g} w={entry['w']:<5g}")
        try:
            print(f"  {label}  {jid}  {job_state(service, jid)}")
        except Exception as e:
            print(f"  {label}  {jid}  ERROR: {e}")


def job_state(service, job_id):
    return str(service.job(job_id).status())