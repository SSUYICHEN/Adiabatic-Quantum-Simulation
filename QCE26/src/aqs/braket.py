"""
braket.py - Run the SSH-Hubbard experiments on AWS Braket (IonQ trapped-ion QPUs).

Relationship to ibm.py / azure.py
---------------------------------
Same pattern, third cloud. The physics circuit and ALL post-processing are
reused unchanged from aqs.core / aqs.observables; only submission differs.
IonQ QPUs (Aria, Forte) are trapped-ion machines with ALL-TO-ALL connectivity,
so -- like Quantinuum and unlike IBM -- the Hubbard CP gates (qubit i <-> i+2N)
and the PBC boundary term need NO SWAP network.

Cost warning
------------
IonQ on Braket bills REAL money (per-task fee + per-shot fee), not credits. A
polarization sweep at 8192 shots can run into hundreds/thousands of USD. This
module therefore:
  * defaults shots LOW (200);
  * always runs on the SV1 cloud simulator first unless you pass a real device;
  * prints an explicit cost estimate and, for a QPU, REQUIRES --i-accept-cost
    before it will submit.

Backends (names via provider.get_backend)
  "SV1"            -- cloud statevector simulator, cheap; validate here first.
  "Aria 1"        -- IonQ Aria-1 QPU (25 qubits).
  "Forte 1"       -- IonQ Forte-1 QPU.
(Exact names come from `check_braket.py`; pass whatever that lists.)

Jobs are async: job_id() is the task ARN. fetch_results() retrieves by ARN via
backend.retrieve_job(), so results can be collected later from the manifest.

Auth: standard AWS credentials (~/.aws/credentials) + region us-east-1, plus
the AmazonBraketFullAccess IAM policy. No code-level auth here.
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

SIMULATOR = "SV1"          # cloud simulator; validate here before any QPU
DEFAULT_QPU = "Forte 1"     # IonQ Forte 1 (adjust to whatever check_braket lists)

# Rough IonQ-on-Braket price points (USD). These change; treat as ballpark and
# confirm against current AWS pricing before spending real money.
_PRICE = {
    "per_task_usd": 0.30,
    "per_shot_usd": 0.03,   # Aria ~0.03; Forte differs
}


# =====================================================================
# Connection
# =====================================================================
_PROVIDER = None


def get_provider():
    global _PROVIDER
    if _PROVIDER is not None:
        return _PROVIDER
    try:
        from qiskit_braket_provider import BraketProvider
    except ImportError as e:
        raise RuntimeError(
            "qiskit-braket-provider not installed:\n"
            "  pip install qiskit-braket-provider amazon-braket-sdk"
        ) from e
    _PROVIDER = BraketProvider()
    return _PROVIDER


def get_backend(name):
    return get_provider().get_backend(name)


def list_targets():
    for b in get_provider().backends():
        try:
            print(f"  {b.name:24s} qubits={getattr(b, 'num_qubits', '?')}")
        except Exception as e:
            print(f"  <error: {e}>")


# =====================================================================
# Circuit construction (identical physics to ibm.py / azure.py)
# =====================================================================
def build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                           T_A=1.0, steps=20, PBC=True):
    model = SSHHModel.from_total_electrons(N_cells, v, w, electrons, PBC=PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B,
                                 T_A=T_A, steps=steps, ramp_U=True)
    qc.measure_all()
    return qc


def _transpile(qc, backend):
    """Transpile for Braket, then ZERO the global phase.

    The decomposed circuit carries a global phase (~pi). Braket cloud devices
    reject the resulting GPhase instruction ("GPhase is not supported"), and
    the adapter's alternative -- silently dropping it -- corrupts the state.
    Global phase has no effect on measurement probabilities, so removing it
    makes the Braket circuit both valid and correct.
    """
    from qiskit import transpile
    isa = transpile(qc, backend=backend, optimization_level=1)
    isa = transpile(isa, basis_gates=["rx", "ry", "rz", "cx"],
                    optimization_level=1)
    isa.global_phase = 0
    return isa


# =====================================================================
# Cost estimate (local, no charge)
# =====================================================================
def estimate_cost_usd(n_tasks, shots):
    return n_tasks * (_PRICE["per_task_usd"] + shots * _PRICE["per_shot_usd"])


def _is_qpu(target):
    return "sv1" not in target.lower() and "sim" not in target.lower()


# =====================================================================
# Result extraction
# =====================================================================
def _counts_of(job, circuit=None):
    result = job.result()
    try:
        return result.get_counts(circuit) if circuit is not None else result.get_counts()
    except Exception:
        return result.get_counts()


# =====================================================================
# Submission (shared core)
# =====================================================================
def _submit(experiment, param_rows, common, target, shots, out_dir,
            estimate_only, accept_cost):
    os.makedirs(out_dir, exist_ok=True)
    n = len(param_rows)
    cost = estimate_cost_usd(n, shots)
    qpu = _is_qpu(target)

    print(f"[{experiment}] target={target}  shots={shots}  {n} tasks")
    print(f"        {common}")
    print(f"        estimated cost ~ ${cost:.2f} USD"
          f"{'  (QPU - REAL money)' if qpu else '  (simulator)'}")

    if estimate_only:
        print("  (estimate only; nothing submitted)")
        return None

    if qpu and not accept_cost:
        print("\n  REFUSING to submit to a QPU without --i-accept-cost.")
        print(f"  Re-run with --i-accept-cost once you've confirmed the ~${cost:.2f} "
              f"budget with your advisor.")
        return None

    backend = get_backend(target)
    manifest = {
        "experiment": experiment, "target": target, "cloud": "aws-braket-ionq",
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
        "estimated_cost_usd": cost,
        "params": {**common, "shots": shots, "target": target},
        "jobs": [],
    }

    for row in param_rows:
        qc = row.pop("_circuit")
        #isa = _transpile(qc, backend)
        job = backend.run(qc, shots=shots)
        entry = dict(row)
        entry["job_id"] = job.job_id()   # task ARN
        manifest["jobs"].append(entry)
        print(f"  {row}  task={job.job_id()}")

    tag = f"{experiment}_{target.replace(' ', '')}_" \
          f"{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
    path = os.path.join(out_dir, f"manifest_{tag}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  manifest -> {path}")
    print(f"  Later:  python aqs_braket.py fetch --manifest {path}")
    return manifest


def submit_polarization(N_cells=6, electrons=14, v=0.5, w=1.5,
                        U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                        T_A=1.0, steps=20, shots=200, binary=1,
                        target=SIMULATOR, out_dir="braket_results/polarization",
                        estimate_only=False, accept_cost=False):
    rows = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                        T_A, steps, PBC=False)
        rows.append({"delta_U": float(delta_U), "U_A": float(U_A),
                     "U_B": float(U_B), "v": float(v), "w": float(w),
                     "_circuit": qc})
    common = {"N_cells": N_cells, "electrons": electrons, "v": v, "w": w,
              "U_A": U_A, "T_A": T_A, "steps": steps, "PBC": False, "binary": binary}
    return _submit("polarization", rows, common, target, shots, out_dir,
                   estimate_only, accept_cost)


def submit_berry(N_cells=6, electrons=12, v=1.0,
                 w_values=(0, 0.5, 0.99, 1.01, 1.5, 2.0),
                 U_A=0.01, delta_U_values=(0.0, 0.01, 0.1, 1.0),
                 T_A=1.0, steps=20, shots=200, binary=1,
                 target=SIMULATOR, out_dir="braket_results/berry",
                 estimate_only=False, accept_cost=False):
    rows = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for w in w_values:
            qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                        T_A, steps, PBC=True)
            rows.append({"delta_U": float(delta_U), "U_A": float(U_A),
                         "U_B": float(U_B), "v": float(v), "w": float(w),
                         "_circuit": qc})
    common = {"N_cells": N_cells, "electrons": electrons, "v": v,
              "U_A": U_A, "T_A": T_A, "steps": steps, "PBC": True, "binary": binary}
    return _submit("berry", rows, common, target, shots, out_dir,
                   estimate_only, accept_cost)


# =====================================================================
# Retrieval + post-processing (paper Eqs. 47-52) -- same schema as ibm/azure
# =====================================================================
def fetch_results(manifest_path, out_dir=None, save_counts=False):
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    exp = manifest["experiment"]
    P = manifest["params"]
    N_cells = P["N_cells"]
    L_sites = 2 * N_cells
    backend = get_backend(manifest["target"])
    out_dir = out_dir or os.path.join(os.path.dirname(manifest_path), "data")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[fetch] {exp}  {len(manifest['jobs'])} tasks  target={manifest['target']}")
    records, pending = [], 0

    for entry in manifest["jobs"]:
        arn = entry.get("job_id")
        if not arn:
            continue
        try:
            job = backend.retrieve_job(arn)
            counts = _counts_of(job)
            if not counts:
                raise ValueError("empty counts")
        except Exception as e:
            print(f"  {arn.split('/')[-1]}  not ready ({e})")
            pending += 1
            continue

        shots = int(sum(counts.values()))
        sigma = shot_noise_estimate(counts)

        if exp == "berry":
            layout = spin_block_layout(N_cells, mode="up_spin")
            z_N = twist_invariant_from_counts(counts, layout)
            data = {
                "L_sites": L_sites, "N_cells": N_cells,
                "U_A": entry["U_A"], "U_B": entry["U_B"],
                "delta_U": entry["delta_U"], "v": entry["v"], "w": entry["w"],
                "t1": entry["v"], "t2": entry["w"],
                "Z_N_real": float(np.real(z_N)), "Z_N_imag": float(np.imag(z_N)),
                "Twist_Amplitude": float(np.abs(z_N)),
                "Berry_Phase_pi": float(np.mod(np.imag(np.log(z_N)), 2 * np.pi) / np.pi),
                "Berry_Phase_pi_wrapped": float(wrap_berry_phase(z_N)),
                "shots": shots, "shot_noise": sigma,
                "job_id": arn, "backend": manifest["target"], "cloud": "aws-ionq",
            }
            if save_counts:
                data["counts"] = counts
            fname = (f"UA_{entry['U_A']:.2f}_UB_{entry['U_B']:.4f}"
                     f"_w_{entry['w']:.2f}_{manifest['target'].replace(' ', '')}.json")
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump({"data": data}, f, indent=2)
            print(f"  dU={entry['delta_U']:<5g} w={entry['w']:<5g} "
                  f"|z|={abs(z_N):.4f} gamma={data['Berry_Phase_pi_wrapped']:+.4f} pi")
            records.append(data)
        else:
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
                             "T_A": P["T_A"], "steps": P["steps"], "PBC": False,
                             "backend": manifest["target"], "cloud": "aws-ionq",
                             "shots": shots, "shot_noise": sigma, "job_id": arn},
                "data": {"unit_cell_j": list(range(N_cells)),
                         "n_A": n_A, "n_B": n_B, "n_A_minus_n_B": diff},
            }
            if save_counts:
                rec["counts"] = counts
            fname = (f"Polarization_UA_{entry['U_A']:.1f}_UB_{entry['U_B']:.4f}"
                     f"_{manifest['target'].replace(' ', '')}.json")
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump(rec, f, indent=2)
            print(f"  dU={entry['delta_U']:<5g} edge(cell 1)={diff[0]:+.4f}")
            records.append(rec)

    print(f"\n  {len(records)} results -> {out_dir}")
    if pending:
        print(f"  {pending} task(s) not ready; re-run fetch later.")
    return records, out_dir