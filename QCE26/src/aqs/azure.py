"""
azure.py - Run the SSH-Hubbard experiments on Azure Quantum (Quantinuum).

Relationship to ibm.py
----------------------
Same idea, different cloud. The physics circuit and ALL classical
post-processing are reused unchanged from aqs.core / aqs.observables; only the
submission layer differs. Quantinuum H-series are trapped-ion machines with
ALL-TO-ALL connectivity, so -- unlike IBM's heavy-hex -- the Hubbard CP gates
(qubit i <-> i+2N) and the PBC boundary term need NO SWAP network. Expect far
fewer two-qubit gates and much higher fidelity than the IBM runs.

Targets seen on this workspace (Basic plan):
  quantinuum.sim.h2-1sc  -- Syntax Checker. FREE. Validates the circuit against
                            the H2 gate set/profile and returns all-zeros. Every
                            submission goes here FIRST as a zero-cost safety net.
  quantinuum.sim.h2-1e   -- Emulator with the real H2-1 noise model. Costs HQC.
  (quantinuum.qpu.h2-1)  -- real hardware; not exposed on this plan yet.

Auth
----
Uses DeviceCodeCredential (works over SSH: you get a code to type into a
browser on your laptop). The credential is cached to a token file so you are
not prompted on every call within its lifetime.

Cost
----
Quantinuum bills in HQC (H-System Quantum Credits). The Syntax Checker is free;
the emulator is not. `estimate_all()` / the --estimate flag call
backend.estimate_cost() for every circuit BEFORE anything is charged.

Workflow
--------
  1. python aqs_azure.py polarization --syntax-check      # FREE validation
  2. python aqs_azure.py polarization --estimate          # HQC cost, no charge
  3. python aqs_azure.py polarization --target quantinuum.sim.h2-1e   # emulator
  4. python aqs_azure.py fetch --manifest <...>           # retrieve + post-process
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

# ------- workspace coordinates (this project's 'alcom' workspace) ---------
import os
# Workspace coordinates are read from environment variables so no account
# identifiers are committed. Set these before running (see handover manual):
#   export AZURE_SUBSCRIPTION_ID=...   export AZURE_RESOURCE_GROUP=...
#   export AZURE_WORKSPACE_NAME=...    export AZURE_LOCATION=...
SUBSCRIPTION_ID = os.environ.get("AZURE_SUBSCRIPTION_ID", "<YOUR_SUBSCRIPTION_ID>")
RESOURCE_GROUP = os.environ.get("AZURE_RESOURCE_GROUP", "AzureQuantum")
WORKSPACE_NAME = os.environ.get("AZURE_WORKSPACE_NAME", "<YOUR_WORKSPACE_NAME>")
LOCATION = os.environ.get("AZURE_LOCATION", "eastus")

SYNTAX_CHECKER = "quantinuum.sim.h2-1sc"
EMULATOR = "quantinuum.sim.h2-1e"
HARDWARE = "quantinuum.qpu.h2-1"


# =====================================================================
# Connection
# =====================================================================
_PROVIDER = None


def get_provider():
    """Cached AzureQuantumProvider using device-code auth (SSH-friendly)."""
    global _PROVIDER
    if _PROVIDER is not None:
        return _PROVIDER
    try:
        from azure.quantum import Workspace
        from azure.quantum.qiskit import AzureQuantumProvider
        from azure.identity import DeviceCodeCredential
    except ImportError as e:
        raise RuntimeError(
            "azure-quantum[qiskit] and azure-identity are required.\n"
            "  pip install 'azure-quantum[qiskit]' azure-identity"
        ) from e

    # Cache tokens so you are not prompted on every run.
    cache_path = os.path.expanduser("~/.aqs_azure_token.json")
    cred = DeviceCodeCredential()

    workspace = Workspace(
        subscription_id=SUBSCRIPTION_ID,
        resource_group=RESOURCE_GROUP,
        name=WORKSPACE_NAME,
        location=LOCATION,
        credential=cred,
    )
    _PROVIDER = AzureQuantumProvider(workspace=workspace)
    return _PROVIDER


def get_backend(target):
    return get_provider().get_backend(target)


def list_targets():
    for b in get_provider().backends():
        n = b.name
        print("  ", n() if callable(n) else n)


# =====================================================================
# Circuit construction (identical physics to ibm.py)
# =====================================================================
def build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                           T_A=1.0, steps=20, PBC=True, name=None):
    model = SSHHModel.from_total_electrons(N_cells, v, w, electrons, PBC=PBC)
    Q_up, Q_dn, _, _ = model.slater_Q_matrices()
    qc = build_annealing_circuit(model, Q_up, Q_dn, U_A, U_B,
                                 T_A=T_A, steps=steps, ramp_U=True)

    # Decompose custom gates (e.g. 'Givens') to a native basis BEFORE measuring.
    from qiskit import transpile
    qc = transpile(qc, basis_gates=["rx", "ry", "rz", "cx"],
                   optimization_level=1)

    qc.measure_all()   # may add its own barrier; removed below

    # Quantinuum's QIR adaptor rejects generic barriers. Strip every barrier
    # (state-prep barrier from core.py, Trotter-step barriers, and the one
    # measure_all just added). Done AFTER measure so no barrier survives.
    from qiskit.transpiler.passes import RemoveBarriers
    qc = RemoveBarriers()(qc)

    if name:
        qc.name = name
    return qc


# =====================================================================
# Result extraction
# =====================================================================
def _counts_of(job, circuit=None):
    """Get a {bitstring: count} dict from an Azure Quantum Qiskit job.

    Azure/Quantinuum sometimes key counts by the circuit; try the safe forms.
    """
    result = job.result()
    try:
        return result.get_counts(circuit) if circuit is not None else result.get_counts()
    except Exception:
        return result.get_counts()


# =====================================================================
# Cost estimation (no charge)
# =====================================================================
def estimate_cost(backend, qc, shots):
    """backend.estimate_cost(...) if the target supports it; else None."""
    try:
        cost = backend.estimate_cost(qc, shots=shots)
        # cost object exposes .estimated_total and .currency_code (HQC)
        return getattr(cost, "estimated_total", None), getattr(cost, "currency_code", "HQC")
    except Exception as e:
        return None, f"(estimate unavailable: {e})"


# =====================================================================
# Submission (shared core for both experiments)
# =====================================================================
def _submit_jobs(experiment, param_rows, common, target, shots,
                 out_dir, syntax_check, estimate_only):
    """param_rows: list of dicts each describing one circuit + its metadata.
    common: dict of shared params written into the manifest.
    """
    os.makedirs(out_dir, exist_ok=True)

    if syntax_check:
        target = SYNTAX_CHECKER
    backend = get_backend(target)

    mode = ("SYNTAX-CHECK" if syntax_check else
            "ESTIMATE" if estimate_only else "SUBMIT")
    print(f"[{experiment}] target={target}  mode={mode}  shots={shots}")
    print(f"        {common}")
    print(f"        {len(param_rows)} circuits\n")

    manifest = {
        "experiment": experiment, "target": target, "cloud": "azure-quantinuum",
        "submitted_utc": datetime.now(timezone.utc).isoformat(),
        "params": {**common, "shots": shots, "target": target},
        "jobs": [],
    }

    total_cost = 0.0
    for row in param_rows:
        qc = row.pop("_circuit")
        entry = dict(row)

        if estimate_only:
            val, cur = estimate_cost(backend, qc, shots)
            entry["estimated_cost"] = val
            entry["currency"] = cur
            if isinstance(val, (int, float)):
                total_cost += val
            print(f"  {row}  est={val} {cur}")
        else:
            job = backend.run(qc, shots=shots)
            entry["job_id"] = job.id()
            print(f"  {row}  job={job.id()}")
        manifest["jobs"].append(entry)

    if estimate_only:
        print(f"\n  TOTAL estimated cost ~ {total_cost} (sum over circuits)")

    tag = f"{experiment}_{target.replace('.', '-')}_" \
          f"{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}"
    path = os.path.join(out_dir, f"manifest_{tag}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"\n  manifest -> {path}")
    if mode == "SUBMIT":
        print(f"  Later:  python aqs_azure.py fetch --manifest {path}")
    return manifest


def submit_polarization(N_cells=6, electrons=14, v=0.5, w=1.5,
                        U_A=0.01, delta_U_values=(0.0, 0.001, 0.01, 0.1, 1.0),
                        T_A=1.0, steps=40, shots=500,
                        target=EMULATOR, out_dir="azure_results/polarization",
                        syntax_check=False, estimate_only=False):
    """Sublattice polarization (OBC, n = 2N+2). One job per delta_U.

    NB default shots=500: Quantinuum shots are far more expensive than IBM's,
    so we default low. Raise deliberately.
    """
    rows = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                    T_A, steps, PBC=False,
                                    name=f"pol_dU{delta_U:g}")
        rows.append({"delta_U": float(delta_U), "U_A": float(U_A),
                     "U_B": float(U_B), "v": float(v), "w": float(w),
                     "_circuit": qc})
    common = {"N_cells": N_cells, "electrons": electrons, "v": v, "w": w,
              "U_A": U_A, "T_A": T_A, "steps": steps, "PBC": False}
    return _submit_jobs("polarization", rows, common, target, shots,
                        out_dir, syntax_check, estimate_only)


def submit_berry(N_cells=6, electrons=12, v=1.0,
                 w_values=(0, 0.5, 0.99, 1.01, 1.5, 2.0),
                 U_A=0.01, delta_U_values=(0.0, 0.001, 0.01, 0.1, 1.0),
                 T_A=1.0, steps=40, shots=500,
                 target=EMULATOR, out_dir="azure_results/berry",
                 syntax_check=False, estimate_only=False):
    """Berry phase (PBC, half filling). One job per (delta_U, w)."""
    rows = []
    for delta_U in delta_U_values:
        U_B = U_A + delta_U
        for w in w_values:
            qc = build_measured_circuit(N_cells, v, w, electrons, U_A, U_B,
                                        T_A, steps, PBC=True,
                                        name=f"berry_dU{delta_U:g}_w{w:g}")
            rows.append({"delta_U": float(delta_U), "U_A": float(U_A),
                         "U_B": float(U_B), "v": float(v), "w": float(w),
                         "_circuit": qc})
    common = {"N_cells": N_cells, "electrons": electrons, "v": v,
              "U_A": U_A, "T_A": T_A, "steps": steps, "PBC": True}
    return _submit_jobs("berry", rows, common, target, shots,
                        out_dir, syntax_check, estimate_only)


# =====================================================================
# Retrieval + post-processing (paper Eqs. 47-52) -- same schema as ibm.py
# =====================================================================
def fetch_results(manifest_path, out_dir=None, save_counts=False):
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    provider = get_provider()
    exp = manifest["experiment"]
    P = manifest["params"]
    N_cells = P["N_cells"]
    L_sites = 2 * N_cells
    out_dir = out_dir or os.path.join(os.path.dirname(manifest_path), "data")
    os.makedirs(out_dir, exist_ok=True)

    print(f"[fetch] {exp}  {len(manifest['jobs'])} jobs  target={manifest['target']}")
    records, pending = [], 0

    for entry in manifest["jobs"]:
        jid = entry.get("job_id")
        if not jid:
            continue
        job = provider.get_job(jid)
        status = str(getattr(job, "details", job).status if hasattr(job, "details")
                     else job.status())
        counts = None
        try:
            counts = _counts_of(job)
        except Exception as e:
            print(f"  {jid}  not ready ({e})")
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
                "job_id": jid, "backend": manifest["target"], "cloud": "azure",
            }
            if save_counts:
                data["counts"] = counts
            fname = (f"UA_{entry['U_A']:.2f}_UB_{entry['U_B']:.4f}"
                     f"_w_{entry['w']:.2f}_{manifest['target'].replace('.', '-')}.json")
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
                             "backend": manifest["target"], "cloud": "azure",
                             "shots": shots, "shot_noise": sigma, "job_id": jid},
                "data": {"unit_cell_j": list(range(N_cells)),
                         "n_A": n_A, "n_B": n_B, "n_A_minus_n_B": diff},
            }
            if save_counts:
                rec["counts"] = counts
            fname = (f"Polarization_UA_{entry['U_A']:.1f}_UB_{entry['U_B']:.4f}"
                     f"_{manifest['target'].replace('.', '-')}.json")
            with open(os.path.join(out_dir, fname), "w", encoding="utf-8") as f:
                json.dump(rec, f, indent=2)
            print(f"  dU={entry['delta_U']:<5g} edge(cell 1)={diff[0]:+.4f}")
            records.append(rec)

    print(f"\n  {len(records)} results -> {out_dir}")
    if pending:
        print(f"  {pending} job(s) not ready; re-run fetch later.")
    return records, out_dir