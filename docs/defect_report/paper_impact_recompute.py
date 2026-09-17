"""Recompute the paper's Fig. 3 / Fig. 4 data (N=6, T=1, L=40) with the pre-fix
Givens (legacy) and the fixed Givens, plus the (-U_A, -dU) equivalence check.
Writes paper_impact.json next to this script.
Run: uv run --no-sync python docs/defect_report/paper_impact_recompute.py  (then paper_impact_plot.py)"""
import os, sys, json, time, contextlib, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from qiskit import QuantumCircuit
import aqs.core as core
from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend
from aqs.observables import twist_invariant, wrap_berry_phase, density_profile

OUT = pathlib.Path(__file__).with_name("paper_impact.json")
try:
    bk = get_backend("cudaq", target="nvidia-fp64"); print("backend: cudaq nvidia-fp64")
except Exception as e:
    print("fp64 target unavailable:", e); bk = get_backend("cudaq"); print("backend: cudaq", bk.target if hasattr(bk, "target") else "")

@contextlib.contextmanager
def patched(module, name, repl):
    old = getattr(module, name); setattr(module, name, repl)
    try: yield
    finally: setattr(module, name, old)

def legacy_givens(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1); gv.rz(-phi, 0); gv.cx(0, 1); gv.cry(-2.0 * theta, 1, 0); gv.cx(0, 1)
    return gv.to_instruction()

def run(N, v, w, Ne, pbc, U_A, U_B, legacy):
    m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    ctx = patched(core, "_givens_instruction", legacy_givens) if legacy else contextlib.nullcontext()
    with ctx:
        qc = build_annealing_circuit(m, Qu, Qd, U_A, U_B, T_A=1.0, steps=40, ramp_U=True)
    sv = bk.statevector(qc)
    out = {}
    if pbc:
        for mode in ("up_spin", "total"):
            z = twist_invariant(sv, 2 * N, mode=mode)
            out[mode] = {"gamma_pi_wrapped": float(wrap_berry_phase(z)), "abs_z": float(abs(z))}
    dens = density_profile(sv, 2 * N)
    out["pol"] = [float(dens[2 * j] - dens[2 * j + 1]) for j in range(N)]
    return out

N = 6
res = {"meta": {"N": N, "T": 1.0, "L": 40, "backend": bk.name}}
t0 = time.time()

# ---- Fig 3 top: v=1, U_A=0.01, w sweep, dU set
W = [0, 0.25, 0.5, 0.75, 0.99, 1.01, 1.25, 1.5, 1.75, 2.0]
DU3 = [0.0, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
res["fig3_top"] = []
for dU in DU3:
    for w in W:
        rec = {"w": w, "delta_U": dU}
        for tag, leg in (("before", True), ("after", False)):
            rec[tag] = run(N, 1.0, w, 12, True, 0.01, 0.01 + dU, leg)
        res["fig3_top"].append(rec); print(f"fig3 top dU={dU} w={w} done {time.time()-t0:.0f}s", flush=True)

# ---- Fig 3 bottom: w=2.0 (and 1.5), dU extended to 3
DU3b = [0.0, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
res["fig3_bottom"] = []
for w in (1.5, 2.0):
    for dU in DU3b:
        rec = {"w": w, "delta_U": dU}
        for tag, leg in (("before", True), ("after", False)):
            rec[tag] = run(N, 1.0, w, 12, True, 0.01, 0.01 + dU, leg)
        # equivalence: fixed code with attractive (-U_A, -U_B) should match legacy
        rec["after_negU"] = run(N, 1.0, w, 12, True, -0.01, -(0.01 + dU), False)
        res["fig3_bottom"].append(rec); print(f"fig3 bottom w={w} dU={dU} done {time.time()-t0:.0f}s", flush=True)

# ---- Fig 4: v=0.1, w=1.0, OBC, 14 electrons, U_A=0.01
DU4 = [0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0]
res["fig4"] = []
for dU in DU4:
    rec = {"delta_U": dU}
    for tag, leg in (("before", True), ("after", False)):
        rec[tag] = run(N, 0.1, 1.0, 14, False, 0.01, 0.01 + dU, leg)
    rec["after_negU"] = run(N, 0.1, 1.0, 14, False, -0.01, -(0.01 + dU), False)
    res["fig4"].append(rec); print(f"fig4 dU={dU} done {time.time()-t0:.0f}s", flush=True)

json.dump(res, open(OUT, "w"), indent=1)
print("wrote", OUT, f"total {time.time()-t0:.0f}s")
