"""Recompute the paper's Fig. 2 (Trotter validation, N=6, v=0.5, w=1.5, 12 e, PBC)
with the pre-fix Givens (legacy) and the fixed Givens. Writes paper_fig2.json next to
this script. Run: uv run --no-sync python docs/defect_report/paper_fig2_recompute.py"""
import os, json, time, contextlib, pathlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from qiskit import QuantumCircuit
import aqs.core as core
from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend

OUT = pathlib.Path(__file__).with_name("paper_fig2.json")
try:
    bk = get_backend("cudaq", target="nvidia-fp64"); print("backend: cudaq nvidia-fp64")
except Exception as e:
    print("fp64 unavailable:", e); bk = get_backend("cudaq")

@contextlib.contextmanager
def patched(module, name, repl):
    old = getattr(module, name); setattr(module, name, repl)
    try: yield
    finally: setattr(module, name, old)

def legacy_givens(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1); gv.rz(-phi, 0); gv.cx(0, 1); gv.cry(-2.0 * theta, 1, 0); gv.cx(0, 1)
    return gv.to_instruction()

def fid(a, b): return float(abs(np.vdot(a, b)) ** 2)

N, v, w, Ne = 6, 0.5, 1.5, 12
T_LIST = (1, 5, 15, 80)
L_LIST = (1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150)
res = {"meta": {"N": N, "v": v, "w": w, "electrons": Ne, "PBC": True, "T": list(T_LIST), "L": list(L_LIST), "backend": bk.name}}
t0 = time.time()
for U, tag in ((0.0, "U0"), (1.0, "U1")):
    res[tag] = {}
    for version, legacy in (("before", True), ("after", False)):
        ctx = patched(core, "_givens_instruction", legacy_givens) if legacy else contextlib.nullcontext()
        with ctx:
            m = SSHHModel.from_total_electrons(N, v, w, Ne, True)
            Qu, Qd, _, _ = m.slater_Q_matrices()
            sv0 = bk.statevector(build_annealing_circuit(m, Qu, Qd, U, U, T_A=0, steps=0, ramp_U=False))
            curves = {}
            for T in T_LIST:
                f_ground, f_adj, prev = [], [], sv0
                for L in L_LIST:
                    sv = bk.statevector(build_annealing_circuit(m, Qu, Qd, U, U, T_A=T, steps=L, ramp_U=False))
                    f_ground.append(fid(sv0, sv)); f_adj.append(fid(prev, sv)); prev = sv
                    print(f"{tag} {version} T={T} L={L} ground={f_ground[-1]:.6f} adjacent={f_adj[-1]:.6f} {time.time()-t0:.0f}s", flush=True)
                curves[str(T)] = {"ground": f_ground, "adjacent": f_adj}
        res[tag][version] = curves
        json.dump(res, open(OUT, "w"), indent=1)
print("wrote", OUT, f"{time.time()-t0:.0f}s")
