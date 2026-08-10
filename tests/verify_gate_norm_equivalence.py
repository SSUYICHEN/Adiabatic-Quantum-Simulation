"""Baseline capture for the gate-primitive normalisation refactor.

Records both the physics outputs AND the raw circuit instruction stream, so the
refactor can be shown to be exactly behaviour-preserving rather than merely
numerically close.
"""
import json, sys
import numpy as np
from qiskit import transpile

from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend
from aqs.experiments import berry_sweep, polarization_sweep, fidelity_scan

bk = get_backend("qiskit")
out = {}

# --- 1. raw circuit instruction stream (gate name + qubits + angles) ---
def circuit_signature(N, v, w, Ne, pbc, UA, UB, TA, steps):
    m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, UA, UB, T_A=TA, steps=steps, ramp_U=True)
    # decompose so the *underlying* rxx/ryy/cp angles are exposed, not the wrapper label
    d = qc.decompose()
    sig = []
    for instr in d.data:
        sig.append([
            instr.operation.name,
            [d.find_bit(b).index for b in instr.qubits],
            [round(float(p), 12) for p in instr.operation.params],
        ])
    return sig

out["circuit_signature"] = {
    f"N{N}_{'PBC' if p else 'OBC'}_U{UA}-{UB}": circuit_signature(N, 1.0, 1.5, 2 * N, p, UA, UB, 1.0, 4)
    for (N, p, UA, UB) in [(2, True, 0.0, 0.0), (2, True, 0.5, 1.3), (3, False, 1.0, 1.0)]
}

# --- 2. statevectors ---
def sv_of(N, v, w, Ne, pbc, UA, UB, TA, steps):
    m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, UA, UB, T_A=TA, steps=steps, ramp_U=True)
    sv = np.asarray(bk.statevector(qc))
    return sv

out["statevector_hashes"] = {}
out["statevector_samples"] = {}
for tag, args in {
    "N2_PBC_U0": (2, 1.0, 1.5, 4, True, 0.0, 0.0, 1.0, 20),
    "N2_PBC_Uneq": (2, 1.0, 1.5, 4, True, 0.5, 1.3, 1.0, 20),
    "N3_OBC_U1": (3, 0.5, 1.5, 8, False, 1.0, 2.0, 1.0, 20),
}.items():
    sv = sv_of(*args)
    out["statevector_hashes"][tag] = float(np.sum(np.abs(sv) ** 2 * np.arange(len(sv))))
    out["statevector_samples"][tag] = [[float(x.real), float(x.imag)] for x in sv[:8]]

# --- 3. physics observables (the published quantities) ---
out["berry"] = [
    {"w": d["w"], "delta_U": d["delta_U"],
     "gamma": d["Berry_Phase_pi_wrapped"], "z": d["Twist_Amplitude"]}
    for d in berry_sweep(N_cells=3, number_of_electrons=6, v=1.0,
                         w_values=(0.5, 1.0, 1.5, 2.0), U_A=0.01,
                         delta_U_values=(0.0, 0.1), T_A=1.0, steps=40, backend=bk)
]
out["polarization"] = [
    {"delta_U": r["metadata"]["delta_U"], "n_A_minus_n_B": r["data"]["n_A_minus_n_B"]}
    for r in polarization_sweep(N_cells=3, total_electrons=8, v=0.5, w=1.5, U_A=1.0,
                                delta_U_values=(0.0, 0.1, 1.0, 3.0), T_A=1.0,
                                steps=40, backend=bk)
]
res = fidelity_scan(N_cells=3, v=0.5, w=1.5, number_of_electrons=6, U=0.0,
                    T_A_values=(1, 5, 15), steps_list=(10, 40, 100), backend=bk)
out["fidelity"] = {str(k): v for k, v in res["curves"].items()}

json.dump(out, open(sys.argv[1], "w"), indent=2)
print("wrote", sys.argv[1])
