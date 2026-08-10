"""Full before/after matrix for every fix on this branch.

Self-contained so it runs unchanged on main and on HEAD:

    git checkout main && python tests/verify_consistency_vs_main.py before.json
    git checkout -    && python tests/verify_consistency_vs_main.py after.json
    python tests/verify_consistency_vs_main.py --compare before.json after.json
"""
import json
import sys

import numpy as np


def _collect():
    import openfermion as of
    from aqs.core import SSHHModel, build_annealing_circuit
    from aqs.backends import get_backend, selftest
    from aqs.experiments import berry_sweep, polarization_sweep, fidelity_scan
    from aqs.hamiltonian import measure_hamiltonian

    bk = get_backend("qiskit")
    out = {}

    # ---- backend selftests -------------------------------------------
    out["selftest"] = {}
    for name in ("qiskit", "cudaq"):
        try:
            out["selftest"][name] = float(selftest(name)[1])
        except Exception as e:
            out["selftest"][name] = f"unavailable: {type(e).__name__}"

    # ---- state preparation energy ------------------------------------
    def bitrev(M, n):
        idx = np.arange(1 << n, dtype=np.int64)
        p = np.zeros_like(idx)
        for b in range(n):
            p |= ((idx >> b) & 1) << (n - 1 - b)
        return M[np.ix_(p, p)]

    out["prep_energy"] = {}
    for N, v, w, Ne, pbc in [(3, 1.0, 1.5, 6, True), (3, 0.5, 1.5, 8, False)]:
        m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
        Qu, Qd, nu, nd = m.slater_Q_matrices()
        sv = bk.statevector(
            build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0, ramp_U=False))
        L = m.L
        ev = np.sort(np.linalg.eigvalsh(m.H[:L, :L]))
        op = of.FermionOperator()
        for i in range(2 * L):
            for j in range(2 * L):
                if abs(m.H[i, j]) > 1e-12:
                    op += of.FermionOperator(f"{i}^ {j}", complex(m.H[i, j]))
        Hs = bitrev(of.get_sparse_operator(op, n_qubits=2 * L).toarray(), 2 * L)
        out["prep_energy"][f"N{N}_{'PBC' if pbc else 'OBC'}"] = {
            "E_prep": float(np.real(np.vdot(sv, Hs @ sv))),
            "E_ground": float(ev[:nu].sum() + ev[:nd].sum()),
        }

    # ---- circuit-level observables -----------------------------------
    out["berry"] = [
        {"w": d["w"], "dU": d["delta_U"], "gamma": d["Berry_Phase_pi_wrapped"],
         "z": d["Twist_Amplitude"]}
        for d in berry_sweep(N_cells=3, number_of_electrons=6, v=1.0,
                             w_values=(0.5, 1.0, 1.5, 2.0), U_A=0.01,
                             delta_U_values=(0.0, 0.1), T_A=1.0, steps=40,
                             backend=bk)]
    out["polarization"] = [
        {"dU": r["metadata"]["delta_U"], "prof": r["data"]["n_A_minus_n_B"]}
        for r in polarization_sweep(N_cells=3, total_electrons=8, v=0.5, w=1.5,
                                    U_A=1.0, delta_U_values=(0.0, 0.1, 1.0, 3.0),
                                    T_A=1.0, steps=40, backend=bk)]
    out["fidelity"] = {
        str(k): v for k, v in
        fidelity_scan(N_cells=3, v=0.5, w=1.5, number_of_electrons=6, U=0.0,
                      T_A_values=(1, 5, 15), steps_list=(10, 40, 100),
                      backend=bk)["curves"].items()}

    # ---- exact-diagonalisation path ----------------------------------
    out["measure"] = {}
    import os
    for name in ("ssh_spinless_N3_topological", "ssh_spinless_N3_trivial"):
        p = os.path.join("examples", f"{name}.json")
        if not os.path.exists(p):
            continue
        r = measure_hamiltonian(p, ["berry", "polarization"], bk)["results"]
        out["measure"][name] = {
            "gamma": r["berry"]["Berry_Phase_pi_wrapped"],
            "z": r["berry"]["Twist_Amplitude"],
            "prof": r["polarization"]["n_A_minus_n_B"],
        }

    # ---- an ASYMMETRIC operator file (the shipped ones cannot detect
    #      the openfermion ordering bug) --------------------------------
    import tempfile
    from aqs.hamiltonian import load_hamiltonian
    from aqs.observables import density_profile_qubits
    n = 4
    terms = [{"coeff": [-8.0, 0.0], "ops": "0^ 0"}]
    for i in range(n - 1):
        terms += [{"coeff": [-0.2, 0.0], "ops": f"{i}^ {i+1}"},
                  {"coeff": [-0.2, 0.0], "ops": f"{i+1}^ {i}"}]
    spec = {"type": "hamiltonian", "format": "fermion_operator",
            "metadata": {"n_qubits": n, "n_cells": 2, "layout": "spinless"},
            "terms": terms}
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "asym.json")
        json.dump(spec, open(p, "w"))
        H, lay, meta = load_hamiltonian(p)
        psi = np.asarray(bk.ground_state(H)).ravel()
        psi = psi / np.linalg.norm(psi)
        dens = density_profile_qubits(psi, n)
        out["asymmetric_measure"] = {
            "density": [float(x) for x in dens],
            "argmax_qubit": int(np.argmax(dens)),
            "ground_energy": float(np.real(np.vdot(psi, H @ psi))),
        }
    return out


def _cmp(label, a, b, tol=0.0):
    """Return (status, detail). status in {same, rounding, CHANGED}."""
    fa, fb = np.atleast_1d(np.array(a, dtype=float)), np.atleast_1d(np.array(b, dtype=float))
    if fa.shape != fb.shape:
        return "CHANGED", "shape"
    d = float(np.max(np.abs(fa - fb)))
    if d == 0.0:
        return "same", "bit-identical"
    if d <= 1e-12:
        return "rounding", f"max|d|={d:.1e}"
    return "CHANGED", f"max|d|={d:.4g}"


def compare(before, after):
    b, a = json.load(open(before)), json.load(open(after))
    rows = []
    for k in ("qiskit", "cudaq"):
        vb, va = b["selftest"].get(k), a["selftest"].get(k)
        if isinstance(vb, float) and isinstance(va, float):
            rows.append((f"selftest[{k}]", *_cmp(k, vb, va), f"{vb:.8f} -> {va:.8f}"))
        else:
            rows.append((f"selftest[{k}]", "n/a", f"{vb} -> {va}", ""))
    for k in b["prep_energy"]:
        vb, va = b["prep_energy"][k], a["prep_energy"][k]
        s, det = _cmp(k, vb["E_prep"], va["E_prep"])
        rows.append((f"prep_energy[{k}]", s, det,
                     f"{vb['E_prep']:+.6f} -> {va['E_prep']:+.6f} (ground {va['E_ground']:+.6f})"))
    rows.append(("berry gamma (dU=0)", *_cmp("g",
        [x["gamma"] for x in b["berry"] if x["dU"] == 0],
        [x["gamma"] for x in a["berry"] if x["dU"] == 0]), ""))
    rows.append(("berry gamma (dU!=0)", *_cmp("g",
        [x["gamma"] for x in b["berry"] if x["dU"] != 0],
        [x["gamma"] for x in a["berry"] if x["dU"] != 0]), ""))
    rows.append(("polarization profiles", *_cmp("p",
        [v for x in b["polarization"] for v in x["prof"]],
        [v for x in a["polarization"] for v in x["prof"]]), ""))
    rows.append(("fidelity scan", *_cmp("f",
        [v for k in b["fidelity"] for v in b["fidelity"][k]],
        [v for k in a["fidelity"] for v in a["fidelity"][k]]), ""))
    for name in b.get("measure", {}):
        rows.append((f"measure[{name}] gamma", *_cmp("g",
            b["measure"][name]["gamma"], a["measure"][name]["gamma"]), ""))
        rows.append((f"measure[{name}] prof", *_cmp("p",
            b["measure"][name]["prof"], a["measure"][name]["prof"]), ""))
    ab, aa = b["asymmetric_measure"], a["asymmetric_measure"]
    rows.append(("asymmetric H: argmax qubit",
                 "same" if ab["argmax_qubit"] == aa["argmax_qubit"] else "CHANGED",
                 f"{ab['argmax_qubit']} -> {aa['argmax_qubit']}", "(trap on mode 0; want 0)"))
    rows.append(("asymmetric H: ground energy",
                 *_cmp("e", ab["ground_energy"], aa["ground_energy"]),
                 "spectrum must be preserved"))

    w = max(len(r[0]) for r in rows) + 2
    print(f"{'quantity':<{w}} {'status':<10} {'detail':<24} notes")
    print("-" * (w + 60))
    for name, status, detail, note in rows:
        print(f"{name:<{w}} {status:<10} {detail:<24} {note}")


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        json.dump(_collect(), open(sys.argv[1], "w"), indent=2)
        print("wrote", sys.argv[1])
