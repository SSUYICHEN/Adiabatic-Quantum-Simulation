"""Quantify the impact of the _givens_instruction state-preparation fix.

Self-contained: the openfermion -> qiskit bit reversal is done inline rather
than imported, so this runs on any revision of the tree.

    python tests/verify_givens_impact.py before.json     # pre-fix revision
    python tests/verify_givens_impact.py after.json      # post-fix revision
"""
import json, sys
import numpy as np
import openfermion as of

from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend
from aqs.experiments import berry_sweep, polarization_sweep, fidelity_scan

bk = get_backend("qiskit")
out = {}


def _bit_reverse(H, n):
    """openfermion places qubit q at bit n-1-q; qiskit places it at bit q."""
    idx = np.arange(1 << n, dtype=np.int64)
    perm = np.zeros_like(idx)
    for p in range(n):
        perm |= ((idx >> p) & 1) << (n - 1 - p)
    return H[np.ix_(perm, perm)]


def prep_energy(N, v, w, Ne, pbc):
    """<H_SSH> of the prepared state, against the lowest/highest fillings.

    A correct Slater-determinant preparation must land exactly on the
    lowest-fill energy. On a chiral-symmetric SSH chain the highest fill is its
    exact negative, so a sign error in the Givens rotation shows up as +E.
    """
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
    Hs = _bit_reverse(
        of.get_sparse_operator(op, n_qubits=2 * L).toarray(), 2 * L)
    return {
        "E_prep": float(np.real(np.vdot(sv, Hs @ sv))),
        "E_lowest_fill": float(ev[:nu].sum() + ev[:nd].sum()),
        "E_highest_fill": float(ev[::-1][:nu].sum() + ev[::-1][:nd].sum()),
    }


out["prep_energy"] = {
    f"N{N}_v{v}_w{w}_Ne{Ne}_{'PBC' if p else 'OBC'}": prep_energy(N, v, w, Ne, p)
    for (N, v, w, Ne, p) in [(3, 1.0, 1.5, 6, True),
                             (3, 0.5, 1.5, 8, False),
                             (4, 1.0, 0.5, 8, True)]
}

out["berry"] = {}
for N, Ne in [(3, 6), (4, 8)]:
    out["berry"][f"N{N}"] = [
        {"w": d["w"], "delta_U": d["delta_U"],
         "gamma": d["Berry_Phase_pi_wrapped"], "z": d["Twist_Amplitude"]}
        for d in berry_sweep(N_cells=N, number_of_electrons=Ne, v=1.0,
                             w_values=(0.5, 1.0, 1.5, 2.0), U_A=0.01,
                             delta_U_values=(0.0, 0.1), T_A=1.0, steps=40,
                             backend=bk)
    ]

out["polarization"] = {}
for N in (3, 4):
    out["polarization"][f"N{N}"] = [
        {"delta_U": r["metadata"]["delta_U"],
         "n_A_minus_n_B": r["data"]["n_A_minus_n_B"]}
        for r in polarization_sweep(N_cells=N, total_electrons=2 * N + 2,
                                    v=0.5, w=1.5, U_A=1.0,
                                    delta_U_values=(0.0, 0.1, 1.0, 3.0),
                                    T_A=1.0, steps=40, backend=bk)
    ]

res = fidelity_scan(N_cells=3, v=0.5, w=1.5, number_of_electrons=6, U=0.0,
                    T_A_values=(1, 5, 15), steps_list=(10, 40, 100), backend=bk)
out["fidelity_N3"] = {str(k): v for k, v in res["curves"].items()}

json.dump(out, open(sys.argv[1], "w"), indent=2)
print("wrote", sys.argv[1])
