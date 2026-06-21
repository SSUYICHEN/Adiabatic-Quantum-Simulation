"""
Generate example Hamiltonian files for `aqs measure`.

Produces a spinless SSH chain as an OpenFermion fermion-operator JSON, in the
self-describing format that aqs.hamiltonian understands. Run:

    python examples/make_example_hamiltonian.py

then measure it:

    aqs measure --hamiltonian examples/ssh_spinless_N3_topological.json \
                --property berry,polarization
"""

import json
import os


def ssh_spinless(N_cells, v, w, pbc=False):
    """Return a Hamiltonian spec dict for a spinless SSH chain.

    Sites 0..L-1 (L=2*N_cells). Bond i<->i+1 has amplitude v if i even else w.
    Layout 'spinless': cell j = qubits (2j, 2j+1), A=2j, B=2j+1.
    """
    L = 2 * N_cells
    terms = []
    for i in range(L):
        if not pbc and i == L - 1:
            continue
        j = (i + 1) % L
        t = v if i % 2 == 0 else w
        # -t (c_i^dag c_j + c_j^dag c_i)
        terms.append({"coeff": [-float(t), 0.0], "ops": f"{i}^ {j}"})
        terms.append({"coeff": [-float(t), 0.0], "ops": f"{j}^ {i}"})
    return {
        "type": "hamiltonian",
        "format": "fermion_operator",
        "metadata": {
            "n_qubits": L,
            "n_cells": N_cells,
            "layout": "spinless",
            "boundary": "PBC" if pbc else "OBC",
            "v": v, "w": w,
        },
        "terms": terms,
    }


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    cases = {
        "ssh_spinless_N3_topological.json": ssh_spinless(3, v=0.5, w=1.5, pbc=False),
        "ssh_spinless_N3_trivial.json":     ssh_spinless(3, v=1.5, w=0.5, pbc=False),
    }
    for name, spec in cases.items():
        path = os.path.join(here, name)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(spec, f, indent=2)
        print("wrote", path)
