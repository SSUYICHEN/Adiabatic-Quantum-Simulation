"""
hamiltonian.py - Load an arbitrary Hamiltonian from a file and measure
topological properties of its ground state.

File format (JSON, self-describing)
-----------------------------------
{
  "type": "hamiltonian",
  "format": "fermion_operator" | "qubit_operator" | "dense_matrix" | "sparse_matrix",

  "metadata": {
     "n_qubits": 8,
     "n_cells": 4,
     "layout": "spin_block" | "spinless" | "custom",
     "mode": "up_spin" | "total",          # spin_block only
     "cell_qubits": [[0,1],[2,3],...],      # custom only
     "A_qubits": [0,2,...], "B_qubits": [1,3,...]   # optional, for polarization
  },

  # ---- one of the following, depending on "format" ----
  "terms": [ {"coeff": [re, im], "ops": "0^ 1"}, ... ],   # operator formats
  "operator_string": "-1.0 [0^ 1] + -1.0 [1^ 0]",         # operator formats (alt)
  "matrix_file": "H.npy",                                  # matrix formats
  "data": [[...], ...]                                     # matrix formats (inline, alt)
}

`fermion_operator` terms are mapped to qubits via Jordan-Wigner (OpenFermion).
`qubit_operator` terms are Pauli strings like "X0 Z1".
Ground state is obtained by exact diagonalisation through the chosen backend
(CPU scipy or GPU cupy).
"""

from __future__ import annotations

import json
import os

import numpy as np

from .observables import spin_block_layout, PROPERTIES


# --------------------------------------------------------------------- layout
def _build_layout(meta):
    kind = meta.get("layout", "spin_block")
    if kind == "spin_block":
        return spin_block_layout(meta["n_cells"], meta.get("mode", "up_spin"))
    n_cells = meta["n_cells"]
    if kind == "spinless":
        cell_qubits = [[2 * j, 2 * j + 1] for j in range(n_cells)]
        A = [2 * j for j in range(n_cells)]
        B = [2 * j + 1 for j in range(n_cells)]
    elif kind == "custom":
        cell_qubits = meta["cell_qubits"]
        A = meta.get("A_qubits", [c[0] for c in cell_qubits])
        B = meta.get("B_qubits", [c[1] for c in cell_qubits if len(c) > 1])
    else:
        raise ValueError(f"unknown layout '{kind}'")
    return {
        "n_qubits": meta["n_qubits"], "n_cells": n_cells,
        "cell_qubits": cell_qubits, "A_qubits": A, "B_qubits": B,
    }


# --------------------------------------------------------------------- operators
def _coeff(c):
    if isinstance(c, (list, tuple)):
        return complex(c[0], c[1] if len(c) > 1 else 0.0)
    return complex(c)


def _build_sparse(spec, base_dir):
    fmt = spec["format"]
    meta = spec["metadata"]
    n_qubits = meta["n_qubits"]

    if fmt in ("fermion_operator", "qubit_operator"):
        from openfermion import FermionOperator, QubitOperator, get_sparse_operator
        Op = FermionOperator if fmt == "fermion_operator" else QubitOperator
        op = Op()
        if "operator_string" in spec:
            op = Op(spec["operator_string"])
        else:
            for term in spec["terms"]:
                op += Op(term.get("ops", ""), _coeff(term["coeff"]))
        return get_sparse_operator(op, n_qubits=n_qubits)

    if fmt in ("dense_matrix", "sparse_matrix"):
        import scipy.sparse as sp
        if "matrix_file" in spec:
            path = os.path.join(base_dir, spec["matrix_file"])
            if path.endswith(".npz"):
                M = sp.load_npz(path)
            else:
                M = np.load(path)
        else:
            M = np.asarray(spec["data"], dtype=complex)
        return sp.csr_matrix(M)

    raise ValueError(f"unknown format '{fmt}'")


# --------------------------------------------------------------------- public API
def load_hamiltonian(path):
    """Read a Hamiltonian JSON file. Returns (sparse_H, layout, metadata)."""
    with open(path, "r", encoding="utf-8") as f:
        spec = json.load(f)
    if spec.get("type") != "hamiltonian":
        raise ValueError("file is not a Hamiltonian spec (\"type\": \"hamiltonian\")")
    base_dir = os.path.dirname(os.path.abspath(path))
    H = _build_sparse(spec, base_dir)
    layout = _build_layout(spec["metadata"])
    return H, layout, spec["metadata"]


def measure_hamiltonian(path, properties, backend):
    """Load H, get its ground state via `backend`, and measure `properties`
    (a list of names from observables.PROPERTIES). Returns a result dict."""
    H, layout, meta = load_hamiltonian(path)
    psi = backend.ground_state(H)
    psi = np.asarray(psi).ravel()
    psi = psi / np.linalg.norm(psi)

    results = {}
    for name in properties:
        if name not in PROPERTIES:
            raise ValueError(f"unknown property '{name}' "
                             f"(available: {sorted(PROPERTIES)})")
        results[name] = PROPERTIES[name](psi, layout)
    return {
        "source": "hamiltonian_file", "file": os.path.basename(path),
        "backend": backend.name, "metadata": meta, "results": results,
    }
