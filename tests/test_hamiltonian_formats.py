"""Tests for the Hamiltonian file loader: all four formats, layouts, errors."""
import json
import os
import tempfile

import numpy as np
import pytest
import scipy.sparse as sp

from aqs.backends import get_backend
from aqs.hamiltonian import load_hamiltonian, measure_hamiltonian

_BK = get_backend("qiskit")


def _write(spec, d, name="h.json"):
    path = os.path.join(d, name)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(spec, f)
    return path


def _meta(n_qubits=4, n_cells=2, layout="spinless", **kw):
    m = {"n_qubits": n_qubits, "n_cells": n_cells, "layout": layout}
    m.update(kw)
    return m


# ------------------------------------------------------------- operator forms
def test_terms_and_operator_string_agree():
    """The two ways of spelling the same FermionOperator must load identically."""
    with tempfile.TemporaryDirectory() as d:
        a = _write({"type": "hamiltonian", "format": "fermion_operator",
                    "metadata": _meta(),
                    "terms": [{"coeff": [-1.0, 0.0], "ops": "0^ 1"},
                              {"coeff": [-1.0, 0.0], "ops": "1^ 0"}]}, d, "a.json")
        b = _write({"type": "hamiltonian", "format": "fermion_operator",
                    "metadata": _meta(),
                    "operator_string": "-1.0 [0^ 1] + -1.0 [1^ 0]"}, d, "b.json")
        Ha, _, _ = load_hamiltonian(a)
        Hb, _, _ = load_hamiltonian(b)
        assert np.allclose(Ha.toarray(), Hb.toarray())


def test_complex_coefficients_round_trip():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "fermion_operator",
                    "metadata": _meta(),
                    "terms": [{"coeff": [-0.5, 0.25], "ops": "0^ 1"},
                              {"coeff": [-0.5, -0.25], "ops": "1^ 0"}]}, d)
        H, _, _ = load_hamiltonian(p)
        A = H.toarray()
        assert np.allclose(A, A.conj().T), "loaded operator must stay Hermitian"
        assert np.abs(A.imag).max() > 1e-9, "imaginary part was dropped"


def test_scalar_coefficient_is_accepted():
    """_coeff accepts a bare number as well as a [re, im] pair."""
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "qubit_operator",
                    "metadata": _meta(),
                    "terms": [{"coeff": 0.5, "ops": "Z0"}]}, d)
        H, _, _ = load_hamiltonian(p)
        assert np.allclose(np.diag(H.toarray()).real[[0, 1]], [0.5, -0.5])


def test_qubit_operator_pauli_strings():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "qubit_operator",
                    "metadata": _meta(n_qubits=2, n_cells=1),
                    "terms": [{"coeff": [1.0, 0.0], "ops": "X0 X1"}]}, d)
        H, _, _ = load_hamiltonian(p)
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        assert np.allclose(H.toarray(), np.kron(X, X))


# -------------------------------------------------------------- matrix forms
def test_inline_dense_matrix():
    M = np.diag([0.0, 1.0, 2.0, 3.0])
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": _meta(), "data": M.tolist()}, d)
        H, _, _ = load_hamiltonian(p)
        assert np.allclose(H.toarray(), M)


def test_dense_matrix_from_npy_file():
    M = np.diag([0.0, 1.0, 2.0, 3.0]).astype(complex)
    with tempfile.TemporaryDirectory() as d:
        np.save(os.path.join(d, "H.npy"), M)
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": _meta(), "matrix_file": "H.npy"}, d)
        H, _, _ = load_hamiltonian(p)
        assert np.allclose(H.toarray(), M)


def test_sparse_matrix_from_npz_file():
    M = sp.csr_matrix(np.diag([0.0, 1.0, 2.0, 3.0]).astype(complex))
    with tempfile.TemporaryDirectory() as d:
        sp.save_npz(os.path.join(d, "H.npz"), M)
        p = _write({"type": "hamiltonian", "format": "sparse_matrix",
                    "metadata": _meta(), "matrix_file": "H.npz"}, d)
        H, _, _ = load_hamiltonian(p)
        assert np.allclose(H.toarray(), M.toarray())


def test_matrix_formats_are_not_bit_reversed():
    """Only operator formats go through the openfermion conversion; a
    user-supplied matrix must be returned exactly as given."""
    M = np.diag([0.0, 1.0, 2.0, 3.0])
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": _meta(), "data": M.tolist()}, d)
        H, _, _ = load_hamiltonian(p)
        assert np.allclose(np.diag(H.toarray()).real, [0.0, 1.0, 2.0, 3.0])


# ------------------------------------------------------------------ layouts
def test_spin_block_layout_from_metadata():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": {"n_qubits": 8, "n_cells": 2,
                                 "layout": "spin_block", "mode": "up_spin"},
                    "data": np.eye(8).tolist()}, d)
        _, lay, _ = load_hamiltonian(p)
        assert lay["cell_qubits"] == [[0, 1], [2, 3]]


def test_custom_layout_uses_supplied_qubit_lists():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": {"n_qubits": 4, "n_cells": 2, "layout": "custom",
                                 "cell_qubits": [[3, 1], [2, 0]],
                                 "A_qubits": [3, 2], "B_qubits": [1, 0]},
                    "data": np.eye(4).tolist()}, d)
        _, lay, _ = load_hamiltonian(p)
        assert lay["cell_qubits"] == [[3, 1], [2, 0]]
        assert lay["A_qubits"] == [3, 2] and lay["B_qubits"] == [1, 0]


def test_custom_layout_infers_sublattices_when_omitted():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": {"n_qubits": 4, "n_cells": 2, "layout": "custom",
                                 "cell_qubits": [[0, 1], [2, 3]]},
                    "data": np.eye(4).tolist()}, d)
        _, lay, _ = load_hamiltonian(p)
        assert lay["A_qubits"] == [0, 2] and lay["B_qubits"] == [1, 3]


# ------------------------------------------------------------------- errors
def test_rejects_a_non_hamiltonian_file():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "something_else"}, d)
        with pytest.raises(ValueError, match="not a Hamiltonian spec"):
            load_hamiltonian(p)


def test_rejects_an_unknown_format():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "nope",
                    "metadata": _meta(), "terms": []}, d)
        with pytest.raises(ValueError, match="unknown format"):
            load_hamiltonian(p)


def test_rejects_an_unknown_layout():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": {"n_qubits": 4, "n_cells": 2, "layout": "nope"},
                    "data": np.eye(4).tolist()}, d)
        with pytest.raises(ValueError, match="unknown layout"):
            load_hamiltonian(p)


def test_measure_rejects_an_unknown_property():
    with tempfile.TemporaryDirectory() as d:
        p = _write({"type": "hamiltonian", "format": "dense_matrix",
                    "metadata": _meta(), "data": np.eye(4).tolist()}, d)
        with pytest.raises(ValueError, match="unknown property"):
            measure_hamiltonian(p, ["not_a_property"], _BK)


# -------------------------------------------------------------- measure API
def test_measure_returns_the_documented_envelope():
    path = os.path.join("examples", "ssh_spinless_N3_topological.json")
    if not os.path.exists(path):
        pytest.skip("example file not present")
    res = measure_hamiltonian(path, ["berry", "polarization"], _BK)
    assert res["source"] == "hamiltonian_file"
    assert res["file"] == "ssh_spinless_N3_topological.json"
    assert res["backend"] == "qiskit"
    assert set(res["results"]) == {"berry", "polarization"}
    assert res["metadata"]["n_cells"] == 3


def test_measure_normalises_the_ground_state():
    path = os.path.join("examples", "ssh_spinless_N3_trivial.json")
    if not os.path.exists(path):
        pytest.skip("example file not present")
    res = measure_hamiltonian(path, ["polarization"], _BK)
    n = res["results"]["polarization"]
    total = sum(n["n_A"]) + sum(n["n_B"])
    assert total == pytest.approx(3.0, abs=1e-6), total
