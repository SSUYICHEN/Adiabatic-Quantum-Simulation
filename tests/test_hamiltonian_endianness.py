"""Guard tests for the OpenFermion -> qiskit qubit-ordering conversion.

openfermion.get_sparse_operator places qubit q at bit position n-1-q
(big-endian), for both FermionOperator and QubitOperator. qiskit's Statevector
-- and therefore every function in observables.PROPERTIES -- places it at bit q.
load_hamiltonian must reconcile the two, or every observable read off a
fermion_operator / qubit_operator file comes back bit-reversed.

Crucially, the Hamiltonians shipped in examples/ CANNOT detect this: their
density profiles are uniform (0.5 on every site), so bit reversal is the
identity on them. Every test below therefore uses a deliberately asymmetric
Hamiltonian.

Runs under pytest, or standalone:  python tests/test_hamiltonian_endianness.py
"""
import json
import os
import tempfile

import numpy as np
import openfermion as of

from aqs.hamiltonian import (_bit_reversal_permutation, _to_little_endian,
                             load_hamiltonian, measure_hamiltonian)
from aqs.backends import get_backend
from aqs.observables import density_profile_qubits

_BK = get_backend("qiskit")


# ------------------------------------------------------------ the permutation
def test_bit_reversal_is_an_involution():
    for n in (1, 2, 3, 4, 6, 8):
        p = _bit_reversal_permutation(n)
        assert np.array_equal(p[p], np.arange(1 << n)), n


def test_bit_reversal_maps_single_bits_as_expected():
    n = 4
    p = _bit_reversal_permutation(n)
    for q in range(n):
        assert p[1 << q] == 1 << (n - 1 - q), q


def test_similarity_transform_preserves_the_spectrum():
    """P H P^T is a permutation similarity, so eigenvalues must be untouched."""
    rng = np.random.default_rng(7)
    for n in (2, 3, 4):
        d = 1 << n
        A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        H = A + A.conj().T
        H_l = _to_little_endian(H, n)
        assert np.allclose(np.linalg.eigvalsh(H), np.linalg.eigvalsh(H_l), atol=1e-10)


# --------------------------------------------------- openfermion's convention
def test_openfermion_is_big_endian():
    """Documents the upstream convention this conversion exists to absorb.

    If openfermion ever changes this, _to_little_endian must be removed rather
    than silently double-reversing.
    """
    n = 3
    for m in range(n):
        H = of.get_sparse_operator(
            of.FermionOperator(f"{m}^ {m}", 1.0), n_qubits=n).toarray()
        occ = np.nonzero(np.real(np.diag(H)) > 0.5)[0]
        bit = [p for p in range(n) if all((int(i) >> p) & 1 for i in occ)]
        assert bit == [n - 1 - m], f"mode {m} landed at bit {bit}, expected {n-1-m}"


# ----------------------------------------- end to end through load_hamiltonian
def _trapped_chain_spec(n, trap_mode, path):
    """Chain with a deep well on one mode -> ground state localised there."""
    terms = [{"coeff": [-8.0, 0.0], "ops": f"{trap_mode}^ {trap_mode}"}]
    for i in range(n - 1):
        terms += [{"coeff": [-0.2, 0.0], "ops": f"{i}^ {i+1}"},
                  {"coeff": [-0.2, 0.0], "ops": f"{i+1}^ {i}"}]
    spec = {
        "type": "hamiltonian", "format": "fermion_operator",
        "metadata": {"n_qubits": n, "n_cells": n // 2, "layout": "spinless"},
        "terms": terms,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(spec, f)
    return path


def test_loaded_operator_is_indexed_the_way_properties_read_it():
    """The decisive check: a particle trapped on mode m must read out at qubit m."""
    n = 4
    with tempfile.TemporaryDirectory() as d:
        for trap in range(n):
            path = _trapped_chain_spec(n, trap, os.path.join(d, "h.json"))
            H, layout, meta = load_hamiltonian(path)
            psi = np.asarray(_BK.ground_state(H)).ravel()
            psi = psi / np.linalg.norm(psi)
            dens = density_profile_qubits(psi, meta["n_qubits"])
            got = int(np.argmax(dens))
            assert got == trap, (
                f"trapped mode {trap} read out at qubit {got}; "
                f"density={np.round(dens, 4)}")


def test_measure_hamiltonian_polarization_has_the_right_sign():
    """Asymmetric case where a bit reversal would flip the reported profile."""
    n = 4
    with tempfile.TemporaryDirectory() as d:
        path = _trapped_chain_spec(n, 0, os.path.join(d, "h.json"))
        res = measure_hamiltonian(path, ["polarization"], _BK)
        prof = res["results"]["polarization"]["n_A_minus_n_B"]
        # trap on mode 0 = A_0, so cell 0 must be strongly A-polarised
        assert prof[0] > 0.5, f"cell 0 polarization {prof[0]:+.4f}, expected > 0.5"


def test_shipped_examples_are_unchanged_by_the_conversion():
    """These are uniform-density, so the fix must be a no-op on them.

    Recorded so that a future change which *does* move these numbers is caught
    and questioned rather than accepted.
    """
    for name, gamma, z in [("ssh_spinless_N3_topological", 0.0, 0.030338),
                           ("ssh_spinless_N3_trivial", 0.0, 0.958801)]:
        path = os.path.join("examples", f"{name}.json")
        if not os.path.exists(path):
            continue
        res = measure_hamiltonian(path, ["berry", "polarization"], _BK)
        b = res["results"]["berry"]
        assert abs(b["Twist_Amplitude"] - z) < 1e-5, (name, b["Twist_Amplitude"])
        assert abs(b["Berry_Phase_pi_wrapped"] - gamma) < 1e-6, name
        assert all(abs(x) < 1e-9
                   for x in res["results"]["polarization"]["n_A_minus_n_B"]), name


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except Exception as e:
            failed += 1
            print(f"  FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
