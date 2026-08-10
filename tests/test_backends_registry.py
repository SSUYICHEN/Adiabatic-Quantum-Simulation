"""Tests for the backend registry, the base class contract and ground_state."""
import numpy as np
import pytest
import scipy.sparse as sp

from aqs.backends import Backend, QiskitBackend, get_backend, selftest


def test_registry_returns_the_right_class_and_defaults_to_qiskit():
    assert isinstance(get_backend("qiskit"), QiskitBackend)
    assert isinstance(get_backend(), QiskitBackend)
    assert isinstance(get_backend(None), QiskitBackend)
    assert isinstance(get_backend("QISKIT"), QiskitBackend)      # case-insensitive


def test_registry_rejects_an_unknown_backend():
    with pytest.raises(ValueError, match="unknown backend"):
        get_backend("does_not_exist")


def test_base_class_methods_are_abstract():
    b = Backend()
    with pytest.raises(NotImplementedError):
        b.statevector(None)
    with pytest.raises(NotImplementedError):
        b.ground_state(None)


def test_selftest_passes_on_qiskit():
    ok, fid = selftest("qiskit")
    assert ok and fid == pytest.approx(1.0, abs=1e-12)


# ------------------------------------------------------------- ground_state
def _check_ground_state(H, want_energy, dense=False):
    bk = QiskitBackend()
    v = np.asarray(bk.ground_state(H)).ravel()
    v = v / np.linalg.norm(v)
    A = H if dense else H.toarray()
    e = float(np.real(np.vdot(v, A @ v)))
    assert e == pytest.approx(want_energy, abs=1e-8), (e, want_energy)


def test_ground_state_sparse_path():
    M = sp.csr_matrix(np.diag([3.0, -2.0, 5.0, 1.0]).astype(complex))
    _check_ground_state(M, -2.0)


def test_ground_state_dense_path():
    M = np.diag([3.0, -2.0, 5.0, 1.0]).astype(complex)
    _check_ground_state(M, -2.0, dense=True)


def test_ground_state_tiny_sparse_falls_back_to_dense():
    """eigsh needs k < n-1, so the 2x2 case must take the dense branch."""
    M = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, -4.0]], dtype=complex))
    _check_ground_state(M, -4.0)


def test_ground_state_matches_numpy_on_a_random_hermitian():
    rng = np.random.default_rng(11)
    d = 16
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = A + A.conj().T
    want = float(np.linalg.eigvalsh(H)[0])
    _check_ground_state(sp.csr_matrix(H), want)


def test_statevector_matches_qiskit_reference():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(3)
    qc.h(0); qc.cx(0, 1); qc.ry(0.4, 2)
    got = QiskitBackend().statevector(qc)
    assert np.allclose(got, np.asarray(Statevector(qc).data))


def test_cudaq_import_guard_is_informative():
    """Without CUDA-Q installed the error must name the package, not TypeError."""
    try:
        import cudaq  # noqa: F401
        pytest.skip("CUDA-Q is installed; the guard path is unreachable here")
    except ImportError:
        pass
    with pytest.raises(RuntimeError, match="CUDA-Q"):
        get_backend("cudaq")


# ------------------------------------------------- GPU-only paths (cudaq)
def _have_cudaq():
    try:
        import cudaq  # noqa: F401
        get_backend("cudaq")
        return True
    except Exception:
        return False


HAVE_CUDAQ = _have_cudaq()
gpu = pytest.mark.skipif(not HAVE_CUDAQ, reason="CUDA-Q unavailable")


@gpu
def test_cudaq_ground_state_sparse_matches_qiskit():
    """CudaqBackend.ground_state routes through cupy; it must agree with scipy."""
    rng = np.random.default_rng(5)
    d = 32
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    H = sp.csr_matrix(A + A.conj().T)
    want = float(np.linalg.eigvalsh(H.toarray())[0])
    v = np.asarray(get_backend("cudaq").ground_state(H)).ravel()
    v = v / np.linalg.norm(v)
    got = float(np.real(np.vdot(v, H.toarray() @ v)))
    assert got == pytest.approx(want, abs=1e-6), (got, want)


@gpu
def test_cudaq_ground_state_dense_path():
    M = np.diag([3.0, -2.0, 5.0, 1.0]).astype(complex)
    v = np.asarray(get_backend("cudaq").ground_state(M)).ravel()
    v = v / np.linalg.norm(v)
    assert float(np.real(np.vdot(v, M @ v))) == pytest.approx(-2.0, abs=1e-8)


@gpu
def test_cudaq_ground_state_tiny_sparse_falls_back_to_dense():
    M = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, -4.0]], dtype=complex))
    v = np.asarray(get_backend("cudaq").ground_state(M)).ravel()
    v = v / np.linalg.norm(v)
    assert float(np.real(np.vdot(v, M.toarray() @ v))) == pytest.approx(-4.0, abs=1e-8)


@gpu
def test_cudaq_rejects_a_gate_it_cannot_translate():
    """The translator handles rx/ry/rz/cx only; anything else must say so loudly."""
    from qiskit import QuantumCircuit
    import qiskit

    qc = QuantumCircuit(2)
    qc.h(0)
    bk = get_backend("cudaq")
    # bypass the transpile step so an untranslatable gate reaches the dispatcher
    orig = qiskit.transpile
    try:
        qiskit.transpile = lambda circuit, **kw: circuit
        with pytest.raises(NotImplementedError, match="not handled by cudaq"):
            bk.statevector(qc)
    finally:
        qiskit.transpile = orig


@gpu
def test_cudaq_falls_back_when_the_target_is_unavailable():
    """An unknown target must degrade to qpp-cpu rather than raising."""
    from aqs.backends import CudaqBackend
    bk = CudaqBackend(target="definitely-not-a-target")
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1)
    got = bk.statevector(qc)
    ref = np.asarray(Statevector(qc).data)
    f = abs(np.vdot(ref, got)) ** 2 / (np.vdot(ref, ref) * np.vdot(got, got)).real
    assert float(f) == pytest.approx(1.0, abs=1e-6)
