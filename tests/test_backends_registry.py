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


# ------------------------------------------- Apple-silicon GPU paths (mlx)
def _have_mlx():
    try:
        get_backend("mlx")
        return True
    except Exception:
        return False


HAVE_MLX = _have_mlx()
metal = pytest.mark.skipif(not HAVE_MLX, reason="MLX / Metal GPU unavailable")


def _fidelity(a, b):
    return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a) * np.vdot(b, b)).real)


def test_mlx_import_guard_is_informative():
    """Without MLX/Metal the error must name the package, not ImportError."""
    if HAVE_MLX:
        pytest.skip("MLX with a Metal GPU is available; the guard path is unreachable here")
    with pytest.raises(RuntimeError, match="MLX"):
        get_backend("mlx")


@metal
def test_mlx_selftest_passes():
    ok, fid = selftest("mlx")
    assert ok and fid == pytest.approx(1.0, abs=1e-6)


@metal
def test_mlx_matches_qiskit_with_global_phase_and_a_wide_gate():
    """Global phase is restored and >2-qubit gates expand through their definition."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(4, global_phase=0.7)
    qc.h(0); qc.ry(0.3, 3); qc.ccx(0, 3, 1); qc.cp(0.9, 1, 2); qc.rxx(0.4, 2, 0)
    ref = np.asarray(Statevector(qc).data)
    got = get_backend("mlx").statevector(qc)
    assert np.allclose(got, ref, atol=1e-6)        # phase-sensitive, not just |<.|.>|


@metal
def test_mlx_single_qubit_circuit():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(1); qc.h(0); qc.s(0)
    assert np.allclose(get_backend("mlx").statevector(qc), Statevector(qc).data, atol=1e-6)


@metal
def test_mlx_rejects_measurement():
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(1, 1); qc.h(0); qc.measure(0, 0)
    with pytest.raises(NotImplementedError, match="not unitary"):
        get_backend("mlx").statevector(qc)


@metal
def test_mlx_hubbard_evolution_matches_qiskit():
    """The full annealing circuit (Givens prep + R/CP Trotter layers) in complex64."""
    from aqs.models.hubbard import HubbardModel, HubbardSim
    sim = HubbardSim(HubbardModel(n_sites=4, t=1.0, U=4.0, PBC=False), T=2.0, L=20)
    ref = sim.evolve()
    got = sim.evolve(backend="mlx")
    assert np.linalg.norm(got) == pytest.approx(1.0, abs=1e-12)
    assert _fidelity(ref, got) == pytest.approx(1.0, abs=1e-9)
    assert sim.energy(got) == pytest.approx(sim.energy(ref), abs=1e-5)


@metal
def test_mlx_ground_state_uses_the_cpu_solver():
    M = sp.csr_matrix(np.diag([3.0, -2.0, 5.0, 1.0]).astype(complex))
    v = np.asarray(get_backend("mlx").ground_state(M)).ravel()
    v = v / np.linalg.norm(v)
    assert float(np.real(np.vdot(v, M.toarray() @ v))) == pytest.approx(-2.0, abs=1e-8)


# ------------------------------------------------ NVIDIA GPU statevector (cupy)
def _have_cupy_gpu():
    try:
        get_backend("cupy")
        return True
    except Exception:
        return False


HAVE_CUPY = _have_cupy_gpu()
cuda = pytest.mark.skipif(not HAVE_CUPY, reason="cupy / CUDA GPU unavailable")


def test_cupy_import_guard_is_informative():
    """Without cupy (or a GPU) the error must name the package, not ImportError."""
    if HAVE_CUPY:
        pytest.skip("cupy with a CUDA GPU is available; the guard path is unreachable here")
    with pytest.raises(RuntimeError, match="cupy"):
        get_backend("cupy")


@cuda
def test_cupy_selftest_passes():
    ok, fid = selftest("cupy")
    assert ok and fid == pytest.approx(1.0, abs=1e-12)


@cuda
def test_cupy_matches_qiskit_with_global_phase_and_a_wide_gate():
    """Element-wise (phase-sensitive) agreement in double precision on an
    asymmetric circuit; >2-qubit gates expand through their definition.
    Killed by: swapping the two qubit indices of a 2-qubit gate in the kernel
    (0.18-type endianness bug) or dropping the global phase."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(4, global_phase=0.7)
    qc.h(0); qc.ry(0.3, 3); qc.ccx(0, 3, 1); qc.cp(0.9, 1, 2); qc.rxx(0.4, 2, 0)
    ref = np.asarray(Statevector(qc).data)
    got = get_backend("cupy").statevector(qc)
    assert np.allclose(got, ref, atol=1e-12)


@cuda
def test_cupy_single_qubit_circuit():
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(1); qc.h(0); qc.s(0)
    assert np.allclose(get_backend("cupy").statevector(qc), Statevector(qc).data, atol=1e-12)


@cuda
def test_cupy_rejects_measurement():
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(1, 1); qc.h(0); qc.measure(0, 0)
    with pytest.raises(NotImplementedError, match="not unitary"):
        get_backend("cupy").statevector(qc)


@cuda
def test_cupy_spin_ssh_annealing_matches_qiskit_in_double_precision():
    """The full annealing circuit (X + Givens prep, R/G/CP layers) with complex
    hoppings, so the G gates and every 2-qubit phase are exercised."""
    from aqs.models.spin_ssh import SpinSSHModel, SpinSSHHSim
    sim = SpinSSHHSim(SpinSSHModel(4, 0.5 + 0.2j, 1.5 - 0.3j, U=2.0, PBC=True), T=1.0, L=20)
    qc = sim.full_circuit()
    ref = get_backend("qiskit").statevector(qc)
    got = get_backend("cupy").statevector(qc)
    assert np.linalg.norm(got) == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(got, ref, atol=1e-11)


@cuda
def test_cupy_ground_state():
    M = sp.csr_matrix(np.diag([3.0, -2.0, 5.0, 1.0]).astype(complex))
    v = np.asarray(get_backend("cupy").ground_state(M)).ravel()
    v = v / np.linalg.norm(v)
    assert float(np.real(np.vdot(v, M.toarray() @ v))) == pytest.approx(-2.0, abs=1e-8)
