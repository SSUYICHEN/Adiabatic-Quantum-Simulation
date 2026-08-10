"""Guard tests for the qiskit <-> CUDA-Q statevector convention.

CUDA-Q is internally inconsistent about qubit ordering, which is what made the
original defect so easy to introduce:

    get_state()  amplitude index   -> qubit q at bit q      (little-endian)
    sample()     bitstring         -> qubit 0 LEFTMOST      (opposite of qiskit)

qiskit's counts put qubit 0 rightmost ('001' for X on qubit 0) while CUDA-Q's
put it leftmost ('100'). Reasoning from the *string* convention to the *index*
convention is exactly the mistake that produced the spurious q[n-1-i] mapping.

GPU tests skip cleanly when CUDA-Q is unavailable.

Runs under pytest, or standalone:  python tests/test_backend_endianness.py
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

from aqs.backends import get_backend, selftest


def _cudaq_available():
    try:
        import cudaq  # noqa: F401
        get_backend("cudaq")
        return True
    except Exception:
        return False


HAVE_CUDAQ = _cudaq_available()


def _skip_without_gpu():
    if not HAVE_CUDAQ:
        try:
            import pytest
            pytest.skip("CUDA-Q unavailable", allow_module_level=False)
        except ImportError:
            pass
        return True
    return False


def _fidelity(a, b):
    a = np.asarray(a).ravel(); b = np.asarray(b).ravel()
    return float(np.abs(np.vdot(a, b)) ** 2 /
                 (np.vdot(a, a) * np.vdot(b, b)).real)


# ------------------------------------------------ qiskit reference convention
def test_qiskit_statevector_is_little_endian():
    """X on qubit q must put the amplitude at index 2**q."""
    n = 3
    for q in range(n):
        qc = QuantumCircuit(n); qc.x(q)
        assert int(np.argmax(np.abs(Statevector(qc).data))) == (1 << q)


# ------------------------------------------------------ CUDA-Q raw convention
def test_cudaq_get_state_is_little_endian_like_qiskit():
    """The convention the translation in CudaqBackend.statevector relies on."""
    if _skip_without_gpu():
        return
    import cudaq
    n = 3
    for target in ("nvidia", "qpp-cpu"):
        cudaq.set_target(target)
        for q in range(n):
            k = cudaq.make_kernel(); qb = k.qalloc(n); k.x(qb[q])
            sv = np.array(cudaq.get_state(k), copy=True).ravel()
            got = int(np.argmax(np.abs(sv)))
            assert got == (1 << q), (
                f"target={target} qubit={q}: index {got}, expected {1 << q}. "
                "If CUDA-Q ever changes this, CudaqBackend.statevector's "
                "identity qubit map must change with it.")


def test_cudaq_bitstrings_are_opposite_to_qiskit_counts():
    """Documents the trap: strings and indices disagree. Not a bug -- a warning.

    If this test ever fails, CUDA-Q changed its bitstring convention, and the
    comment in CudaqBackend.statevector should be revisited.
    """
    if _skip_without_gpu():
        return
    import cudaq
    cudaq.set_target("qpp-cpu")
    n = 3
    k = cudaq.make_kernel(); qb = k.qalloc(n); k.x(qb[0]); k.mz(qb)
    key = list(cudaq.sample(k, shots_count=8))[0]
    assert key == "100", f"CUDA-Q bitstring for X@q0 was {key!r}, expected '100'"


# ------------------------------------------------------- backend agreement
def test_cudaq_backend_matches_qiskit_on_the_selftest_circuit():
    if _skip_without_gpu():
        return
    ok, fid = selftest("cudaq")
    assert ok, f"cudaq selftest fidelity {fid:.8f} (pre-fix value was 0.18124364)"


def test_cudaq_backend_matches_qiskit_on_asymmetric_circuits():
    """Randomised circuits deliberately NOT symmetric under qubit reversal.

    A symmetric circuit cannot distinguish the two conventions, which is how the
    original defect slipped through casual checks.
    """
    if _skip_without_gpu():
        return
    bk = get_backend("cudaq")
    rng = np.random.default_rng(20240810)
    worst = 1.0
    for _ in range(12):
        n = int(rng.integers(2, 6))
        qc = QuantumCircuit(n)
        for _ in range(12):
            g = int(rng.integers(0, 4))
            if g == 3 and n >= 2:
                a, b = rng.choice(n, size=2, replace=False)
                qc.cx(int(a), int(b))
            else:
                t = int(rng.integers(0, n))
                [qc.rx, qc.ry, qc.rz][g](float(rng.uniform(0, 2 * np.pi)), t)
        worst = min(worst, _fidelity(np.asarray(Statevector(qc).data),
                                     bk.statevector(qc)))
    assert worst > 1 - 1e-6, f"worst fidelity {worst:.10f}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = skipped = 0
    for fn in fns:
        gpu_only = "cudaq" in fn.__name__
        if gpu_only and not HAVE_CUDAQ:
            skipped += 1
            print(f"  SKIP  {fn.__name__} (CUDA-Q unavailable)")
            continue
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except Exception as e:
            failed += 1
            print(f"  FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed - skipped}/{len(fns) - skipped} passed"
          + (f", {skipped} skipped" if skipped else ""))
    raise SystemExit(1 if failed else 0)
