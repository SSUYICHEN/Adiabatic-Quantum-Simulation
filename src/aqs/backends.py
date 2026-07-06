"""
backends.py - Simulation backends.

A backend provides two capabilities:

  * statevector(circuit)  -> np.ndarray
        Final statevector amplitudes of a qiskit QuantumCircuit, in qiskit's
        little-endian convention (qubit k = bit k of the index). The observables
        rely on this convention.

  * ground_state(H)       -> np.ndarray
        Lowest-eigenvalue eigenvector of a (sparse) Hermitian Hamiltonian, used
        for the arbitrary-Hamiltonian path (exact diagonalisation).

Two backends:
  * 'qiskit' : CPU. qiskit Statevector + scipy sparse eigsh.       (default, tested)
  * 'cudaq'  : GPU. NVIDIA CUDA-Q statevector + cupy GPU eigsh.    (Linux+GPU only)

The cudaq backend is import-guarded: importing it without CUDA-Q installed
raises a clear error. Because CUDA-Q is Linux/GPU-only it cannot be exercised on
a Windows/CPU box; run `aqs selftest --backend cudaq` on the GPU machine to
verify it (including the qiskit<->cudaq endianness mapping) against qiskit.
"""

from __future__ import annotations

import numpy as np


# =====================================================================
# Base
# =====================================================================
class Backend:
    name = "base"

    def statevector(self, circuit) -> np.ndarray:
        raise NotImplementedError

    def ground_state(self, H) -> np.ndarray:
        raise NotImplementedError


# =====================================================================
# qiskit (CPU)
# =====================================================================
class QiskitBackend(Backend):
    name = "qiskit"

    def statevector(self, circuit) -> np.ndarray:
        from qiskit.quantum_info import Statevector
        return np.asarray(Statevector(circuit).data)

    def ground_state(self, H) -> np.ndarray:
        import scipy.sparse as sp
        from scipy.sparse.linalg import eigsh
        if sp.issparse(H):
            if H.shape[0] <= 2:  # eigsh needs k < n-1; tiny case -> dense
                vals, vecs = np.linalg.eigh(H.toarray())
                return vecs[:, int(np.argmin(vals))]
            _, vecs = eigsh(H.astype(complex), k=1, which="SA")
            return vecs[:, 0]
        vals, vecs = np.linalg.eigh(np.asarray(H))
        return vecs[:, int(np.argmin(vals))]


# =====================================================================
# cudaq (GPU)  -- Linux + NVIDIA GPU only
# =====================================================================
class CudaqBackend(Backend):
    name = "cudaq"

    def __init__(self, target="nvidia"):
        try:
            import cudaq  # noqa: F401
        except Exception as e:  # pragma: no cover - GPU only
            raise RuntimeError(
                "CUDA-Q (cudaq) is not available in this environment. The cudaq "
                "backend requires a Linux machine with an NVIDIA GPU and the "
                "`cudaq`/`cupy` packages (see requirements-gpu.txt). On Windows, "
                "use --backend qiskit, or run inside WSL with the GPU stack."
            ) from e
        import cudaq
        self._cudaq = cudaq
        try:
            cudaq.set_target(target)
        except Exception:
            cudaq.set_target("qpp-cpu")  # fall back to CUDA-Q's CPU simulator

    # -- circuit statevector via a CUDA-Q kernel -----------------------
    def statevector(self, circuit) -> np.ndarray:
        cudaq = self._cudaq
        from qiskit import transpile
        # Decompose to a minimal gate set CUDA-Q's builder supports directly.
        tqc = transpile(circuit, basis_gates=["rx", "ry", "rz", "cx"],
                        optimization_level=0)
        n = tqc.num_qubits

        kernel = cudaq.make_kernel()
        q = kernel.qalloc(n)

        # qiskit is little-endian (qubit 0 = LSB); CUDA-Q orders qubit 0 as the
        # most-significant bit. Map qiskit qubit i -> CUDA-Q qubit (n-1-i) so the
        # returned amplitude index matches qiskit's convention bit-for-bit.
        def cq(qiskit_index):
            return q[n - 1 - qiskit_index]

        for instr in tqc.data:
            name = instr.operation.name
            qubits = [tqc.find_bit(b).index for b in instr.qubits]
            params = [float(p) for p in instr.operation.params]
            if name == "rx":
                kernel.rx(params[0], cq(qubits[0]))
            elif name == "ry":
                kernel.ry(params[0], cq(qubits[0]))
            elif name == "rz":
                kernel.rz(params[0], cq(qubits[0]))
            elif name == "cx":
                kernel.cx(cq(qubits[0]), cq(qubits[1]))
            elif name in ("barrier", "id"):
                continue
            else:  # pragma: no cover
                raise NotImplementedError(
                    f"gate '{name}' not handled by cudaq translation; "
                    "extend backends.CudaqBackend.statevector")
        state = cudaq.get_state(kernel)
        return np.array(state, copy=True).astype(complex).ravel()

    # -- ground state via cupy GPU diagonalisation ---------------------
    def ground_state(self, H) -> np.ndarray:
        try:
            import cupy as cp
            import cupyx.scipy.sparse as csp
            from cupyx.scipy.sparse.linalg import eigsh as cu_eigsh
        except Exception as e:  # pragma: no cover - GPU only
            raise RuntimeError(
                "cupy is required for GPU diagonalisation in the cudaq backend."
            ) from e
        import scipy.sparse as sp
        if sp.issparse(H):
            Hg = csp.csr_matrix(H.astype(complex))
            if Hg.shape[0] <= 2:
                vals, vecs = cp.linalg.eigh(Hg.toarray())
                return cp.asnumpy(vecs[:, int(cp.argmin(vals))])
            _, vecs = cu_eigsh(Hg, k=1, which="SA")
            return cp.asnumpy(vecs[:, 0])
        Hg = cp.asarray(np.asarray(H))
        vals, vecs = cp.linalg.eigh(Hg)
        return cp.asnumpy(vecs[:, int(cp.argmin(vals))])


# =====================================================================
# registry + self-test
# =====================================================================
def get_backend(name="qiskit", **kwargs) -> Backend:
    name = (name or "qiskit").lower()
    if name == "qiskit":
        return QiskitBackend()
    if name == "cudaq":
        return CudaqBackend(**kwargs)
    raise ValueError(f"unknown backend '{name}' (choose 'qiskit' or 'cudaq')")


def selftest(backend_name="qiskit"):
    """Build a small entangling circuit and compare the backend's statevector to
    qiskit's reference (up to global phase). Returns (ok, fidelity)."""
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Statevector
    qc = QuantumCircuit(4)
    qc.h(0); qc.cx(0, 1); qc.ry(0.7, 2); qc.cx(2, 3)
    qc.rz(0.3, 1); qc.cx(1, 2); qc.rx(1.1, 0)
    ref = np.asarray(Statevector(qc).data)
    got = get_backend(backend_name).statevector(qc)
    fid = float(np.abs(np.vdot(ref, got)) ** 2 / (np.vdot(ref, ref) * np.vdot(got, got)).real)
    return (fid > 1 - 1e-6), fid
