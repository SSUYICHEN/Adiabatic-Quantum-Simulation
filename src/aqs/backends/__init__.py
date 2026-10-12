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

Four backends:
  * 'qiskit' : CPU. qiskit Statevector + scipy sparse eigsh.       (default, tested)
  * 'cudaq'  : GPU. NVIDIA CUDA-Q statevector + cupy GPU eigsh.    (Linux+GPU only)
  * 'cupy'   : GPU. custom CUDA kernel statevector (complex128)
               + cupy GPU eigsh.                                   (Linux+GPU only)
  * 'mlx'    : GPU. Apple MLX custom Metal kernel statevector
               (complex64) + scipy eigsh on the CPU.               (Apple silicon only)

The GPU backends are import-guarded: constructing one without its package
installed raises a clear error. Because CUDA-Q is Linux/GPU-only it cannot be
exercised on a Windows/CPU box; run `aqs selftest --backend cudaq` (or `cupy`)
on the GPU machine to verify it (including the qiskit<->cudaq endianness
mapping) against qiskit. Likewise `aqs selftest --backend mlx` on a Mac.

Why both 'cudaq' and 'cupy' on NVIDIA: the cudaq path transpiles every gate to
{rx, ry, rz, cx} (about 14 basis ops per two-qubit gate) and hands CUDA-Q one
kernel with all of them, whose compilation grows superlinearly with length: a
16-qubit annealing circuit of 320 layers (54k basis ops) takes minutes while the
GPU work is well under a second. The cupy backend applies each gate's exact
unitary with one kernel launch (the MLX design in CUDA), so the time is the
memory traffic of the statevector: the same circuit runs in seconds.
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
                "use --backend qiskit, or run inside WSL with the GPU stack; "
                "on an Apple-silicon Mac use --backend mlx."
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

        # CUDA-Q's get_state() indexes amplitudes with qubit j at bit position j,
        # i.e. the same little-endian convention qiskit's Statevector uses.
        # Verified on the nvidia, nvidia-fp64 and qpp-cpu targets: X on qubit
        # 0/1/2 of a 3-qubit kernel yields amplitude index 1/2/4.
        #
        # Do NOT "fix" this into a reversal. CUDA-Q is internally inconsistent:
        # its sample() BITSTRINGS put qubit 0 leftmost (X on qubit 0 -> '100'),
        # the opposite of qiskit's counts ('001'). That string convention is
        # what the previous comment here described, but it does not apply to
        # the statevector. `aqs selftest --backend cudaq` guards this.
        def cq(qiskit_index):
            return q[qiskit_index]

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
        return _cupy_ground_state(H, "the cudaq backend")


def _cupy_ground_state(H, who) -> np.ndarray:
    """Lowest eigenvector of H with cupy on the GPU (sparse eigsh or dense eigh)."""
    try:
        import cupy as cp
        import cupyx.scipy.sparse as csp
        from cupyx.scipy.sparse.linalg import eigsh as cu_eigsh
    except Exception as e:  # pragma: no cover - GPU only
        raise RuntimeError(f"cupy is required for GPU diagonalisation in {who}.") from e
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


def _flatten_gates(circuit, qubit_map=None):
    """([(unitary, [qubit, ...]), ...], global phase) of a circuit.

    Every instruction on <= 2 qubits becomes its exact unitary (qiskit Operator,
    little-endian: bit b of the matrix index <-> qubits[b]); wider instructions
    are expanded through their definitions; non-unitary ones raise.
    """
    from qiskit.quantum_info import Operator
    phase = float(circuit.global_phase)
    gates = []
    for instr in circuit.data:
        op = instr.operation
        qubits = [circuit.find_bit(b).index for b in instr.qubits]
        if qubit_map is not None:
            qubits = [qubit_map[q] for q in qubits]
        if op.name in ("barrier", "id"):
            continue
        if instr.clbits or op.name in ("measure", "reset"):
            raise NotImplementedError(
                f"'{op.name}' is not unitary; this backend only simulates unitary circuits")
        if len(qubits) <= 2:
            gates.append((np.asarray(Operator(op).data), qubits))
        elif op.definition is not None:
            sub, sub_phase = _flatten_gates(op.definition, qubits)
            gates.extend(sub)
            phase += sub_phase
        else:  # pragma: no cover
            raise NotImplementedError(
                f"{len(qubits)}-qubit gate '{op.name}' has no definition to expand")
    return gates, phase


# =====================================================================
# cupy (GPU)  -- Linux + NVIDIA GPU only
# =====================================================================
# Same decomposition as the Metal kernels below, in CUDA and complex128: one
# thread per block of 2^k amplitudes that a k-qubit gate mixes; the block is
# read and written by that thread only, so the update is in place. U is the
# gate's little-endian matrix (bit b of its index <-> qubit qs[b]).
_CUPY_SOURCE = r"""
#include <cupy/complex.cuh>
extern "C" {
__global__ void apply_1q(complex<double>* psi, const complex<double>* U, int q0,
                         unsigned long long n_blocks) {
    unsigned long long t = blockIdx.x * (unsigned long long)blockDim.x + threadIdx.x;
    if (t >= n_blocks) return;
    unsigned long long i = ((t >> q0) << (q0 + 1)) | (t & ((1ULL << q0) - 1ULL));
    unsigned long long idx[2] = {i, i | (1ULL << q0)};
    complex<double> a[2];
    for (int j = 0; j < 2; ++j) a[j] = psi[idx[j]];
    for (int k = 0; k < 2; ++k) {
        complex<double> s(0.0, 0.0);
        for (int j = 0; j < 2; ++j) s += U[2 * k + j] * a[j];
        psi[idx[k]] = s;
    }
}
__global__ void apply_2q(complex<double>* psi, const complex<double>* U, int q0, int q1,
                         unsigned long long n_blocks) {
    unsigned long long t = blockIdx.x * (unsigned long long)blockDim.x + threadIdx.x;
    if (t >= n_blocks) return;
    int lo = q0 < q1 ? q0 : q1, hi = q0 < q1 ? q1 : q0;
    unsigned long long i = ((t >> lo) << (lo + 1)) | (t & ((1ULL << lo) - 1ULL));
    i = ((i >> hi) << (hi + 1)) | (i & ((1ULL << hi) - 1ULL));
    unsigned long long idx[4] = {i, i | (1ULL << q0), i | (1ULL << q1),
                                 i | (1ULL << q0) | (1ULL << q1)};
    complex<double> a[4];
    for (int j = 0; j < 4; ++j) a[j] = psi[idx[j]];
    for (int k = 0; k < 4; ++k) {
        complex<double> s(0.0, 0.0);
        for (int j = 0; j < 4; ++j) s += U[4 * k + j] * a[j];
        psi[idx[k]] = s;
    }
}
}
"""


class CupyBackend(Backend):
    """Statevector on an NVIDIA GPU through a custom CUDA kernel (cupy RawModule).

    Every instruction acting on <= 2 qubits is applied as its exact unitary in
    complex128, so the result agrees with qiskit to round-off (1 - F ~ 1e-13 on
    the annealing circuits) and the cost is one kernel launch and one pass over
    the statevector per gate. ground_state uses cupy's sparse eigsh on the GPU.
    """
    name = "cupy"

    def __init__(self):
        try:
            import cupy as cp
            if cp.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("no CUDA device")
            self._module = cp.RawModule(code=_CUPY_SOURCE)
            self._k1 = self._module.get_function("apply_1q")
            self._k2 = self._module.get_function("apply_2q")
        except Exception as e:
            raise RuntimeError(
                "cupy with a CUDA GPU is not available in this environment. The cupy "
                "backend requires a Linux machine with an NVIDIA GPU and the `cupy` "
                "package (`uv sync --extra cu12` or `cu13`). Elsewhere use --backend "
                "qiskit (CPU) or, on an Apple-silicon Mac, --backend mlx."
            ) from e
        self._cp = cp

    def statevector(self, circuit) -> np.ndarray:
        cp = self._cp
        n = circuit.num_qubits
        dim = 1 << n
        gates, phase = _flatten_gates(circuit)
        psi = cp.zeros(dim, dtype=np.complex128)
        psi[0] = 1.0
        # all gate matrices go up in one transfer; each launch reads its slice
        U1 = [U for U, q in gates if len(q) == 1]
        U2 = [U for U, q in gates if len(q) == 2]
        U1g = cp.asarray(np.array(U1, dtype=np.complex128).reshape(-1, 4)) if U1 else None
        U2g = cp.asarray(np.array(U2, dtype=np.complex128).reshape(-1, 16)) if U2 else None
        i1 = i2 = 0
        block = 256
        for U, qubits in gates:
            if len(qubits) == 1:
                n_blocks = dim >> 1
                grid = (max(1, (n_blocks + block - 1) // block),)
                self._k1(grid, (block,), (psi, U1g[i1], np.int32(qubits[0]), np.uint64(n_blocks)))
                i1 += 1
            else:
                n_blocks = dim >> 2
                grid = (max(1, (n_blocks + block - 1) // block),)
                self._k2(grid, (block,), (psi, U2g[i2], np.int32(qubits[0]), np.int32(qubits[1]),
                                          np.uint64(n_blocks)))
                i2 += 1
        return cp.asnumpy(psi) * np.exp(1j * phase)

    def ground_state(self, H) -> np.ndarray:
        return _cupy_ground_state(H, "the cupy backend")


# =====================================================================
# mlx (GPU)  -- Apple silicon (M-series) only
# =====================================================================
# One thread per block of 2^k amplitudes that a k-qubit gate mixes. The state
# is interleaved float32 (re, im) because Metal has no complex type; U is the
# gate's little-endian matrix flattened the same way, and qs its qubits in
# qiskit order, so amplitude offset bit b of the block <-> qubit qs[b].
_MLX_2Q_SOURCE = """
    uint t = thread_position_in_grid.x;
    uint q0 = qs[0], q1 = qs[1];
    uint lo = min(q0, q1), hi = max(q0, q1);
    uint i = ((t >> lo) << (lo + 1)) | (t & ((1u << lo) - 1u));
    i = ((i >> hi) << (hi + 1)) | (i & ((1u << hi) - 1u));
    uint idx[4] = {i, i | (1u << q0), i | (1u << q1), i | (1u << q0) | (1u << q1)};
    float2 a[4];
    for (int j = 0; j < 4; ++j) a[j] = float2(psi[2 * idx[j]], psi[2 * idx[j] + 1]);
    for (int k = 0; k < 4; ++k) {
        float re = 0.0f, im = 0.0f;
        for (int j = 0; j < 4; ++j) {
            float ur = U[2 * (4 * k + j)], ui = U[2 * (4 * k + j) + 1];
            re += ur * a[j].x - ui * a[j].y;
            im += ur * a[j].y + ui * a[j].x;
        }
        out[2 * idx[k]] = re;
        out[2 * idx[k] + 1] = im;
    }
"""

_MLX_1Q_SOURCE = """
    uint t = thread_position_in_grid.x;
    uint q = qs[0];
    uint i = ((t >> q) << (q + 1)) | (t & ((1u << q) - 1u));
    uint idx[2] = {i, i | (1u << q)};
    float2 a[2];
    for (int j = 0; j < 2; ++j) a[j] = float2(psi[2 * idx[j]], psi[2 * idx[j] + 1]);
    for (int k = 0; k < 2; ++k) {
        float re = 0.0f, im = 0.0f;
        for (int j = 0; j < 2; ++j) {
            float ur = U[2 * (2 * k + j)], ui = U[2 * (2 * k + j) + 1];
            re += ur * a[j].x - ui * a[j].y;
            im += ur * a[j].y + ui * a[j].x;
        }
        out[2 * idx[k]] = re;
        out[2 * idx[k] + 1] = im;
    }
"""


class MlxBackend(Backend):
    """Statevector on the Apple-silicon GPU through a custom MLX Metal kernel.

    Every instruction acting on <= 2 qubits is applied as its exact unitary
    (qiskit Operator, little-endian); wider instructions are expanded through
    their definitions. Amplitudes are complex64 on the GPU -- MLX has no
    float64 on the GPU -- which keeps 1 - F ~ 1e-12 on the annealing circuits
    and is returned as complex128 for the observables.

    ground_state stays on the CPU (scipy): MLX's eigh is CPU-only and float32.
    """
    name = "mlx"

    def __init__(self, sync_every=64):
        try:
            import mlx.core as mx
            if not mx.metal.is_available():
                raise RuntimeError("Metal GPU not available")
        except Exception as e:
            raise RuntimeError(
                "MLX with a Metal GPU is not available in this environment. The "
                "mlx backend requires an Apple-silicon (M-series) Mac and the "
                "`mlx` package (installed by `uv sync` on such a Mac). Elsewhere "
                "use --backend qiskit (CPU) or --backend cudaq (NVIDIA GPU)."
            ) from e
        self._mx = mx
        self._sync_every = int(sync_every)  # bound the lazy graph's size
        self._kernels = {
            1: mx.fast.metal_kernel(name="aqs_apply_1q", input_names=["psi", "U", "qs"],
                                    output_names=["out"], source=_MLX_1Q_SOURCE),
            2: mx.fast.metal_kernel(name="aqs_apply_2q", input_names=["psi", "U", "qs"],
                                    output_names=["out"], source=_MLX_2Q_SOURCE),
        }

    _gates = staticmethod(_flatten_gates)   # kept as the historical name

    def statevector(self, circuit) -> np.ndarray:
        mx = self._mx
        n = circuit.num_qubits
        N = 1 << n
        gates, phase = self._gates(circuit)

        psi0 = np.zeros(2 * N, dtype=np.float32)
        psi0[0] = 1.0
        psi = mx.array(psi0)
        for k, (U, qubits) in enumerate(gates):
            width = len(qubits)
            U_ri = np.stack([U.real, U.imag], -1).astype(np.float32).reshape(-1)
            threads = N >> width
            psi = self._kernels[width](
                inputs=[psi, mx.array(U_ri), mx.array(np.asarray(qubits, dtype=np.uint32))],
                grid=(threads, 1, 1), threadgroup=(min(256, threads), 1, 1),
                output_shapes=[(2 * N,)], output_dtypes=[mx.float32])[0]
            if self._sync_every and k % self._sync_every == self._sync_every - 1:
                mx.eval(psi)
        mx.eval(psi)
        sv = np.array(psi).view(np.complex64).astype(np.complex128)
        # The circuit is unitary, so the only norm change is float32 round-off
        # accumulated over the gates (~1e-5 after ~1e3 gates); remove it so
        # expectation values are not scaled by it.
        sv /= np.linalg.norm(sv)
        return sv * np.exp(1j * phase)

    def ground_state(self, H) -> np.ndarray:
        return QiskitBackend().ground_state(H)


# =====================================================================
# registry + self-test
# =====================================================================
def get_backend(name="qiskit", **kwargs) -> Backend:
    name = (name or "qiskit").lower()
    if name == "qiskit":
        return QiskitBackend()
    if name == "cudaq":
        return CudaqBackend(**kwargs)
    if name == "cupy":
        return CupyBackend(**kwargs)
    if name == "mlx":
        return MlxBackend(**kwargs)
    raise ValueError(f"unknown backend '{name}' (choose 'qiskit', 'cudaq', 'cupy' or 'mlx')")


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
