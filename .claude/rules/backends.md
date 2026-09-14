---
paths:
  - "src/aqs/backends.py"
  - "src/aqs/backends/**"
---

# Backend Contract

## The one contract

Every backend returns statevectors indexed in **qiskit little-endian convention**
(qubit q → bit q of the amplitude index). Endianness normalization happens **inside the
backend** — callers, observables, and experiments never compensate for a backend's native
convention.

## CUDA-Q specifics

- `get_state()` is already little-endian like qiskit: the qubit map between the transpiled
  qiskit circuit and the CUDA-Q kernel is the **identity** (`cq(i) = q[i]`). Do not
  introduce a `q[n−1−i]` reversal — the observation "CUDA-Q bitstrings look reversed"
  applies to `sample()` strings only, not to statevector indices (see
  `endianness-and-mapping.md`; the reversal bug cost fidelity 0.181).
- Execution path: transpile to the `{rx, ry, rz, cx}` basis, then rebuild as a
  `cudaq.make_kernel()`. Unsupported gates must raise, not silently skip.
- `cudaq` import stays lazy/optional: a CPU-only environment (plain `uv sync`, no extras)
  must import `aqs` and run the qiskit backend cleanly; GPU tests skip cleanly.
- Default target `nvidia` is fp32 — relevant when comparing against fp64 qiskit
  statevectors; fall back gracefully when the target is unavailable.

## Selftest = endianness gatekeeper

`aqs selftest` compares backends on an **asymmetric** (non-uniform-density,
non-reflection-symmetric) state. A uniform or symmetric selftest state is structurally
blind to bit-order bugs and would be a fake guard. Keep `selftest --backend cudaq` as the
cross-backend gate for any backend or transpilation change.

## Adding a backend

Requires all of: registry entry via `get_backend`; inclusion in selftest; a docstring
stating both of its native conventions ("index convention of `get_state`, string
convention of `sample`") with the normalization applied; skip-clean behavior when its
runtime is absent.
