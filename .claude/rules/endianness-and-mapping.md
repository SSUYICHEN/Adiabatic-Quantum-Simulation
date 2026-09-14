---
paths:
  - "src/aqs/mapping.py"
  - "src/aqs/hamiltonians.py"
  - "src/aqs/hamiltonian.py"
  - "src/aqs/backends.py"
  - "src/aqs/backends/**"
  - "src/aqs/observables.py"
---

# Endianness & Qubit Mapping (the #1 bug source of this repo)

Bit order caused three separate silent defects. This table is **measured fact**
(`docs/experience.md` §1) — never re-derive it from memory or from library docs.

## The endianness table (LAW)

| Interface | Convention | X on q0, 3 qubits |
|---|---|---|
| qiskit `Statevector` (indices) | little-endian (qubit q → bit q) | index `1` |
| qiskit `Operator` / `np.kron` | little-endian (q0 = **rightmost** kron factor) | `kron(I, X)` |
| CUDA-Q `get_state()` | little-endian — **same as qiskit** | index `1` |
| CUDA-Q `sample()` bitstrings | qubit 0 **leftmost** | `'100'` |
| qiskit `counts` / `to_dict()` strings | qubit 0 **rightmost** | `'001'` |
| OpenFermion `get_sparse_operator` | **big-endian** (mode q → bit n−1−q) | bit `2` |

Consequences that must never be "fixed" away:
- CUDA-Q is internally inconsistent: statevector little-endian, sample strings reversed.
  Aligning statevectors between qiskit and CUDA-Q uses the **identity** qubit map — a
  `q[n-1-i]` reversal there is a bug (it produced the fidelity-0.181 selftest failure).
- OpenFermion is the only big-endian interface, and **real-hopping terms cannot expose it**:
  the mode swap conjugates the off-diagonal elements carrying `Im(t)` while leaving
  `XX+YY` terms invariant. Any derivation involving `Im` must convert first.

## Rules

- Every OpenFermion-derived matrix/operator passes `_to_little_endian` (bit-reversal
  similarity transform; spectrum-preserving) **before** any comparison with a statevector
  or a qiskit operator. No exceptions — not even "it's real so it doesn't matter" (then the
  conversion is a cheap no-op; when it does matter you won't notice otherwise).
- Bit extraction from amplitude indices is little-endian: `(idx >> q) & 1` reads qubit q.
- Qubit layouts are defined in exactly one module (`mapping.py` in the target layout;
  currently split across `core.py`/`observables.py` layout helpers). No other module may
  hardcode the layout arithmetic independently:
  - **Spinful SSHH (paper Eq. 36):** 4N qubits; spin-up block `[0, 2N)`, spin-down block
    `[2N, 4N)`; within a block even index = sublattice A, odd = B.
  - **Spinless (spec §1, `reconstruct_goals.md`):** 2N qubits; `q_2j = A_j`, `q_2j+1 = B_j`.
- Matrix-format Hamiltonian inputs (`dense_matrix`/`sparse_matrix`/inline `data`) are
  taken **as supplied** — already little-endian by contract — and must NOT be bit-reversed.
  Only OpenFermion-built operators (`fermion_operator`/`qubit_operator` formats) convert.
- Cross-block operators (e.g. Hubbard `n↑n↓`) are diagonal and JW-string-free by
  construction — never "add" a JW string or parity factor to them.
- Any change in these files requires the asymmetric-input + complex-hopping test discipline
  of `testing-and-verification.md`.
