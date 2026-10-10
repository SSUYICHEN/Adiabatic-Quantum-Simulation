"""mapping.py - The one place where qubit bit order is decided.

Measured conventions (.claude/rules/endianness-and-mapping.md):

    qiskit Statevector / this package : qubit q = bit q of the index   (little-endian)
    OpenFermion get_sparse_operator    : qubit q = bit (n-1-q)          (big-endian)
    CUDA-Q get_state()                 : little-endian like qiskit (sample() strings differ)

Everything OpenFermion-built passes through `_to_little_endian` before it is
compared with a circuit state. Nothing else in the package re-derives bit order.
"""

from __future__ import annotations

import numpy as np


def _bit_reversal_permutation(n_qubits):
    """Index map between OpenFermion's and qiskit's qubit orderings.

    OpenFermion's get_sparse_operator places qubit q at bit position
    (n_qubits - 1 - q) -- big-endian -- for both FermionOperator and
    QubitOperator. qiskit's Statevector, and therefore every function in
    observables.PROPERTIES, places qubit q at bit position q. Reversing the bits
    of a basis index converts between the two.

    Memory: this materialises 2**n_qubits int64 entries (128 MiB at n = 24).
    That is well below the sparse operator it permutes, but it does scale
    exponentially -- if the loader ever needs to go much beyond ~24 qubits,
    permute the eigenvector instead of the matrix.
    """
    idx = np.arange(1 << n_qubits, dtype=np.int64)
    perm = np.zeros_like(idx)
    for p in range(n_qubits):
        perm |= ((idx >> p) & 1) << (n_qubits - 1 - p)
    return perm


def _to_little_endian(H, n_qubits):
    """Relabel a big-endian operator into qiskit's basis: H -> P H P^T.

    P is a permutation matrix, so this is a unitary similarity transform: the
    spectrum is preserved exactly and only the basis labelling changes.
    """
    perm = _bit_reversal_permutation(n_qubits)
    return H[perm][:, perm]
