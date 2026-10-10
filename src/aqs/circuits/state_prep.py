"""state_prep.py - qiskit side of OpenFermion's Slater-determinant preparation.

This is the only module that knows how OpenFermion describes a state-
preparation circuit. OpenFermion supplies the algorithm
(``slater_determinant_preparation_circuit``: an orthonormal n x M matrix Q of
occupied orbitals -> layers of Givens rotations ``(i, j, theta, phi)`` on
neighbouring modes); this module turns that description into qiskit gates built
from the repository's own primitives, and nothing else in ``aqs`` touches the
description format.

Conventions taken from OpenFermion itself
(``openfermion.circuits.primitives.state_preparation``):

* the circuit acts on the reference state with the first n modes occupied;
* each tuple ``(i, j, theta, phi)`` is  ``Ryxxy(theta)`` on ``(i, j)``, then
  ``Z(j) ** (phi / pi)``,  with  ``Ryxxy(t) = exp(-i t (Y_i X_j - X_i Y_j)/2)``.

``Ryxxy(theta)`` is this repository's ``create_G_gate(-theta)`` and
``Z ** (phi/pi)`` is qiskit's phase gate ``p(phi)``; see ``givens_gate``.
Guarded by tests/test_state_preparation.py against OpenFermion's own cirq
realisation and against determinants built with ``jw_slater_determinant``.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import openfermion
from qiskit import QuantumCircuit
from qiskit.circuit import Instruction

from aqs.circuits.gates import create_G_gate


def givens_gate(theta: float, phi: float) -> Instruction:
    """qiskit gate for one OpenFermion Givens tuple ``(i, j, theta, phi)``.

    Equals OpenFermion's ``Ryxxy(theta)`` followed by ``Z(j) ** (phi/pi)``,
    i.e. in the (|00>, |01>, |10>, |11>) occupation basis of qubits (i, j)

        diag-block { 1 ;  [[cos t, sin t], [-e^{i p} sin t, e^{i p} cos t]] ;  e^{i p} } .

    Operand order matters: the gate is placed on ``[i, j]`` as OpenFermion
    lists them (``i < j``).
    """
    qc = QuantumCircuit(2, name="Givens")
    qc.append(create_G_gate(-theta), [0, 1])  # Ryxxy(theta) == G(-theta)
    if phi != 0.0:
        qc.p(phi, 1)  # det G = e^{i phi} on |11>
    return qc.to_instruction()


def prepare_slater_determinant(
    qc: QuantumCircuit, Q: np.ndarray, qubits: Sequence[int]
) -> QuantumCircuit:
    """Append the circuit preparing the Slater determinant of ``Q`` to ``qc``.

    ``Q`` is the n x M matrix whose rows are the occupied orbitals in the mode
    basis; ``qubits`` lists the M qubits that carry those modes, in mode order.
    X gates on the first n of them set the reference occupation, then
    OpenFermion's Givens network rotates it into the determinant. The result
    has exactly n particles on ``qubits`` and leaves all other qubits alone.
    """
    Q = np.asarray(Q)
    n_occupied, n_modes = Q.shape
    if len(qubits) != n_modes:
        raise ValueError(
            f"Q has {n_modes} modes but {len(qubits)} qubits were given")
    for q in qubits[:n_occupied]:
        qc.x(q)
    description = openfermion.circuits.slater_determinant_preparation_circuit(Q)
    for parallel_ops in description:
        for op in parallel_ops:
            if not (isinstance(op, tuple) and len(op) == 4):
                raise NotImplementedError(
                    f"unsupported operation from OpenFermion: {op!r}")
            i, j, theta, phi = op
            qc.append(givens_gate(theta, phi), [qubits[i], qubits[j]])
    return qc
