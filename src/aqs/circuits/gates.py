"""gates.py - The three two-qubit Pauli-rotation primitives (paper Eq. 37-39, spec 3.1).

    R(theta)  = exp(-i theta/2 (X X + Y Y))         real hopping
    G(theta)  = exp(-i theta/2 (X_i Y_j - Y_i X_j)) imaginary hopping; X on the FIRST operand
    CP(phi)   = exp(-i phi n_i n_j)                  density-density / Hubbard

Every angle sign and operand order here is a cited physics decision pinned by
tests/test_gate_primitives.py and the H_eff round-trip tests; changing one is
a physics change, never a refactor side effect (.claude/rules/gates-and-circuits.md).
"""

from __future__ import annotations

from qiskit import QuantumCircuit


def create_R_gate(theta):
    """R(theta) = exp(-i theta/2 (XX + YY)) : real (Hermitian) hopping term.

    Since rxx(p) = exp(-i p/2 XX) and XX commutes with YY, this is exactly
    rxx(theta) . ryy(theta) -- no factor of two anywhere.
    """
    qc = QuantumCircuit(2, name=f"R({theta:.3f})")
    qc.rxx(theta, 0, 1)
    qc.ryy(theta, 0, 1)
    return qc.to_instruction()


def create_G_gate(theta):
    """G(theta) = exp(-i theta/2 (X0 Y1 - Y0 X1)) : imaginary hopping term.

    Note the operator ordering: X on the FIRST qubit passed, Y on the second.
    G is antisymmetric under swapping the two qubits, unlike R.
    """
    qc = QuantumCircuit(2, name=f"G({theta:.3f})")
    qc.sdg(0); qc.rxx(-theta, 0, 1); qc.s(0)
    qc.sdg(1); qc.rxx(theta, 0, 1); qc.s(1)
    return qc.to_instruction()


def create_CP_gate(theta):
    """CP(theta) = exp(-i theta n_i n_j) : density-density / on-site Hubbard term.

    qiskit's cp(l) = diag(1, 1, 1, e^{+i l}), so the argument is negated to
    realise exp(-i theta n_i n_j) = diag(1, 1, 1, e^{-i theta}).
    """
    qc = QuantumCircuit(2, name=f"CP({theta:.3f})")
    qc.cp(-theta, 0, 1)
    return qc.to_instruction()
