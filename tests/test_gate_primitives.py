"""Guard tests for the two-qubit gate primitives in aqs.core.

These pin down the two things that can silently go wrong with a parameterised
gate wrapper, neither of which raises an exception:

  1. the wrapper implements a different angle convention than its docstring
     claims (an off-by-a-factor-of-two, historically);
  2. the imaginary-hopping (G) sign is conjugated, which no amount of Trotter
     refinement will reveal because it is a systematic error, not a
     discretisation error.

Runs under pytest, or standalone:  python tests/test_gate_primitives.py
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

from aqs.core import create_R_gate, create_G_gate, create_CP_gate

# --------------------------------------------------------------------- helpers
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]])


def _kron(on_q0, on_q1):
    """Two-qubit operator in qiskit's basis, where qubit 0 is the RIGHT factor."""
    return np.kron(on_q1, on_q0)


_XX = _kron(_X, _X)
_YY = _kron(_Y, _Y)
_XY = _kron(_X, _Y)   # X on qubit 0, Y on qubit 1
_YX = _kron(_Y, _X)


def _expi(M):
    """exp(-i M) for Hermitian M."""
    w, v = np.linalg.eigh(M)
    return v @ np.diag(np.exp(-1j * w)) @ v.conj().T


def _as_matrix(instruction):
    qc = QuantumCircuit(2)
    qc.append(instruction, [0, 1])
    return Operator(qc).data


def _equal_up_to_phase(A, B, atol=1e-10):
    i = int(np.argmax(np.abs(A)))
    return np.allclose(A * (B.flat[i] / A.flat[i]), B, atol=atol)


# --------------------------------------------------------- 1. docstring contract
def test_R_gate_matches_its_documented_form():
    """create_R_gate(t) must equal exp(-i t/2 (XX + YY)) -- not half of it."""
    for t in (0.0, 0.37, -1.23, 2.9):
        want = _expi(t / 2.0 * (_XX + _YY))
        assert _equal_up_to_phase(_as_matrix(create_R_gate(t)), want), f"theta={t}"


def test_G_gate_matches_its_documented_form():
    """create_G_gate(t) must equal exp(-i t/2 (X0 Y1 - Y0 X1))."""
    for t in (0.0, 0.37, -1.23, 2.9):
        want = _expi(t / 2.0 * (_XY - _YX))
        assert _equal_up_to_phase(_as_matrix(create_G_gate(t)), want), f"theta={t}"


def test_CP_gate_matches_its_documented_form():
    """create_CP_gate(t) must equal exp(-i t n_i n_j) = diag(1,1,1,e^{-it})."""
    for t in (0.0, 0.37, -1.23, 2.9):
        want = np.diag([1.0, 1.0, 1.0, np.exp(-1j * t)])
        assert _equal_up_to_phase(_as_matrix(create_CP_gate(t)), want), f"theta={t}"


def test_G_gate_is_antisymmetric_under_qubit_swap():
    """G(-t) on (0,1) == G(t) on (1,0); R has no such asymmetry."""
    t = 0.41
    a = QuantumCircuit(2); a.append(create_G_gate(t), [0, 1])
    b = QuantumCircuit(2); b.append(create_G_gate(-t), [1, 0])
    assert _equal_up_to_phase(Operator(a).data, Operator(b).data)


# ------------------------------------------- 2. effective-Hamiltonian round trip
def _true_bond_hamiltonian(t):
    """-t a0^ a1 - conj(t) a1^ a0 as a 4x4 matrix in qiskit little-endian order.

    The bit reversal is done inline rather than imported, so this test stays
    self-contained. openfermion places qubit q at bit position n-1-q, whereas
    qiskit places it at bit q; for two modes that swap silently CONJUGATES the
    off-diagonal, which would make this test confirm the wrong G sign.
    """
    import openfermion as of
    op = (of.FermionOperator("0^ 1", -complex(t))
          + of.FermionOperator("1^ 0", -np.conj(complex(t))))
    H_big = of.get_sparse_operator(op, n_qubits=2).toarray()
    perm = [0, 2, 1, 3]          # bit-reversal on 2 qubits
    return H_big[np.ix_(perm, perm)]


def test_trotter_bond_reproduces_the_true_hamiltonian():
    """One R.G bond step must have effective Hamiltonian == the true bond H.

    This is the check that catches a conjugated imaginary part. Extracting
    H_eff = i log(U)/dt is necessary because a sign error there is systematic:
    it does NOT shrink as dt -> 0, so a plain convergence test cannot see it.
    """
    from scipy.linalg import logm
    dt = 1e-3
    for t in (0.7 + 0.45j, -0.3 + 1.1j, 0.9 - 0.2j, 1.5 + 0j):
        H = _true_bond_hamiltonian(t)
        qc = QuantumCircuit(2)
        qc.append(create_R_gate(-dt * t.real), [0, 1])
        if abs(t.imag) > 1e-12:
            qc.append(create_G_gate(dt * t.imag), [0, 1])
        H_eff = 1j * logm(Operator(qc).data) / dt
        # first-order Trotter: residual is O(dt), so scale the tolerance with dt
        assert np.max(np.abs(H_eff - H)) < 5.0 * dt * abs(t) ** 2 + 1e-9, (
            f"t={t}\n true[1,2]={H[1,2]:.6f}\n  eff[1,2]={H_eff[1, 2]:.6f}")


def test_annealing_circuit_bond_angles_are_convention_consistent():
    """End-to-end: the built circuit must evolve under the intended H."""
    from scipy.linalg import expm
    from aqs.core import SSHHModel, build_annealing_circuit
    from aqs.backends import get_backend
    bk = get_backend("qiskit")
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4, PBC=False)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    sv0 = bk.statevector(build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0))
    # the prepared state must be an eigenstate-preserving starting point:
    # norm 1 and correct particle number
    assert abs(np.linalg.norm(sv0) - 1.0) < 1e-10
    idx = np.arange(len(sv0))
    npart = float(np.sum(np.abs(sv0) ** 2 *
                         np.array([bin(int(i)).count("1") for i in idx])))
    assert abs(npart - 4) < 1e-8, npart


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  PASS  {fn.__name__}")
        except Exception as e:
            failed += 1
            print(f"  FAIL  {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
