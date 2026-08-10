"""Guard tests for Slater-determinant state preparation (core._givens_instruction).

The defect these pin down was silent for a long time because fidelity_scan --
the only existing check that touches state preparation -- compares the evolved
state against the PREPARED state. That is self-referential, so it is completely
insensitive to *which* determinant was prepared.

The tests below are deliberately not self-referential: they compare against
independently computed references (the single-particle spectrum, and a Slater
determinant built directly with OpenFermion).

Runs under pytest, or standalone:  python tests/test_state_preparation.py
"""
import numpy as np
import openfermion as of
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

from aqs.core import SSHHModel, build_annealing_circuit, _givens_instruction
from aqs.backends import get_backend

_BK = get_backend("qiskit")


# --------------------------------------------------------------------- helpers
def _bit_reverse(M, n):
    """openfermion places qubit q at bit n-1-q; qiskit places it at bit q."""
    idx = np.arange(1 << n, dtype=np.int64)
    perm = np.zeros_like(idx)
    for p in range(n):
        perm |= ((idx >> p) & 1) << (n - 1 - p)
    return M[np.ix_(perm, perm)] if M.ndim == 2 else M[perm]


def _fidelity(a, b):
    a = np.asarray(a).ravel(); b = np.asarray(b).ravel()
    return float(np.abs(np.vdot(a, b)) ** 2 /
                 (np.vdot(a, a) * np.vdot(b, b)).real)


def _reference_slater(Q, n_qubits):
    """prod_a (sum_i Q[a,i] a_i^dag) |vac>, in qiskit little-endian order."""
    op = None
    for a in range(Q.shape[0]):
        term = of.FermionOperator()
        for i in range(Q.shape[1]):
            term += of.FermionOperator(f"{i}^", complex(Q[a, i]))
        op = term if op is None else op * term
    M = of.get_sparse_operator(op, n_qubits=n_qubits).toarray()
    vac = np.zeros(1 << n_qubits, dtype=complex); vac[0] = 1.0
    psi = _bit_reverse(M @ vac, n_qubits)
    return psi / np.linalg.norm(psi)


# ------------------------------------------------- 1. the Givens gate itself
def _as_matrix(instr):
    qc = QuantumCircuit(2); qc.append(instr, [0, 1])
    return Operator(qc).data


def test_givens_matches_openfermion_convention():
    """G leaves |00>, rotates the single-occupancy block, phases |11> by det G."""
    for theta, phi in [(0.0, 0.0), (0.7, 0.0), (0.7, 0.9), (-1.3, -0.4)]:
        c, s, e = np.cos(theta), np.sin(theta), np.exp(1j * phi)
        want = np.array([
            [1, 0, 0, 0],
            [0, c, s, 0],
            [0, -s * e, c * e, 0],
            [0, 0, 0, e],
        ], dtype=complex)
        got = _as_matrix(_givens_instruction(theta, phi))
        assert np.allclose(got, want, atol=1e-10), f"(theta,phi)=({theta},{phi})"


def test_givens_carries_the_determinant_phase_on_the_doubly_occupied_state():
    """|11> must pick up det G = e^{i phi}; the old decomposition left it at 1."""
    phi = 0.9
    got = _as_matrix(_givens_instruction(0.7, phi))
    assert abs(got[3, 3] - np.exp(1j * phi)) < 1e-10, got[3, 3]


# ------------------------------------------- 2. preparation lands on the GROUND state
def test_prepared_state_has_the_lowest_fill_energy():
    """The decisive check: <H_SSH> must equal the sum of the LOWEST orbitals.

    A negated Givens angle prepares the highest-orbital determinant instead,
    whose energy is the exact negative on a chiral-symmetric chain -- so this
    assertion flips sign rather than drifting, and cannot be missed.
    """
    for N, v, w, Ne, pbc in [(3, 1.0, 1.5, 6, True),
                             (3, 0.5, 1.5, 8, False),
                             (2, 1.0, 0.5, 4, True)]:
        m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
        Qu, Qd, nu, nd = m.slater_Q_matrices()
        sv = _BK.statevector(
            build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0, ramp_U=False))
        L = m.L
        ev = np.sort(np.linalg.eigvalsh(m.H[:L, :L]))
        lowest = ev[:nu].sum() + ev[:nd].sum()

        op = of.FermionOperator()
        for i in range(2 * L):
            for j in range(2 * L):
                if abs(m.H[i, j]) > 1e-12:
                    op += of.FermionOperator(f"{i}^ {j}", complex(m.H[i, j]))
        Hs = _bit_reverse(of.get_sparse_operator(op, n_qubits=2 * L).toarray(), 2 * L)
        E = float(np.real(np.vdot(sv, Hs @ sv)))
        assert abs(E - lowest) < 1e-8, (
            f"N={N} v={v} w={w} PBC={pbc}: prepared E={E:+.6f}, "
            f"lowest fill={lowest:+.6f} (highest fill={-lowest:+.6f})")


def test_prepared_state_equals_the_reference_slater_determinant():
    """Full-state comparison against a determinant built directly in OpenFermion."""
    for N, v, w, Ne, pbc in [(2, 1.0, 1.5, 4, True), (3, 0.5, 1.5, 6, False)]:
        m = SSHHModel.from_total_electrons(N, v, w, Ne, pbc)
        Qu, Qd, nu, nd = m.slater_Q_matrices()
        sv = _BK.statevector(
            build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0, ramp_U=False))
        L = m.L
        # spin-up block occupies qubits [0,L), spin-down [L,2L); build each and
        # take the tensor product in qiskit's little-endian ordering
        up = _reference_slater(Qu, L)
        dn = _reference_slater(Qd, L)
        ref = np.kron(dn, up)          # qiskit: lower qubits are the RIGHT factor
        assert _fidelity(sv, ref) > 1 - 1e-9, f"N={N} PBC={pbc}"


def test_particle_number_is_preserved_by_preparation():
    m = SSHHModel.from_total_electrons(3, 0.5, 1.5, 6, True)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    sv = _BK.statevector(
        build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0, ramp_U=False))
    idx = np.arange(len(sv))
    n = float(np.sum(np.abs(sv) ** 2 *
                     np.array([bin(int(i)).count("1") for i in idx])))
    assert abs(n - 6) < 1e-8, n


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
