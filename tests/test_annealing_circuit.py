"""Tests for the SSH-Hubbard model and its Trotterised annealing circuit.

core.build_annealing_circuit's Trotter loop had NO test coverage, despite being
where the gate-angle conventions are actually consumed. These tests cover the
structural rules (which gates get emitted where, and the PBC parity correction)
and the numerical contract (the circuit evolves under the intended Hamiltonian).
"""
import numpy as np
import pytest
from scipy.linalg import expm

from aqs.core import SSHHModel, build_annealing_circuit
from aqs.backends import get_backend

_BK = get_backend("qiskit")


def _fid(a, b):
    a = np.asarray(a).ravel(); b = np.asarray(b).ravel()
    return float(np.abs(np.vdot(a, b)) ** 2 /
                 (np.vdot(a, a) * np.vdot(b, b)).real)


def _gate_names(qc):
    return [inst.operation.name for inst in qc.data]


def _decomposed(qc, name):
    """Yield (angle, sorted_qubits) for every `name` gate one level down.

    create_*_gate wrap their bodies with to_instruction(), which exposes an
    empty params list -- the angle survives only inside the sub-circuit and in
    the display label. decompose() once to read the real parameters.
    """
    d = qc.decompose()
    for inst in d.data:
        if inst.operation.name == name:
            yield (float(inst.operation.params[0]),
                   sorted(d.find_bit(b).index for b in inst.qubits))


# --------------------------------------------------------------- the model
def test_from_total_electrons_splits_spin_sectors():
    for Ne, want in [(0, (0, 0)), (1, (1, 0)), (4, (2, 2)), (7, (4, 3))]:
        m = SSHHModel.from_total_electrons(3, 0.5, 1.5, Ne)
        assert (m.num_up, m.num_dn) == want, Ne
        assert m.number_of_electrons == Ne


def test_single_particle_hamiltonian_is_hermitian_and_bipartite():
    for pbc in (True, False):
        m = SSHHModel.from_total_electrons(3, 0.5 + 0.2j, 1.5 - 0.1j, 6, PBC=pbc)
        H = m.H
        assert np.allclose(H, H.conj().T), pbc
        # two decoupled spin blocks of equal size
        L = m.L
        assert np.allclose(H[:L, L:], 0) and np.allclose(H[L:, :L], 0)
        assert np.allclose(H[:L, :L], H[L:, L:])


def test_obc_drops_the_wrap_bond_pbc_keeps_it():
    L = 6
    pbc = SSHHModel.from_total_electrons(3, 0.5, 1.5, 6, PBC=True).H
    obc = SSHHModel.from_total_electrons(3, 0.5, 1.5, 6, PBC=False).H
    assert abs(pbc[L - 1, 0]) > 0        # wrap bond present
    assert abs(obc[L - 1, 0]) == 0       # and absent
    # every other matrix element identical
    d = pbc - obc
    d[L - 1, 0] = d[0, L - 1] = d[2 * L - 1, L] = d[L, 2 * L - 1] = 0
    assert np.allclose(d, 0)


def test_slater_q_selects_the_lowest_orbitals():
    m = SSHHModel.from_total_electrons(3, 0.5, 1.5, 6)
    Qu, Qd, nu, nd = m.slater_Q_matrices()
    L = m.L
    ev, evec = np.linalg.eigh(m.H[:L, :L])
    order = np.argsort(ev)
    assert Qu.shape == (nu, L) and Qd.shape == (nd, L)
    # rows span the same subspace as the lowest nu eigenvectors
    want = evec[:, order][:, :nu]
    proj_want = want @ want.conj().T
    proj_got = Qu.conj().T @ Qu
    assert np.allclose(proj_want, proj_got, atol=1e-10)


# ------------------------------------------------- circuit structure
def test_steps_zero_emits_only_state_preparation():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 1.0, 1.0, T_A=1.0, steps=0)
    names = set(_gate_names(qc))
    assert not any(n.startswith(("R(", "G(", "CP(")) for n in names), names


def test_real_hopping_emits_no_G_gates():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=1.0, steps=3)
    assert not any(n.startswith("G(") for n in _gate_names(qc))
    assert any(n.startswith("R(") for n in _gate_names(qc))


def test_complex_hopping_emits_G_gates():
    m = SSHHModel.from_total_electrons(2, 1.0 + 0.4j, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=1.0, steps=3)
    assert any(n.startswith("G(") for n in _gate_names(qc))


def test_zero_U_emits_no_CP_gates():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 0.0, 0.0, T_A=1.0, steps=3)
    assert not any(n.startswith("CP(") for n in _gate_names(qc))


def test_ramp_makes_CP_angles_grow_linearly():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    steps = 4
    qc = build_annealing_circuit(m, Qu, Qd, 1.0, 1.0, T_A=1.0, steps=steps, ramp_U=True)
    angles = sorted({round(a, 12) for a, _ in _decomposed(qc, "cp")})
    # lambda_l = (2l-1)/(2L) for l = 1..L  ->  equally spaced
    assert len(angles) == steps, angles
    diffs = np.diff(angles)
    assert np.allclose(diffs, diffs[0], atol=1e-12), angles


def test_no_ramp_gives_a_single_CP_angle():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 1.0, 1.0, T_A=1.0, steps=4, ramp_U=False)
    angles = {round(a, 12) for a, _ in _decomposed(qc, "cp")}
    assert len(angles) == 1, angles


def test_pbc_wrap_bond_carries_the_parity_factor():
    """The wrap bond's R angle is -(-1)^n times the bulk angle; OBC has no wrap bond."""
    for Ne in (4, 6, 7):
        m = SSHHModel.from_total_electrons(3, 1.0, 1.5, Ne, PBC=True)
        Qu, Qd, nu, nd = m.slater_Q_matrices()
        qc = build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=1.0, steps=1)
        L = m.L
        # real hoppings -> every rxx comes from an R gate, none from G
        rxx = dict((tuple(q), a) for a, q in _decomposed(qc, "rxx"))
        wrap, bulk = rxx.get((0, L - 1)), rxx.get((1, 2))
        assert wrap is not None and bulk is not None, sorted(rxx)
        assert abs(wrap - (-((-1) ** nu) * bulk)) < 1e-12, (Ne, wrap, bulk, nu)


def test_obc_has_no_wrap_bond_gate():
    m = SSHHModel.from_total_electrons(3, 1.0, 1.5, 6, PBC=False)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    qc = build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=1.0, steps=1)
    L = m.L
    pairs = {tuple(q) for _, q in _decomposed(qc, "rxx")}
    assert (0, L - 1) not in pairs, pairs


# ------------------------------------------------- numerical contract
def test_noninteracting_evolution_is_stationary_and_trotter_converges():
    """U = 0: the prepared state is a Slater determinant of H_SSH eigenstates, so
    exact evolution only multiplies it by a phase and |<psi0|psi(T)>| = 1.

    The circuit is Trotterised, so the deficit is pure discretisation error and
    must shrink as steps grows. Asserting the trend, not just a threshold,
    distinguishes 'converging correctly' from 'happens to be close'.
    """
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4, PBC=False)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    sv0 = _BK.statevector(build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0, steps=0))

    deficits = []
    for steps in (20, 80, 320):
        sv = _BK.statevector(
            build_annealing_circuit(m, Qu, Qd, 0, 0, T_A=0.6, steps=steps))
        deficits.append(1.0 - _fid(sv, sv0))
    assert all(b < a for a, b in zip(deficits, deficits[1:])), deficits
    assert deficits[-1] < 1e-5, deficits


def test_particle_number_is_conserved_through_the_ramp():
    m = SSHHModel.from_total_electrons(2, 1.0, 1.5, 4, PBC=True)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    sv = _BK.statevector(
        build_annealing_circuit(m, Qu, Qd, 0.7, 1.4, T_A=1.0, steps=10))
    idx = np.arange(len(sv))
    n = float(np.sum(np.abs(sv) ** 2 *
                     np.array([bin(int(i)).count("1") for i in idx])))
    assert abs(n - 4) < 1e-8, n


def test_circuit_stays_normalised():
    m = SSHHModel.from_total_electrons(2, 0.8 + 0.3j, 1.2, 4, PBC=True)
    Qu, Qd, _, _ = m.slater_Q_matrices()
    for steps in (0, 1, 7):
        sv = _BK.statevector(
            build_annealing_circuit(m, Qu, Qd, 0.5, 0.9, T_A=1.0, steps=steps))
        assert abs(np.linalg.norm(sv) - 1.0) < 1e-10, steps
