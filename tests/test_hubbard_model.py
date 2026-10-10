"""Guard tests for the 1D Hubbard model (SSHH with the A/B distinction off).

The Hubbard chain is the paper's spinful SSH-Hubbard model with v = w = t and
U_A = U_B = U, so the circuit is ``aqs.core.build_annealing_circuit`` unchanged.
Layout: site i spin-up on qubit i, spin-down on qubit n_sites + i.

Every physics assertion compares against an independent reference:

* the second-quantized Hamiltonian built here with OpenFermion and converted
  to little-endian (never the module's own Pauli construction alone),
* the analytic tight-binding spectrum -2|t| cos(k + arg t),
* exact diagonalization restricted to the fixed (n_up, n_dn) sector,
* a Trotter-free integration of the same U ramp.

Inputs are asymmetric (n_up != n_dn, hand-placed product states) and use a
complex hopping where a symmetric or real input would be blind to the bug.

Reverse verification (2026-10-05): each mutation below was applied to
src/aqs/models/hubbard.py and the suite confirmed to fail, then reverted:

* conjugate the bulk Im sign in the Pauli builder          -> 2 tests fail
* drop the Jordan-Wigner string on the wrap bond           -> 2 tests fail
* fill the determinant from the highest orbitals           -> 4 tests fail
* pair (i, i+1) instead of (i, n_sites+i) in the occupations -> 4 tests fail
* S^z from n_up + n_dn instead of n_up - n_dn              -> 2 tests fail
* pass conj(t) to the circuit builder                      -> 2 tests fail
"""

from __future__ import annotations

import numpy as np
import openfermion
import pytest
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from scipy.linalg import expm

from aqs.mapping import _to_little_endian
from aqs.models.hubbard import (
    HubbardModel,
    HubbardSim,
    hubbard_pauli_hamiltonian,
    hubbard_hamiltonian,
    hubbard_ground_state_fixed_n,
    sector_indices,
)


# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

def _openfermion_hamiltonian(n_sites, t, U, pbc):
    """-t c_i^dag c_{i+1} + h.c. per spin, + U n_up n_dn; little-endian dense."""
    op = openfermion.FermionOperator()
    bonds = [(i, i + 1) for i in range(n_sites - 1)]
    if pbc:
        bonds.append((n_sites - 1, 0))
    for offset in (0, n_sites):
        for i, j in bonds:
            op += openfermion.FermionOperator(((offset + i, 1), (offset + j, 0)), -t)
            op += openfermion.FermionOperator(((offset + j, 1), (offset + i, 0)), -np.conj(t))
    for i in range(n_sites):
        op += openfermion.FermionOperator(
            ((i, 1), (i, 0), (n_sites + i, 1), (n_sites + i, 0)), U)
    big = openfermion.get_sparse_operator(op, n_qubits=2 * n_sites).tocsr()
    return _to_little_endian(big, 2 * n_sites).toarray()


def _fid(a, b):
    return float(abs(np.vdot(a, b)) ** 2)


def _exact_ramp(model, T, psi0, steps=4000):
    """Trotter-free evolution under H(U=0) + (t/T) (H(U) - H(U=0))."""
    idx = sector_indices(model)
    H0 = _openfermion_hamiltonian(model.n_sites, model.t, 0.0, model.PBC)[np.ix_(idx, idx)]
    HU = _openfermion_hamiltonian(model.n_sites, model.t, model.U, model.PBC)[np.ix_(idx, idx)]
    hu = np.real(np.diag(HU - H0))
    dt = T / steps
    U0 = expm(-1j * dt * H0)
    v = psi0[idx].copy()
    for k in range(steps):
        half = np.exp(-0.5j * dt * hu * (k + 0.5) / steps)
        v = half * (U0 @ (half * v))
    out = np.zeros_like(psi0)
    out[idx] = v
    return out


def _product_state(n_sites, up_sites, dn_sites):
    qc = QuantumCircuit(2 * n_sites)
    for i in up_sites:
        qc.x(i)
    for i in dn_sites:
        qc.x(n_sites + i)
    return Statevector(qc).data


# ---------------------------------------------------------------------------
# Model definition
# ---------------------------------------------------------------------------

class TestModelDefinition:
    def test_is_the_sshh_model_with_equal_hoppings(self):
        model = HubbardModel(n_sites=6, t=0.8, U=2.0)
        sshh = model.sshh()
        assert sshh.v == sshh.w == 0.8
        assert sshh.L == 6 and sshh.N_cells == 3
        assert (sshh.num_up, sshh.num_dn) == (3, 3)

    def test_defaults_to_half_filling(self):
        model = HubbardModel(n_sites=6)
        assert (model.num_up, model.num_dn) == (3, 3)

    def test_odd_site_count_is_rejected(self):
        with pytest.raises(ValueError, match="even"):
            HubbardModel(n_sites=5)

    def test_single_particle_spectrum_is_the_tight_binding_band(self):
        """Ring levels are -2|t| cos(k + arg t), k = 2 pi m / n; complex t on purpose."""
        t = 0.9 * np.exp(0.3j)
        model = HubbardModel(n_sites=6, t=t, PBC=True)
        got = np.sort(np.linalg.eigvalsh(model.single_particle_hamiltonian()))
        k = 2 * np.pi * np.arange(6) / 6
        want = np.sort(-2 * abs(t) * np.cos(k + np.angle(t)))
        assert np.allclose(got, want, atol=1e-12)

    def test_shell_gap_flags_the_open_shell_ring(self):
        assert HubbardModel(n_sites=6, PBC=True).shell_gap() == pytest.approx(2.0)
        assert HubbardModel(n_sites=8, PBC=True).shell_gap() == pytest.approx(0.0, abs=1e-12)
        assert HubbardModel(n_sites=8, PBC=False).shell_gap() > 0.5


# ---------------------------------------------------------------------------
# Reference Hamiltonian (little-endian) versus OpenFermion
# ---------------------------------------------------------------------------

class TestReferenceHamiltonian:
    @pytest.mark.parametrize("pbc", [True, False])
    @pytest.mark.parametrize("t", [1.0, 0.7 + 0.4j])
    def test_pauli_hamiltonian_matches_openfermion(self, pbc, t):
        model = HubbardModel(n_sites=4, t=t, U=1.7, PBC=pbc, num_up=2, num_dn=1)
        got = hubbard_pauli_hamiltonian(model).to_matrix()
        want = _openfermion_hamiltonian(4, t, 1.7, pbc)
        assert np.abs(got - want).max() < 1e-12

    def test_interaction_counts_doubly_occupied_sites(self):
        model = HubbardModel(n_sites=4, t=1.0, U=3.0, PBC=False)
        H = hubbard_hamiltonian(model).toarray()
        idx = int(np.argmax(np.abs(_product_state(4, up_sites=[0, 2], dn_sites=[2, 3]))))
        assert H[idx, idx].real == pytest.approx(3.0)      # only site 2 is doubly occupied

    def test_sector_indices_fix_both_spin_counts(self):
        model = HubbardModel(n_sites=4, num_up=2, num_dn=1)
        idx = sector_indices(model)
        assert len(idx) == 6 * 4
        for i in idx:
            assert bin(int(i) & 0b1111).count("1") == 2
            assert bin(int(i) >> 4).count("1") == 1


# ---------------------------------------------------------------------------
# Circuits
# ---------------------------------------------------------------------------

class TestCircuits:
    @pytest.mark.parametrize("pbc,n_up,n_dn", [(False, 2, 1), (True, 1, 3), (True, 3, 1)])
    def test_prepared_state_is_the_lowest_fill_determinant(self, pbc, n_up, n_dn):
        """Energy = sum of the lowest levels and fidelity 1 with fixed-sector ED.

        The historical failure prepares the highest band, at the opposite energy.
        """
        model = HubbardModel(n_sites=4, t=1.0, U=0.0, PBC=pbc, num_up=n_up, num_dn=n_dn)
        sim = HubbardSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        levels = np.sort(np.linalg.eigvalsh(model.single_particle_hamiltonian()))
        want = levels[:n_up].sum() + levels[:n_dn].sum()
        assert sim.energy(psi) == pytest.approx(want, abs=1e-10)
        assert _fid(hubbard_ground_state_fixed_n(model), psi) == pytest.approx(1.0, abs=1e-10)

    def test_open_shell_ring_is_refused(self):
        """An 8-site ring at half filling has a degenerate U = 0 ground state."""
        with pytest.raises(ValueError, match="open shell"):
            HubbardSim(HubbardModel(n_sites=8, PBC=True), T=1.0, L=4)

    def test_evolution_conserves_both_spin_counts(self):
        model = HubbardModel(n_sites=4, t=1.0, U=2.5, PBC=False, num_up=2, num_dn=1)
        psi = HubbardSim(model, T=1.0, L=6).evolve()
        idx = sector_indices(model)
        assert np.vdot(psi[idx], psi[idx]).real == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.parametrize("pbc,n_up,n_dn", [(False, 2, 1), (True, 3, 1)])
    def test_trotter_error_vanishes_with_L_for_complex_hopping(self, pbc, n_up, n_dn):
        """Against a Trotter-free ramp the infidelity must fall as 1/L^2.

        A conjugated hopping or wrong wrap sign evolves under a different
        Hamiltonian; that error does not shrink with L.
        """
        model = HubbardModel(n_sites=4, t=0.8 * np.exp(0.35j), U=2.0, PBC=pbc,
                             num_up=n_up, num_dn=n_dn)
        psi0 = Statevector(HubbardSim(model, T=1.0, L=1).preparation_circuit()).data
        ref = _exact_ramp(model, 1.0, psi0)
        err = {L: 1 - _fid(ref, HubbardSim(model, T=1.0, L=L).evolve()) for L in (20, 40)}
        assert err[40] < 2e-3
        assert err[20] / err[40] == pytest.approx(4.0, rel=0.15)

    def test_slow_ramp_reaches_the_interacting_ground_state(self):
        model = HubbardModel(n_sites=4, t=1.0, U=2.0, PBC=False)
        gs = hubbard_ground_state_fixed_n(model)
        start = Statevector(HubbardSim(model, T=1.0, L=1).preparation_circuit()).data
        end = HubbardSim(model, T=8.0, L=320).evolve()
        assert _fid(gs, end) > 0.995
        assert _fid(gs, end) > _fid(gs, start)


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------

class TestObservables:
    def _sim(self, n_sites=4):
        return HubbardSim(HubbardModel(n_sites=n_sites, PBC=False), T=1.0, L=1)

    def test_double_occupancy_of_a_hand_placed_state(self):
        """Up on sites 0, 2; down on sites 2, 3: one doubly occupied site of four."""
        psi = _product_state(4, up_sites=[0, 2], dn_sites=[2, 3])
        assert self._sim().double_occupancy(psi) == pytest.approx(0.25)

    def test_spin_correlation_of_a_hand_placed_state(self):
        """S^z = (+1/2, 0, 0, -1/2) for up on 0, 2 and down on 2, 3."""
        psi = _product_state(4, up_sites=[0, 2], dn_sites=[2, 3])
        C = self._sim().spin_correlation(psi)
        sz = np.array([0.5, 0.0, 0.0, -0.5])
        assert np.allclose(C, np.outer(sz, sz), atol=1e-12)

    def test_structure_factor_of_the_neel_state(self):
        """Up on even sites, down on odd sites: S(pi) = n_sites / 4."""
        psi = _product_state(4, up_sites=[0, 2], dn_sites=[1, 3])
        sim = self._sim()
        assert sim.spin_structure_factor(psi, q=np.pi) == pytest.approx(1.0)
        assert sim.spin_structure_factor(psi, q=0.0) == pytest.approx(0.0, abs=1e-12)

    def test_free_fermion_double_occupancy_is_the_density_product(self):
        """U = 0 on a closed-shell ring: uniform density 1/2 per spin, D = 1/4."""
        model = HubbardModel(n_sites=6, PBC=True)
        sim = HubbardSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        assert sim.double_occupancy(psi) == pytest.approx(0.25, abs=1e-10)

    def test_repulsion_suppresses_double_occupancy_in_the_exact_ground_state(self):
        sim = self._sim()
        d = [sim.double_occupancy(hubbard_ground_state_fixed_n(
                HubbardModel(n_sites=4, U=U, PBC=False))) for U in (0.0, 4.0, 8.0)]
        assert d[0] > d[1] > d[2]
        assert d[2] < 0.06

    def test_energy_matches_the_reference_hamiltonian(self):
        model = HubbardModel(n_sites=4, t=1.0, U=2.0, PBC=False, num_up=2, num_dn=1)
        sim = HubbardSim(model, T=1.0, L=6)
        psi = sim.evolve()
        H = _openfermion_hamiltonian(4, 1.0, 2.0, False)
        assert sim.energy(psi) == pytest.approx(np.vdot(psi, H @ psi).real, abs=1e-10)
