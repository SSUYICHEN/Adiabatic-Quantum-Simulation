"""Physics smoke tests: paper signatures realized in the spinless model.

Each test asserts one *physical characteristic* from the paper
(.claude/specs/paper_theory/SSHH.md) in its spinless-model form
(.claude/specs/reconstruct_goals.md). These are qualitative gates -- fast,
statevector-exact, thresholds calibrated with wide margins against measured
values -- meant to catch "the physics silently changed", not to re-verify
conventions (that is tests/test_spinless_model.py's job).

Paper section map:
  chiral spectrum / gap closing ......... Sec. III (Eq. 24)
  Berry quantization + pi step .......... Sec. III-B, Fig. 3 (Eq. 19)
  robustness / breakdown under V ........ Sec. III-B, V-B
  edge-localized polarization ........... Sec. III-B, Fig. 4 (Eq. 22)
  edge-occupation staircase ............. Sec. III-B (Eq. 21)
  Trotter/adiabatic convergence ......... Sec. V-A, Fig. 2

Note on interactions: the paper breaks chiral symmetry with on-site
U_A != U_B. The spinless model's nearest-neighbor V_v, V_w live on bonds and
preserve inversion symmetry, so the Berry phase stays *quantized* for any
(V_v, V_w); the topological signature instead breaks down through a
first-order-like transition at strong coupling (measured: the topological
gamma flips by pi between V=2 and V=4 at v=0.5, w=1.5). The robustness and
breakdown tests below assert exactly that behavior.
"""

from __future__ import annotations

import numpy as np
import pytest
from qiskit.quantum_info import Statevector

from aqs.models.spinless import SpinlessSSHModel, SpinlessSSHHSim
from aqs.observables import density_profile_qubits


def _prepared_state(model: SpinlessSSHModel) -> np.ndarray:
    """Exact SSH ground-state determinant (no evolution needed)."""
    sim = SpinlessSSHHSim(model, T=1.0, L=1)
    return Statevector(sim.preparation_circuit()).data


def _berry(model: SpinlessSSHModel, psi: np.ndarray) -> float:
    return SpinlessSSHHSim(model, T=1.0, L=1).berry_phase(psi)


def _pi_distance(angle: float) -> float:
    """Distance of an angle from the nearest multiple of pi."""
    return float(abs((angle + np.pi / 2) % np.pi - np.pi / 2))


# ---------------------------------------------------------------------------
# Band structure (paper Sec. III, Eq. 24)
# ---------------------------------------------------------------------------

class TestBandStructure:
    def test_chiral_symmetry_makes_the_spectrum_symmetric(self):
        """A<->B chiral symmetry: single-particle levels come in +/-E pairs."""
        for pbc in (True, False):
            model = SpinlessSSHModel(N_cells=6, v=0.7, w=1.3, PBC=pbc)
            eigs = np.sort(
                np.linalg.eigvalsh(model.single_particle_hamiltonian())
            )
            np.testing.assert_allclose(eigs, -eigs[::-1], atol=1e-12)

    def test_gap_follows_eq24_and_closes_at_the_transition(self):
        """PBC gap = 2 min_k |v + w e^{ik}|; exactly zero at v = w (even N)."""
        def gap(v, w):
            model = SpinlessSSHModel(N_cells=6, v=v, w=w, PBC=True)
            eigs = np.linalg.eigvalsh(model.single_particle_hamiltonian())
            return 2 * np.min(np.abs(eigs))

        # away from the transition the finite-size gap is bounded below by
        # the thermodynamic-limit formula Delta = 2 min{|v+w|, |v-w|}
        for v, w in ((1.5, 0.5), (0.5, 1.5), (1.0, 0.25)):
            assert gap(v, w) >= 2 * min(abs(v + w), abs(v - w)) - 1e-12
        # N=6 contains k=pi, so the gap closes exactly at the critical point
        assert gap(1.0, 1.0) < 1e-12
        assert gap(1.0, 1.0) < gap(0.99, 1.01) < gap(0.5, 1.5)


# ---------------------------------------------------------------------------
# Many-body Berry phase (paper Sec. III-B / Fig. 3)
# ---------------------------------------------------------------------------

class TestBerryPhase:
    @pytest.mark.parametrize("v,w", [(1.5, 0.5), (0.5, 1.5), (1.0, 0.25), (0.25, 1.0)])
    def test_berry_phase_is_quantized_to_multiples_of_pi(self, v, w):
        """Inversion symmetry pins gamma to 0 or pi (mod 2 pi) in both phases."""
        model = SpinlessSSHModel(N_cells=6, v=v, w=w, PBC=True)
        gamma = _berry(model, _prepared_state(model))
        assert _pi_distance(gamma) < 1e-6

    def test_berry_phase_steps_by_pi_across_the_transition(self):
        """The Fig. 3 step: gamma jumps by pi between w < v and w > v.

        Only the jump is asserted; absolute values carry a constant
        filling-dependent offset (rules: models-and-observables.md).
        """
        gammas = []
        for v, w in ((1.0, 0.75), (0.75, 1.0)):  # just either side of v = w
            model = SpinlessSSHModel(N_cells=6, v=v, w=w, PBC=True)
            gammas.append(_berry(model, _prepared_state(model)))
        jump = abs(gammas[1] - gammas[0]) % (2 * np.pi)
        assert abs(jump - np.pi) < 1e-6

    def test_berry_phase_is_robust_under_weak_interaction(self):
        """Weak V (inside the Eq. 25-analog window): gamma unchanged.

        v=0.5, w=1.5 -> min{|v+w|, |v-w|} = 1, so V = 0.5 is 'weak'.
        Adiabatic evolution at the spec operating point T=1, L=40.
        """
        model0 = SpinlessSSHModel(N_cells=6, v=0.5, w=1.5, PBC=True)
        gamma0 = _berry(model0, _prepared_state(model0))
        sim = SpinlessSSHHSim(model0, T=1.0, L=40)
        for V in (0.1, 0.5):
            (rec,) = sim.interaction_sweep(V_points=((V, V),))
            diff = abs(rec["berry_phase"] - gamma0) % (2 * np.pi)
            assert min(diff, 2 * np.pi - diff) < 1e-3

    def test_topological_signature_breaks_down_at_strong_interaction(self):
        """Strong V: the topological gamma flips by pi (measured at V=4).

        The bond interaction preserves inversion, so the breakdown appears as
        a quantized pi flip (CDW-like transition) rather than the smooth
        de-quantization the paper's on-site Delta-U produces.
        """
        model0 = SpinlessSSHModel(N_cells=6, v=0.5, w=1.5, PBC=True)
        gamma0 = _berry(model0, _prepared_state(model0))
        sim = SpinlessSSHHSim(model0, T=1.0, L=40)
        (rec,) = sim.interaction_sweep(V_points=((4.0, 4.0),))
        diff = abs(rec["berry_phase"] - gamma0) % (2 * np.pi)
        assert abs(diff - np.pi) < 1e-2
        assert _pi_distance(rec["berry_phase"]) < 1e-2  # still quantized


# ---------------------------------------------------------------------------
# Sublattice polarization under OBC (paper Sec. III-B / Fig. 4)
# ---------------------------------------------------------------------------

class TestPolarization:
    def test_edge_localized_in_the_topological_phase(self):
        """w > v, n = N+1: polarization lives at the edges, not the bulk."""
        model = SpinlessSSHModel(N_cells=6, v=0.1, w=1.0, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        profile = sim.polarization_profile(_prepared_state(model))
        assert abs(profile[0]) > 0.4
        assert max(abs(p) for p in profile[1:-1]) < 0.05
        # inversion antisymmetry of the profile
        np.testing.assert_allclose(profile, -profile[::-1], atol=1e-9)

    def test_no_edge_polarization_in_the_trivial_phase(self):
        """v > w: no edge modes (bulk-boundary correspondence, w=0 case)."""
        model = SpinlessSSHModel(
            N_cells=6, v=1.0, w=0.1, PBC=False, n_particles=6
        )
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        profile = sim.polarization_profile(_prepared_state(model))
        assert max(abs(p) for p in profile) < 0.05

    def test_edge_polarization_survives_weak_interaction(self):
        """Fig. 4 bottom analog: weak V keeps the edge feature intact."""
        model = SpinlessSSHModel(N_cells=6, v=0.1, w=1.0, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=40)
        recs = sim.interaction_sweep(V_points=((0.0, 0.0), (0.1, 0.1)))
        assert abs(recs[1]["edge_polarization"] - recs[0]["edge_polarization"]) < 0.05
        assert abs(recs[1]["edge_polarization"]) > 0.4


# ---------------------------------------------------------------------------
# Edge-occupation staircase (paper Eq. 21)
# ---------------------------------------------------------------------------

class TestEdgeOccupation:
    def test_edge_occupation_jumps_when_the_edge_modes_fill(self):
        """<n_edge> vs n is a staircase: the two near-zero edge modes of the
        topological OBC chain fill at n = N and n = N+1 (spinless counting),
        producing increments far larger than any bulk-level increment.

        Measured at N=6, v=0.1, w=1.0: edge steps contribute 0.99 each while
        every bulk step is <= 0.34 (N=4 is too short -- its bulk levels carry
        up to 0.50 of edge weight because the edge cells are half the chain).
        """
        N = 6
        edge_qubits = (0, 1, 2 * N - 2, 2 * N - 1)  # cells 0 and N-1
        occupations = []
        for n in range(1, 2 * N + 1):
            model = SpinlessSSHModel(
                N_cells=N, v=0.1, w=1.0, PBC=False, n_particles=n
            )
            density = density_profile_qubits(_prepared_state(model), 2 * N)
            occupations.append(sum(density[q] for q in edge_qubits))
        increments = np.diff([0.0] + occupations)
        edge_steps = {N, N + 1}  # n at which the edge pair fills
        for n in range(1, 2 * N + 1):
            if n in edge_steps:
                assert increments[n - 1] > 0.8
            else:
                assert increments[n - 1] < 0.5
        assert min(increments[n - 1] for n in edge_steps) > 2 * max(
            increments[n - 1] for n in range(1, 2 * N + 1) if n not in edge_steps
        )


# ---------------------------------------------------------------------------
# Trotter / adiabatic convergence (paper Sec. V-A / Fig. 2)
# ---------------------------------------------------------------------------

class TestConvergencePattern:
    def test_fidelity_improves_with_L_and_degrades_with_T(self):
        """Fig. 2 pattern: at fixed L the fidelity is worse for longer T;
        at fixed T it improves monotonically with L."""
        model = SpinlessSSHModel(N_cells=3, v=0.5, w=1.5, PBC=True)
        fid = {}
        for T in (1.0, 15.0):
            sim = SpinlessSSHHSim(model, T=T, L=1)
            records = sim.trotter_benchmark(L_values=(10, 40))
            fid[T] = [r["fidelity"] for r in records]
            assert fid[T][1] > fid[T][0]
        assert fid[1.0][1] > 0.999
        assert fid[1.0][1] > fid[15.0][1]
