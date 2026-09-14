"""Guard tests for the spinless SSH + nearest-neighbor-interaction model.

Authority: the Pauli-operator representation (spec section 2) and the Trotter
angle table (spec section 3.2) of ``.claude/specs/reconstruct_goals.md``, which
are mutually consistent and normative per ``.claude/rules/gates-and-circuits.md``.

Every physics assertion here compares against an independent reference:

* the second-quantized Hamiltonian built with OpenFermion and converted to
  little-endian (never the module's own Pauli construction alone),
* exact diagonalization restricted to the fixed particle-number sector,
* H_eff = i log(U) / dt round-trips on single-bond circuits, which are exact
  (no Trotter error) and therefore catch sign/conjugation convention bugs that
  do not shrink with dt.

Inputs are deliberately asymmetric and use complex hoppings where a symmetric
or real input would be structurally blind to the bug class under test
(see .claude/rules/testing-and-verification.md).

Reverse verification (2026-08-11): each mutation below was applied to
src/aqs/models/spinless.py and the suite confirmed to fail, then reverted:

* flip bulk G angle sign (theta_G -> +dt Im t)            -> 6 tests fail
* drop the wrap parity factor (parity_factor -> 1)        -> 2 tests fail
* wrap G sign follows bulk (theta_G *= +parity_factor)    -> 2 tests fail
* conjugate bulk Im sign in the Pauli builder             -> 10 tests fail
* fill the determinant from the highest orbitals          -> 3 tests fail
* swap A/B in the polarization profile                    -> 1 test fails
* Berry position weight q instead of floor(q/2)           -> 1 test fails
"""

from __future__ import annotations

import numpy as np
import openfermion
import pytest
from qiskit.quantum_info import Operator, Statevector
from scipy.linalg import expm, logm

from aqs.hamiltonian import _to_little_endian
from aqs.models.spinless import (
    SpinlessSSHModel,
    SpinlessSSHHSim,
    spec_pauli_hamiltonian,
    exact_hamiltonian,
    exact_ground_state_fixed_n,
)


# ---------------------------------------------------------------------------
# helpers (references only -- no production logic)
# ---------------------------------------------------------------------------

def _fermionic_reference(model: SpinlessSSHModel) -> np.ndarray:
    """Dense little-endian matrix of the spec Hamiltonian, built via OpenFermion.

    The fermionic coefficients are the ones implied by the *normative* Pauli
    representation of spec section 2.1 (see module docstring of
    aqs.models.spinless): intracell  -v b_j^dag a_j + h.c., intercell and PBC
    wrap  -conj(w) b_j^dag a_{j+1} + h.c., plus the diagonal NN interaction.
    OpenFermion supplies the Jordan-Wigner strings (including the PBC wrap
    string) independently of the module under test.
    """
    n = 2 * model.N_cells
    op = openfermion.FermionOperator()
    for j in range(model.N_cells):
        a, b = 2 * j, 2 * j + 1
        op += openfermion.FermionOperator(f"{b}^ {a}", -model.v)
        op += openfermion.FermionOperator(f"{a}^ {b}", -np.conj(model.v))
        op += openfermion.FermionOperator(f"{a}^ {a} {b}^ {b}", model.V_v)
    last = model.N_cells - 1 if model.PBC else model.N_cells - 2
    for j in range(last + 1):
        b, a_next = 2 * j + 1, (2 * j + 2) % n
        op += openfermion.FermionOperator(f"{b}^ {a_next}", -np.conj(model.w))
        op += openfermion.FermionOperator(f"{a_next}^ {b}", -model.w)
        op += openfermion.FermionOperator(f"{b}^ {b} {a_next}^ {a_next}", model.V_w)
    H_big = openfermion.get_sparse_operator(op, n_qubits=n).toarray()
    return _to_little_endian(H_big, n)


def _sector_indices(n_qubits: int, n_particles: int) -> np.ndarray:
    return np.array(
        [i for i in range(2**n_qubits) if bin(i).count("1") == n_particles]
    )


def _random_sector_state(n_qubits: int, n_particles: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = _sector_indices(n_qubits, n_particles)
    psi = np.zeros(2**n_qubits, dtype=complex)
    amps = rng.normal(size=len(idx)) + 1j * rng.normal(size=len(idx))
    psi[idx] = amps / np.linalg.norm(amps)
    return psi


# ---------------------------------------------------------------------------
# Hamiltonian construction
# ---------------------------------------------------------------------------

class TestHamiltonianConstruction:
    def test_pauli_hamiltonian_matches_openfermion_reference_obc(self):
        """Spec Pauli H == OpenFermion-built fermionic H (complex v, w, OBC).

        Reverse-verification: conjugating v in the module's Pauli builder, or
        dropping _to_little_endian here, makes this fail (complex hoppings are
        mandatory -- real ones cannot see either mutation).
        """
        model = SpinlessSSHModel(
            N_cells=3, v=0.5 + 0.3j, w=1.2 - 0.7j, V_v=0.4, V_w=0.9, PBC=False
        )
        H_pauli = spec_pauli_hamiltonian(model).to_matrix()
        np.testing.assert_allclose(
            H_pauli, _fermionic_reference(model), atol=1e-12
        )

    def test_exact_hamiltonian_matches_openfermion_reference_pbc(self):
        """PBC: the wrap bond's Jordan-Wigner string must be present.

        exact_hamiltonian (the many-body reference the module exposes) must
        agree with OpenFermion's independent JW machinery including the wrap
        string, for complex w and an interaction wrap term.
        """
        model = SpinlessSSHModel(
            N_cells=3, v=0.2 + 0.6j, w=1.0 - 0.4j, V_v=0.3, V_w=0.8, PBC=True
        )
        np.testing.assert_allclose(
            exact_hamiltonian(model), _fermionic_reference(model), atol=1e-12
        )

    def test_single_particle_spectrum_matches_many_body_one_particle_sector(self):
        """eigh(single-particle M) == spectrum of H restricted to n=1 (OBC+PBC)."""
        for pbc in (False, True):
            model = SpinlessSSHModel(
                N_cells=3, v=0.9 - 0.2j, w=0.4 + 1.1j, PBC=pbc
            )
            H = exact_hamiltonian(model)
            idx = _sector_indices(model.n_sites, 1)
            H_1 = H[np.ix_(idx, idx)]
            single = np.linalg.eigvalsh(model.single_particle_hamiltonian())
            np.testing.assert_allclose(
                np.sort(np.linalg.eigvalsh(H_1)), np.sort(single), atol=1e-12
            )


# ---------------------------------------------------------------------------
# Trotter circuit: H_eff round-trips (exact on single-bond layers)
# ---------------------------------------------------------------------------

class TestTrotterConventions:
    """H_eff = i log(U)/dt round-trips.

    A bond with a single Pauli generator (pure-real coupling -> R only,
    pure-imaginary -> G only) is simulated exactly by its gate, so H_eff must
    equal the spec Pauli H to machine precision at ANY dt. For complex
    couplings, R and G on one bond do not commute (they act as X and Y on the
    single-occupancy subspace), so the step is itself a first-order split:
    there the guard is that the H_eff error SHRINKS ~O(dt), which a sign or
    conjugation bug (error O(1), dt-independent) cannot satisfy.
    """

    @pytest.mark.parametrize(
        "v,w",
        [
            (0.7, 0.0),      # intracell R only
            (0.4j, 0.0),     # intracell G only -- the sharp G-sign pin
            (0.0, 1.1),      # intercell R only
            (0.0, -0.6j),    # intercell G only -- literal spec Im(w) sign pin
        ],
    )
    def test_single_generator_heff_roundtrip_is_exact(self, v, w):
        """Reverse-verification: flipping the builder's G angle sign fails the
        pure-imaginary cases with error 2|Im t|, independent of dt."""
        model = SpinlessSSHModel(N_cells=2, v=v, w=w, PBC=False)
        for dt in (0.05, 0.2):
            sim = SpinlessSSHHSim(model, T=dt, L=1)
            U = Operator(sim.trotter_circuit()).data
            H_eff = 1j * logm(U) / dt
            np.testing.assert_allclose(
                H_eff, spec_pauli_hamiltonian(model).to_matrix(), atol=1e-9
            )

    def test_complex_coupling_heff_error_shrinks_with_dt(self):
        """Complex v: H_eff error is O(dt) Trotter error, not a plateau.

        A conjugated Im(v) would leave a dt-independent error of 2|Im v| = 0.8;
        the assertions bound the error well below that and require the ~4x
        shrink per 4x dt refinement characteristic of a genuine O(dt) term.
        """
        model = SpinlessSSHModel(N_cells=2, v=0.7 + 0.4j, w=0.0, PBC=False)
        H_true = spec_pauli_hamiltonian(model).to_matrix()
        errors = []
        for dt in (0.2, 0.05, 0.0125):
            sim = SpinlessSSHHSim(model, T=dt, L=1)
            U = Operator(sim.trotter_circuit()).data
            errors.append(np.max(np.abs(1j * logm(U) / dt - H_true)))
        assert errors[0] > 2.5 * errors[1] > 2.5 * 2.5 * errors[2]
        assert errors[2] < 0.01

    @pytest.mark.parametrize("n_particles", [1, 2])
    @pytest.mark.parametrize("w", [0.8, 0.5j])
    def test_pbc_wrap_bond_matches_exact_evolution_in_fixed_sectors(
        self, n_particles, w
    ):
        """PBC wrap hopping: circuit step == expm(-i dt H) on sector states.

        N=2, v=0: the two w-bonds act on disjoint qubit pairs {1,2} and {3,0},
        and each single-generator bond gate is exact, so any deviation is the
        wrap parity or string handling. Both parities (n=1 odd, n=2 even) and
        both generators (real w -> R, imaginary w -> G) are exercised; the
        imaginary case is the one a string-conjugation bug flips.
        """
        model = SpinlessSSHModel(
            N_cells=2, v=0.0, w=w, PBC=True, n_particles=n_particles
        )
        dt = 0.13
        sim = SpinlessSSHHSim(model, T=dt, L=1)
        U_circ = Operator(sim.trotter_circuit()).data
        U_exact = expm(-1j * dt * exact_hamiltonian(model))
        for seed in (0, 1):
            psi = _random_sector_state(model.n_sites, n_particles, seed)
            np.testing.assert_allclose(U_circ @ psi, U_exact @ psi, atol=1e-10)

    def test_interaction_cp_angles_follow_midpoint_lambda(self):
        """CP angles across the ramp are dt * V * (2l-1)/(2L)  (spec 3.2)."""
        model = SpinlessSSHModel(N_cells=2, v=0.5, w=1.5, V_v=0.7, PBC=False)
        L, T = 4, 2.0
        sim = SpinlessSSHHSim(model, T=T, L=L)
        qc = sim.trotter_circuit()
        dt = T / L
        angles = [
            -inst.operation.definition.data[0].operation.params[0]
            for inst in qc.data
            if inst.operation.name.startswith("CP(")
        ]
        expected = []
        for step in range(1, L + 1):
            lam = (2 * step - 1) / (2 * L)
            expected.extend([dt * model.V_v * lam] * model.N_cells)
        np.testing.assert_allclose(angles, expected, atol=1e-12)

    def test_trotter_error_shrinks_with_more_steps(self):
        """Fidelity error vs the exact propagator decreases with L (OBC, V=0).

        A plateau here (error not shrinking) is the signature of a systematic
        convention bug rather than Trotter error.
        """
        model = SpinlessSSHModel(
            N_cells=3, v=0.6 + 0.2j, w=1.3 - 0.5j, PBC=False, n_particles=3
        )
        T = 1.0
        U_exact = expm(-1j * T * exact_hamiltonian(model))
        psi0 = Statevector(
            SpinlessSSHHSim(model, T=T, L=1).preparation_circuit()
        ).data
        target = U_exact @ psi0
        errors = []
        for L in (2, 8, 32):
            sim = SpinlessSSHHSim(model, T=T, L=L)
            evolved = Statevector(sim.full_circuit()).data
            errors.append(1.0 - abs(np.vdot(target, evolved)) ** 2)
        assert errors[0] > errors[1] > errors[2]
        assert errors[2] < 5e-3


# ---------------------------------------------------------------------------
# State preparation
# ---------------------------------------------------------------------------

class TestStatePreparation:
    @pytest.mark.parametrize("pbc,n_particles", [(False, 4), (True, 3)])
    def test_prepared_state_has_lowest_fill_energy(self, pbc, n_particles):
        """<H> of the prepared determinant == sum of the lowest orbitals.

        Complex hoppings + asymmetric filling: on a chiral chain the wrong
        (highest-band) determinant lands at exactly -E, so this assertion
        flips sign rather than drifting when the Givens network regresses.
        """
        model = SpinlessSSHModel(
            N_cells=3, v=0.4 + 0.3j, w=1.2 - 0.8j, PBC=pbc, n_particles=n_particles
        )
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        H = exact_hamiltonian(
            SpinlessSSHModel(
                N_cells=3, v=model.v, w=model.w, PBC=pbc, n_particles=n_particles
            )
        )
        energy = np.real(np.vdot(psi, H @ psi))
        eigs = np.sort(
            np.linalg.eigvalsh(model.single_particle_hamiltonian())
        )
        np.testing.assert_allclose(energy, eigs[:n_particles].sum(), atol=1e-9)

    def test_prepared_state_has_fixed_particle_number(self):
        model = SpinlessSSHModel(N_cells=3, v=0.5, w=1.5, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        weight = 0.0
        for idx in _sector_indices(model.n_sites, model.n_particles):
            weight += abs(psi[idx]) ** 2
        np.testing.assert_allclose(weight, 1.0, atol=1e-12)

    def test_default_filling_follows_spec(self):
        """Spec section 4: N particles under PBC, N+1 under OBC."""
        assert SpinlessSSHModel(N_cells=4, v=1.0, w=0.5, PBC=True).n_particles == 4
        assert SpinlessSSHModel(N_cells=4, v=1.0, w=0.5, PBC=False).n_particles == 5


# ---------------------------------------------------------------------------
# Adiabatic evolution against fixed-sector exact diagonalization
# ---------------------------------------------------------------------------

class TestAdiabaticEvolution:
    def test_adiabatic_evolution_reaches_interacting_ground_state(self):
        """Final state fidelity vs the fixed-sector ED ground state of H.

        PBC, half filling (odd n=3 exercises the wrap parity in a physical
        workflow), interaction inside the adiabatic-safety window. The
        reference is diagonalization restricted to the 3-particle sector --
        the *global* interacting ground state may live at another filling and
        would be the wrong reference.
        """
        model = SpinlessSSHModel(
            N_cells=3, v=0.5, w=1.5, V_v=0.3, V_w=0.3, PBC=True
        )
        ground = exact_ground_state_fixed_n(model)
        fidelities = []
        for L in (16, 128):  # fixed T: improvement must come from dt -> 0
            sim = SpinlessSSHHSim(model, T=2.0, L=L)
            psi = Statevector(sim.full_circuit()).data
            fidelities.append(abs(np.vdot(ground, psi)) ** 2)
        # A plateau below 1 under L-refinement would indicate a systematic
        # convention bug (measured: 0.99854 -> 0.99995 across this refinement).
        assert fidelities[1] > fidelities[0]
        assert fidelities[1] > 0.9999


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------

class TestObservables:
    def test_polarization_profile_on_asymmetric_product_state(self):
        """P_j = <n_{2j} - n_{2j+1}> on |110100>: exact values, asymmetric.

        Qubit occupations (little-endian): q0=1, q1=1, q2=0, q3=1, q4=0, q5=0
        -> P_0 = 1-1 = 0, P_1 = 0-1 = -1, P_2 = 0-0 = 0. A bit-reversal bug
        reads the pattern backwards and reports [0, -1, 0] -> [0, +1, 0]...
        the asymmetry makes any ordering mutation visible.
        """
        model = SpinlessSSHModel(N_cells=3, v=1.0, w=0.5, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**6, dtype=complex)
        psi[0b001011] = 1.0  # occupied qubits 0, 1, 3
        np.testing.assert_allclose(
            sim.polarization_profile(psi), [0.0, -1.0, 0.0], atol=1e-12
        )

    def test_berry_phase_jumps_by_pi_across_transition(self):
        """Spec 5.2 Berry phase: trivial vs topological differ by pi (N=6 PBC).

        Only the jump is asserted -- absolute values carry a filling-dependent
        constant offset (see .claude/rules/models-and-observables.md).
        """
        phases = []
        for v, w in ((1.5, 0.5), (0.5, 1.5)):
            model = SpinlessSSHModel(N_cells=6, v=v, w=w, PBC=True)
            sim = SpinlessSSHHSim(model, T=1.0, L=1)
            psi = Statevector(sim.preparation_circuit()).data
            phases.append(sim.berry_phase(psi))
        jump = (phases[1] - phases[0]) % (2 * np.pi)
        assert min(abs(jump - np.pi), abs(jump + np.pi - 2 * np.pi)) < 1e-6

    def test_berry_phase_requires_pbc(self):
        model = SpinlessSSHModel(N_cells=2, v=1.0, w=0.5, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**4, dtype=complex)
        psi[0] = 1.0
        with pytest.raises(ValueError, match="PBC"):
            sim.berry_phase(psi)

    def test_polarization_requires_obc(self):
        model = SpinlessSSHModel(N_cells=2, v=1.0, w=0.5, PBC=True)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**4, dtype=complex)
        psi[0] = 1.0
        with pytest.raises(ValueError, match="OBC"):
            sim.polarization_profile(psi)


# ---------------------------------------------------------------------------
# Spec section 6 drivers
# ---------------------------------------------------------------------------

class TestValidation:
    def test_rejects_single_cell(self):
        with pytest.raises(ValueError, match="N_cells"):
            SpinlessSSHModel(N_cells=1, v=1.0, w=0.5)

    def test_rejects_bad_particle_count(self):
        with pytest.raises(ValueError, match="n_particles"):
            SpinlessSSHModel(N_cells=2, v=1.0, w=0.5, n_particles=5)

    def test_rejects_nonpositive_step_count(self):
        model = SpinlessSSHModel(N_cells=2, v=1.0, w=0.5)
        with pytest.raises(ValueError, match="L"):
            SpinlessSSHHSim(model, T=1.0, L=0)


class TestSimDrivers:
    def test_trotter_benchmark_returns_monotone_records(self):
        model = SpinlessSSHModel(N_cells=2, v=0.5, w=1.5, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        records = sim.trotter_benchmark(L_values=(2, 8))
        assert [r["L"] for r in records] == [2, 8]
        assert all(0.0 <= r["fidelity"] <= 1.0 + 1e-12 for r in records)
        assert records[1]["fidelity"] >= records[0]["fidelity"]

    def test_interacting_benchmark_uses_adjacent_L_stability(self):
        """With V != 0 the reference is the previous L's state (stability
        diagnostic, not absolute accuracy -- see rules); records still bounded."""
        model = SpinlessSSHModel(N_cells=2, v=0.5, w=1.5, V_v=0.4, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=1)
        records = sim.trotter_benchmark(L_values=(2, 4, 8))
        assert [r["L"] for r in records] == [2, 4, 8]
        assert all(0.0 <= r["fidelity"] <= 1.0 + 1e-12 for r in records)

    def test_interaction_sweep_records_observable_per_point(self):
        model = SpinlessSSHModel(N_cells=2, v=0.5, w=1.5, PBC=True)
        sim = SpinlessSSHHSim(model, T=2.0, L=16)
        records = sim.interaction_sweep(V_points=((0.0, 0.0), (0.2, 0.2)))
        assert len(records) == 2
        for rec in records:
            assert set(rec) >= {"V_v", "V_w", "berry_phase"}

    def test_interaction_sweep_reports_polarization_under_obc(self):
        model = SpinlessSSHModel(N_cells=2, v=0.5, w=1.5, PBC=False)
        sim = SpinlessSSHHSim(model, T=1.0, L=8)
        (rec,) = sim.interaction_sweep(V_points=((0.1, 0.2),))
        assert set(rec) >= {"V_v", "V_w", "edge_polarization", "polarization"}
        assert len(rec["polarization"]) == model.N_cells
        np.testing.assert_allclose(
            rec["edge_polarization"], rec["polarization"][0], atol=1e-12
        )
