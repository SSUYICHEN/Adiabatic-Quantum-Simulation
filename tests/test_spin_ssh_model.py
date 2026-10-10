"""Guard tests for the SSH-Hubbard chain with spin as the sublattice.

Authority: ``.claude/specs/reconstruct_goals.md`` (revision 2026-10-06): the
fermionic Hamiltonian of spec section 2.1 is normative, and section 2.2 is its
Jordan-Wigner image. Modes: q_{2j} = (j, up), q_{2j+1} = (j, down).

Every physics assertion compares against an independent reference:

* the second-quantized Hamiltonian of spec 2.1 built with OpenFermion and
  converted to little-endian (never the module's own Pauli construction alone),
* exact diagonalization restricted to the fixed particle-number sector,
* H_eff = i log(U) / dt round-trips on single-bond circuits, which are exact
  (no Trotter error) and therefore catch sign/conjugation convention bugs that
  do not shrink with dt.

Inputs are deliberately asymmetric and use complex hoppings where a symmetric
or real input would be structurally blind to the bug class under test.

Reverse verification (2026-10-06, after the rename and sign-convention change):
each mutation below was applied to src/aqs/models/spin_ssh.py and the suite
confirmed to fail, then reverted:

* flip the bulk G angle sign (theta_G -> +dt Im t)        -> 6 tests fail
* flip the bulk R angle sign (theta_R -> -dt Re t)        -> 8 tests fail
* drop the wrap parity factor (parity_factor -> 1)        -> 3 tests fail
* wrap G sign follows bulk (theta_G *= +parity_factor)    -> 2 tests fail
* conjugate the bulk Im sign in the Pauli builder         -> 11 tests fail
* fill the determinant from the highest orbitals          -> 4 tests fail
* swap up/down in the magnetisation profile               -> 1 test fails
* Berry position weight q instead of floor(q/2)           -> 2 tests fail
(counts over test_spin_ssh_model.py + test_spin_ssh_physics_smoke.py)
"""

from __future__ import annotations

import numpy as np
import openfermion
import pytest
from qiskit.quantum_info import Operator, Statevector
from scipy.linalg import expm, logm

from aqs.mapping import _to_little_endian
from aqs.models.spin_ssh import (
    SpinSSHModel,
    SpinSSHHSim,
    pauli_hamiltonian,
    exact_hamiltonian,
    exact_ground_state_fixed_n,
)


# ---------------------------------------------------------------------------
# helpers (references only -- no production logic)
# ---------------------------------------------------------------------------

def _fermionic_reference(model: SpinSSHModel) -> np.ndarray:
    """Dense little-endian matrix of the spec-2.1 Hamiltonian, built via OpenFermion.

    v c_{j up}^dag c_{j dn} + h.c.,  w c_{j dn}^dag c_{j+1 up} + h.c.  (PBC: the
    wrap term creates on the LAST mode and annihilates on mode 0),  U n_up n_dn.
    OpenFermion supplies the Jordan-Wigner strings (including the PBC wrap
    string) independently of the module under test.
    """
    n = 2 * model.N
    op = openfermion.FermionOperator()
    for j in range(model.N):
        up, dn = 2 * j, 2 * j + 1
        op += openfermion.FermionOperator(f"{up}^ {dn}", model.v)
        op += openfermion.FermionOperator(f"{dn}^ {up}", np.conj(model.v))
        op += openfermion.FermionOperator(f"{up}^ {up} {dn}^ {dn}", model.U)
    for j in range(model.N - 1):
        dn, up_next = 2 * j + 1, 2 * j + 2
        op += openfermion.FermionOperator(f"{dn}^ {up_next}", model.w)
        op += openfermion.FermionOperator(f"{up_next}^ {dn}", np.conj(model.w))
    if model.PBC:
        op += openfermion.FermionOperator(f"{n - 1}^ 0", model.w)
        op += openfermion.FermionOperator(f"0^ {n - 1}", np.conj(model.w))
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
        """Module Pauli H == OpenFermion-built spec-2.1 H (complex v, w, OBC).

        Reverse-verification: conjugating v in the module's Pauli builder, or
        dropping _to_little_endian here, makes this fail (complex hoppings are
        mandatory -- real ones cannot see either mutation).
        """
        model = SpinSSHModel(
            N=3, v=0.5 + 0.3j, w=1.2 - 0.7j, U=0.4, PBC=False
        )
        H_pauli = pauli_hamiltonian(model).to_matrix()
        np.testing.assert_allclose(
            H_pauli, _fermionic_reference(model), atol=1e-12
        )

    def test_exact_hamiltonian_matches_openfermion_reference_pbc(self):
        """PBC: the wrap bond's Jordan-Wigner string must be present.

        exact_hamiltonian (the many-body reference the module exposes) must
        agree with OpenFermion's independent JW machinery including the wrap
        string, for complex w.
        """
        model = SpinSSHModel(
            N=3, v=0.2 + 0.6j, w=1.0 - 0.4j, U=0.3, PBC=True
        )
        np.testing.assert_allclose(
            exact_hamiltonian(model), _fermionic_reference(model), atol=1e-12
        )

    def test_single_particle_spectrum_matches_many_body_one_particle_sector(self):
        """eigh(single-particle M) == spectrum of H restricted to n=1 (OBC+PBC)."""
        for pbc in (False, True):
            model = SpinSSHModel(
                N=3, v=0.9 - 0.2j, w=0.4 + 1.1j, PBC=pbc
            )
            H = exact_hamiltonian(model)
            idx = _sector_indices(model.n_modes, 1)
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
            (0.7, 0.0),      # spin flip, R only
            (0.4j, 0.0),     # spin flip, G only -- the sharp G-sign pin
            (0.0, 1.1),      # intersite R only
            (0.0, -0.6j),    # intersite G only -- spec Im(w) sign pin
        ],
    )
    def test_single_generator_heff_roundtrip_is_exact(self, v, w):
        """Reverse-verification: flipping the builder's G angle sign fails the
        pure-imaginary cases with error 2|Im t|, independent of dt."""
        model = SpinSSHModel(N=2, v=v, w=w, PBC=False)
        for dt in (0.05, 0.2):
            sim = SpinSSHHSim(model, T=dt, L=1)
            U = Operator(sim.trotter_circuit()).data
            H_eff = 1j * logm(U) / dt
            np.testing.assert_allclose(
                H_eff, pauli_hamiltonian(model).to_matrix(), atol=1e-9
            )

    def test_complex_coupling_heff_error_shrinks_with_dt(self):
        """Complex v: H_eff error is O(dt) Trotter error, not a plateau.

        A conjugated Im(v) would leave a dt-independent error of 2|Im v| = 0.8;
        the assertions bound the error well below that and require the ~4x
        shrink per 4x dt refinement characteristic of a genuine O(dt) term.
        """
        model = SpinSSHModel(N=2, v=0.7 + 0.4j, w=0.0, PBC=False)
        H_true = pauli_hamiltonian(model).to_matrix()
        errors = []
        for dt in (0.2, 0.05, 0.0125):
            sim = SpinSSHHSim(model, T=dt, L=1)
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
        model = SpinSSHModel(
            N=2, v=0.0, w=w, PBC=True, n_particles=n_particles
        )
        dt = 0.13
        sim = SpinSSHHSim(model, T=dt, L=1)
        U_circ = Operator(sim.trotter_circuit()).data
        U_exact = expm(-1j * dt * exact_hamiltonian(model))
        for seed in (0, 1):
            psi = _random_sector_state(model.n_modes, n_particles, seed)
            np.testing.assert_allclose(U_circ @ psi, U_exact @ psi, atol=1e-10)

    def test_interaction_cp_angles_follow_midpoint_lambda(self):
        """CP angles across the ramp are dt * U * (2l-1)/(2L)  (spec 3.2)."""
        model = SpinSSHModel(N=2, v=0.5, w=1.5, U=0.7, PBC=False)
        L, T = 4, 2.0
        sim = SpinSSHHSim(model, T=T, L=L)
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
            expected.extend([dt * model.U * lam] * model.N)
        np.testing.assert_allclose(angles, expected, atol=1e-12)

    def test_trotter_error_shrinks_with_more_steps(self):
        """Fidelity error vs the exact propagator decreases with L (OBC, U=0).

        A plateau here (error not shrinking) is the signature of a systematic
        convention bug rather than Trotter error.
        """
        model = SpinSSHModel(
            N=3, v=0.6 + 0.2j, w=1.3 - 0.5j, PBC=False, n_particles=3
        )
        T = 1.0
        U_exact = expm(-1j * T * exact_hamiltonian(model))
        psi0 = Statevector(
            SpinSSHHSim(model, T=T, L=1).preparation_circuit()
        ).data
        target = U_exact @ psi0
        errors = []
        for L in (2, 8, 32):
            sim = SpinSSHHSim(model, T=T, L=L)
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
        model = SpinSSHModel(
            N=3, v=0.4 + 0.3j, w=1.2 - 0.8j, PBC=pbc, n_particles=n_particles
        )
        sim = SpinSSHHSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        H = exact_hamiltonian(
            SpinSSHModel(
                N=3, v=model.v, w=model.w, PBC=pbc, n_particles=n_particles
            )
        )
        energy = np.real(np.vdot(psi, H @ psi))
        eigs = np.sort(
            np.linalg.eigvalsh(model.single_particle_hamiltonian())
        )
        np.testing.assert_allclose(energy, eigs[:n_particles].sum(), atol=1e-9)

    def test_prepared_state_has_fixed_particle_number(self):
        model = SpinSSHModel(N=3, v=0.5, w=1.5, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        psi = Statevector(sim.preparation_circuit()).data
        weight = 0.0
        for idx in _sector_indices(model.n_modes, model.n_particles):
            weight += abs(psi[idx]) ** 2
        np.testing.assert_allclose(weight, 1.0, atol=1e-12)

    def test_default_filling_follows_spec(self):
        """Spec section 4: N particles under PBC, N+1 under OBC."""
        assert SpinSSHModel(N=4, v=1.0, w=0.5, PBC=True).n_particles == 4
        assert SpinSSHModel(N=4, v=1.0, w=0.5, PBC=False).n_particles == 5


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
        model = SpinSSHModel(
            N=3, v=0.5, w=1.5, U=0.3, PBC=True
        )
        ground = exact_ground_state_fixed_n(model)
        fidelities = []
        for L in (16, 128):  # fixed T: improvement must come from dt -> 0
            sim = SpinSSHHSim(model, T=2.0, L=L)
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
    def test_magnetisation_profile_on_asymmetric_product_state(self):
        """m_j = <n_{j up} - n_{j dn}> on |110100>: exact values, asymmetric.

        Qubit occupations (little-endian): q0=1, q1=1, q2=0, q3=1, q4=0, q5=0
        -> m_0 = 1-1 = 0, m_1 = 0-1 = -1, m_2 = 0-0 = 0. A bit-reversal bug
        reads the pattern backwards and reports [0, -1, 0] -> [0, +1, 0]...
        the asymmetry makes any ordering mutation visible.
        """
        model = SpinSSHModel(N=3, v=1.0, w=0.5, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**6, dtype=complex)
        psi[0b001011] = 1.0  # occupied qubits 0, 1, 3
        np.testing.assert_allclose(
            sim.magnetisation_profile(psi), [0.0, -1.0, 0.0], atol=1e-12
        )

    def test_berry_phase_jumps_by_pi_across_transition(self):
        """Spec 5.2 Berry phase: trivial vs topological differ by pi (N=6 PBC).

        Only the jump is asserted -- absolute values carry a filling-dependent
        constant offset (see .claude/rules/models-and-observables.md).
        """
        phases = []
        for v, w in ((1.5, 0.5), (0.5, 1.5)):
            model = SpinSSHModel(N=6, v=v, w=w, PBC=True)
            sim = SpinSSHHSim(model, T=1.0, L=1)
            psi = Statevector(sim.preparation_circuit()).data
            phases.append(sim.berry_phase(psi))
        jump = (phases[1] - phases[0]) % (2 * np.pi)
        assert min(abs(jump - np.pi), abs(jump + np.pi - 2 * np.pi)) < 1e-6

    def test_berry_phase_requires_pbc(self):
        model = SpinSSHModel(N=2, v=1.0, w=0.5, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**4, dtype=complex)
        psi[0] = 1.0
        with pytest.raises(ValueError, match="PBC"):
            sim.berry_phase(psi)

    def test_magnetisation_requires_obc(self):
        model = SpinSSHModel(N=2, v=1.0, w=0.5, PBC=True)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        psi = np.zeros(2**4, dtype=complex)
        psi[0] = 1.0
        with pytest.raises(ValueError, match="OBC"):
            sim.magnetisation_profile(psi)


# ---------------------------------------------------------------------------
# Spec section 6 drivers
# ---------------------------------------------------------------------------

class TestValidation:
    def test_rejects_single_cell(self):
        with pytest.raises(ValueError, match="N"):
            SpinSSHModel(N=1, v=1.0, w=0.5)

    def test_rejects_bad_particle_count(self):
        with pytest.raises(ValueError, match="n_particles"):
            SpinSSHModel(N=2, v=1.0, w=0.5, n_particles=5)

    def test_rejects_nonpositive_step_count(self):
        model = SpinSSHModel(N=2, v=1.0, w=0.5)
        with pytest.raises(ValueError, match="L"):
            SpinSSHHSim(model, T=1.0, L=0)


class TestSimDrivers:
    def test_trotter_benchmark_returns_monotone_records(self):
        model = SpinSSHModel(N=2, v=0.5, w=1.5, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        records = sim.trotter_benchmark(L_values=(2, 8))
        assert [r["L"] for r in records] == [2, 8]
        assert all(0.0 <= r["fidelity"] <= 1.0 + 1e-12 for r in records)
        assert records[1]["fidelity"] >= records[0]["fidelity"]

    def test_interacting_benchmark_uses_adjacent_L_stability(self):
        """With U != 0 the reference is the previous L's state (stability
        diagnostic, not absolute accuracy -- see rules); records still bounded."""
        model = SpinSSHModel(N=2, v=0.5, w=1.5, U=0.4, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=1)
        records = sim.trotter_benchmark(L_values=(2, 4, 8))
        assert [r["L"] for r in records] == [2, 4, 8]
        assert all(0.0 <= r["fidelity"] <= 1.0 + 1e-12 for r in records)

    def test_interaction_sweep_records_observable_per_point(self):
        model = SpinSSHModel(N=2, v=0.5, w=1.5, PBC=True)
        sim = SpinSSHHSim(model, T=2.0, L=16)
        records = sim.interaction_sweep(U_values=(0.0, 0.2))
        assert len(records) == 2
        for rec in records:
            assert set(rec) >= {"U", "berry_phase"}

    def test_interaction_sweep_reports_magnetisation_under_obc(self):
        model = SpinSSHModel(N=2, v=0.5, w=1.5, PBC=False)
        sim = SpinSSHHSim(model, T=1.0, L=8)
        (rec,) = sim.interaction_sweep(U_values=(0.1,))
        assert set(rec) >= {"U", "edge_magnetisation", "magnetisation"}
        assert len(rec["magnetisation"]) == model.N
        np.testing.assert_allclose(
            rec["edge_magnetisation"], rec["magnetisation"][0], atol=1e-12
        )
