"""Tests for observables.py: layouts, densities, twist invariant, registry.

The twist invariant is checked against closed-form limits rather than recorded
numbers, so the tests state what the physics should be instead of freezing
whatever the code currently emits.
"""
import numpy as np
import pytest

from aqs.observables import (PROPERTIES, density_profile, density_profile_qubits,
                             spin_block_layout, twist_invariant,
                             twist_invariant_layout, wrap_berry_phase)


def _basis_state(n_qubits, occupied):
    """Little-endian computational basis state with `occupied` qubits set."""
    idx = 0
    for q in occupied:
        idx |= 1 << q
    sv = np.zeros(1 << n_qubits, dtype=complex)
    sv[idx] = 1.0
    return sv


# --------------------------------------------------------------- layouts
def test_spin_block_layout_up_spin():
    lay = spin_block_layout(3, mode="up_spin")
    assert lay["n_qubits"] == 12 and lay["n_cells"] == 3 and lay["L_sites"] == 6
    assert lay["cell_qubits"] == [[0, 1], [2, 3], [4, 5]]
    assert lay["A_qubits"] == [0, 2, 4] and lay["B_qubits"] == [1, 3, 5]


def test_spin_block_layout_total_includes_the_down_block():
    lay = spin_block_layout(3, mode="total")
    assert lay["cell_qubits"] == [[0, 1, 6, 7], [2, 3, 8, 9], [4, 5, 10, 11]]
    assert lay["A_qubits"] == [0, 6, 2, 8, 4, 10]
    assert lay["B_qubits"] == [1, 7, 3, 9, 5, 11]


# --------------------------------------------------------------- densities
def test_density_profile_qubits_reads_little_endian():
    n = 4
    for q in range(n):
        d = density_profile_qubits(_basis_state(n, [q]), n)
        assert d[q] == pytest.approx(1.0), q
        assert sum(d) == pytest.approx(1.0)


def test_density_profile_sums_both_spin_blocks():
    L = 3                       # 3 spatial sites -> 6 qubits
    sv = _basis_state(2 * L, [1, L + 1])         # site 1 doubly occupied
    d = density_profile(sv, L)
    assert d[1] == pytest.approx(2.0)
    assert d[0] == pytest.approx(0.0) and d[2] == pytest.approx(0.0)


def test_density_of_a_superposition_is_the_weighted_average():
    n = 2
    sv = (_basis_state(n, [0]) + _basis_state(n, [1])) / np.sqrt(2)
    d = density_profile_qubits(sv, n)
    assert d == pytest.approx([0.5, 0.5])


def test_empty_state_returns_zero_density():
    assert density_profile_qubits(np.zeros(8, dtype=complex), 3) == pytest.approx([0, 0, 0])
    assert density_profile(np.zeros(64, dtype=complex), 3) == pytest.approx([0, 0, 0])


# --------------------------------------------------------- twist invariant
def test_twist_invariant_of_the_vacuum_is_one():
    lay = spin_block_layout(3, mode="up_spin")
    z = twist_invariant_layout(_basis_state(12, []), lay)
    assert z == pytest.approx(1.0 + 0j)


def test_twist_invariant_of_a_number_eigenstate_is_a_pure_phase():
    """X = sum_j (j+1) n_j is diagonal, so a basis state gives |z| = 1 exactly."""
    lay = spin_block_layout(3, mode="up_spin")
    N = 3
    for occ in ([0], [2], [0, 2, 4], [1, 3]):
        z = twist_invariant_layout(_basis_state(12, occ), lay)
        assert abs(z) == pytest.approx(1.0)
        X = sum(q // 2 + 1 for q in occ)
        assert z == pytest.approx(np.exp(1j * 2 * np.pi / N * X))


def test_twist_invariant_zero_for_a_null_vector():
    lay = spin_block_layout(3, mode="up_spin")
    assert twist_invariant_layout(np.zeros(4096, dtype=complex), lay) == 0


def test_twist_invariant_wrapper_matches_the_layout_version():
    sv = _basis_state(12, [0, 3, 4])
    assert twist_invariant(sv, 6, mode="up_spin") == pytest.approx(
        twist_invariant_layout(sv, spin_block_layout(3, "up_spin")))


def test_uniform_one_particle_per_cell_gives_the_parity_offset():
    """Fully dimerised limit: every cell holds exactly one particle, so
    z_N = exp(i pi (N+1)) -- the documented filling/parity offset."""
    for N in (3, 4, 5, 6):
        lay = spin_block_layout(N, mode="up_spin")
        occ = [2 * j for j in range(N)]            # one particle per cell
        z = twist_invariant_layout(_basis_state(4 * N, occ), lay)
        assert z == pytest.approx(np.exp(1j * np.pi * (N + 1)), abs=1e-9), N


# ----------------------------------------------------------- berry wrapping
def test_wrap_berry_phase_range_and_values():
    assert wrap_berry_phase(1.0 + 0j) == pytest.approx(0.0)
    assert abs(wrap_berry_phase(-1.0 + 0j)) == pytest.approx(1.0)
    # np.mod gives [0, 2pi), so mathematically the range is [-1/2, 3/2); allow a
    # tolerance because rounding can land exactly on either boundary
    for ang in np.linspace(-np.pi, np.pi, 401):
        g = wrap_berry_phase(np.exp(1j * ang))
        assert -0.5 - 1e-12 <= g <= 1.5 + 1e-12, (ang, g)


def test_wrap_berry_phase_is_periodic():
    """Adding 2*pi to the phase must not change the wrapped result."""
    for ang in (0.3, -1.1, 2.7):
        a = wrap_berry_phase(np.exp(1j * ang))
        b = wrap_berry_phase(np.exp(1j * (ang + 2 * np.pi)))
        assert a == pytest.approx(b, abs=1e-9), ang


# ------------------------------------------------------------ the registry
def test_registry_exposes_the_documented_properties():
    assert set(PROPERTIES) == {"berry", "polarization"}


def test_berry_property_fields_are_self_consistent():
    lay = spin_block_layout(3, mode="up_spin")
    res = PROPERTIES["berry"](_basis_state(12, [0, 2, 4]), lay)
    z = complex(res["Z_N_real"], res["Z_N_imag"])
    assert res["Twist_Amplitude"] == pytest.approx(abs(z))
    assert res["Berry_Phase_pi_wrapped"] == pytest.approx(wrap_berry_phase(z))
    assert 0 <= res["Berry_Phase_pi"] < 2


def test_polarization_property_matches_the_density_difference():
    lay = spin_block_layout(3, mode="up_spin")
    sv = _basis_state(12, [0, 3])              # A_0 and B_1 occupied
    res = PROPERTIES["polarization"](sv, lay)
    assert res["n_A"] == pytest.approx([1.0, 0.0, 0.0])
    assert res["n_B"] == pytest.approx([0.0, 1.0, 0.0])
    assert res["n_A_minus_n_B"] == pytest.approx([1.0, -1.0, 0.0])
    assert res["unit_cell_j"] == [0, 1, 2]


def test_properties_depend_only_on_probabilities():
    """Both registered properties are diagonal, so a global rephasing of each
    amplitude must leave them unchanged. This is what makes shot-based
    estimation a drop-in replacement."""
    lay = spin_block_layout(3, mode="up_spin")
    rng = np.random.default_rng(3)
    sv = rng.normal(size=4096) + 1j * rng.normal(size=4096)
    sv /= np.linalg.norm(sv)
    phased = sv * np.exp(1j * rng.uniform(0, 2 * np.pi, size=sv.shape))
    for name in PROPERTIES:
        a, b = PROPERTIES[name](sv, lay), PROPERTIES[name](phased, lay)
        for k in a:
            if isinstance(a[k], float):
                assert a[k] == pytest.approx(b[k], abs=1e-9), (name, k)


def test_polarization_rejects_a_ragged_custom_layout_clearly():
    """Unequal cell sizes cannot define A/B sublattices; the error must say so
    rather than surfacing as a bare IndexError from inside the loop."""
    lay = {"n_qubits": 3, "n_cells": 2, "cell_qubits": [[0, 1], [2]],
           "A_qubits": [0, 2], "B_qubits": [1]}
    sv = _basis_state(3, [0])
    with pytest.raises(ValueError, match="tile into 2 cells"):
        PROPERTIES["polarization"](sv, lay)


def test_berry_still_works_on_a_ragged_layout():
    """berry only needs cell_qubits, so it must not be blocked by the check above."""
    lay = {"n_qubits": 3, "n_cells": 2, "cell_qubits": [[0, 1], [2]],
           "A_qubits": [0, 2], "B_qubits": [1]}
    res = PROPERTIES["berry"](_basis_state(3, [0]), lay)
    assert res["Twist_Amplitude"] == pytest.approx(1.0)
