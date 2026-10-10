"""Model definitions (reconstruction target layout).

`spinless`  -- extended spinless SSH + nearest-neighbor interaction
               (spec: .claude/specs/reconstruct_goals.md).
`hubbard`   -- 1D Hubbard chain: the spinful SSHH model with v = w and
               U_A = U_B (no sublattice), reusing the ``aqs.core`` circuit.
The spinful SSHH model of the paper still lives in ``aqs.core`` until the
full restructure lands.
"""

from aqs.models.spinless import (
    SpinlessSSHModel,
    SpinlessSSHHSim,
    spec_pauli_hamiltonian,
    exact_hamiltonian,
    exact_ground_state_fixed_n,
)
from aqs.models.hubbard import (
    HubbardModel,
    HubbardSim,
    hubbard_pauli_hamiltonian,
    hubbard_hamiltonian,
    hubbard_ground_state_fixed_n,
)

__all__ = [
    "HubbardModel",
    "HubbardSim",
    "hubbard_pauli_hamiltonian",
    "hubbard_hamiltonian",
    "hubbard_ground_state_fixed_n",
    "SpinlessSSHModel",
    "SpinlessSSHHSim",
    "spec_pauli_hamiltonian",
    "exact_hamiltonian",
    "exact_ground_state_fixed_n",
]
