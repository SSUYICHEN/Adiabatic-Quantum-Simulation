"""Model definitions (reconstruction target layout).

`spinless`  -- extended spinless SSH + nearest-neighbor interaction
               (spec: .claude/specs/reconstruct_goals.md).
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

__all__ = [
    "SpinlessSSHModel",
    "SpinlessSSHHSim",
    "spec_pauli_hamiltonian",
    "exact_hamiltonian",
    "exact_ground_state_fixed_n",
]
