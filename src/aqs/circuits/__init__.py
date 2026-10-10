"""Circuit building blocks: gate primitives, state preparation, Trotter evolution."""

from aqs.circuits.gates import create_R_gate, create_G_gate, create_CP_gate
from aqs.circuits.state_prep import givens_gate, prepare_slater_determinant
from aqs.circuits.trotter import build_annealing_circuit

__all__ = [
    "create_R_gate", "create_G_gate", "create_CP_gate",
    "givens_gate", "prepare_slater_determinant",
    "build_annealing_circuit",
]
