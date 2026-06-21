"""
aqs - Adiabatic Quantum Simulation toolkit for topological properties.
======================================================================

A practical CLI for measuring topological properties on a quantum simulator.

Two model sources:
  * built-in SSH-Hubbard model    (parametrised by v, w, U, N, filling)
  * an arbitrary Hamiltonian file  (`aqs measure --hamiltonian H.json`)

Two simulation backends:
  * qiskit  - CPU statevector (default)
  * cudaq   - NVIDIA CUDA-Q GPU statevector / cupy GPU diagonalisation

Two topological properties (extensible registry, see observables.py):
  * berry         - Berry phase (twist invariant, PBC)
  * polarization  - electron polarization <n_A> - <n_B> (edge response, OBC)

This was previously the SSH-only `sshh_qa` package; it is now generalised.
"""

__version__ = "2.0.0"

from .core import (
    SSHHModel,
    build_annealing_circuit,
    create_R_gate,
    create_G_gate,
    create_CP_gate,
)
from .observables import (
    twist_invariant,
    twist_invariant_layout,
    density_profile,
    wrap_berry_phase,
    PROPERTIES,
)

__all__ = [
    "SSHHModel",
    "build_annealing_circuit",
    "create_R_gate",
    "create_G_gate",
    "create_CP_gate",
    "twist_invariant",
    "twist_invariant_layout",
    "density_profile",
    "wrap_berry_phase",
    "PROPERTIES",
    "__version__",
]
