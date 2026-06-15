"""
sshh_qa - SSH-Hubbard Adiabatic Quantum Simulation Toolkit
==========================================================

A unified, practical tool that consolidates the research workflow for the
Su-Schrieffer-Heeger-Hubbard (SSH-H) model into a single CLI:

  1. fidelity      - verify Trotter convergence to choose T_A and steps (L)
  2. berry         - measure the Berry phase (bulk topological invariant, PBC)
  3. polarization  - measure the electron polarization / edge response (OBC)
  4. plot          - render publication-quality figures from saved JSON

The physics core (model, gates, annealing circuit, observables) is shared by
every command so results are identical to the original standalone scripts.
"""

__version__ = "1.0.0"

from .core import (
    SSHHModel,
    build_annealing_circuit,
    create_R_gate,
    create_G_gate,
    create_CP_gate,
)
from .observables import twist_invariant, density_profile

__all__ = [
    "SSHHModel",
    "build_annealing_circuit",
    "create_R_gate",
    "create_G_gate",
    "create_CP_gate",
    "twist_invariant",
    "density_profile",
    "__version__",
]
