"""core.py - Compatibility shim.

The contents moved on 2026-10-10 to the layout of core-rules.md:

    SSHHModel                     -> aqs.models.sshh
    create_R_gate / G / CP        -> aqs.circuits.gates
    build_annealing_circuit       -> aqs.circuits.trotter
    (Slater preparation)          -> aqs.circuits.state_prep

Import from those modules in new code; this file only re-exports them so that
existing scripts and the frozen public API keep working.
"""

from aqs.models.sshh import SSHHModel
from aqs.circuits.gates import create_R_gate, create_G_gate, create_CP_gate
from aqs.circuits.trotter import build_annealing_circuit

__all__ = ["SSHHModel", "create_R_gate", "create_G_gate", "create_CP_gate", "build_annealing_circuit"]
