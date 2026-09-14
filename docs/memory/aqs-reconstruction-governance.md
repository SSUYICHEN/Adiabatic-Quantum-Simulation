---
name: aqs-reconstruction-governance
description: "Governance layer + key decisions for the AQS codebase reconstruction (two models, two spec authorities, gate-convention hazard)"
metadata: 
  node_type: memory
  type: project
  originSessionId: a823726b-dc71-47bc-8c1f-06edf54d5fb4
  modified: 2026-08-11T07:21:16.610Z
---

As of 2026-08-11, the AQS repo has a governance layer for a planned full reconstruction:
`CLAUDE.md` + seven rule files in `.claude/rules/` + spec `.claude/specs/reconstruct_goals.md`.

User decisions (2026-08-10/11 session):
- `.claude/specs/reconstruct_goals.md` is the user's spec, saved **verbatim**: extended
  **spinless** SSH + nearest-neighbor interaction, 2N qubits, static H, λ_ℓ=(2ℓ−1)/2L ramp
  on V only, deliverable class `SpinlessSSHHSim`, benchmark N=6/12 qubits, T=1, L=40.
- **Both models coexist** in the reconstruction: spinful SSHH (paper
  `.claude/specs/paper_theory/SSHH.md`) + spinless (spec), sharing gates/observables/backends.
- **Keep both backends** (qiskit primary, CUDA-Q kept).
- Rules target a **new restructured layout** (`src/aqs/models/`, `circuits/`, `backends/`,
  `mapping.py`) with current flat paths as extra globs; the restructure itself has NOT
  happened yet — src/aqs is still flat.

Critical hazard to re-check each session: the spec's fermionic H is `+v b†a + h.c.` while
`core.py` builds `H[i,j]=−t` (conjugate convention) — never port angle signs between the
two models by analogy; settle G-sign questions only via H_eff = i·log(U)/δt round-trip
against the spec's Pauli H (little-endian) with complex v,w. The other branch's
`spinless.py` uses the old gate conventions and must not be copied.
