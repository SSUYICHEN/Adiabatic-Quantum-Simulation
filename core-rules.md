# Core rules — Adiabatic Quantum Simulation (AQS)

Gate-based adiabatic quantum simulation of 1D topological fermion models. Two authorities
govern all physics in this repo:

1. **Paper theory** — `.claude/specs/paper_theory/SSHH.md`: spinful SSH–Hubbard (SSHH)
   model, 4N qubits (paper equations cited throughout the rules).
2. **Reconstruction spec** — `.claude/specs/reconstruct_goals.md`: extended **spinless**
   SSH model with nearest-neighbor interaction, 2N qubits, static H with λ-ramped
   interaction; deliverable class `SpinlessSSHHSim`, benchmark N=6, T=1.0, L=40.

Both models coexist and share gates, observables, and backends. Prior debugging history and
its lessons: `docs/experience.md` (read it before touching any gate, mapping, or state prep).

## Architecture map (theory ↔ code)

Target layout for the reconstruction; current flat file in parentheses until the
restructure lands. As of 2026-09-14, `src/aqs/models/spinless.py` exists in the
working tree; the remaining layout is still largely flat. Check the tree and Git
status before treating this target map or historical branch notes as current state.

| Module | Responsibility | Theory |
|---|---|---|
| `models/sshh.py` (`core.py`) | Spinful SSHH model, single-particle H, Slater Q matrices | paper Eq. 9–15 |
| `models/spinless.py` (present in working tree) | Spinless SSH+NN model, `SpinlessSSHHSim` | spec §1–2 |
| `mapping.py` (`hamiltonian.py` + layout helpers) | JW conventions, qubit layouts, `_to_little_endian` — the only place bit order is decided | paper Eq. 26–29, 36; spec §1 |
| `circuits/gates.py` (`core.py`) | R / G / CP primitives, Givens rotation | paper Eq. 37–39 ≡ spec §3.1 |
| `circuits/state_prep.py` (`core.py`) | X pattern + OpenFermion Givens network | paper Eq. 42–45; spec §4 |
| `circuits/trotter.py` (`core.py`) | Trotterized adiabatic evolution, midpoint ramp, PBC parity | paper Eq. 33–35, 40–41; spec §3.2 |
| `observables.py` | Twist invariant z_N / Berry phase, densities, polarization, `PROPERTIES` registry | paper Eq. 16–20, 46–52; spec §5 |
| `hamiltonians.py` (`hamiltonian.py`) | JSON Hamiltonian I/O + exact references (fixed-sector ED) | extension |
| `backends/` (`backends.py`) | qiskit + CUDA-Q behind one endianness contract; `selftest` | — |
| `experiments.py` | `fidelity_scan` (Fig 2), `berry_sweep` (Fig 3, PBC), `polarization_sweep` (Fig 4, OBC), spinless sweeps | paper Sec. V; spec §6 |
| `plotting.py` / `cli.py` | Figure rendering; frozen CLI surface | — |

## Absolute non-negotiables

1. **Never trust "runs without error."** Every historical defect here was a silent
   convention mismatch with plausible output. Physics changes are validated only against
   independent references (fixed-sector ED, OpenFermion-built determinants, lowest-fill
   energy, H_eff round-trip) — see `.claude/rules/testing-and-verification.md`.
2. **All endianness goes through the conversion/layout module.** OpenFermion is big-endian;
   everything OpenFermion-built passes `_to_little_endian` before comparison. The measured
   endianness table in `.claude/rules/endianness-and-mapping.md` is law.
3. **Convention choices are cited decisions, never refactor side effects.** Gate angle
   signs, operand order, conjugation, parity factors — each has a docstring citation and a
   guard test; changing one requires an H_eff round-trip proof.
4. **Never mix the two gate-convention eras.** This branch's primitives match their
   docstrings; `main` and `fix/cudaq-bit-order-and-uv-migration` compensate at call sites.
   Blending them is a silent semantic conflict git will not flag.
5. **Guard tests are reverse-verified** (historical bug reintroduced ⇒ test fails) and use
   **asymmetric inputs + complex hoppings** — symmetric/uniform inputs are structurally
   blind to the bug classes of this repo.
6. **GPU machine:** always `uv sync --extra cu12`; never `git stash -u`; `git add` explicit
   filenames only; check current disk space with `df -h` before large operations.

## Core invariants (details live in the rules)

- Qubit layouts — spinful: 4N, up `[0,2N)` / down `[2N,4N)`, even=A odd=B; spinless: 2N,
  `q_2j=A_j`, `q_2j+1=B_j`. → `endianness-and-mapping.md`
- PBC wrap: hopping bond carries parity `−(−1)^{n_s}`; interaction terms carry **none**. → `gates-and-circuits.md`
- Interaction ramp at layer midpoint `(2ℓ−1)/(2L)`; hopping angles constant. → `gates-and-circuits.md`
- Polarization ≡ 0 under PBC (edge physics needs OBC); Berry phase only meaningful as the
  π jump, absolute values carry a constant offset. → `models-and-observables.md`
- Prepared state energy = sum of lowest occupied orbitals (wrong determinant gives exactly
  −E, not a drift). → `gates-and-circuits.md`
- References always in the fixed particle-number sector. → `models-and-observables.md`
- Frozen contracts: package name `aqs`, CLI entry + subcommands, JSON Hamiltonian schema,
  `PROPERTIES` keys. → `experiments-and-cli.md`

## Build & test commands

```bash
uv sync --extra cu12                 # GPU machine (plain `uv sync` removes GPU wheels)
uv run --no-sync pytest tests/ -q              # full suite
uv run --no-sync coverage run --source=src/aqs -m pytest tests/ && uv run --no-sync coverage report -m
uv run --no-sync aqs selftest --backend qiskit
uv run --no-sync aqs selftest --backend cudaq  # cross-backend endianness gatekeeper
uv run --no-sync python tests/verify_consistency_vs_main.py out.json   # A/B vs another revision
```

Cross-revision comparison: `git worktree` + `PYTHONPATH` (recipe in
`.claude/rules/environment.md`), never stash.

## Self-verification ladder

After any physics change, run in order — each rung catches what the previous can't:

1. Unit guards (`uv run --no-sync pytest tests/ -q`) — the four convention suites.
2. **H_eff = i·log(U)/δt round-trip** at multiple δt — catches systematic sign/conjugation
   errors that do not shrink with δt.
3. **L-convergence scan** — fidelity must improve with L; a plateau below 1 means a
   convention bug, not Trotter error (fix one of T/L, scan the other).
4. **Fixed-sector ED cross-check** of prepared/evolved states and observables.
5. **Figure reproduction** at the reference operating point (N=6, T=1, L=40).
6. **Dual-backend selftest** (qiskit vs CUDA-Q on an asymmetric state).

## Rules index (`.claude/rules/`)

- `testing-and-verification.md` — verification methodology; what counts as evidence.
- `endianness-and-mapping.md` — the endianness table, JW/layout single-source rules.
- `gates-and-circuits.md` — gate conventions, both Trotter layers, state prep.
- `models-and-observables.md` — the two models' authorities, observable invariants.
- `backends.md` — backend statevector contract, CUDA-Q specifics, selftest.
- `experiments-and-cli.md` — figure contracts, frozen CLI/JSON surfaces, examples.
- `environment.md` — uv/GPU/disk/git/typst constraints.

## Shared instruction loading and durable memory

This file is the single source for project-wide core rules. `AGENTS.md` and
`CLAUDE.md` are thin entry points; edit shared policy here instead of copying it
into either entry point.

At the start of each task, read this file and all seven files in `.claude/rules/`
listed above. Both agents must explicitly read those files: do not rely on a
platform automatically discovering that directory or interpreting its `paths`
frontmatter. Treat those paths as topic hints; physics verification applies to
physics changes even when no test file has been edited. Relative paths in this
file are relative to the repository root. The detailed rules remain in their
existing locations so code and theory citations keep resolving.

Read `docs/memory/README.md` for persistent project decisions. Record new durable
project decisions, unresolved issues, and verification evidence in repository
Markdown files; platform memory must not be the only copy. Historical notes are
not proof of the current branch, hardware, package versions, or test results.
Inspect current evidence before acting on them. Never store credentials in these
files.

The build/test commands above assume the environment has already been synced
with the appropriate extra. `--no-sync` avoids changing installed dependencies
while running checks. On this GPU workstation, preserve `cu12`; on a different
machine, inspect the driver and Python before choosing an extra.

No project skills exist at this audit. Future portable skills should have one
canonical `SKILL.md` per skill under `.agents/skills/<name>/`, with a relative
link at `.claude/skills/<name>` pointing to `../../.agents/skills/<name>`.
Keep platform-specific permissions, hooks, and tool configuration out of shared
skill instructions. Verify discovery in both clients after adding a skill.
