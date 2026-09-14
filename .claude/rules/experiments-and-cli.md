---
paths:
  - "src/aqs/experiments.py"
  - "src/aqs/plotting.py"
  - "src/aqs/cli.py"
  - "examples/**"
---

# Experiments, Plotting, CLI & Examples

## Figure contracts (paper Sec. V)

| Runner | Reproduces | Boundary condition |
|---|---|---|
| `fidelity_scan` | Fig. 2 — Trotter/adiabatic validation | caller-selectable |
| `berry_sweep` | Fig. 3 — Berry phase vs w, ΔU | **PBC forced in the runner** |
| `polarization_sweep` | Fig. 4 — polarization profile vs ΔU | **OBC forced in the runner** |

Boundary-condition forcing lives in the experiment runner, not the CLI — physics
constraints (Berry needs PBC, polarization needs OBC) must be impossible to override from
the command line by accident.

- Reference operating points: spinful N=6 (24 qubits), spinless N=6 (12 qubits); T=1,
  L=40 (paper Sec. V-A; spec §6). Sweeps must accept these and reproduce the qualitative
  signatures: fidelity → 1 with growing L (a plateau below 1 is a bug — see
  `testing-and-verification.md`), π Berry jump across |v|=|w|, edge-localized polarization
  in the topological phase.
- Fidelity references: without interaction, compare against the exact final state; with
  interaction, adjacent-L comparison is a **stability** diagnostic only, not a proof of
  proximity to the true state — never present it as accuracy.
- Spinless sweeps (spec §6.4): scan `V_v, V_w` to locate the topological-breakdown
  threshold; same discipline.

## CLI & format freezes

- Frozen public contracts: entry point `aqs = aqs.cli:main` (+ `python -m aqs`),
  subcommand names (`fidelity`, `berry`, `polarization`, `measure`, `selftest`, `plot`),
  the JSON Hamiltonian schema (`type`/`format`/`metadata`/`terms`/`operator_string`/
  `matrix_file`/`data`), and `PROPERTIES` registry keys. Extensions are additive.
- Any CLI default that diverges from the paper's parameters is documented in `--help`
  (e.g. polarization defaults v=0.5, w=1.5 vs paper's v=0.1, w=1.0) — divergence is a
  recorded decision, never implicit.
- `--backend` is available on every measurement command; `selftest` exits non-zero on
  failure.

## Examples

- Shipped examples must include at least one **asymmetric-density** case. The uniform
  spinless examples (`ssh_spinless_N3_*.json`, density 0.5 everywhere) are format demos
  only — they are structurally blind to ordering bugs and must never be the only smoke
  test cited as evidence of correctness.
