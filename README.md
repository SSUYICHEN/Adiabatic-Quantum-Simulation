# SSH-Hubbard QA Toolkit (`sshh_qa`)

A practical command-line tool that turns the research scripts for the
**Su-Schrieffer-Heeger-Hubbard (SSH-H)** model into one reproducible workflow:

1. **`fidelity`** — verify Trotter convergence and choose the experiment
   parameters **T (annealing time)** and **steps / L**.
2. **`berry`** — measure the **Berry phase** (bulk topological invariant, PBC)
   over a parameter sweep, saving one JSON per point and a phase-diagram figure.
3. **`polarization`** — measure the **electron polarization**
   `⟨n_A⟩ − ⟨n_B⟩` (edge response, OBC half-filling) vs `ΔU`, saving JSON + figure.
4. **`plot`** — (re)draw any figure from saved JSON.

All commands share one physics core (`sshh_qa/core.py`), so the numerics are
identical to the original standalone scripts (`TrotterforPBC.py`, `deltaU.py`,
`halffillingedgewithdeltaU.py`, …).

---

## Install

The dependencies adapt to your Python version automatically (Python 3.9 gets
conservative pins; Python 3.10+ gets current releases), so the same package
installs everywhere.

**Option A — install the pre-built wheel (recommended for sharing):**

```powershell
py -m venv .venv
.venv\Scripts\python -m pip install dist\sshh_qa-1.0.0-py3-none-any.whl
```

**Option B — install from this source folder:**

```powershell
py -m venv .venv
.venv\Scripts\python -m pip install .
```

**Option C — just the dependencies, run as a module (no install):**

```powershell
py -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

A ready-made `.venv` already exists in this folder, and a `sshh-qa.bat` launcher
is provided that uses it automatically.

> **Build a fresh wheel to give to others:**
> `.venv\Scripts\python -m build --wheel` → `dist\sshh_qa-1.0.0-py3-none-any.whl`.
> Do **not** copy the `.venv\` folder to another machine — let them recreate it.

## Run

After Option A/B the `sshh-qa` command exists; with Option C use `python -m sshh_qa`.
Any of these work (the `.bat` picks up the venv for you):

```powershell
sshh-qa fidelity --N 6 --t1 0.5 --t2 1.5 --electrons 12 --U 0
.venv\Scripts\sshh-qa fidelity --N 6 --t1 0.5 --t2 1.5 --electrons 12 --U 0
.venv\Scripts\python -m sshh_qa fidelity --N 6 --t1 0.5 --t2 1.5 --electrons 12 --U 0
```

---

## The intended workflow

### Step 0 — decide T and L with the fidelity check

Non-interacting Trotter error (`U = 0`, fidelity vs the prepared ground state):

```powershell
sshh-qa fidelity --N 6 --t1 0.5 --t2 1.5 --electrons 12 --U 0 --reference ground
```

Interacting self-convergence (`U ≠ 0`, fidelity between adjacent step counts):

```powershell
sshh-qa fidelity --N 6 --t1 0.5 --t2 1.5 --electrons 12 --U 1 --reference adjacent
```

The console prints a fidelity table and the first `(T_A, steps)` that reaches
fidelity ≥ 0.999 — use that as the operating point for the measurements below.

### Step 1 — Berry phase (topological property #1, PBC)

```powershell
sshh-qa berry --N 6 --electrons 12 --t1 1.0 ^
    --t2 0,0.25,0.5,0.75,0.99,1.01,1.25,1.5,1.75,2.0 ^
    --UA 0.01 --delta-U 0,0.01,0.1,0.3,1,3 --TA 1 --steps 40
```

Produces one JSON per `(ΔU, w)` point and `BerryPhase_PhaseDiagram.png`
(solid line = topology intact, dotted = gap closed by the interaction).

### Step 2 — Electron polarization (topological property #2, OBC half-filling)

```powershell
sshh-qa polarization --N 6 --t1 0.5 --t2 1.5 --UA 1.0 ^
    --delta-U 0,0.0001,0.001,0.01,0.1,1,3 --TA 1 --steps 40
```

Produces one JSON per `ΔU` and `Polarization_vs_DeltaU.png`.

### Re-plot later from saved JSON

```powershell
sshh-qa plot --kind berry        --data QA_results\berry\N6_TA1_steps40_UA0.01
sshh-qa plot --kind polarization --data QA_results\polarization\N6_TA1_steps40_t10.5
sshh-qa plot --kind fidelity     --data QA_results\fidelity\fidelity_N6_U0_ground.json
```

---

## Output layout

```
QA_results/
  fidelity/      fidelity_N6_U0_ground.json   + .png
  berry/<run>/   UA_..._UB_..._t2_....json    + BerryPhase_PhaseDiagram.png
  polarization/<run>/  Polarization_UA_..._UB_....json + Polarization_vs_DeltaU.png
```

## Key options (all commands)

| Option | Meaning |
|--------|---------|
| `--N` | number of unit cells (qubits used = `2 * 2N`) |
| `--t1` / `--t2` | intra-cell `v` / inter-cell `w` hopping |
| `--electrons` | total electron count |
| `--UA` / `--delta-U` | A-sublattice U and the staggering `ΔU = U_B − U_A` |
| `--TA` / `--steps` | annealing time and Trotter steps (from the fidelity check) |
| `--no-plot` | compute + save JSON only, skip the figure |

Run `sshh-qa <command> -h` for the full option list.

## Module map

| File | Role |
|------|------|
| `sshh_qa/core.py` | model, gates, parity-corrected annealing circuit |
| `sshh_qa/observables.py` | twist invariant (Berry phase), density profile |
| `sshh_qa/experiments.py` | `fidelity_scan`, `berry_sweep`, `polarization_sweep` |
| `sshh_qa/plotting.py` | the three figure types |
| `sshh_qa/cli.py` | argparse command-line interface |
