# Adiabatic Quantum Simulation

Adiabatic Quantum Simulation of topological properties on a quantum simulator — for the built-in
**SSH-Hubbard** model *or* **arbitrary Hamiltonian**, on a
**CPU (qiskit)** or **NVIDIA GPU (CUDA-Q)** backend.
This repository contains the official implementation of the paper:
[![arXiv](https://img.shields.io/badge/arXiv-2605.11823-b31b1b.svg)](https://arxiv.org/abs/2605.11823)

## What it does

| Command | Purpose |
|---------|---------|
| `aqs fidelity` | Trotter-convergence check → choose annealing time **T** and **steps/L** |
| `aqs berry` | Berry phase (twist invariant, **PBC**) vs inter-cell hopping `w` |
| `aqs polarization` | Electron polarization `⟨n_A⟩−⟨n_B⟩` (edge response, **OBC**) vs `ΔU` |
| `aqs measure` | Topological properties of an **arbitrary Hamiltonian file** |
| `aqs selftest` | Verify a backend's statevector against qiskit (e.g. check CUDA-Q) |
| `aqs plot` | (Re)draw any figure from saved JSON |

## Install with uv (recommended)

The project is managed with [uv](https://docs.astral.sh/uv/); `uv.lock` pins the
whole environment, and `aqs` is always installed in **editable** mode, so edits
under `src/aqs/` take effect immediately with no reinstall.

```bash
git clone https://github.com/SSUYICHEN/Adiabatic-Quantum-Simulation.git
cd Adiabatic-Quantum-Simulation

uv sync                 # CPU only (qiskit backend) - works on any platform
uv sync --extra cu12    # + NVIDIA GPU, CUDA 12.x driver
uv sync --extra cu13    # + NVIDIA GPU, CUDA 13.x driver

uv run aqs selftest --backend cudaq   # verify the GPU backend
uv run aqs berry --N 6 ...            # run any command
```

Pick the extra matching your **driver's** CUDA version, shown top-right in
`nvidia-smi`: a `CUDA Version: 12.x` driver needs `cu12`, `13.x` needs `cu13`.
They are mutually exclusive (declared under `[tool.uv] conflicts`) - the cu13
wheels fail at import against a CUDA 12 driver, so there is no single `gpu`
extra that guesses for you. `cu13` additionally requires Python >= 3.11.

> **Note:** `uv sync` *without* an extra removes the GPU packages (~4 GB) from
> the environment, since uv makes the venv exactly match what you asked for.
> On a GPU box always pass `--extra cu12` (or `cu13`).

`.python-version` pins Python 3.10, the interpreter this project's GPU results
were verified on (CUDA-Q 0.12.0). Newer CUDA-Q (0.15.x) requires Python >= 3.11;
to use it run `uv sync -p 3.12 --extra cu12`.

CUDA-Q is Linux/GPU-only; on Windows use `--backend qiskit` (or run in WSL).

### Install with pip (alternative)

```powershell
py -m venv .venv
pip install -e .                 # CPU
pip install -e ".[cu12]"         # + GPU, CUDA 12.x driver
pip install -r requirements-dev.txt   # tests
```

## Tests

`pytest` and `coverage` are in the `dev` dependency group, which `uv sync`
installs by default — no separate step.

```bash
uv run pytest tests/ -q

uv run coverage run --source=src/aqs -m pytest tests/
uv run coverage report -m
```

The suite covers the conventions that this codebase gets silently wrong: qubit
bit-ordering (qiskit / CUDA-Q / OpenFermion each differ), gate-angle
conventions, and Slater-determinant state preparation. The GPU tests skip
cleanly without CUDA-Q. See `docs/experience.md` for what these guard against
and why the shipped example Hamiltonians cannot detect ordering bugs.

One-off verification scripts (not part of the suite) compare two revisions:

```bash
uv run python tests/verify_consistency_vs_main.py out.json
uv run python tests/verify_consistency_vs_main.py --compare before.json after.json
```

## Workflow (built-in SSH-Hubbard model)

```powershell
# Step 0: choose T and steps
aqs fidelity --N 6 --v 0.5 --w 1.5 --electrons 12 --U 0

# Step 1: Berry phase (PBC)
aqs berry --N 6 --electrons 12 --v 1.0 --w 0,0.25,0.5,0.75,0.99,1.01,1.25,1.5,1.75,2.0 --UA 0.01 --delta-U 0,0.1,1 --TA 1 --steps 40

# Step 2: electron polarization (OBC, default:half-filling)
aqs polarization --N 6 --v 0.5 --w 1.5 --UA 1.0 --delta-U 0,0.01,0.1,1,3 --TA 1 --steps 40

# Run any of these on GPU:
aqs berry --N 6 ... --backend cudaq
```

## Arbitrary Hamiltonian

```powershell
aqs measure --hamiltonian examples\ssh_spinless_N3_topological.json --property berry,polarization --backend qiskit
```
The ground state is obtained by exact diagonalisation (CPU `scipy` or GPU `cupy`), then the requested properties are measured.

### Hamiltonian file format (JSON, self-describing)

```json
{
  "type": "hamiltonian",
  "format": "fermion_operator | qubit_operator | dense_matrix | sparse_matrix",
  "metadata": {
    "n_qubits": 6, "n_cells": 3,
    "layout": "spin_block | spinless | custom",
    "mode": "up_spin | total",
    "cell_qubits": [[0,1],[2,3],[4,5]],
    "A_qubits": [0,2,4], "B_qubits": [1,3,5]
  },
  "terms": [ {"coeff": [-0.5, 0.0], "ops": "0^ 1"}, ... ],
  "operator_string": "-1.0 [0^ 1] + -1.0 [1^ 0]",
  "matrix_file": "H.npy"
}
```
- `fermion_operator` → mapped to qubits via Jordan-Wigner (OpenFermion).
- `qubit_operator` → Pauli strings (e.g. `"X0 Z1"`).
- `dense_matrix` / `sparse_matrix` → `matrix_file` (`.npy`/`.npz`) or inline `"data"`.

See `examples/make_example_hamiltonian.py` for a generator.

### Adding a new property

Register a function `f(statevector, layout) -> dict` in
`src/observables.py::PROPERTIES`. It is then available everywhere via
`--property <name>`.

## Backends

| `--backend` | Hardware | Circuit statevector | Ground state (measure) |
|-------------|----------|---------------------|------------------------|
| `qiskit` (default) | CPU | qiskit `Statevector` | `scipy` sparse eigsh |
| `cudaq` | NVIDIA GPU | CUDA-Q kernel + `get_state` | `cupy` GPU eigsh |

`aqs selftest --backend cudaq` builds a small circuit on both backends and
checks the statevectors agree (also validates the qiskit↔CUDA-Q endianness map).

## Output layout

```
aqs_results/
  fidelity/      fidelity_N6_U0_ground.json + .png
  berry/<run>/   UA_..._w_....json          + BerryPhase_PhaseDiagram.png
  polarization/<run>/  Polarization_*.json  + Polarization_vs_DeltaU.png
  measure/       <hamiltonian>_<backend>.json (+ _polarization.png)
```

## Module map

| File | Role |
|------|------|
| `src/aqs/core.py` | SSH-Hubbard model, gates, parity-corrected annealing circuit |
| `src/aqs/observables.py` | twist invariant, density, **property registry** |
| `src/aqs/backends.py` | qiskit (CPU) / cudaq (GPU) backends + selftest |
| `src/aqs/hamiltonian.py` | load arbitrary Hamiltonian files → ground state + measure |
| `src/aqs/experiments.py` | `fidelity_scan`, `berry_sweep`, `polarization_sweep` |
| `src/aqs/plotting.py` | the three figure types |
| `src/aqs/cli.py` | the `aqs` command-line interface |

Run `aqs <command> -h` for the full option list of each command.
