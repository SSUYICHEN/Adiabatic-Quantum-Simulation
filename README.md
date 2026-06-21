# AQS Toolkit (`aqs`)

**Adiabatic Quantum Simulation of topological properties.**
Measure topological invariants on a quantum simulator — for the built-in
**SSH-Hubbard** model *or* for an **arbitrary Hamiltonian** you supply, on a
**CPU (qiskit)** or **NVIDIA GPU (CUDA-Q)** backend.

## What it does

| Command | Purpose |
|---------|---------|
| `aqs fidelity` | Trotter-convergence check → choose annealing time **T** and **steps/L** |
| `aqs berry` | Berry phase (twist invariant, **PBC**) vs inter-cell hopping `w` |
| `aqs polarization` | Electron polarization `⟨n_A⟩−⟨n_B⟩` (edge response, **OBC**) vs `ΔU` |
| `aqs measure` | Topological properties of an **arbitrary Hamiltonian file** |
| `aqs selftest` | Verify a backend's statevector against qiskit (e.g. check CUDA-Q) |
| `aqs plot` | (Re)draw any figure from saved JSON |

## Install (CPU)

```powershell
py -m venv .venv
.venv\Scripts\python -m pip install dist\aqs-2.0.0-py3-none-any.whl
# or from source:  .venv\Scripts\python -m pip install .
```

## Install (GPU, NVIDIA CUDA-Q)

```bash
python -m pip install aqs-2.0.0-py3-none-any.whl
python -m pip install -r requirements-gpu.txt      # cudaq + cupy
aqs selftest --backend cudaq                       # verify the GPU backend
```
CUDA-Q is Linux/GPU-only; on Windows use `--backend qiskit` (or run in WSL).

## Workflow (built-in SSH-Hubbard model)

```powershell
# Step 0: choose T and steps
aqs fidelity --N 6 --v 0.5 --w 1.5 --electrons 12 --U 0

# Step 1: Berry phase (PBC)
aqs berry --N 6 --electrons 12 --v 1.0 ^
    --w 0,0.25,0.5,0.75,0.99,1.01,1.25,1.5,1.75,2.0 ^
    --UA 0.01 --delta-U 0,0.1,1 --TA 1 --steps 40

# Step 2: electron polarization (OBC half-filling)
aqs polarization --N 6 --v 0.5 --w 1.5 --UA 1.0 ^
    --delta-U 0,0.01,0.1,1,3 --TA 1 --steps 40

# Run any of these on GPU:
aqs berry --N 6 ... --backend cudaq
```

## Arbitrary Hamiltonian

```powershell
aqs measure --hamiltonian examples\ssh_spinless_N3_topological.json ^
            --property berry,polarization --backend qiskit
```
The ground state is obtained by exact diagonalisation (CPU `scipy` or GPU
`cupy`), then the requested properties are measured.

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
| `src/core.py` | SSH-Hubbard model, gates, parity-corrected annealing circuit |
| `src/observables.py` | twist invariant, density, **property registry** |
| `src/backends.py` | qiskit (CPU) / cudaq (GPU) backends + selftest |
| `src/hamiltonian.py` | load arbitrary Hamiltonian files → ground state + measure |
| `src/experiments.py` | `fidelity_scan`, `berry_sweep`, `polarization_sweep` |
| `src/plotting.py` | the three figure types |
| `src/cli.py` | the `aqs` command-line interface |

Run `aqs <command> -h` for the full option list of each command.
