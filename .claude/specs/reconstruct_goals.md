# Quantum Circuit Simulation Specification: Extended Spinless SSH Model with Nearest-Neighbor Interaction (Static Hamiltonian)

## 1. System Overview & Problem Formulation
This project implements the Adiabatic Quantum Simulation (AQS) of the 1D **Extended Spinless Su-Schrieffer-Heeger (SSH) Model with Nearest-Neighbor Interaction** using Qiskit.

The system Hamiltonian is **time-independent**:
$$H = H_{\text{SSH}} + H_{\text{NN}}$$

The simulation adiabatically transforms the ground state of $H_{\text{SSH}}$ into the target ground state of $H$ by slowly ramping up the interaction strength parameter $\lambda(t) = \frac{t}{T}$ during Trotterized time evolution.

### System Parameters & Mapping
- **System Size:** $N$ unit cells ($2N$ total spatial lattice sites / qubits).
- **Boundary Conditions:** Periodic Boundary Conditions (PBC) or Open Boundary Conditions (OBC).
- **Total Qubit Count:** $2N$ qubits (1 qubit per spatial sublattice site $A_j, B_j$).
- **Qubit Indexing:**
  - $q_{2j} \rightarrow \text{Site } A_j$ (Unit cell $j$)
  - $q_{2j+1} \rightarrow \text{Site } B_j$ (Unit cell $j$)

---

## 2. Target Static Hamiltonian ($H = H_{\text{SSH}} + H_{\text{NN}}$)

### 2.1 Single-Particle SSH Hopping Term ($H_{\text{SSH}}$)
Fermionic Real-Space Representation:
$$H_{\text{SSH}} = \sum_{j=0}^{N-1} \left( v b_{j}^\dagger a_{j} + v^* a_{j}^\dagger b_{j} \right) + \sum_{j=0}^{N-2} \left( w b_{j}^\dagger a_{j+1} + w^* a_{j+1}^\dagger b_{j} \right) + \text{PBC terms}$$

Pauli Operator Representation (via Jordan-Wigner Transformation):
- **Intracell Hopping ($v$):**
  $$H_{\text{intra}} = -\frac{\text{Re}(v)}{2} \sum_{j=0}^{N-1} (X_{2j} X_{2j+1} + Y_{2j} Y_{2j+1}) - \frac{\text{Im}(v)}{2} \sum_{j=0}^{N-1} (X_{2j} Y_{2j+1} - Y_{2j} X_{2j+1})$$
- **Intercell Hopping ($w$):**
  $$H_{\text{inter}} = -\frac{\text{Re}(w)}{2} \sum_{j=0}^{N-2} (X_{2j+1} X_{2j+2} + Y_{2j+1} Y_{2j+2}) - \frac{\text{Im}(w)}{2} \sum_{j=0}^{N-2} (X_{2j+1} Y_{2j+2} - Y_{2j+1} X_{2j+2})$$

### 2.2 Nearest-Neighbor Interaction Term ($H_{\text{NN}}$)
Fermionic Real-Space Representation:
$$H_{\text{NN}} = \sum_{j=0}^{N-1} V_v n_{Aj} n_{Bj} + \sum_{j=0}^{N-2} V_w n_{Bj} n_{A, j+1}$$
where $n_{Aj} = a_j^\dagger a_j = \frac{I - Z_{2j}}{2}$ and $n_{Bj} = b_j^\dagger b_j = \frac{I - Z_{2j+1}}{2}$.

Pauli Operator Representation (Ising / Controlled-Phase Form):
- **Intracell Interaction ($V_v$):**
  $$H_{\text{NN, intra}} = \frac{V_v}{4} \sum_{j=0}^{N-1} (I_{2j} - Z_{2j})(I_{2j+1} - Z_{2j+1})$$
- **Intercell Interaction ($V_w$):**
  $$H_{\text{NN, inter}} = \frac{V_w}{4} \sum_{j=0}^{N-2} (I_{2j+1} - Z_{2j+1})(I_{2j+2} - Z_{2j+2})$$

---

## 3. Quantum Circuit Gate Primitives & Trotter Layer Construction

### 3.1 Two-Qubit Gate Primitives
Define the following parameterized gates:
1. **$R_{i,j}(\theta)$ Gate** (XX + YY rotation):
   $$R_{i,j}(\theta) = \exp\left( -i \frac{\theta}{2} (X_i X_j + Y_i Y_j) \right)$$
2. **$G_{i,j}(\theta)$ Gate** (XY - YX rotation):
   $$G_{i,j}(\theta) = \exp\left( -i \frac{\theta}{2} (X_i Y_j - Y_i X_j) \right)$$
3. **$CP_{i,j}(\phi)$ Gate** (Controlled-Phase rotation):
   $$CP_{i,j}(\phi) = \exp\left( -i \frac{\phi}{4} (I_i - Z_i)(I_j - Z_j) \right) = \text{CPhase}(\phi) \text{ on qubits } (i, j)$$

### 3.2 Adiabatic Trotter Layer $U_l$ ($l = 1, \dots, L$)
For total evolution time $T$ divided into $L$ intervals with step size $\delta t = T / L$:

The dynamic interpolation weight at step $l$ is set to $\lambda_l = \frac{2l - 1}{2L}$.

Each Trotter Layer $U_l$ is constructed in sequence:
1. **Apply $H_{\text{SSH}}$ Propagators (Constant Parameters per step):**
   - Intracell: $R_{2j, 2j+1}(-\delta t \text{Re}(v))$ and $G_{2j, 2j+1}(-\delta t \text{Im}(v))$ for $j=0, \dots, N-1$.
   - Intercell: $R_{2j+1, 2j+2}(-\delta t \text{Re}(w))$ and $G_{2j+1, 2j+2}(-\delta t \text{Im}(w))$ for $j=0, \dots, N-2$.
2. **Apply $H_{\text{NN}}$ Propagators (Scaled by $\lambda_l$):**
   - Intracell CP gates: $CP_{2j, 2j+1}(\phi_{v, l})$ where $\phi_{v, l} = \delta t \cdot V_v \cdot \lambda_l$.
   - Intercell CP gates: $CP_{2j+1, 2j+2}(\phi_{w, l})$ where $\phi_{w, l} = \delta t \cdot V_w \cdot \lambda_l$.

---

## 4. Initial State Preparation
- Prepare the non-interacting Slater determinant ground state $|\Psi_0\rangle$ of $H_{\text{SSH}}$ on $2N$ qubits.
- Fill $N$ particles (half-filling for PBC) or $N+1$ particles (for OBC edge modes) using OpenFermion Givens-rotation state preparation circuits.

---

## 5. Observables & Measurement Protocols

### 5.1 Sublattice Polarization (OBC)
Evaluate spatial charge imbalance from Z-basis bitstrings:
$$P_j^e = \langle n_{Aj} - n_{Bj} \rangle = \frac{\langle Z_{2j+1} - Z_{2j} \rangle}{2}$$

### 5.2 Many-Body Berry Phase (PBC)
Construct position operator phase factor $\hat{X} = \sum_{q=0}^{2N-1} \lfloor q/2 \rfloor \cdot n_q$:
1. Measure $M$ shots in the computational basis.
2. Estimate complex phase $\bar{z}_N = \frac{1}{M} \sum_{m=1}^{M} \exp\left( \frac{i 2\pi}{N} \langle b^{(m)} | \hat{X} | b^{(m)} \rangle \right)$.
3. Extracted Berry Phase: $\gamma = \text{Im}(\ln \bar{z}_N) \pmod{2\pi}$.

---

## 6. Execution Tasks for Agent
1. Build a modular Python class `SpinlessSSHHSim` using Qiskit.
2. Setup a 12-qubit circuit ($N=6$ unit cells).
3. Benchmark Trotter performance ($L=40, T=1.0$).
4. Sweep interaction parameters ($V_v, V_w$) to detect topological phase breakdown.
