介紹 OpenFermion 裡面的 slater determinant preparation circuit
他可以使用 Slater Determinant State Preparation Algorithm
---
Input : nu times N 的軌域係數矩陣Q, Q^\dagger Q = I_M

where N : 軌域總數
\nu 電子數量

Output : { (j_1, k_1, theta_1, phi_1), (j_2, k_2, theta_2, phi_2), ...}
---
之後再重建的時候使用 Given matrix 來製作initial state

每一個 Givens 旋轉在兩個 qubit 的佔據基底 $\{|00\rangle, |01\rangle, |10\rangle, |11\rangle\}$ 上應該是

$$
G(\theta,\phi)=
\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & \cos\theta & \sin\theta & 0\\
0 & -e^{i\phi}\sin\theta & e^{i\phi}\cos\theta & 0\\
0 & 0 & 0 & e^{i\phi}
\end{pmatrix}.
$$


可以使用 @/home/b11202015/sda/Adiabatic-Quantum-Simulation/docs/defect_report/slater_determinant_circuit_example.ipynb 作為範例

我們使用了
```
# ---- 修正前：逐字抄自 main 0fe9d1d 的 core._givens_instruction ----
def legacy_givens(theta, phi):
    gv = QuantumCircuit(2, name="Givens")
    gv.rz(phi, 1)
    gv.rz(-phi, 0)
    gv.cx(0, 1)
    gv.cry(-2.0 * theta, 1, 0)      # <-- 角度反號
    gv.cx(0, 1)
    return gv.to_instruction()
```

也就是 main branch 的程式，進行測試比較的是使用

```
c, s = np.cos(theta), np.sin(theta)
    e = np.exp(1j * phi)
    U = np.zeros((4, 4), dtype=complex)
    U[0, 0] = 1.0
    U[1, 1] = c
    U[2, 2] = c * e
    U[1, 2] = s
    U[2, 1] = -s * e
    U[3, 3] = e
    return UnitaryGate(U, label="Givens")
```

使用最檢單的氫分子模型，發現兩者製備出來的並不相同。

因此我使用數學檢測，可以知道

Qubit 0 --(R_z(-phi))--(R_y(-2theta)) ---.---
                             |           |   
Qubit 1 --(R_z(phi))---------.-----------x---

$$U_{rz} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & e^{-i\phi} & 0 & 0 \\ 0 & 0 & e^{i\phi} & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

$$M = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & \cos\theta & \sin\theta \\ 0 & 0 & -\sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \end{pmatrix}$$

$$U_{\text{mid}} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 \end{pmatrix} \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & \sin\theta & \cos\theta & 0 \\ 0 & \cos\theta & -\sin\theta & 0 \end{pmatrix}$$

$$U = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos\theta & -\sin\theta & 0 \\ 0 & \sin\theta & \cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & e^{-i\phi} & 0 & 0 \\ 0 & 0 & e^{i\phi} & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & e^{-i\phi}\cos\theta & -e^{i\phi}\sin\theta & 0 \\ 0 & e^{-i\phi}\sin\theta & e^{i\phi}\cos\theta & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

與Given確實不相同

$$
G(\theta,\phi)=
\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & \cos\theta & \sin\theta & 0\\
0 & -e^{i\phi}\sin\theta & e^{i\phi}\cos\theta & 0\\
0 & 0 & 0 & e^{i\phi}
\end{pmatrix}.
$$

因此我們檢測這個不同之處對於學姊論文產生的影響。

/home/b11202015/sda/Adiabatic-Quantum-Simulation/docs/defect_report/paper_impact_recompute.py

以及
/home/b11202015/sda/Adiabatic-Quantum-Simulation/docs/defect_report/paper_impact.md

