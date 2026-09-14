# Adiabatic Quantum Simulation of the Topological Su–Schrieffer–Heeger–Hubbard Model

Ssu-Yi Chen\*<sup>द</sup>, Bo-Hung Chen\*<sup>‡§∥</sup>, Dah-Wei Chiou\*<sup>‡§\*\*</sup>, and Jie-Hong Roland Jiang\*<sup>†‡§††</sup>
\*Graduate Institute of Electronics Engineering, National Taiwan University, Taipei, Taiwan

† Department of Electrical Engineering, National Taiwan University, Taipei, Taiwan

<sup>‡</sup> Center for Quantum Science and Engineering, National Taiwan University, Taipei, Taiwan

§ Physics Division, National Center for Theoretical Sciences, Taipei, Taiwan

Email: ¶r12943161@ntu.edu.tw, ∥kenny81778189@gmail.com, \*\*dwchiou@gmail.com, ††jhjiang@ntu.edu.tw

Abstract—We develop an adiabatic quantum simulation framework on gate-based quantum computers to probe topological signatures of the one-dimensional fermionic Su-Schrieffer-Heeger-Hubbard (SSHH) model. We present explicit quantum-circuit constructions for initial-state preparation and time evolution, together with a practical measurement protocol and classical post-processing procedure for extracting the many-body Berry phase and the spatial profile of the sublattice polarization. Using classical simulations of the proposed circuits, we demonstratefor the first time within a genuine many-body framework—that the topological characteristics of the SSH model remain robust against weak Hubbard interactions but eventually break down as the chiral-symmetry-breaking component of the interaction exceeds a threshold. The required qubit number, gate complexity, measurement shots, and classical pre- and post-processing costs all scale polynomially with system size. Our results provide a proof-of-concept framework for probing topological properties of interacting many-body systems via adiabatic quantum simulation on future large-scale quantum computers.

Index Terms—Adiabatic quantum simulation, Hamiltonian simulation, Topological invariants, Su-Schrieffer-Heeger model, Hubbard model.

# I. INTRODUCTION

Simulating interacting many-body systems is widely regarded as one of the most promising near-term applications of quantum computation in the noisy intermediate-scale quantum (NISQ) era [1], [2]. Classical numerical methods often become intractable because the Hilbert space grows exponentially with system size [3]–[5]. Quantum computers, by contrast, provide a natural platform for representing, evolving, and measuring quantum states [6], potentially enabling the study of interacting quantum systems beyond the reach of classical approaches [7].

A representative NISQ algorithm for such tasks is the *variational quantum eigensolver* (VQE) [8], [9], a hybrid quantum–classical method that estimates the ground-state energy of a target quantum system by iteratively optimizing the parameters of a parameterized quantum circuit [10], [11]. The method benefits from the variational principle: if the approximated state deviates from the exact target state by  $O(\varepsilon)$ , the resulting energy error is typically  $O(\varepsilon^2)$  [12]. However, an accurate energy does not necessarily imply an accurate wavefunction. Unless the *ansatz* (i.e., VQE circuit pattern) incorporates the relevant physical structure, the optimized state may fail to capture global, symmetry-protected, or topological properties

of the true ground state [13]–[15]. Consequently, for problems concerned with broader physical characteristics rather than energy estimation alone, VQE may not reliably reproduce the desired state.

An alternative approach is *adiabatic quantum simulation*, in which the target state is obtained by slowly deforming a known ground state of a simple Hamiltonian into that of the target system [16], [17]. Provided that the evolution is sufficiently slow and the energy gap to excited states remains finite, the final state remains close to the desired ground state.

Adiabatic dynamics can be implemented either on analog quantum annealers or on gate-based quantum computers. Analog quantum annealers, such as D-Wave systems [18], naturally realize annealing dynamics but are typically limited to transverse-field Ising Hamiltonians [19], [20]. This restriction limits their direct applicability to interacting fermionic models [21], whose qubit representations generally involve noncommuting and often nonlocal Pauli operators beyond the Ising form [22], [23]. Gate-based quantum computers, by contrast, can in principle simulate arbitrary quantum dynamics [24], [25]. In particular, adiabatic evolution can be realized within this framework, and may even be implemented digitally through fault-tolerant quantum computation [26].

In this work, we investigate the one-dimensional fermionic Su-Schrieffer-Heeger-Hubbard (SSHH) model [27] and demonstrate that nontrivial topological properties of many-body systems can be probed through adiabatic quantum simulation on gate-based quantum computers.

The SSHH model extends the Su–Schrieffer–Heeger (SSH) model [28]—a paradigmatic model for one-dimensional topological insulators—by including the Hubbard interaction that accounts for on-site electron–electron Coulomb repulsion [29]. In the absence of many-body interactions, the SSH model admits a simple single-particle description based on its band structure. The inclusion of the Hubbard term introduces particle–particle interactions and turns the system into a genuine many-body problem. The fermionic SSHH model therefore provides a natural platform for demonstrating the advantages of adiabatic quantum simulation and for exploring its limitations

The SSHH model is also of broad interest in condensed matter physics, as it can exhibit a variety of phenomena depending on the model parameters, including antiferromagnetic order, topological and semimetallic phases, and superconductor-Mott-insulator transitions [30], [31].

The SSH model has been extensively studied, including generalizations with arbitrary long-range hopping [32]. The effects of many-body interactions have also been explored from single-particle or approximate perspectives [33], and quantum simulations of SSH-type many-body systems have been demonstrated on gate-based quantum computers [34]. These studies suggest that the topological signatures of the SSH model are robust against weak interactions. However, this robustness has not yet been demonstrated from a genuinely many-body perspective.

Here we address this question by studying the SSHH model using adiabatic quantum simulation on gate-based quantum computers. Starting from a many-body ground state of the SSH model, we adiabatically introduce the Hubbard interaction. The initial state can be systematically prepared in quantum circuits based on pre-solved single-particle states, and the adiabatic evolution can likewise be implemented through a sequence of quantum gates. By extracting the many-body Berry phase and the spatial profile of the sublattice polarization through a practical measurement protocol and classical post-processing, we demonstrate-within a fully many-body framework—that the topological characteristics of the SSH model remain robust against weak Hubbard interactions but eventually break down when the chiral-symmetry-breaking component of the interaction becomes sufficiently strong.

Both the gate count and circuit complexity scale polynomially with system size for the initial-state preparation and the adiabatic evolution. The number of measurement shots required by the protocol and the runtime of the post-processing procedure also scale polynomially. Although our results are obtained from classical simulations for small systems, they provide a proof of concept that the same framework can be implemented on future large-scale quantum computers to investigate topological properties of larger interacting manybody systems.

#### II. ADIABATIC QUANTUM SIMULATION FRAMEWORK

We consider a many-body quantum system governed by the Hamiltonian  $H = H_0 + H_1$ , where  $H_0$  represents the one-body contribution and  $H_1$  the many-body interactions. The spectrum or ground state of  $H_0$ , or at least the corresponding properties, can usually be obtained efficiently by classical methods and, in some cases, even analytically. Once  $H_1$  is included, however, determining the properties of the full Hamiltonian H rapidly becomes intractable beyond modest system sizes, because the number of degrees of freedom of many-body systems grows exponentially with the system size.

Adiabatic quantum simulation offers a promising approach to overcome this difficulty in accessing ground-state properties. The central idea is to vary the Hamiltonian slowly from  $H_0$  to H, thereby transforming the ground state of  $H_0$  into that of H. This can be implemented by initializing the system in the ground state of  $H_0$  and evolving it under a time-dependent Hamiltonian, for example,

<span id="page-1-0"></span>
$$H(t) = H_0 + \frac{t}{T}H_1, \qquad t \in [0, T],$$
 (1)

which continuously interpolates between  $H_0$  at t=0 and H at t = T. If the interpolation is sufficiently slow (i.e., if T is large enough), the adiabatic theorem ensures that the evolving state remains close to the instantaneous ground state of H(t) throughout the evolution. Consequently, the final state at t = T approximates the ground state of H, allowing physical observables of the interacting many-body system to be extracted from measurements on the final state.

Let U(t,t') denote the time-evolution operator from t' to t generated by H(t). It satisfies the Schrödinger equation

$$i\hbar \frac{d}{dt}U(t,t') = H(t)U(t,t'), \tag{2}$$

with the initial condition  $U(t',t')=\mathbb{1}$ . Formally, the solution can be written as the time-ordered  $(\mathcal{T})$  exponential

$$U(T,0) = \mathcal{T} \exp\left(-\frac{i}{\hbar} \int_0^T H(t) dt\right). \tag{3}$$

Although this time-ordered exponential is generally difficult to evaluate exactly, a useful approximation can be obtained by considering the infinitesimal propagator over a short interval  $[t, t + \delta t]$ . To first order in  $\delta t$ , one finds

$$U(t + \delta t, t) = \mathbb{1} - \frac{i}{\hbar} \int_{t}^{t + \delta t} H(t') dt' + O(\delta t^{2})$$

$$= \mathbb{1} - \frac{i \delta t}{\hbar} \left[ H_{0} + \frac{t'}{T} H_{1} \right]_{t' = t + \delta t/2} + O(\delta t^{2})$$

$$= \exp \left[ -\frac{i \delta t}{\hbar} \left( \frac{t + \delta t/2}{T} \right) H_{1} \right] \exp \left[ -\frac{i \delta t}{\hbar} H_{0} \right] + O(\delta t^{2}).$$
This corresponds to the first order Teetter decomposition in

This corresponds to the first-order Trotter decomposition, in which the Hamiltonian is treated as constant within each short interval and evaluated at the midpoint.

The full propagator U(T,0) can therefore be approximated by a product of many short-time propagators. Setting

$$\delta t = T/L, \tag{5}$$

we obtain

$$U(T,0) = U(T,T-\delta t) \cdots U(\delta t,0) + O(L\delta t^2)$$
  

$$\equiv U_L U_{L-1} \cdots U_1 + O(T^2/L), \tag{6}$$

<span id="page-1-1"></span>

where the propagator for the 
$$\ell$$
-th interval is
$$U_{\ell} = \exp\left(-\frac{i\delta t}{\hbar} \frac{2\ell - 1}{2L} H_1\right) \exp\left(-\frac{i\delta t}{\hbar} H_0\right). \tag{7}$$

The requirement that the interpolation be sufficiently slow is quantified by the adiabatic condition, which demands that transitions between different instantaneous eigenstates of H(t)remain strongly suppressed throughout the evolution. Let  $|n(t)\rangle$  denote the instantaneous eigenstates of H(t) with corresponding eigenvalues  $E_n(t)$ . Using time-dependent perturbation theory to estimate transition amplitudes yields the adiabatic condition

$$\max_{t \in [0,T]} \max_{m \neq n} \frac{\left| \langle m(t) | \dot{H}(t) | n(t) \rangle \right|}{\left| E_m(t) - E_n(t) \right|^2} \ll \frac{1}{\hbar}.$$
 (8)

For the interpolation in (1), the Hamiltonian varies at a constant rate,  $\dot{H}(t) = H_1/T$ . The requirement that T be sufficiently large can therefore be expressed by the worst-case estimate  $T \gg \hbar \|H_1\|/\Delta_{\min}^2$ , where  $\|\cdot\|$  denotes the operator norm and  $\Delta_{\min} = \min_{t \in [0,T]} \min_{m \neq n} |E_m(t) - E_n(t)|$  is the minimum instantaneous energy gap encountered along the interpolation path.

In practical applications this bound can often be tightened because the initial state is prepared in a specific form. For instance, if the system is initialized in the ground state, the relevant index pairs (m, n) include only those with ncorresponding to the ground state, which substantially reduces the bound.

It should be noted that this nondegenerate adiabatic analysis breaks down if the relevant instantaneous energy gap closes at some point during the evolution, for example because of an energy level crossing or a degeneracy enforced by symmetry. In such cases the standard adiabatic bound becomes singular, and the evolution must instead be analyzed within a degenerate adiabatic framework or along a modified interpolation path that maintains a finite gap. In many physical systems, however, the closing of the gap signals a genuine quantum phase transition. In these situations the singularity cannot be removed by such reformulations, as it reflects an intrinsic physical property of the system rather than a limitation of the analysis. This issue will be examined more concretely in the context of the SSHH model discussed below.

#### III. FERMIONIC SSHH MODEL

The Hamiltonian of the fermionic SSHH model can be expressed as

$$H = H_{\text{SSH},\uparrow} + H_{\text{SSH},\downarrow} + H_{\text{Hubbard}},$$
 (9)

where  $H_{\rm SSH,\uparrow}$  and  $H_{\rm SSH,\downarrow}$  describe the SSH model for the two spin sectors ( $\uparrow$  and  $\downarrow$ ), and  $H_{\text{Hubbard}}$  represents the onsite Hubbard interaction.

Consider a system of electrons in a one-dimensional lattice with N unit cells, each consisting of two sublattice sites, denoted A and B. The lattice is assumed to satisfy either periodic boundary conditions (PBC) or open boundary conditions

Let  $a_{is}^{\dagger}(a_{js})$  and  $b_{is}^{\dagger}(b_{js})$  denote the fermionic creation (annihilation) operators on the A and B sublattices, respectively, in the j-th unit cell, where  $s \in \{\uparrow, \downarrow\}$  labels the spin. These operators satisfy the canonical anticommutation relations

$$\{a_{is},a_{js'}\}=\{a_{is}^{\dagger},a_{js'}^{\dagger}\}=\{b_{is},b_{js'}\}=\{b_{is}^{\dagger},b_{js'}^{\dagger}\}=0,\\ \{a_{is},a_{js'}^{\dagger}\}=\{b_{is},b_{js'}^{\dagger}\}=\delta_{ij}\delta_{ss'},\qquad \{a_{is},b_{js'}^{\dagger}\}=0, \ \ (10)$$
 which encode the Fermi–Dirac statistics of electrons. The corresponding number operators are

$$n_{Ajs}=a_{js}^{\dagger}a_{js}, \qquad n_{Bjs}=b_{js}^{\dagger}b_{js}.$$
 (11)  
The real-space Hamiltonian of the spinful SSH model is

$$H_{\text{SSH},s} = \sum_{j=1}^{N} v \, b_{js}^{\dagger} a_{js} + v^* \, a_{js}^{\dagger} b_{js}$$

$$+ \sum_{j=1}^{N-1} w \, b_{js}^{\dagger} a_{j+1,s} + w^* \, a_{j+1,s}^{\dagger} b_{js}$$
+ periodic-boundary terms, (12)

where v and w denote the intracell and intercell hopping amplitudes, respectively. The periodic-boundary contribution is  $w b_{Ns}^{\dagger} a_{1s} + w^* a_{1s}^{\dagger} b_{Ns}$  for PBC, accounting for the hopping between the rightmost and leftmost cells, and 0 for OBC. Equation (12) is invariant under the exchange of A and Bsublattices, a property known as *chiral symmetry*, which plays a crucial role in protecting the topological features of the

Since the spinful SSH model contains no many-body interactions, it admits a single-particle description. Its phases are characterized by the topological invariant w, known as the winding number. The system is a trivial insulator with w = 0when |v| > |w|, and a topological insulator with w = 1when |v| < |w|. At the critical point |v| = |w|, the band gap closes and the system becomes metallic; this point marks the topological phase transition, where the winding number is ill-defined.

The Hubbard interaction is described by the real-space Hamiltonian

<span id="page-2-3"></span>
$$H_{\text{Hubbard}} = \sum_{j=1}^{N} \left( U_A \, n_{Aj\uparrow} n_{Aj\downarrow} + U_B \, n_{Bj\uparrow} n_{Bj\downarrow} \right), \quad (13)$$

where  $U_A \geq 0$  and  $U_B \geq 0$  represent the on-site Coulomb repulsion on the A and B sublattices, respectively. In the conventional SSHH model one typically assumes  $U_A = U_B \equiv$ U. Here we also consider the more general unbalanced case  $U_A \neq U_B$ , which breaks the chiral symmetry and allows us to investigate how Hubbard interactions modify the topological features of the original SSH system. In the special case  $U_A = U_B \equiv U$  and v = w, the SSHH model reduces to the standard one-dimensional Hubbard model.

In the strong-coupling limit  $U \gg |v|, |w|$ , charge fluctuations are strongly suppressed and the system is effectively described by localized magnetic moments. When U/|v| is not large, the competition between hopping and Coulomb repulsion leads to nontrivial correlation effects.

In this work we focus on the transition between the trivial and topological insulating phases of the SSHH model. In particular, we investigate how robust the topological signatures of the original SSH model remain in the presence of particleparticle interactions, especially the unbalanced interactions that explicitly break chiral symmetry. This question is of interest because topological properties are generally expected to be robust against weak perturbations, yet their stability against many-body interactions has not been fully established.

### <span id="page-2-2"></span>A. Adiabatic evolution of the SSHH system

Treating  $H_{\text{SSH},\uparrow} + H_{\text{SSH},\downarrow}$  as  $H_0$  and  $H_{\text{Hubbard}}$  as  $H_1$  in (1), we can implement adiabatic quantum simulation on a quantum circuit. By repeatedly executing the circuit and measuring the final state, one can extract topological information about the SSHH model.

<span id="page-2-1"></span>Note that H(t) does not induce spin flips. More precisely,

<span id="page-2-0"></span>
$$[H, \hat{n}_s] = [H_{\text{SSH},s'}, \hat{n}_s] = [H_{\text{Hubbard}}, \hat{n}_s] = 0,$$
 (14)  
where the spin number operators are  $\hat{n}_s \equiv \sum_{\alpha \in \{A,B\}} \sum_j n_{\alpha j s}$ , with eigenvalues denoted by  $n_{\uparrow}$  and  $n_{\downarrow}$ .

Consequently, the Hilbert space of the many-body system decomposes into subspaces with fixed spin populations, corresponding to the simultaneous eigenspaces of  $\hat{n}_{\uparrow}$  and  $\hat{n}_{\downarrow}$ . These subspaces are invariant under H(t), and therefore the spin populations remain constant throughout the evolution governed by (1).

In realistic physical systems, however, spin-flip interactions are rarely completely absent. When the Hubbard interaction becomes sufficiently strong, the system may lower its energy by redistributing the spin populations, since spin flips can reduce on-site repulsion. From a mathematical perspective, each fixed-spin-population subspace labeled by  $(n_{\uparrow}, n_{\downarrow})$  therefore forms a dynamically disconnected sector of H(t). Different sectors, however, need not be separated by a nonzero energy gap.

This leads to an important caveat: the evolution described by (1) may represent an adiabatic process only in an artificial sense, namely when the dynamics are restricted to a fixed  $(n_{\uparrow}, n_{\downarrow})$  sector. In general, starting from the ground state of  $H_0$ , the final state obtained through (1) need not coincide with the true ground state of the full Hamiltonian H(T).

Fortunately, in the situations considered in this paper the relevant ground states correspond to configurations in which the lower energy levels are fully occupied and separated from higher levels by a finite excitation gap. Because the occupied lower levels are completely filled, there is no available phase space for spin flipping within them. Consequently, provided that the Hubbard interaction remains sufficiently weak relative to the band gap and that T is sufficiently large, the evolution governed by (1) represents a physically meaningful adiabatic process, and the final state faithfully approximates the ground state of the target Hamiltonian.

The initial state is taken to be the ground state of the spinful SSH model, obtained by filling the single-particle energy eigenstates in ascending order of energy. The corresponding many-body wavefunction is the Slater determinant of the occupied eigenstates, which equivalently can be written as

<span id="page-3-0"></span>
$$|\Psi_{n_{\uparrow},n_{\downarrow}}\rangle = a_{\epsilon_{1,\uparrow}}^{\dagger} \dots a_{\epsilon_{n_{\uparrow},\uparrow}}^{\dagger} a_{\epsilon_{1,\downarrow}}^{\dagger} \dots a_{\epsilon_{n_{\perp},\downarrow}}^{\dagger} |0\rangle , \qquad (15)$$

where  $a_{\epsilon_{j,s}}^{\dagger}$  is the creation operator for the j-th lowest single-particle eigenstate of  $H_{\mathrm{SSH},s}$ , and  $n=n_{\uparrow}+n_{\downarrow}$  is the total number of electrons.

Because of (14), the spectrum of H(t) is invariant under exchanging  $\uparrow$  and  $\downarrow$ , resulting in a two-fold degeneracy whenever  $n_{\uparrow} \neq n_{\downarrow}$ . The lattice with N unit cells has 4N single-particle degrees of freedom, corresponding to spin-up and spin-down electrons on the A and B sublattice sites of each cell.

Under PBC, the SSH model exhibits the spectrum of an insulator: the lower band, consisting of 2N energy levels below 0, is separated from the upper band of 2N levels above 0 by a nonzero excitation gap. At half filling (n=2N), the lower band is completely occupied, and the total many-body energy is separated from higher-energy states by a finite gap. Topological signatures can then be extracted from the many-body Berry phase of the n=2N system.

Under OBC, the SSH model with w = 0 exhibits a similar

band structure: the lower band of 2N levels below 0 remains separated from the upper band by a nonzero gap, so the n=2N state is again gapped from above.

In contrast, under OBC with w=1, the SSH model exhibits edge modes arising from its nontrivial topology. The bulk spectrum consists of 2N-2 levels below 0 and 2N-2 above 0, separated by a nonzero gap, together with four states near zero energy whose wavefunctions are localized at the boundaries (edges).

For the w=1 case under OBC, if we consider the n=2N+2 state, both the lower band and the four edge modes are fully occupied. The total energy of the system is then separated from higher-energy states by a finite gap, and topological signatures can be extracted from the edge localization of the sublattice polarization as the edge modes are all occupied.

For both the n=2N and n=2N+2 states considered above, we have  $n_{\uparrow}=n_{\downarrow}$ , so the spin-exchange degeneracy is lifted. Moreover, the total many-body energy is well separated from higher-energy states. Starting from these states, we can therefore safely perform the adiabatic simulation governed by (1); the resulting final state faithfully approximates the ground state of the target Hamiltonian, from which the topological signatures can be extracted.

### B. Topological characterization of the SSHH model

Topological properties of the SSHH model can be characterized from several many-body perspectives, including the many-body Berry phase, the charge center, bulk charge polarization, and edge-state localization. These features are expected to be robust against weak perturbations, whether arising from one-body or many-body effects. Consequently, the SSHH model should retain the topological characteristics of the SSH model provided that the Hubbard interaction remains sufficiently weak.

For a many-body state  $|\Psi\rangle$ , the total particle number is

<span id="page-3-1"></span>
$$n := \langle \Psi | \, \hat{n} \, | \Psi \rangle \equiv \langle \Psi | \, (\hat{n}_{\uparrow} + \hat{n}_{\perp}) \, | \Psi \rangle \,. \tag{16}$$

If  $|\Psi\rangle$  belongs to an eigenspace of  $\hat{n}$ , it describes an n-particle system. For a lattice with N unit cells, the filling ratio is defined as n/N. In our setting n is an integer and  $0 \le n/N \le 4$ . Writing the filling ratio in lowest terms as  $\tilde{n}/\tilde{N}$  (with  $\tilde{n}$  and  $\tilde{N}$  coprime), we define the position operator

$$\hat{X} := \sum_{j=1}^{N} \sum_{s \in \{\uparrow, \downarrow\}} j\left(a_{j,s}^{\dagger} a_{j,s} + b_{j,s}^{\dagger} b_{j,s}\right), \tag{17a}$$

$$\equiv \sum_{j=1}^{N} j \left( n_{Aj\uparrow} + n_{Aj\downarrow} + n_{Bj\uparrow} + n_{Bj\downarrow} \right), \tag{17b}$$

and the quantity

<span id="page-3-2"></span>
$$z_N[\tilde{n}/\tilde{N}] = \langle \Psi | \exp\left(\frac{i2\pi\tilde{N}}{N}\hat{X}\right) | \Psi \rangle, \qquad (18)$$

which depends only on the filling ratio in lowest terms [35], [36].

From  $z_N[\tilde{n}/\tilde{N}]$  one defines the many-body (more precisely, n-body) Berry phase [35]–[37]

<span id="page-3-3"></span>
$$\gamma_N = \operatorname{Im} \ln z_N[\tilde{n}/\tilde{N}],\tag{19}$$

which reduces to the Zak phase when n = 1, i.e., the Berry phase obtained from integrating the Berry connection over the Brillouin zone of a single-particle band [38].

Under PBC, the expectation value of the position operator is related to  $z_N$  by

$$\langle X_c \rangle = \frac{N}{2\pi \tilde{N}} \operatorname{Im} \ln z_N [\tilde{n}/\tilde{N}],$$
 (20)

which is commonly referred to as the charge center [35]–[37].

For the single-particle SSH model, the Zak phase is related to the winding number via  $\gamma_N = w\pi$  [39], and therefore serves as an indicator of nontrivial topology. For a general n-electron state occupying different energy levels of the same singleparticle SSH model, however, the n-body Berry phase defined above is not necessarily quantized. Nevertheless, when  $|\Psi\rangle$ takes the form of (15) and corresponds to the half-filled ground state occupying all negative-energy levels, the n-body Berry phase again becomes quantized as  $\gamma_N = w\pi$ . This occurs because the occupied states completely fill the Brillouin zone of the negative-energy band.

In the thermodynamic limit  $n, N \to \infty$  with fixed filling ratio, the many-body Berry phase is directly related to the charge polarization (dipole moment per unit cell),  $P_{\rm e}=e\langle X_c\rangle=$  $\frac{e}{2\pi}\gamma_N$ . The many-body Berry phase therefore corresponds to a physical observable and is, in principle, experimentally measurable.

Another hallmark of topological systems is the presence of robust edge states under OBC. In the single-particle SSH model, bulk-boundary correspondence implies that the trivial phase (w = 0) hosts no edge states, whereas the topological phase (w = 1) supports 2w zero-energy edge modes [32].

Edge modes can also be characterized from a many-body perspective. Under OBC we define the edge-occupation operator

$$\hat{n}_{\text{edge}} = \sum_{\alpha \in \{A,B\}} \sum_{s \in \{\uparrow,\downarrow\}} (\hat{n}_{\alpha 1s} + \hat{n}_{\alpha Ns}). \tag{21}$$

If single-particle energy levels are filled sequentially as in (15), the expectation value  $\langle \hat{n}_{\rm edge} \rangle$  behaves as a staircase function of n when edge modes are present: it varies smoothly with nexcept at certain critical values where it jumps abruptly as an edge mode becomes occupied.

Since edge states themselves are not directly measurable, it is useful to consider experimentally accessible quantities. For this purpose we introduce the sublattice polarization operator

$$P_j^{\mathrm{e}} = \sum_{s \in \{\uparrow, \downarrow\}} (n_{Ajs} - n_{Bjs}). \tag{22}$$

Although  $\langle P_i^{\rm e} \rangle$  vanishes for the half-filled ground state (n =2N), it becomes nonzero at the boundary cells for the ground state with two additional electrons (n = 2N + 2) if edge modes are present. This occurs because the edge states are fully occupied in this configuration. The spatial profile of  $\langle P_i^{\rm e} \rangle$ under OBC therefore provides an experimentally meaningful indicator of the system's topological character.

These topological signatures of the SSH model are expected to persist under weak perturbations, including the Hubbard interaction introduced in the SSHH model. When the interaction strength is sufficiently small compared with the hopping amplitudes, the single-particle description remains a good approximation. In this regime the Hubbard term can be treated within a mean-field approximation. For the state (15), the effective Hubbard Hamiltonian becomes

$$H_{\rm Hubbard}^{\rm eff} \qquad (23)$$

$$= \sum_{j=1}^{N} \left[ U_A \left( \langle n \rangle_{Aj\uparrow} \, n_{Aj\downarrow} + n_{Aj\uparrow} \, \langle n \rangle_{Aj\downarrow} + \langle n \rangle_{Aj\uparrow} \, \langle n \rangle_{Aj\downarrow} \right) + U_B \left( \langle n \rangle_{Bj\uparrow} \, n_{Bj\downarrow} + n_{Bj\uparrow} \, \langle n \rangle_{Bj\downarrow} + \langle n \rangle_{Bj\uparrow} \, \langle n \rangle_{Bj\downarrow} \right). \right]$$
which effectively adds on-site one-particle potentials to

which effectively adds on-site one-particle potentials to  $H_{\text{SSH},\uparrow}$  and  $H_{\text{SSH},\downarrow}$ .

Therefore, as long as the difference  $\Delta U := U_B - U_A$ remains sufficiently small, the many-body Berry phase and the edge sublattice polarization inherited from the SSH model remain intact even though chiral symmetry is explicitly broken when  $\Delta U \neq 0$ . By examining these indicators for the n=2Nand n = 2N + 2 ground states, we can demonstrate that the nontrivial topology persists under weak perturbations but eventually breaks down when interactions become sufficiently strong.

As discussed in Sec. III-A, the adiabatic evolution governed by (1) may fail to produce the true ground state of the target Hamiltonian because the dynamics are restricted to a fixedspin-population sector  $(n_{\uparrow}, n_{\perp})$ . To ensure that the evolution corresponds to a physically meaningful process, the Hubbard interaction must remain weak enough so that spin flips do not occur during the evolution.

For the n = 2N and n = 2N + 2 ground states, this requirement translates into a constraint on the energy cost of flipping a spin, characterized by  $\max\{U_A, U_B, |\Delta U|\}$ . The relevant energy scale is the single-particle band gap of the

$$\Delta_{\rm gap} = 2 \min_k |v + e^{-ik}w| = 2 \min\{|v + w|, |v - w|\}.$$
 (24) For the  $n = 2N$  ground state, the energy cost of a spin flip must be smaller than  $\Delta_{\rm gap}$ , whereas for the  $n = 2N+2$  ground state it must be smaller than  $\Delta_{\rm gap}$ .

state it must be smaller than  $\Delta_{\rm gap}/2$ . A sufficient condition is therefore estimated as

<span id="page-4-1"></span>
$$\max\{U_A, U_B, |\Delta U|\} < \min\{|v + w|, |v - w|\}. \tag{25}$$

In Sec. V, we perform numerical simulations of the adiabatic evolution in this parameter regime and examine the resulting topological signatures of the SSHH model.<sup>1</sup>

# IV. QUANTUM-CIRCUIT IMPLEMENTATION OF ADIABATIC SIMULATION

To implement adiabatic quantum simulation on a gatebased quantum computer, the continuous evolution generated by H(t) must be approximated by a sequence of quantum gates. For the simulation to be efficient, the total qubit number

<span id="page-4-0"></span><sup>1</sup>When  $\Delta U = 0$ , evolving the half-filled (n = 2N) ground state of  $H_0$ under (1) preserves the same topological characteristics in both the manybody Berry phase and the sublattice polarization profile regardless of the value of  $U_A = U_B$ , since the Hubbard interaction merely adds a uniform onsite potential. However, beyond the condition in (25), the final state obtained from the adiabatic evolution can deviate considerably from the true ground state of the target Hamiltonian.

and gate complexity should scale at most polynomially with the size of the simulated system. Likewise, the measurement protocol must remain scalable, requiring only polynomially many shots, and the classical post-processing used to estimate observable expectation values from measurement outcomes must also be computationally efficient. In this section, we develop such a framework for the SSHH model.

We begin by applying the Jordan–Wigner transformation to express the second-quantized fermionic Hamiltonian in terms of Pauli operators. Using the resulting qubit representation, we then construct explicit quantum circuits for preparing the initial state and implementing the time evolution generated by H(t), both with polynomial gate depth. Finally, we describe a practical measurement protocol and post-processing procedure for extracting the many-body Berry phase and the spatial profile of the sublattice polarization.

# A. Jordan-Wigner transformation

We employ the Jordan–Wigner transformation to express the fermionic algebra in terms of Pauli matrices. We denote the Pauli-Z eigenstates as  $|Z+\rangle \equiv |0\rangle$  and  $|Z-\rangle \equiv |1\rangle$ , corresponding to an empty and an occupied fermionic state,

More explicitly, for a system consisting of N unit cells, the Jordan–Wigner transformation maps the fermionic creation and annihilation operators according to

$$a_{j\uparrow} \to \frac{1}{2} (X_{Aj\uparrow} + iY_{Aj\uparrow}) \prod_{k=1}^{j-1} Z_{Ak\uparrow} Z_{Bk\uparrow},$$
 (26a)

$$b_{j\uparrow} \to \frac{1}{2} (X_{Bj\uparrow} + iY_{Bj\uparrow}) Z_{Aj\uparrow} \prod_{k=1}^{j-1} Z_{Ak\uparrow} Z_{Bk\uparrow},$$
 (26b)

$$b_{j\uparrow} \to \frac{1}{2} (X_{Bj\uparrow} + iY_{Bj\uparrow}) Z_{Aj\uparrow} \prod_{k=1}^{j-1} Z_{Ak\uparrow} Z_{Bk\uparrow}, \qquad (26b)$$

$$a_{j\downarrow} \to \frac{1}{2} (X_{Aj\downarrow} + iY_{Aj\downarrow}) P_{\uparrow}^{z} \prod_{k=1}^{j-1} Z_{Ak\downarrow} Z_{Bk\downarrow}, \qquad (26c)$$

$$b_{j\downarrow} \to \frac{1}{2} (X_{Bj\downarrow} + iY_{Bj\downarrow}) P_{\uparrow}^z Z_{Aj\downarrow} \prod_{k=1}^{j-1} Z_{Ak\downarrow} Z_{Bk\downarrow},$$
 (26d)

where the product is understood to be the identity for j = 1. The operator  $P_s^z$  represent the fermionic parity of the entire spin-s sector, defined as

<span id="page-5-1"></span>
$$P_s^z := \prod_{k=1}^N Z_{Aks} Z_{Bks}.$$
 (27)

Under this mapping, the local number operators take the simple diagonal form

$$n_{Ajs} \equiv a_{js}^{\dagger} a_{js} = (I_{Ajs} - Z_{Ajs})/2,$$
 (28a)

<span id="page-5-0"></span>
$$n_{Bis} \equiv \vec{b}_{is}^{\dagger} b_{is} = (I_{Bis} - Z_{Bis})/2,$$
 (28b)

whereas the intracell, intercell, and boundary hopping terms are mapped as

$$b_{js}^{\dagger}a_{js} = (X_{Bjs} - iY_{Bjs})(X_{Ajs} + iY_{Ajs})/4,$$
 (29a)

$$b_{js}^{\dagger} a_{j+1s} = (X_{Bjs} - iY_{Bjs})(X_{Aj+1s} + iY_{Aj+1s})/4,$$
 (29b)

$$b_{Ns}^{\dagger} a_{1s} = (X_{BNs} - iY_{BNs})(X_{A1s} + iY_{A1s})P_s^z/4.$$
 (29c)

Substituting (29) into (12), we obtain the Pauli-operator representation of the spinful SSH Hamiltonian,

$$H_{\text{SSH},s} = -\frac{1}{2} \Big\{ \sum_{i=1}^{N} \Big[ \text{Re}(v) \big( X_{Ajs} X_{Bjs} + Y_{Ajs} Y_{Bjs} \big) \Big] \Big\}$$

<span id="page-5-2"></span>
$$+\operatorname{Im}(v)\left(X_{Ajs}Y_{Bjs}-Y_{Ajs}X_{Bjs}\right)\Big]$$

$$+\sum_{j=1}^{N-1}\Big[\operatorname{Re}(w)\left(X_{Aj+1s}X_{Bjs}+Y_{Aj+1s}Y_{Bjs}\right)$$

$$+\operatorname{Im}(w)\left(X_{Aj+1s}Y_{Bjs}-Y_{Aj+1s}X_{Bjs}\right)\Big]$$
+ periodic-boundary terms \(\begin{array}{c}\), (30)

where the periodic-boundary contribution is 0 for OBC and for PBC given by

$$\left[\operatorname{Re}(w)\left(X_{A1s}X_{BNs} + Y_{A1s}Y_{BNs}\right) + \operatorname{Im}(w)\left(X_{A1s}Y_{BNs} - Y_{A1s}X_{BNs}\right)\right]P_s^z, \tag{31}$$

accounting for the hopping between the 1st and N-th cells.

Similarly, substituting (28) into (13), we obtain the Paulioperator representation of the Hubbard Hamiltonian,

$$H_{\text{Hubbard}} = \frac{U}{4} \sum_{s} \sum_{j=1}^{N} \left[ (I_{Aj\uparrow} - Z_{Aj\uparrow})(I_{Aj\downarrow} - Z_{Aj\downarrow}) + (I_{Bj\uparrow} - Z_{Bj\uparrow})(I_{Bj\downarrow} - Z_{Bj\downarrow}) \right]. \quad (32)$$

In the same manner, we consider the propagator on the  $\ell$ -th interval defined in (7),

<span id="page-5-4"></span><span id="page-5-3"></span>
$$U_{\ell} = U_{\ell}^{\text{Hubbard}} U_{\ell}^{\text{SSH},\uparrow} U_{\ell}^{\text{SSH},\downarrow}.$$
 (33)

<span id="page-5-5"></span>Substituting (30) into (33) and applying the first-order Trotter decomposition, the spinful SSH propagator in the Paulioperator representation is given by

$$U_{\ell}^{SSH,s}$$

$$= \exp\left(\frac{\delta t \operatorname{Re}(w)}{-2i} (X_{A1s} X_{BLs} + Y_{A1s} Y_{BLs}) P_{s}^{z}\right)$$

$$\times \exp\left(\frac{\delta t \operatorname{Re}(w)}{-2i} \sum_{i=1}^{N-1} (X_{Ai+1s} X_{Bis} + Y_{Ai+1s} Y_{Bis})\right)$$

$$\times \exp\left(\frac{\delta t \operatorname{Re}(v)}{-2i} \sum_{i=1}^{N} (X_{Ais} X_{Bis} + Y_{Ais} Y_{Bis})\right)$$

$$\times \exp\left(\frac{\delta t \operatorname{Im}(w)}{-2i} (X_{A1s} Y_{BLs} - Y_{A1s} X_{BLs}) P_{s}^{z}\right)$$

$$\times \exp\left(\frac{\delta t \operatorname{Im}(w)}{-2i} \sum_{i=1}^{N-1} (X_{Ai+1s} Y_{Bis} - Y_{Ai+1s} X_{Bis})\right)$$

$$\times \exp\left(\frac{\delta t \operatorname{Im}(v)}{-2i} \sum_{i=1}^{N} (X_{Ais} Y_{Bis} - Y_{Ais} X_{Bis})\right) + O(\delta t^{2})$$

where the error term represents the approximation inaccuracy of the first-order Trotter decomposition.

Within a fixed  $(n_{\uparrow}, n_{\downarrow})$  sector, the whole system has a definite parity. Consequently,  $P_s^z$  acts as a c-number and can be replaced by  $P_s^z \to (-1)^{n_s}$ . All the terms involving  $P_s^z$ appear only in PBC and are absent in OBC.

Similarly, substituting (32) into (33), we have

$$= \exp\left(\frac{\delta t U_A}{4i} \left(\frac{2\ell - 1}{2L}\right) \sum_{j=1}^{N} (I_{Aj\uparrow} - Z_{Aj\uparrow}) (I_{Aj\downarrow} - Z_{Aj\downarrow})\right) \times \exp\left(\frac{\delta t U_B}{4i} \left(\frac{2\ell - 1}{2L}\right) \sum_{j=1}^{N} (I_{Bj\uparrow} - Z_{Bj\uparrow}) (I_{Bj\downarrow} - Z_{Bj\downarrow})\right) + O(\delta t^2).$$
(35)

where the error term again arises from the first-order Trotter decomposition.

Under this mapping, each fermionic mode  $(\alpha, j, s)$  is associated with a qubit whose computational basis encodes its occupation number. Consequently, the many-body wavefunction on a lattice with N unit cells can be represented using 4N qubits, enumerated as  $q \in \{0, \dots, 4N-1\}$ . The presence (absence) of a fermion in the mode  $(\alpha, j, s)$  is mapped to the qubit state  $|1\rangle$  ( $|0\rangle$ ) of the corresponding qubit, for example, according to the explicit prescription

$$q(A, j, \uparrow) = 2j - 2,$$
  $q(A, j, \downarrow) = 2N + 2j - 2,$  (36a)  
 $q(B, j, \uparrow) = 2j - 1,$   $q(B, j, \downarrow) = 2N + 2j - 1,$  (36b)

where the 4N qubits are partitioned into two segments corre-

sponding to the spin-up and spin-down sectors. For later use, we define three types of two-qubit Pauli-

rotation gates acting on the qubit pair (i, j):

$$R_{i,j}(\theta) = \exp\left(\left(X_i X_j + Y_i Y_j\right) \theta / 2i\right),\tag{37}$$

$$R_{i,j}(\theta) = \exp\left(\left(X_i X_j + Y_i Y_j\right) \theta/2i\right), \tag{37}$$

$$G_{i,j}(\theta) = \exp\left(\left(X_i Y_j - Y_i X_j\right) \theta/2i\right), \tag{38}$$

$$CP_{i,j}(\theta) = \exp\left(\left(I_i - Z_i\right)\left(I_j - Z_j\right)\theta/4i\right). \tag{39}$$

In principle, each of these gates can be decomposed into two CNOT gates and at most four single-qubit rotation gates.

With these definitions, the gate-level propagator for the  $\ell$ -th time interval is implemented as

$$U_{\ell} = \prod_{j=0}^{N-1} CP_{2j,2j+2N}(\phi_{A,\ell}) CP_{2j+1,2j+2N+1}(\phi_{B,\ell})$$

$$\times \prod_{j=0,2N} R_{j+2N-1,j}((-1)^{n_s}\theta_w) G_{j+2N-1,j}((-1)^{n_s}\vartheta_w)$$

$$\times \prod_{j=0}^{2N-2} R_{j,j+1}(\theta_w) R_{j+2N,j+2N+1}(\theta_w)$$

$$\times \prod_{j=1}^{2N-2} G_{j,j+1}(\vartheta_w) G_{j+2N,j+2N+1}(\vartheta_w)$$

$$\times \prod_{j=0}^{4N-1} R_{j,j+1}(\theta_v) G_{j,j+1}(\vartheta_v), \tag{40}$$

where

$$\theta_v = -\delta t \operatorname{Re}(v), \qquad \qquad \theta_v = -\delta t \operatorname{Im}(v), \qquad (41a)$$

$$\theta_w = -\delta t \operatorname{Re}(w), \qquad \qquad \theta_w = -\delta t \operatorname{Im}(w), \qquad (41b)$$

$$\theta_{v} = -\delta t \operatorname{Re}(v), \qquad \theta_{v} = -\delta t \operatorname{Im}(v), \qquad (41a)$$

$$\theta_{w} = -\delta t \operatorname{Re}(w), \qquad \theta_{w} = -\delta t \operatorname{Im}(w), \qquad (41b)$$

$$\phi_{A,\ell} = \delta t U_{A} \left(\frac{2\ell - 1}{2L}\right), \quad \phi_{B,\ell} = \delta t U_{B} \left(\frac{2\ell - 1}{2L}\right). \quad (41c)$$

The circuit representation of  $U_{\ell}$  is shown in Fig. 1. Under PBC, each propagator contains  $4N R_{ij}$  gates,  $4N G_{ij}$  gates, and 2N  $CP_{ij}$  gates, giving a total of 10N two-qubit Paulirotation gates. Under OBC, the corresponding counts are 4N-2, 4N-2, and 2N, respectively, giving a total of 10N-4 gates. Since each such gate decomposes into two CNOT gates and four single-qubit rotation gates, implementing L propagators requires 2L(10N) CNOT gates and 4L(10N)single-qubit rotation gates for PBC, and 2L(10N-4) CNOT gates and 4L(10N-4) single-qubit rotation gates for OBC. Therefore, the total gate count scales as O(NL).

<span id="page-6-0"></span>![](_page_6_Figure_16.jpeg)

<span id="page-6-2"></span>Fig. 1: Quantum circuit implementation of a single Trotter step for the model with N=2 unit cells. Qubits  $q_0$  to  $q_3$  and  $q_4$  to  $q_7$  encode the spin-up and spin-down orbitals, respectively. The circuit is partitioned into three main stages: (1) The red shaded part represents the intra-cell hopping operations, utilizing R and G gates to realizing the real and imaginary components. (2) The blue shaded part represents the inter-cell hopping operations. The R and G gates in the gray boxes represent the periodic-boundary hopping terms, present only in PBC. (3) The green shaded part represents the Hubbard interaction, implemented via Controlled-Phase gates, which entangle the corresponding spin-up and spin-down orbitals at each site.

### B. Initial state preparation

To prepare the fermionic many-body ground state of the SSH model defined in (15) in the lattice-site basis, we employ the state-preparation framework introduced in Ref. [40] and implemented in the Python library OpenFermion [41], [42]. Within this framework, the target state can be prepared by a quantum circuit whose depth scales polynomially with the number of qubits. Below, we explain how this preparation routine is integrated into our circuit construction.

Suppose that the target state, obtained either analytically or numerically, takes the form given in (15). Our goal is to prepare this state in the lattice-site basis defined by the creation operators  $a_{js}^{\dagger}$  and  $b_{js}^{\dagger}$  on quantum circuits. For this purpose, we express the single-particle creation operator  $a_{\epsilon s}^{\dagger}$  in terms of the lattice-site operators as

$$a_{\epsilon_{\jmath,s}}^{\dagger} = \sum_{i} \left( \langle A, j, s | \epsilon_{\jmath,s} \rangle \, a_{js}^{\dagger} + \langle B, j, s | \epsilon_{\jmath,s} \rangle \, b_{js}^{\dagger} \right), \quad (42)$$

where  $|\epsilon_s\rangle$  denotes the single-particle energy eigenstate in the spin-s sector. Substituting this transformation into (15), the many-body state becomes

$$|\Psi_{n_{\uparrow},n_{\downarrow}}\rangle = \prod_{s \in \{\uparrow,\downarrow\}} \prod_{\epsilon=\epsilon_{1}}^{\epsilon_{n_{s}}} \left[ \sum_{j} \left( \langle A, j, s | \epsilon_{\jmath,s} \rangle \, a_{js}^{\dagger} + \langle B, j, s | \epsilon_{\jmath,s} \rangle \, b_{js}^{\dagger} \right) \right] |0\rangle \,. \tag{43}$$

In the logical-qubit representation, the index  $(\epsilon_{i,s})$  appearing in (15) is mapped to the qubit index  $p \equiv p(\epsilon_{\eta,s})$  as

<span id="page-6-1"></span>
$$p(\epsilon_{j,\uparrow}) = j - 1, \qquad p(\epsilon_{j,\downarrow}) = 2N + j - 1, \tag{44}$$

following the same indexing convention used by Open-Fermion.

Although a direct preparation of the state in (43) appears prohibitively complicated for large systems, the method of Ref. [40] provides an efficient construction by exploiting the fermionic anticommutation relations. The key observation is that the antisymmetric structure imposed by the Pauli exclusion principle allows the state-preparation problem to be reduced to a sequence of structured matrix decompositions. As a result, both the classical preprocessing cost and the quantum circuit construction depend only on the relevant rectangular transformation matrix

$$Q_{pq} = \langle \epsilon_{j,s} | \alpha, j, s \rangle, \qquad p \in S_p, \quad q \in S_q, \tag{45}$$

where  $S_p$  indexes the occupied single-particle modes and  $S_q$  indexes the lattice-site modes used in the qubit representation. Accordingly,  $|S_p|$  is the number of fermions, whereas  $|S_q|$  is the number of lattice-site modes included in the mapping, which is typically the full single-particle Hilbert-space dimension of the system.

This state-preparation procedure first applies n X gates to prepare the occupation pattern in (15). It then implements the basis transformation from (15) to (43) using  $|S_p|(|S_q|-|S_p|)$   $G_{ij}$  gates. For a system with  $N_f$  single-particle degrees of freedom, the half-filled case  $|S_p| = N_f/2$  and  $|S_q| = N_f$  therefore requires  $N_f^2/4$   $G_{ij}$  gates. Since each  $G_{ij}$  gate is decomposed into two CNOT gates and four single-qubit rotation gates, the corresponding gate counts are  $N_f^2/2$  CNOT gates and  $N_f^2$  single-qubit rotation gates, in addition to the  $N_f/2$  initial X gates. Equivalently, if the initial X gates are included among the single-qubit gates, the total number of single-qubit gates is  $N_f(N_f+1/2)$ .

In the SSH model considered here, the two spin sectors are independent in the initial-state preparation. The transformations in the spin-up and spin-down sectors can therefore be implemented in parallel. Although the full system contains 4N single-particle degrees of freedom, each spin sector contains only 2N lattice-site modes. At half filling,  $|S_p|_{\uparrow}=n_{\uparrow}=N,$   $|S_p|_{\downarrow}=n_{\downarrow}=N,$  and  $|S_q|_{\uparrow}=|S_q|_{\downarrow}=2N.$  Thus, each spin sector requires  $N^2$   $G_{ij}$  gates, corresponding to  $2N^2$  CNOT gates and  $4N^2$  single-qubit rotation gates, plus N initial X gates. Since the two spin sectors can be executed in parallel, the circuit depth is governed by the cost of a single sector rather than by the sum over both sectors. Consequently, the initial-state preparation remains polynomial in the number of unit cells.

# C. Measurement protocol

Although the Berry phase is not itself a direct physical observable, it can be inferred from the expectation value of the periodic-boundary-condition-resolved position operator defined in (17). By applying the Jordan–Wigner transformation in (26) together with the qubit-index mapping in (36), this operator becomes

$$\hat{X} = \sum_{j=0}^{N-1} (j+1) \left( n_{2j} + n_{2j+1} + n_{2N+2j} + n_{2N+2j+1} \right), \tag{46}$$

where  $n_j$  denotes the number operator on the j-th logical qubit. Since  $\hat{X}$  is diagonal in the computational basis, all of its components commute with the standard bit-string measurement. Therefore, the quantity required for the Berry phase

can be estimated from a single set of computational-basis measurements.

Let  $|b\rangle \equiv |b_0 \dots b_{4N-1}\rangle$  denote a computational-basis bit string, where  $b_j \in \{0,1\}$ , and define  $p(b) = |\langle b|\psi\rangle|^2$ . Using the spectral decomposition in the computational basis, the quantity  $z_N[\tilde{n}/\tilde{N}]$  in (18) can be written as

$$z_{N}[\tilde{n}/\tilde{N}] = \sum_{bb'} \langle \psi|b\rangle\langle b| e^{\frac{i2\pi\tilde{N}}{N}\hat{X}} |b'\rangle\langle b'|\psi\rangle$$

$$= \sum_{b} p(b)e^{\frac{i2\pi\tilde{N}}{N}\langle b|\hat{X}|b\rangle}, \tag{47}$$

with  $\langle b|n_j|b\rangle=b_j$ . The Berry phase is then obtained from the complex phase of  $z_N[\tilde{n}/\tilde{N}]$ , as specified in (19).

In practice, however, it is unnecessary, and generally inefficient for large systems, to store all sampled bit strings and reconstruct the full empirical distribution p(b). Instead,  $z_N[\tilde{n}/\tilde{N}]$  can be estimated directly by averaging the phase factor associated with each measurement outcome. Let  $|b^{(m)}\rangle$  denote the m-th bit string obtained from the quantum device. For each measurement outcome, we compute

$$z_N^{(m)}[\tilde{n}/\tilde{N}] = e^{\frac{i2\pi\tilde{N}}{N}\langle b^{(m)}|\hat{X}|b^{(m)}\rangle}.$$
 (48)

After M measurements, the estimator of  $z_N[\tilde{n}/\tilde{N}]$  is given by

$$\bar{z}_N[\tilde{n}/\tilde{N}] = \frac{1}{M} \sum_{m=1}^M z_N^{(m)}[\tilde{n}/\tilde{N}].$$
 (49)

Accordingly, the Berry phase is estimated as

$$\bar{\gamma} = \operatorname{Im} \ln \bar{z}_N [\tilde{n}/\tilde{N}],$$
 (50)

where the result is understood modulo  $2\pi$ .

The electron polarization can be estimated from the same computational-basis measurements. For each sampled bit string  $|b^{(m)}\rangle$ , the local sublattice-resolved charge imbalance at unit cell j is evaluated as

$$p_j^{e,(m)} = \langle b^{(m)} | n_{2j} + n_{2N+2j} - n_{2j+1} - n_{2N+2j+1} | b^{(m)} \rangle,$$
(51)

for j = 0, ..., N - 1. The corresponding estimator is then

$$\bar{p}_{j}^{e} = \frac{1}{M} \sum_{m=1}^{M} p_{j}^{e,(m)}.$$
 (52)

Thus, both the Berry phase and the electron polarization can be obtained from standard computational-basis measurements. The Berry phase is extracted by classically post-processing the sampled bit strings into the complex estimator  $\bar{z}_N[\tilde{n}/\tilde{N}]$ , whereas the electron polarization is obtained by directly averaging the corresponding occupation-number imbalance. In contrast to generic observables such as the ground-state energy, whose Pauli-string decomposition generally requires measurements in multiple bases, the present quantities are diagonal in the computational basis after the Jordan–Wigner transformation. Consequently, their estimation requires neither additional basis rotations nor separate measurement settings, and therefore introduces no extra circuit depth beyond that required for state preparation and adiabatic evolution.

### V. NUMERICAL SIMULATIONS

<span id="page-7-0"></span>We implement and simulate the quantum adiabatic simulation for the one-dimensional SSHH model within the Oiskit

framework [\[43\]](#page-11-2). To resolve the relevant topological signatures within the limits of classical circuit simulation, we use N = 6 unit cells, corresponding to 24-qubit circuits. The Berry phase is evaluated for the half-filled periodic SSHH ground state with 12 electrons, whereas the sublattice polarization is computed for the open-boundary ground state with 14 electrons. The results show that the SSH topological signatures persist under weak Hubbard interactions but break down under sufficiently strong chiral-symmetry-breaking interactions.

# *A. Validation of the Trotter decomposition*

Before evaluating the topological observables, we first examine the validation of the first-order Trotter decomposition with respect to the total evolution time T and the number of propagator intervals L. This analysis fixes the simulation setting used in the subsequent Berry-phase and sublatticepolarization calculations.

We assess convergence using two fidelity criteria. Without Hubbard interactions, the circuit-simulated final state is compared with the analytically predicted final state. With Hubbard interactions, where no analytical reference state is readily available, we instead compare final states obtained from adjacent values of the propagator interval number L.

For both assessments, we fix the hopping amplitudes at v = 0.5 and w = 1.5. The total evolution time is varied over T = 1, 5, 15, and 80, while the number of intervals is chosen from

$$L \in \{1, 10, 20, 30, 40, 50, 60, 80, 100, 120, 150\}.$$

Note that, for the interacting case, the fidelity obtained at L = 1 is computed with the initial state, severing as a reference baseline.

The results for the noninteracting case U = 0 and the interacting case U = 1 are shown in the top and bottom panels of Fig. [2,](#page-8-0) respectively. For short evolution times, T = 1 and T = 5, the final-state fidelity rapidly approaches unity, with near-ideal values already achieved at relatively small L. For the intermediate evolution time T = 15, convergence with respect to L is slower; nevertheless, the fidelity increases over most of the sampled range and becomes close to unity at L = 150. In contrast, for the long evolution time T = 80, the chosen range of interval numbers is clearly insufficient, as the fidelity remains substantially below unity even at the largest value of L considered.

These results indicate that longer evolution times require a larger number of intervals to maintain the accuracy of the adiabatic simulation. We also observe that the interacting case appears to approach unity with fewer intervals than the noninteracting case. However, this behavior should not be interpreted as evidence that the interacting simulation is more accurate in an absolute sense. In the interacting case, the fidelity is evaluated between final states obtained at two adjacent entries in the chosen sequence of interval numbers, rather than between the simulated state and an exact reference state. Thus, this quantity only diagnoses the stability of the final state under increasing L; a value close to unity indicates that the result is insensitive to further refinement of L, but

<span id="page-8-0"></span>![](_page_8_Figure_9.jpeg)

Fig. 2: Final-state fidelity for the case with U = 0 (top) and the case with U = 1 (bottom).

does not by itself guarantee proximity to the ideal adiabatically evolved state.

Based on these observations, we choose T = 1 and L = 40 as the simulation setting for the subsequent numerical experiments. This choice provides high fidelity in the noninteracting benchmark and stable behavior in the interacting case, while keeping the circuit depth within a tractable range.

# *B. Many-body Berry phase*

Next, we investigate the Berry phase of the SSHH model for different values of ∆U. We fix U<sup>A</sup> = 0.01 and the intracell hopping amplitude v = 1, while varying the intercell hopping amplitude w from 0 to 2 in increments of 0.25. In addition, we include two sampling points near the transition, w = 0.99 and w = 1.01, while omitting the critical point w = 1, where the Berry phase is ill-defined because the gap closes.

The results are shown in the top panel of Fig. [3.](#page-9-0) For ∆U = 0, the model reduces to the SSH model, and the Berry phase exhibits a sharp step-like change across w = 1. As ∆U is increased from 0.0003 to 0.3, this step-like signature is progressively diminished due to the breaking of chiral symmetry. Moreover, in the both topologically trivial and nontrivial regime, the Berry phase is no longer pinned to the quantized value γ = 1 and γ = 0, respectively, in our phase convention. This behavior indicates that the topological phase breaks down when the symmetry protection responsible for Berry-phase quantization is removed.

To further assess the robustness of the Berry phase in the weak ∆U regime, we fix w = 1.5 and extend the set of ∆U toward smaller values. The results are shown in the bottom panel of Fig. [3.](#page-9-0) In the present setting, the Berry phase remains essentially unchanged for ∆U ≲ 10<sup>−</sup><sup>2</sup> . This indicates that the interaction imbalance acts only as a weak perturbation in this regime, so that the topological signature is preserved. In contrast, for ∆U ≳ 10<sup>−</sup><sup>1</sup> , the Berry phase is shifted appreciably away from the quantized value γ = 0, indicating that the symmetry-breaking perturbation is no longer negligible. This

<span id="page-9-0"></span>![](_page_9_Figure_0.jpeg)

Fig. 3: Many-body Berry phase for different values of ∆U.

trend agrees with the theoretical expectation that Berry-phase quantization is robust only when the unbalanced Hubbard interaction remains sufficiently weak.

# *C. Spatial profile of sublattice polarization*

We finally examine the spatial profile of the sublattice polarization in the open-boundary SSHH model and its response to ∆U. Following the Berry-phase analysis, we use the same sets of U<sup>A</sup> and ∆U values, but set v = 0.1 and w = 1.0 to enhance the dimerization so as to better separate edge and bulk polarization features in 6-unit-cell systems.

As shown in the top panel of Fig. [4,](#page-9-1) for ∆U = 0, the spatial profile of the sublattice polarization is localized at the two edges, while the bulk remains essentially unpolarized. When ∆U ̸= 0, this profile is progressively distorted, and a finite polarization offset develops in the bulk.

To resolve the small-imbalance regime more clearly, we further examine smaller values of ∆U, focusing on the response at the first unit cell, namely the edge on the A sublattice. The results are shown in the bottom panel of Fig. [4.](#page-9-1) For ∆U ≲ 0.01, the edge-localized feature is retained, with the corresponding sublattice polarization remaining close to unity. For larger values, ∆U ≳ 0.5, the sublattice polarization gradually degrades as ∆U increase.

These results indicate that the spatial profile of the sublattice polarization retains its edge-localized character as long as the interaction imbalance ∆U is sufficiently weak.

# VI. SUMMARY

We have investigated the topological properties of the SSHH model in detail, focusing on the many-body Berry phase under PBC and the spatial profile of sublattice polarization under OBC. We also examined the feasibility and limitations of obtaining the true many-body ground state of the SSHH model via quantum adiabatic simulation, starting from an initial state constructed from pre-solved one-particle states of the SSH model. Subtle issues such as spin-population conservation and chiral-symmetry breaking were analyzed and clarified.

<span id="page-9-1"></span>![](_page_9_Figure_10.jpeg)

Fig. 4: Sublattice polarization for different values of ∆U.

Using the Jordan–Wigner transformation, the many-body wavefunction of the SSHH model on a lattice with N unit cells can be mapped to a quantum register of 4N qubits. Based on a first-order Trotter decomposition, we constructed explicit quantum circuits implementing the adiabatic time evolution, together with a systematic circuit-based preparation of the required initial states. We also proposed practical measurement protocols for extracting the many-body Berry phase and the spatial profile of sublattice polarization.

We performed classical numerical simulations of the proposed circuits using 24 qubits for a 6-unit-cell SSHH model. Under PBC, the many-body Berry phase retains the expected step-like behavior, pinned to 0 and π, as long as the imbalanced Hubbard interaction ∆U remains weak. As ∆U exceeds a threshold, this topological signature begins to break down as the Berry phase deviates from these quantized values. Under OBC, the spatial profile of the sublattice polarization remains edge-localized for weak ∆U, but this localization gradually degrades once ∆U crosses the threshold.

Taken together, these results demonstrate that adiabatic quantum simulation provides a practical route for probing physically meaningful observables in interacting many-body systems beyond ground-state energy estimation. Although the present results are obtained from classical simulations of small systems (N = 6), the required qubit number (4N), gate complexity, measurement shots, and classical pre- and post-processing costs all scale polynomially with system size. Our work therefore establishes a proof-of-concept framework for studying nontrivial topological and spatial properties of interacting many-body systems via adiabatic quantum simulation, pointing toward potential implementations for large-size systems on future large-scale quantum computers.

# AI TOOL USAGE STATEMENT

ChatGPT and Gemini were used for text refinement and Python code support. All content, code, and results were reviewed and validated by the authors, who take full responsibility for the final work.

# REFERENCES

- <span id="page-10-0"></span>[1] J. Preskill, "Quantum computing in the NISQ era and beyond," *Quantum*, vol. 2, p. 79, 2018.
- <span id="page-10-1"></span>[2] K. Bharti, A. Cervera-Lierta, T. H. Kyaw, T. Haug, S. Alperin-Lea, A. Anand, M. Degroote, H. Heimonen, J. S. Kottmann, T. Menke, W.-K. Mok, S. Sim, L.-C. Kwek, and A. Aspuru-Guzik, "Noisy intermediatescale quantum algorithms," *Reviews of Modern Physics*, vol. 94, p. 015004, 2022.
- <span id="page-10-2"></span>[3] U. Schollwock, "The density-matrix renormalization group," ¨ *Reviews of Modern Physics*, vol. 77, pp. 259–315, 2005.
- [4] M. Troyer and U.-J. Wiese, "Computational complexity and fundamental limitations to fermionic quantum monte carlo simulations," *Physical Review Letters*, vol. 94, p. 170201, 2005.
- <span id="page-10-3"></span>[5] J. Kempe, A. Kitaev, and O. Regev, "The complexity of the local hamiltonian problem," *SIAM Journal on Computing*, vol. 35, no. 5, pp. 1070–1097, 2006.
- <span id="page-10-4"></span>[6] R. P. Feynman, "Simulating physics with computers," *International Journal of Theoretical Physics*, vol. 21, no. 6–7, pp. 467–488, 1982.
- <span id="page-10-5"></span>[7] S. Lloyd, "Universal quantum simulators," *Science*, vol. 273, no. 5278, pp. 1073–1078, 1996.
- <span id="page-10-6"></span>[8] A. Peruzzo, J. McClean, P. Shadbolt, M.-H. Yung, X.-Q. Zhou, P. J. Love, A. Aspuru-Guzik, and J. L. O'Brien, "A variational eigenvalue solver on a photonic quantum processor," *Nature Communications*, vol. 5, no. 1, July 2014. [Online]. Available: [http://dx.doi.org/10.1038/](http://dx.doi.org/10.1038/ncomms5213) [ncomms5213](http://dx.doi.org/10.1038/ncomms5213)
- <span id="page-10-7"></span>[9] A. Kandala, A. Mezzacapo, K. Temme, M. Takita, M. Brink, J. M. Chow, and J. M. Gambetta, "Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets," *Nature*, vol. 549, pp. 242–246, 2017.
- <span id="page-10-8"></span>[10] M. Cerezo, A. Arrasmith, R. Babbush, S. C. Benjamin, S. Endo, K. Fujii, J. R. McClean, K. Mitarai, X. Yuan, L. Cincio *et al.*, "Variational quantum algorithms," *Nature Reviews Physics*, vol. 3, no. 9, pp. 625– 644, 2021.
- <span id="page-10-9"></span>[11] J. Tilly, H. Chen, S. Cao, D. Picozzi, K. Setia, Y. Li, E. Grant, L. Wossnig, I. Rungger, G. H. Booth, and J. Tennyson, "The variational quantum eigensolver: A review of methods and best practices," *Physics Reports*, vol. 986, p. 1–128, Nov. 2022. [Online]. Available: <http://dx.doi.org/10.1016/j.physrep.2022.08.003>
- <span id="page-10-10"></span>[12] J. R. McClean, J. Romero, R. Babbush, and A. Aspuru-Guzik, "The theory of variational hybrid quantum-classical algorithms," *New Journal of Physics*, vol. 18, no. 2, p. 023023, feb 2016. [Online]. Available: <https://doi.org/10.1088/1367-2630/18/2/023023>
- <span id="page-10-11"></span>[13] N. Vaquero-Sabater, A. Carreras, R. Orus, N. J. Mayhall, and ´ D. Casanova, "Physically motivated improvements of variational quantum eigensolvers," *Journal of Chemical Theory and Computation*, vol. 20, no. 12, pp. 5133–5144, 2024.
- [14] R.-Y. Sun, T. Shirakawa, and S. Yunoki, "Efficient variational quantum circuit structure for correlated topological phases," *Physical Review B*, vol. 108, p. 075127, 2023.
- <span id="page-10-12"></span>[15] G. B. Yoo, J. Kim, and Y. Kwon, "Symmetry-protecting ansatz for the variational quantum eigensolver," *Physical Review A*, vol. 111, p. 032615, 2025.
- <span id="page-10-13"></span>[16] E. Farhi, J. Goldstone, S. Gutmann, J. Lapan, A. Lundgren, and D. Preda, "A quantum adiabatic evolution algorithm applied to random instances of an np-complete problem," *Science*, vol. 292, no. 5516, pp. 472–475, 2001.
- <span id="page-10-14"></span>[17] T. Albash and D. A. Lidar, "Adiabatic quantum computation," *Reviews of Modern Physics*, vol. 90, p. 015002, 2018.
- <span id="page-10-15"></span>[18] M. W. Johnson, M. H. S. Amin, S. Gildert, T. Lanting, F. Hamze, N. Dickson, R. Harris, A. J. Berkley, J. Johansson, P. Bunyk, E. M. Chapple, C. Enderud, J. P. Hilton, K. Karimi, E. Ladizinsky, N. Ladizinsky, T. Oh, I. Perminov, C. Rich, M. C. Thom, E. Tolkacheva, C. J. S. Truncik, S. Uchaikin, J. Wang, B. Wilson, and G. Rose, "Quantum annealing with manufactured spins," *Nature*, vol. 473, pp. 194–198, 2011.
- <span id="page-10-16"></span>[19] T. Kadowaki and H. Nishimori, "Quantum annealing in the transverse Ising model," *Physical Review E*, vol. 58, pp. 5355–5363, 1998.
- <span id="page-10-17"></span>[20] P. Hauke, H. G. Katzgraber, W. Lechner, H. Nishimori, and W. D. Oliver, "Perspectives of quantum annealing: methods and implementations," *Reports on Progress in Physics*, vol. 83, no. 5, p. 054401, may 2020. [Online]. Available:<https://doi.org/10.1088/1361-6633/ab85b8>

- <span id="page-10-18"></span>[21] R. Levy, Z. G. Izquierdo, Z. Wang, J. Marshall, J. Barreto, L. Fry-Bouriaux, D. T. O'Connor, P. A. Warburton, N. Wiebe, E. Rieffel *et al.*, "Towards solving the fermi-hubbard model via tailored quantum annealers," *arXiv preprint arXiv:2207.14374*, 2022.
- <span id="page-10-19"></span>[22] S. B. Bravyi and A. Y. Kitaev, "Fermionic quantum computation," *Annals of Physics*, vol. 298, no. 1, pp. 210–226, 2002.
- <span id="page-10-20"></span>[23] J. T. Seeley, M. J. Richard, and P. J. Love, "The Bravyi-Kitaev transformation for quantum computation of electronic structure," *The Journal of Chemical Physics*, vol. 137, no. 22, p. 224109, 2012.
- <span id="page-10-21"></span>[24] F. Tacchino, A. Chiesa, S. Carretta, and D. Gerace, "Quantum computers as universal quantum simulators: State-of-the-art and perspectives," *Advanced Quantum Technologies*, vol. 3, no. 3, p. 1900052, 2020. [Online]. Available: [https://advanced.onlinelibrary.wiley.com/doi/abs/10.](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/qute.201900052) [1002/qute.201900052](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/qute.201900052)
- <span id="page-10-22"></span>[25] B. Fauseweh, "Quantum many-body simulations on digital quantum computers: State-of-the-art and future challenges," *Nature Communications*, vol. 15, no. 1, p. 2123, 2024.
- <span id="page-10-23"></span>[26] D. Aharonov, W. van Dam, J. Kempe, Z. Landau, S. Lloyd, and O. Regev, "Adiabatic quantum computation is equivalent to standard quantum computation," *SIAM Review*, vol. 50, no. 4, pp. 755–787, 2008.
- <span id="page-10-24"></span>[27] S. Kivelson and D. E. Heim, "Hubbard versus peierls and the su-schrieffer-heeger model of polyacetylene," *Phys. Rev. B*, vol. 26, pp. 4278–4292, Oct 1982. [Online]. Available: [https:](https://link.aps.org/doi/10.1103/PhysRevB.26.4278) [//link.aps.org/doi/10.1103/PhysRevB.26.4278](https://link.aps.org/doi/10.1103/PhysRevB.26.4278)
- <span id="page-10-25"></span>[28] W. P. Su, J. R. Schrieffer, and A. J. Heeger, "Solitons in polyacetylene," *Phys. Rev. Lett.*, vol. 42, pp. 1698–1701, Jun 1979. [Online]. Available: <https://link.aps.org/doi/10.1103/PhysRevLett.42.1698>
- <span id="page-10-26"></span>[29] J. Hubbard, "Electron correlations in narrow energy bands," *Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences*, vol. 276, no. 1365, pp. 238–257, 11 1963. [Online]. Available: <https://doi.org/10.1098/rspa.1963.0204>
- <span id="page-10-27"></span>[30] C. Feng, B. Xing, D. Poletti, R. Scalettar, and G. Batrouni, "Phase diagram of the su-schrieffer-heeger-hubbard model on a square lattice," *Physical Review B*, vol. 106, no. 8, Aug. 2022. [Online]. Available: <http://dx.doi.org/10.1103/PhysRevB.106.L081114>
- <span id="page-10-28"></span>[31] N. H. Le, A. J. Fisher, N. J. Curson, and E. Ginossar, "Topological phases of a dimerized fermi–hubbard model for semiconductor nanolattices," *npj Quantum Information*, vol. 6, no. 1, Feb. 2020. [Online]. Available:<http://dx.doi.org/10.1038/s41534-020-0253-9>
- <span id="page-10-29"></span>[32] B.-H. Chen and D.-W. Chiou, "An elementary rigorous proof of bulk-boundary correspondence in the generalized su-schrieffer-heeger model," *Physics Letters A*, vol. 384, no. 7, p. 126168, Mar. 2020. [Online]. Available:<http://dx.doi.org/10.1016/j.physleta.2019.126168>
- <span id="page-10-30"></span>[33] E. Di Salvo, A. Moustaj, C. Xu, L. Fritz, A. K. Mitchell, C. M. Smith, and D. Schuricht, "Topological phases of the interacting su-schrieffer-heeger model: An analytical study," *Physical Review B*, vol. 110, no. 16, Oct. 2024. [Online]. Available: <http://dx.doi.org/10.1103/PhysRevB.110.165145>
- <span id="page-10-31"></span>[34] H.-C. Chang, H.-C. Hsu, and Y.-C. Lin, "Probing entanglement dynamics and topological transitions on noisy intermediate-scale quantum computers," *Physical Review Research*, vol. 7, no. 1, Jan. 2025. [Online]. Available: [http://dx.doi.org/10.1103/PhysRevResearch.](http://dx.doi.org/10.1103/PhysRevResearch.7.013043) [7.013043](http://dx.doi.org/10.1103/PhysRevResearch.7.013043)
- <span id="page-10-32"></span>[35] A. A. Aligia and G. Ortiz, "Quantum mechanical position operator and localization in extended systems," *Physical Review Letters*, vol. 82, no. 12, p. 2560–2563, Mar. 1999. [Online]. Available: <http://dx.doi.org/10.1103/PhysRevLett.82.2560>
- <span id="page-10-33"></span>[36] R. Resta and S. Sorella, "Electron localization in the insulating state," *Physical Review Letters*, vol. 82, no. 2, p. 370–373, Jan. 1999. [Online]. Available:<http://dx.doi.org/10.1103/PhysRevLett.82.370>
- <span id="page-10-34"></span>[37] R. Resta, "Quantum-mechanical position operator in extended systems," *Physical Review Letters*, vol. 80, no. 9, p. 1800–1803, Mar. 1998. [Online]. Available:<http://dx.doi.org/10.1103/PhysRevLett.80.1800>
- <span id="page-10-35"></span>[38] J. Zak, "Berry's phase for energy bands in solids," *Phys. Rev. Lett.*, vol. 62, pp. 2747–2750, Jun 1989. [Online]. Available: <https://link.aps.org/doi/10.1103/PhysRevLett.62.2747>
- <span id="page-10-36"></span>[39] J. K. Asboth, L. Oroszl ´ any, and A. P ´ alyi, ´ *A Short Course on Topological Insulators*. Springer International Publishing, 2016. [Online]. Available:<http://dx.doi.org/10.1007/978-3-319-25607-8>
- <span id="page-10-37"></span>[40] Z. Jiang, K. J. Sung, K. Kechedzhi, V. N. Smelyanskiy, and S. Boixo, "Quantum algorithms to simulate many-body physics of correlated fermions," *Physical Review Applied*, vol. 9, no. 4, Apr. 2018. [Online]. Available:<http://dx.doi.org/10.1103/PhysRevApplied.9.044036>

- <span id="page-11-0"></span>[41] J. R. McClean, K. J. Sung, I. D. Kivlichan, Y. Cao, C. Dai, E. S. Fried, C. Gidney, B. Gimby, P. Gokhale, T. Haner, T. Hardikar, V. Havl ¨ ´ıcek, ˇ O. Higgott, C. Huang, J. Izaac, Z. Jiang, X. Liu, S. McArdle, M. Neeley, T. O'Brien, B. O'Gorman, I. Ozfidan, M. D. Radin, J. Romero, N. Rubin, N. P. D. Sawaya, K. Setia, S. Sim, D. S. Steiger, M. Steudtner, Q. Sun, W. Sun, D. Wang, F. Zhang, and R. Babbush, "Openfermion: The electronic structure package for quantum computers," 2019. [Online]. Available:<https://arxiv.org/abs/1710.07629>
- <span id="page-11-1"></span>[42] G. Q. AI. Openfermion. [Online]. Available: [https://quantumai.google/](https://quantumai.google/openfermion) [openfermion](https://quantumai.google/openfermion)
- <span id="page-11-2"></span>[43] A. Javadi-Abhari, M. Treinish, K. Krsulich, C. J. Wood, J. Lishman, J. Gacon, S. Martiel, P. D. Nation, L. S. Bishop, A. W. Cross, B. R. Johnson, and J. M. Gambetta, "Quantum computing with Qiskit," 2024.