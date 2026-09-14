---
paths:
  - "src/aqs/models/**"
  - "src/aqs/core.py"
  - "src/aqs/observables.py"
---

# Model Definitions & Observables

## Two models, two authorities

| Model | Authority | Qubits | Ramped term |
|---|---|---|---|
| Spinful SSHH | paper `.claude/specs/paper_theory/SSHH.md` (Eq. 9–13) | 4N | on-site Hubbard `U_A, U_B` |
| Spinless SSH+NN | spec `.claude/specs/reconstruct_goals.md` §1–2 | 2N | NN interaction `V_v, V_w` |

For the spinless model, the spec's **Pauli-operator representation (§2)** is normative —
its fermionic H is `+v b†a + h.c.`. The spinful module uses `H[i,j] = −t` (Hermitian
conjugate of the paper's `v b†a` convention), so **code v,w = conj(paper v,w)** there.
Moot while hoppings are real; a hard blocker to resolve before any complex-hopping
feature. Never port sign expressions between the two models by analogy.

## Model invariants

- v = intracell hopping (even bonds), w = intercell hopping (odd bonds) — standard SSH
  naming in both models.
- Chiral symmetry: with no interaction imbalance, the single-particle matrix is invariant
  under A↔B exchange — keep a structural test for it (it also underlies the ±E trap in
  state prep).
- Boundary conditions are a model parameter; the single-particle H builder and the circuit
  builder must agree on it (OBC skips the wrap bond in both).
- Adiabatic-safety bound (paper Eq. 25): `max{U_A, U_B, |ΔU|} < min{|v+w|, |v−w|}` —
  experiments outside it may not track the true ground state; plotting already renders
  such points distinctly (dotted).

## Observable invariants

- **Berry phase (PBC only).** Twist invariant `z_N = ⟨exp(i 2π/N · X̂)⟩` estimated from
  computational-basis samples (paper Eq. 47–50; spec §5.2). Never assert absolute phase
  values: the paper's `X̂ = Σ (j+1) n_j` carries a constant `π(N+1)` filling/parity offset
  (for even n the trivial phase sits at γ = π and the topological at 0 — inverted from
  naive labels); the spec's `X̂ = Σ ⌊q/2⌋ n_q` is 0-based and offsets differently. The
  physical signal is the **π jump across the transition**, nothing else.
- **Sublattice polarization (OBC only).** `P_j = ⟨n_Aj − n_Bj⟩` (both spins summed in the
  spinful model, paper Eq. 51–52). Under PBC polarization is identically zero — a nonzero
  PBC polarization is a caller bug, and a "polarization result" computed under PBC is
  meaningless. Force/assert OBC in every polarization path.
- **Edge signature.** Edge response requires OBC with the edge modes occupied (spinful
  n = 2N+2; spinless N+1). Cell 0 / cell N−1 polarization is the edge readout.
- Exact reference states come from diagonalization restricted to the **fixed
  particle-number sector** of the circuit (the circuits conserve particle number; the
  interacting global ground state may live at a different filling — comparing against it
  is wrong, e.g. the N=3 interacting example's global minimum has 2 particles, not 3).
- Bit extraction in observables is little-endian (`(idx >> q) & 1`); layout arithmetic
  comes from the layout helpers, never re-derived locally (see
  `endianness-and-mapping.md`).

## Frozen contracts

- `PROPERTIES` registry keys are part of the CLI/JSON contract — additions are fine,
  renames/removals are breaking changes and need explicit sign-off.
- Recorded decision: `berry` currently defaults to `mode="up_spin"`; the paper's Eq. 46
  sums both spin blocks (`mode="total"`). Whichever default the reconstruction picks,
  keep both modes implemented and tested, and state the choice in `--help`.
