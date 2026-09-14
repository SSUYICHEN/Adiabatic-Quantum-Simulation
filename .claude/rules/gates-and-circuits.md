---
paths:
  - "src/aqs/circuits/**"
  - "src/aqs/core.py"
---

# Gate Primitives & Circuit Construction

## Gate conventions (paper Eq. 37–39 ≡ spec §3.1 — they agree; both are normative)

- `R(θ) = exp(−i θ/2 (XX + YY))` — implemented as `rxx(θ)·ryy(θ)`, **no factor of two
  anywhere**.
- `G(θ) = exp(−i θ/2 (X_i Y_j − Y_i X_j))` — X acts on the **first** operand passed.
  G is antisymmetric under qubit swap (R is symmetric); operand order is load-bearing
  whenever `Im(t) ≠ 0`.
- `CP(φ) = exp(−i φ/4 (I−Z)(I−Z)) = exp(−i φ n_i n_j)` — qiskit's `cp(λ)` gives
  `diag(1,1,1,e^{+iλ})`, so the implementation is `cp(−φ)`.

Docstrings must state the exact matrix/exponent and cite the equation. **Changing an angle
sign or factor is a physics decision requiring the H_eff round-trip proof — never a
refactor side effect.** The `main` branch and `fix/cudaq-bit-order-and-uv-migration` use an
older convention (factors of 2, `_G_SIGN`) with compensating call sites; the two
conventions must never be mixed in one tree.

- G operand order resolves a paper ambiguity (Eq. 30 vs Eq. 40): this codebase follows the
  **Eq. 40 reading — lower qubit index first**, `G(i, i+1)`. Documented, not to be
  "corrected" silently.

## Spinful SSHH Trotter layer (paper Eq. 33–35, 40–41)

- First-order Trotter; per-step order: hopping layers (v then w), then Hubbard CP layer.
- Hopping angles per step: `θ_R = −δt·Re(t)`, `θ_G = +δt·Im(t)` — derived from this code's
  single-particle H (`H[i,j] = −t`, the Hermitian conjugate of the paper's convention, so
  code v,w = conj(paper v,w); see `models-and-observables.md`).
- Hubbard ramp evaluated at the layer midpoint: `s_ℓ = (2ℓ−1)/(2L)`; CP angle
  `δt·s_ℓ·U_{A/B}`, staggered A/B, acting between spin blocks `(i, 2N+i)`.
- **PBC wrap bond:** the hopping term carries the JW parity factor — wrap angle
  `= −(−1)^{n_s} ×` bulk angle (the extra minus beyond the paper's bare `(−1)^{n_s}`
  absorbs the reversed-order wrap bond's JW string sign; pinned by
  `test_annealing_circuit.py::test_pbc_wrap_bond_carries_the_parity_factor`).
  Interaction/density terms are diagonal and carry **no** parity correction. This
  asymmetry is deliberate and extremely easy to get wrong.

## Spinless SSH+NN Trotter layer (spec §3.2, `reconstruct_goals.md`)

- Per layer ℓ: hopping propagators with **constant** angles
  `R(−δt·Re(v))`, `G(−δt·Im(v))` intracell `(2j, 2j+1)`;
  `R(−δt·Re(w))`, `G(−δt·Im(w))` intercell `(2j+1, 2j+2)` —
  then NN interaction CP gates scaled by `λ_ℓ = (2ℓ−1)/(2L)`:
  `CP(δt·V_v·λ_ℓ)` intracell, `CP(δt·V_w·λ_ℓ)` intercell.
- These angles follow from the spec's **own Pauli representation (§2.1)**, which is the
  normative source for the spinless model. Note the spec's fermionic H is `+v b†a + h.c.`
  (opposite hopping sign to the spinful module's `−t` convention), so do not port angle
  expressions between the two models by pattern matching.
- ⚠️ The other branch's `spinless.py` calls the **old** gate conventions
  (`th_R = −2.0·dt·t_re`, `_G_SIGN·2.0·dt·t_im`, `CP(−dt·V·λ)`). Never copy its call
  sites. Any G-sign doubt is settled one way only: H_eff round-trip against the spec's
  Pauli H (little-endian) with **complex** v, w.
- PBC wrap for the spinless chain: hopping wrap bond carries the parity factor
  `−(−1)^{N_f}`; the NN interaction wrap term `n_{B,N−1} n_{A,0}` is diagonal and carries
  **none**.

## State preparation (paper Eq. 42–45, spec §4)

- X gates set the occupation pattern, then OpenFermion
  `slater_determinant_preparation_circuit(Q)` Givens tuples `(j, k, θ, φ)` lowered by
  `_givens_instruction`: the explicit 4×4 `UnitaryGate`
  `diag-block {1; [[c, s],[−s·e^{iφ}, c·e^{iφ}]]; e^{iφ}}` (OpenFermion convention,
  including the `det G = e^{iφ}` phase on `|11⟩`).
- Do not replace the UnitaryGate with a gate decomposition without a unitary-equivalence
  proof (the historical `rz/cx/cry/cx` decomposition was wrong twice over: negated angle →
  prepared the **highest** band at `+E`; missing determinant phase, unfixable by any
  parameter choice).
- Q matrices are the **lowest** occupied orbitals of the single-particle H. After any
  change, the prepared state must pass the lowest-fill-energy check and match an
  OpenFermion-built reference determinant (see `testing-and-verification.md`).
- Fillings: spinful — `num_up`/`num_dn` per sector; spinless — N particles (PBC,
  half-filling) or N+1 (OBC, edge modes occupied).

## Structural invariants (cheap tests, always keep)

- Particle number is conserved by every circuit (prep and evolution).
- OBC circuits contain no wrap-bond gate.
- CP angles grow linearly across the ramp.
- The circuit stays normalized with complex hoppings.
- Any circuit-structure change keeps the H_eff round-trip passing at multiple δt — an
  error that does not shrink with δt is a convention bug, not Trotter error.
