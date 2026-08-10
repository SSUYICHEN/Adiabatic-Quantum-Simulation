# Normalise the two-qubit gate primitives to match their documented form

| | |
|---|---|
| **Source branch** | `fix/gate-primitive-normalization` |
| **Target branch** | `main` (`0fe9d1d`) |
| **Type** | Pure refactor — bit-identical output |
| **Risk** | Low on its own; **see “Interaction with the other branch” before merging both** |

---

## Problem

The three gate wrappers in `src/aqs/core.py` did not implement the operators
their own docstrings describe:

| function | docstring claims | actually implemented | discrepancy |
|---|---|---|---|
| `create_R_gate(t)` | `exp(-i t/2 (XX+YY))` | `R_spec(t/2)` | factor 2 |
| `create_G_gate(t)` | `exp(-i t/2 (XY-YX))` | `G_spec(-t/2)` | factor 2 **and** negated |
| `create_CP_gate(t)` | controlled phase | `CP_spec(-t)` | negated |

This was **not a deliberate alternative convention** — the docstrings state the
spec form, so the code contradicted its own contract. The factors were
incidental artefacts of wrapping qiskit primitives without normalising them:

- `rxx(p) = exp(-i p/2 XX)` and `XX` commutes with `YY`, so
  `R_spec(t) == rxx(t) · ryy(t)` exactly — the stray `/2.0` *was* the whole defect.
- `cp(l) = diag(1,1,1,e^{+il})`, so `CP_spec(t)` is simply `cp(-t)`.

`build_annealing_circuit` compensated at every call site (`-2.0*tau*w`,
`-tau*ramp*U`), which is why the physics was right despite the wrappers being
wrong.

**Why it is worth fixing even though nothing is numerically wrong today:**
anyone calling `create_R_gate(theta)` straight from its docstring silently gets
half the intended angle — no exception, no failing test, just wrong physics.
This repository has already produced several defects with exactly that failure
signature.

## Change

```diff
-    qc.rxx(theta / 2.0, 0, 1)
-    qc.ryy(theta / 2.0, 0, 1)
+    qc.rxx(theta, 0, 1)
+    qc.ryy(theta, 0, 1)
```

```diff
-    qc.sdg(0); qc.rxx(theta / 2.0, 0, 1); qc.s(0)
-    qc.sdg(1); qc.rxx(-theta / 2.0, 0, 1); qc.s(1)
+    qc.sdg(0); qc.rxx(-theta, 0, 1); qc.s(0)
+    qc.sdg(1); qc.rxx(theta, 0, 1); qc.s(1)
```

```diff
-    qc.cp(theta, 0, 1)
+    qc.cp(-theta, 0, 1)
```

The call sites simplify correspondingly — the factor of two disappears:

```diff
-                    theta_R = -2.0 * tau * w_R
-                    theta_I = -2.0 * tau * w_I
+                    theta_R = -tau * w_R
+                    theta_I = tau * w_I
```

```diff
-                        theta_U = -tau * ramp * current_U
+                        theta_U = tau * ramp * current_U
```

Net effect: the emitted `rxx` / `ryy` / `cp` angles are the *same floating-point
values* as before.

## Verification — bit-identical, not merely close

`tests/verify_gate_norm_equivalence.py` captures three levels of output before
and after and diffs them exactly:

```
1. 原始電路指令串（gate 名稱 + qubits + 角度，12 位小數）
  [IDENTICAL] N2_PBC_U0.0-0.0  (109 instructions)
  [IDENTICAL] N2_PBC_U0.5-1.3  (125 instructions)
  [IDENTICAL] N3_OBC_U1.0-1.0  (201 instructions)
2. 態向量
  [IDENTICAL] weighted moments
  [IDENTICAL] leading amplitudes
3. 物理可觀測量
  [IDENTICAL] berry sweep (8 points)
  [IDENTICAL] polarization sweep (4 curves)
  [IDENTICAL] fidelity scan (3x3)

=> BIT-IDENTICAL, 行為完全保留
```

Reproduce:

```bash
git checkout main   && python tests/verify_gate_norm_equivalence.py before.json
git checkout fix/gate-primitive-normalization
                       python tests/verify_gate_norm_equivalence.py after.json
diff <(python -m json.tool before.json) <(python -m json.tool after.json)   # empty
```

No published result changes. Berry phase, sublattice polarization and the
fidelity scan are unaffected.

## New guard tests

`tests/test_gate_primitives.py` (6 tests, runs under `pytest` or standalone)
pins the convention so this cannot silently regress:

- each primitive against its documented operator, up to global phase;
- `G`’s antisymmetry under qubit swap (`R` has no such asymmetry);
- an **effective-Hamiltonian round trip**: reconstruct `H_eff = i·log(U)/dt`
  for one bond and compare against the true little-endian bond Hamiltonian.

That last one earns its place. A conjugated imaginary part is a *systematic*
error — it does **not** shrink as `dt → 0`, so an ordinary Trotter-convergence
test cannot detect it. Only reconstructing the generator exposes it.

Confirmed the tests actually bite: against the pre-refactor `core.py` they fail
4 of 6, with the diagnostic

```
FAIL test_trotter_bond_reproduces_the_true_hamiltonian: t=(0.7+0.45j)
     true[1,2] = -0.700000-0.450000j
      eff[1,2] = -0.350000+0.225000j        <- half magnitude AND conjugated
```

and pass 6 of 6 after.

```bash
python tests/test_gate_primitives.py     # 6/6 passed
```

## A trap documented in the code

Deriving the `G` sign requires care: `openfermion.get_sparse_operator` is
big-endian, so for a two-mode system it swaps modes 0 and 1 and **silently
conjugates** the off-diagonal that carries `Im(t)`. The real part `XX+YY` is
invariant under that swap and therefore never reveals the mistake. A comment at
the call site and a docstring note in the test record this; the test builds its
reference little-endian explicitly rather than importing a helper, so it stays
self-contained.

## Interaction with the other branch — read before merging both

`fix/cudaq-bit-order-and-uv-migration` adds `src/aqs/spinless.py`, which calls
these same primitives with the **old** convention:

```python
th_R = -2.0 * dt * t_re
th_G = _G_SIGN * 2.0 * dt * t_im     # _G_SIGN = -1
qc.append(create_CP_gate(-dt * model.V_v * lam), ...)
```

If both branches land without reconciliation, the spinless model breaks
**silently** — wrong by a factor of two on the hopping angles and by a sign on
the interaction, with no test failure on that branch. Whichever merges second
must apply:

```python
th_R = -dt * t_re
th_G = dt * t_im                      # _G_SIGN no longer needed
qc.append(create_CP_gate(dt * model.V_v * lam), ...)
```

Git will **not** flag this: the two branches touch different files, so there is
no textual conflict — only a semantic one. Recommended order: merge this MR
first, then rebase the other branch and re-run its validation suite (state-prep
fidelity 1.0, adiabatic convergence ≥ 0.9997).

## Checklist

- [x] Bit-identical output verified at instruction / statevector / observable level
- [x] Guard tests added; confirmed failing before and passing after
- [x] No change to any published figure or number
- [ ] Reconcile `spinless.py` call sites when the other branch merges
