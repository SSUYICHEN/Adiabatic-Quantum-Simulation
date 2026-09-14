---
paths:
  - "tests/**"
---

# Testing & Verification Methodology

Every historical defect in this repo was a **silent convention mismatch**: no exception, no
test failure, plausible-looking output, wrong physics (see `docs/experience.md`). The rules
below exist because ordinary testing practice did not catch any of them.

## Prime directive

"Runs without error" and "the numbers look reasonable" are evidence of **nothing**.
A physics result counts as verified only when it matches an **independent reference**.

Acceptable independent references:
- Exact diagonalization restricted to the **fixed particle-number sector** of the prepared
  state (the interacting global ground state may live at a different filling).
- Slater determinants constructed directly with OpenFermion — **endian-converted** with
  `_to_little_endian` before comparison.
- The lowest-fill-energy check: the prepared state's `⟨H⟩` must equal the sum of the
  occupied lowest single-particle eigenvalues. On a chiral-symmetric chain the wrong
  determinant gives exactly `−E`, so this check flips sign rather than drifting.

Forbidden as references:
- The code's own earlier output (`verify_consistency_vs_main.py` is a regression tripwire,
  not a correctness proof).
- The prepared state itself as the fidelity baseline (`fidelity_scan` is structurally blind
  to which determinant was prepared — this is exactly how the Givens defect hid).

## Systematic-error detectors

- **H_eff round-trip:** extract `H_eff = i·log(U)/δt` from the single-step propagator and
  compare element-wise to the target Hamiltonian, at more than one δt. Trotter error shrinks
  with δt; a conjugated imaginary part or wrong sign **does not**. This is the only test
  class that catches "the imaginary part was silently conjugated".
- **L-scaling discrimination:** fidelity improving with L ⇒ Trotter error (fine). Fidelity
  plateauing below 1 ⇒ systematic convention bug — go hunting. Never scan T and L together
  at fixed δt = T/L: the plateau that produces is constant Trotter error, not adiabatic
  failure. Fix one, scan the other.

## Test-input requirements

- Ordering/endianness guards MUST use **asymmetric inputs** (non-uniform density, e.g. a deep
  potential well on one mode). Uniform-density states are invariant under bit reversal and
  structurally cannot detect ordering bugs. Reflection-symmetric circuits likewise cannot
  distinguish the two qubit conventions.
- Sign/conjugation guards MUST include **complex hoppings** (`Im(v), Im(w) ≠ 0`). All the
  real-hopping terms are invariant under the OpenFermion mode swap; only complex terms
  expose big-endian contamination.

## Reverse verification

A guard test only counts once it has been shown to **fail with the historical bug
reintroduced** (mutate the code back, watch it go red, revert). Record the mutation that
kills the test in the test's docstring. A guard that has never failed guards nothing.
All four existing convention-guard suites were validated this way; hold new tests to the
same bar.

## Structural rules

- The four convention-guard suites are load-bearing and must survive any restructure with
  their intent intact: `test_gate_primitives.py`, `test_state_preparation.py`,
  `test_backend_endianness.py`, `test_hamiltonian_endianness.py`.
- Standalone `tests/verify_*.py` A/B harnesses stay runnable as plain scripts (they are
  designed to run against arbitrary revisions via `PYTHONPATH` + `git worktree`).
- Coverage target is ~99%, but guard quality outranks coverage: never delete or weaken a
  guard to simplify a refactor.
- GPU-dependent tests must skip cleanly on CPU-only environments.
