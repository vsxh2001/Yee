# ADR-0053 — Phase 1.3.1.1 step 5.2: dielectric β-extraction fix

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.2 (solver correctness fix)

## Context

Step 5.1's independently-verified slab-loaded-guide reference disagreed
with the step-5 mixed cross-section eigensolver by ~2.9× (solver
ε_eff≈1.35 vs reference 8.17 for a half-ε_r=10.2 fill — physically
impossible). Root-cause derivation (spec §2): the solver forms
`S x = k_c² T_ε x` with an ε_r-weighted mass `T_ε = ∫ε_r N·N`, then
extracts `β² = k₀² − k_c²` with vacuum `k₀`. The physical transverse
equation `∇×(1/μ_r ∇×E_t) = (k₀²ε_r − β²)E_t` gives
`(k₀² T_ε − S) x = β² T_1 x` with **unweighted** RHS mass `T_1 = ∫N·N`
and eigenvalue β² directly. The current `β² = k₀² − k_c²` is equivalent
only when `ε_r ≡ 1`; for any `ε_r ≠ 1` (uniform or inhomogeneous) it
under-counts the dielectric. The ε_r=1 homogeneous canary passed (4e-14)
because the two forms coincide there, masking the bug; the V2′ bracket
gate (158–655) was too wide to catch the wrong loaded β=201.

## Decision

Reformulate the β-extraction to solve `(k₀² T_ε − S) x = β² T_1 x`
(eigenvalue β² directly, **unweighted** RHS mass), for both the
transverse `solve_dense` and the mixed `solve_dense_mixed` paths. Add a
**uniformly-filled-guide** analytic test (`β = √(ε_r k₀² − (π/a)²)`) as
the smoking-gun isolation + primary anchor, then tighten the
inhomogeneous gate to a published-benchmark comparison against
`reference.rs` (≤5%) and close the CLAUDE.md §4 gap.

## Rationale

(1) The bug is structural, not a tuning issue: the dispersion relation
`β² = k₀² − k_c²` is the homogeneous-air special case, hard-coded.
Reformulating to the eigenvalue-is-β² arrangement makes the physical
quantity the direct output and removes the special-case assumption.

(2) The uniform-fill analytic is a *fully independent* anchor (closed
form, no transverse resonance, no FEM) — it isolates the β-extraction
bug from inhomogeneity and from the coupling block (which step 5.1
verified is correct). It is the cheapest possible confirmation and the
strongest published-benchmark anchor for the fix.

(3) The wide-bracket gate gave false confidence (the reviewer flagged
this at step-5 as P1-2). Tightening to the reference / analytic closes
the §4 gap honestly. The stale regression values (180.23 / 201.52) were
wrong and are replaced, not preserved.

## Consequences

* `solve_dense` / `solve_dense_mixed` reformulated; `assemble_*` exposes
  an unweighted mass `T_1` alongside `T_ε`. The ε_r=1 path is preserved
  bit-identically (canary 4e-14).
* The step-5 inhomogeneous β values were **wrong** and are corrected;
  ADR-0051/0052 are superseded on the β-extraction point.
* The inhomogeneous gate upgrades from the V2′ bracket to a published
  benchmark; §4 gap closed.
* Mode selection becomes "largest valid β²" under the indefinite
  `(k₀²T_ε − S)` operator (re-validated).
* Lossy/heterogeneous complex-ε_r β-extraction remains Phase 1.3.1.2.

## As-built amendment (2026-05-23)

Two refinements during implementation; the merged code (`3b6d899`) is
authoritative where this differs from the Decision above.

1. **Mixed path is a hybrid, not pure option A.** Option A (solving the
   β-direct pencil `(k₀²B − A) x = β² B₁ x` directly) was tried and
   **drifts off the physical mode** onto a spurious `E_z ≈ 0` β-direct
   branch that interleaves with the gradient cluster at the top of the
   spectrum (where `(K − σB₁) ≈ −A` is near-singular, so shift-invert
   thrashes). The shipped mixed path therefore **selects** the dominant
   mode on the *cutoff* pencil `A x = k_c² B x` (gradient null-space
   cleanly at `k_c² ≈ 0`) and then **extracts** `β²` as the β-direct
   Rayleigh quotient on that eigenvector. The transverse `solve_dense`
   *does* use the clean option-A form (no E_z block, no drift).

2. **§4 gap is closed for the extraction, not yet for the inhomogeneous
   case.** The uniformly-filled-guide analytic (β=√(ε_r k₀²−(π/a)²),
   rel 1.5e-4) is a closed-form published benchmark and **passes** —
   certifying the β-extraction. The *inhomogeneous* slab-loaded gate is
   **not** at ≤5%: a mesh-stable ~17% residual vs the step-5.1 reference
   remains, and the reconciliation stays a non-failing diagnostic (the
   V2′ bracket is retained as the inhomogeneous floor). So the Decision's
   "§4 gap closed" holds for the extraction anchor only.

3. **The ~17% residual is (a)+(b) inseparable** at first-order elements
   on the high-contrast (ε_r=10.2) interface: (a) discretization, and
   (b) a Rayleigh-quotient eigenvector mismatch — the hybrid evaluates
   `β²` on the *cutoff*-pencil eigenvector, which differs from the true
   *β-direct* eigenvector for inhomogeneous ε_r (they coincide on
   uniform fill, hence the exact uniform anchor; reviewed magnitude
   check: 17% in β ≈ 31% in β² ⇒ ~34° eigenvector angle if (b)
   dominates — plausible at 10.2:1 contrast). **step-5.3 = the DIRECT
   β-direct-pencil solve with a sparse shift-and-invert** that targets
   the physical β² past the spurious cluster — this resolves both (a)
   and (b); a finer mesh with the same hybrid would plateau at another
   biased value.

## References

* Jin, *FEM in EM* 3rd ed. §8.2-8.4.
* Pozar, *Microwave Engineering* 4th ed. §3.3.
* ADR-0051 (step-5 as-built), ADR-0052 (step-5.1 finding).
* Step 5.2 spec + plan.
* `crates/yee-mom/src/eigensolver/{solve,assembly,reference}.rs`.
