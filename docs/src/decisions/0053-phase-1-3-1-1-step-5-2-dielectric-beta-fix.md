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

## References

* Jin, *FEM in EM* 3rd ed. §8.2-8.4.
* Pozar, *Microwave Engineering* 4th ed. §3.3.
* ADR-0051 (step-5 as-built), ADR-0052 (step-5.1 finding).
* Step 5.2 spec + plan.
* `crates/yee-mom/src/eigensolver/{solve,assembly,reference}.rs`.
