# ADR-0057 — Phase 1.3.1.1 step 5.6: p=2-robust mode selection (ε_eff-screened)

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.6 (high-contrast inhomogeneous closure)

## Context

step 5.5 validated p=2 elements but ε_r=10.2 did not close — p=2 is
worse than p1 because `solve_dense_mixed` selects the smallest
transverse-dominated k_c², and at p=2 the curl-free gradient edge
functions enlarge the near-null cluster with spurious modes that pass the
transverse-energy filter, so the selection locks onto a non-dominant
(low-ε_eff) mode (ADR-0056 as-built; reviewer-confirmed). The physical
dominant quasi-TEM mode is distinguished by **high ε_eff** (field
concentrated in the dielectric); the spurious gradient modes spread
uniformly (low ε_eff, below the area-average).

## Decision

Select the dominant mode by the **β-direct Rayleigh quotient / ε_eff
maximum** among transverse-energy-dominated propagating candidates,
rather than the smallest cutoff k_c². The physical dominant mode is the
slowest wave = largest β² = highest ε_eff; this directly rejects the
gradient-cluster contamination at p=2. Optionally seed the β-direct
shift-invert σ₀ from a physics ε_eff estimate. The new rule must reduce
to the current behaviour where it is already correct (p1 / homogeneous /
FR-4). Also: wire `ElementOrder::Second` through `NumericalCrossSection`
(end-to-end p=2) and add the p=2 uniform-fill analytic anchor (step-5.5
review P1-1).

## Rationale

(1) **Targets the precise, reviewer-confirmed blocker.** The diagnosis
(selection picks a low-ε_eff gradient mode at p=2) yields a direct
discriminator: ε_eff. Selecting the highest-ε_eff transverse mode is the
physical definition of the dominant quasi-TEM mode and is order-robust
(it does not depend on the null-space dimension, which is what grew at
p=2).

(2) **High-confidence, in-lane, well-scoped.** Unlike the
breadth-rotation alternatives (grind-risky), this is a precise fix to a
localised cause with a crisp DoD (ε_r=10.2 ≤5%) and strong guards (FR-4
/ homogeneous gates must stay green). step 5.5 proved p=2 elements are
correct, so the only remaining lever is selection.

(3) **Closes the certification gaps too.** The p=2 uniform-fill anchor
(P1-1) closes the ε_r≠1 `assemble_mixed_p2` coverage gap; the ports.rs
wiring makes the validated p=2 path usable end-to-end.

## Consequences

* `solve_dense_mixed` selection changes from smallest-k_c² to
  highest-ε_eff (β-direct RQ) among valid candidates; p1 behaviour
  preserved (the new rule reduces to it, or is scoped p=2-only if needed).
* `ElementOrder::Second` reachable via `NumericalCrossSection`.
* If selection is fixed but ε_r=10.2 still > 5% at the dense-cap mesh, a
  finer `--release` example + documented finding (selection-fix validated
  by ε_eff recovery toward 8.17) — not a failure.
* The dense cutoff selection's `O(n³)` cost (the p=2 CI burden) is a
  perf follow-on (sparse selection); not blocking the closure.

## References

* Jin §9 (spurious modes / edge-element null space). Boffi-Brezzi-
  Demkowicz. ADR-0054/0055/0056. Step 5.6 spec + plan.
* `crates/yee-mom/src/eigensolver/{solve,assembly,reference}.rs`,
  `crates/yee-mom/src/ports.rs`.
