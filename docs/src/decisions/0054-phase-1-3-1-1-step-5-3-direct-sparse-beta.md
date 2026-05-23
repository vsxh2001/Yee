# ADR-0054 — Phase 1.3.1.1 step 5.3: direct β-direct sparse shift-and-invert

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.3 (inhomogeneous §4 closure)

## Context

Step 5.2 fixed the β-extraction (uniform-fill anchor exact) but the
*inhomogeneous* mixed path ships a **hybrid**: it selects the dominant
mode on the cutoff pencil `A x = k_c² B x`, then extracts β² as a
Rayleigh quotient on that eigenvector. Because the cutoff-pencil
eigenvector differs from the true β-direct eigenvector for inhomogeneous
ε_r, this leaves a mesh-stable ~17% bias at ε_r=10.2 (reviewed: an
inseparable mix of (a) discretization and (b) the RQ eigenvector
mismatch). A naive *direct* β-direct solve drifts onto the spurious
`E_z≈0` gradient cluster at `β²≈k₀²⟨ε_r⟩`.

## Decision

Solve the β-direct pencil `(k₀²B − A) x = β² B₁ x` **directly** via a
**faer sparse shift-and-invert** with a physics-informed shift
`σ₀ = (k₀²−k_c²)⟨ε_r⟩` (the hybrid's own β² estimate, within ~17% of
physical — ample to target the physical eigenpair and avoid the spurious
cluster). Inverse-iterate to the eigenpair nearest σ₀, screen with the
transverse-energy filter, take β² as the Rayleigh quotient on the
converged **true** β-direct eigenvector. faer is already a yee-mom
dependency — **no new dependency**. Close the §4 inhomogeneous gap at a
representative contrast (FR-4 ε_r=4.4, ≤5% vs the step-5.1 reference);
ε_r=10.2 is a stretch that, if still short, localises a pure
discretization residual for step-5.4.

## Rationale

(1) **Targeted, not global.** The spurious cluster is what made a global
β-direct eigen-sweep drift. A shift placed near the physical β² isolates
and amplifies the physical eigenpair — the standard remedy for
interior/clustered spectra.

(2) **Recovers the true eigenvector** ⇒ β² is exact for that mode (RQ
stationary at its own eigenvector), eliminating residual source (b)
directly. Mesh refinement (now tractable via sparse LU) then addresses
(a), and the convergence study discriminates the two — the step-5.2
review's open question.

(3) **No new dependency / pure-Rust.** faer's sparse LU (already used by
the yee-fem `build_shifted` shift-invert pattern) handles the indefinite
shifted operator; reuses the dependency yee-mom already links.

(4) **FR-4 as the primary §4 anchor.** ε_r=4.4 is the most common PCB
substrate and a representative contrast; closing it is the meaningful
benchmark. ε_r=10.2 (RT/duroid 6010) is a high-contrast stretch whose
residual, if any, isolates discretization for higher-order work.

## Consequences

* New sparse shift-and-invert path in `solve.rs`; the dense path stays
  as the small-`n` reference / fallback.
* The inhomogeneous reconciliation upgrades from a non-failing
  diagnostic to a failing ≤5% gate at FR-4 (§4 closed there).
* A mesh-convergence study resolves the step-5.2 (a)-vs-(b) question.
* If ε_r=10.2 stays short with the true eigenvector, that residual is
  pure discretization → step-5.4 (higher-order / p-refinement).
* Lossy complex-ε_r sparse solve remains Phase 1.3.1.2.

## As-built outcome (2026-05-23, merge `6eca76a` + hardening `d051ebf`)

The decision held as designed. Confirmed results:

* **§4 closed at FR-4:** ε_r=4.4 horizontal slab β=328.57 vs the verified
  reference 324.05 — rel 1.39% ≤5%, a failing→passing gate.
* **(a)-vs-(b) settled empirically:** the true β-direct eigenvector
  improved on the step-5.2 hybrid by only **~1%** (483→489), so the
  RQ-eigenvector-mismatch (b) is ~1%; the ε_r=10.2 ~16% residual is
  **(a) discretization-dominated** (β plateaus 489→487→486 under
  refinement) → step-5.4 (higher-order elements).
* **Eigenvector finding:** the true β-direct dominant mode is
  near-pure-transverse (‖E_z‖/‖E_t‖≈2e-5); the step-5.2 hybrid's 0.0105
  was a *cutoff-pencil* artifact, not the physical mode. The coupling
  block is nonetheless strongly load-bearing (zeroing `B_tz` shifts β
  ≈49%), so the coupling guard was re-anchored from the E_z-fraction
  onto that β-sensitivity (a stronger test).
* **Hardening (`d051ebf`):** `beta_direct_shift_invert` now errors on
  non-convergence rather than returning an unconverged β² (review P1-1) —
  landed before step-5.4 introduces new geometries.

## References

* Saad, *Numerical Methods for Large Eigenvalue Problems* (shift-invert).
* Jin, *FEM in EM* 3rd ed. §8.4.
* ADR-0051/0052/0053 (the step-5 chain).
* Step 5.3 spec + plan.
* `crates/yee-mom/src/eigensolver/{solve,reference}.rs`,
  `crates/yee-fem/src/solve.rs` (the shift-invert pattern).
