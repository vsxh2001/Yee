# ADR-0058 — Phase 1.3.1.1 step 5.7: sparse cutoff-mode selection

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** 1.3.1.1 step 5.7 (mesh scaling + the last ε_r=10.2 lever)

## Context

step 5.6 fixed mode selection but found p=2 ≈ p1 at the dense-cap mesh
(~6×6, n≈457) — the dense `complex_eigenvalues` cutoff-pencil
eigendecomposition is `O(n³)` and caps practical meshes there. The p=2
higher-order advantage + the ε_r=10.2 high-contrast convergence need
finer meshes (ADR-0057 as-built). The dense cap also limits the *whole*
cross-section eigensolver, not just the ε_r=10.2 stretch.

The obstacle to a sparse selection is the curl-free **gradient null
cluster** at `k_c² ≈ 0` (enlarged at p=2), which dominates a naive sparse
solve near zero; the physical dominant mode sits at a small but strictly
positive `k_c²` above it.

## Decision

Replace the dense cutoff eigendecomposition with a **sparse
shift-and-invert of the cutoff pencil at a small positive shift σ** —
placed above the gradient floor (`k_c² ≈ 0`) and below the physical
dominant cutoff (σ ~ a fraction of the analytic air-cutoff `(π/a)²`) — to
return the low-cutoff *physical* candidates without forming the dense
spectrum and without the gradient cluster. Hand those candidates to the
**unchanged step-5.6 selection** (β-direct shift-invert + converged-
eigenvector transverse screen + highest-β²). Keep the dense path as a
small-`n` fallback/reference; the sparse path is production for large `n`
and must reproduce the dense dominant mode at the validation meshes.
Reuse the step-4 `LobpcgEigen` + step-5.3 faer sparse-LU machinery (no
new dependency).

## Rationale

(1) **Foundational, not stretch-case polish.** The `O(n³)` dense cap
limits every cross-section solve to coarse meshes; sparse selection
unlocks realistic mesh sizes generally. The ε_r=10.2 closure is one
beneficiary; mesh scaling for all contrasts is the broader value.

(2) **Cleanly dispatchable + leverages prior work.** The mechanism
(shift-invert + LOBPCG + faer sparse LU) is exactly what steps 4 and 5.3
already built and validated; this applies it to the cutoff pencil. Clear
DoD (dense-agreement at existing meshes + finer-mesh ε_r=10.2). This is
the decision-rule-favoured choice (value + dispatchability) over the
grind-risky breadth alternatives.

(3) **Gradient cluster handled by a positive shift** — the standard
spectral-transformation remedy; gradient-deflation (discrete-gradient
range projection) is the escape-hatch if the shift alone is insufficient.

## Consequences

* `cutoff_candidates` gains a sparse shift-invert path; the dense path is
  retained as the small-`n` reference (DoF-thresholded).
* If finer meshes close ε_r=10.2 ≤5%, the §4 inhomogeneous gate is closed
  across contrasts; else a deeper-than-discretization finding → step-5.8.
* The whole cross-section eigensolver scales past the ~457 DoF cap.
* The breadth-rotation tracks remain documented-as-grind-risky.

## References

* Saad (shift-invert). Knyazev 2001 (LOBPCG). Jin §9 (Nedelec null
  space / deflation). ADR-0054/0056/0057. Step 5.7 spec + plan.
* `crates/yee-fem/src/solve.rs` (pattern), `crates/yee-mom/src/eigensolver/solve.rs`.
