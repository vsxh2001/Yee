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

## As-built + step-5.8 resolution (2026-05-24, merges `1db2f51`, `157a401`)

**Sparse selection shipped** (merge `1db2f51`): σ-ladder + positive-k_c²
union (no deflation needed), dense-agreement bit-identical at all
validation meshes, 20×20 (n≈1242) in ~9 s, lib suite 78s→17s.

**The ε_r=10.2 gap is a MODE-FAMILY mismatch, not a solver defect — the
cross-section chain closes here.** With the dense cap lifted, p=1 AND
converged-p=2 β both plateau at ≈485 (ε_eff≈5.79), flat across
8×8→24×24 — retiring the discretization hypothesis. step-5.8's
multi-root census (`157a401`) then enumerated the LSM-to-y transcendental
roots: **{582.95, 216, 161, 158}** — the FEM's β≈485 is **not** an
LSM-to-y root (≈17% from the nearest), but sits near the **LSE-to-y**
dominant root (≈465, verified step-5.1). So the FEM's dominant mode is
~LSE-to-y, while the reference (and the "16.6% gap") compared **LSM-to-y**
(583) — ADR-0052's LSM-to-y family assignment was made from the step-5.2
**contaminated** cutoff eigenvector (the one step-5.6 found non-transverse).

**Disposition:** **FR-4 (ε_r=4.4, 1.39%) stands as THE §4 inhomogeneous
closure.** The ε_r=10.2 high-contrast comparison is a **documented
mode-family-identity caveat** (FEM dominant ≈ LSE-to-y 465 vs the
LSM-to-y 583 the reference reports) — not a solver accuracy bug (the FEM
is mesh-converged, p1≡p2, all element matrices / β-extraction / selection
reviewer-validated). The cross-section eigensolver (steps 4→5.8) is
**production-complete**: validated at homogeneous / uniform / FR-4,
p=2 elements + sparse mesh-scaling + Python bindings, with the
high-contrast modal-family identity as the one documented open question
(a future modal-classification study if a high-contrast use case demands
it — not a continuation of this chain). Review P1-1 (n>260 dense-sparse
agreement guard, `#[ignore]`'d) + P1-2 (recover_k0_sq invariant) landed
in `157a401`.

## References

* Saad (shift-invert). Knyazev 2001 (LOBPCG). Jin §9 (Nedelec null
  space / deflation). ADR-0054/0056/0057. Step 5.7 spec + plan.
* `crates/yee-fem/src/solve.rs` (pattern), `crates/yee-mom/src/eigensolver/solve.rs`.
