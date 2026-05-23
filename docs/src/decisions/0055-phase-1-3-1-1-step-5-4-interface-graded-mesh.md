# ADR-0055 — Phase 1.3.1.1 step 5.4: interface-graded mesh before p-refinement

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.4 (high-contrast inhomogeneous accuracy)

## Context

Step 5.3 closed the §4 inhomogeneous gate at FR-4 (ε_r=4.4, 1.39%) and
proved the residual ε_r=10.2 gap (β=489 vs reference 583, ~16%) is
**first-order-element discretization** of the dielectric-interface field
peak — not the β-extraction (uniform anchor exact) nor the eigenvector
mismatch (only ~1%). β plateaus under *uniform* refinement.

To resolve a discretization-limited interface field peak there are two
textbook levers: `h`-refinement (smaller elements where the field varies
fastest) and `p`-refinement (higher polynomial order). `p`-refinement on
curl-conforming Nedelec elements is a large implementation (second-order
edge elements, non-constant curl → embedded quadrature, new DoF
bookkeeping).

## Decision

Try **interface-graded `h`-refinement first** — a `TriMesh2D` builder
clustering rows geometrically toward the dielectric interface, reusing
the existing first-order element matrices and the step-5.3 sparse solve.
Only if graded `h` plateaus short of ≤5% do we escalate to
`p`-refinement (step-5.5), with the `h`-plateau as the evidence that
justifies the larger investment.

## Rationale

(1) **Cheap-first / walking-skeleton.** Graded `h` needs only a
non-uniform mesh generator + the solver that already exists; `p` needs a
new element family. The diagnosed cause (interface field peak
under-resolution) is exactly what graded `h` targets. Spend the small
effort before the large one.

(2) **Leverages step 5.3.** The sparse shift-and-invert made finer
meshes tractable; graded `h` is the natural consumer of that capability.

(3) **Either outcome is progress.** If graded `h` closes ≤5%, §4 is
closed at high contrast cheaply. If it plateaus, the improvement curve is
direct evidence for the `p`-refinement decision (step-5.5) — converting a
vague "needs higher-order" into a quantified one.

## Consequences

* New graded-mesh builder in `mesh.rs`; the uniform builders stay
  (additive — no regression to FR-4/uniform/homogeneous gates).
* High-DoF convergence may run as a `--release` example rather than an
  inline test (the dense cutoff-pencil selection caps inline meshes); a
  sparse cutoff selection is deferred unless the inline gate needs it.
* `h` may not suffice → step-5.5 (p-refinement) queued with evidence.
* `solve.rs` / `assembly.rs` / `reference.rs` untouched (graded `h`
  reuses them); a solver/assembly change would be a separate finding.

## As-built outcome (2026-05-23, merge `95ac64f`)

**Graded `h` PLATEAUS — the cheap path is ruled out.** At ε_r=10.2 the
best graded mesh (r=1.5) reached β=491.71 vs reference 582.95 (rel
15.65%), only ~0.5pp better than the uniform plateau (489.03, 16.11%),
and adding DoF at fixed grading drifts *back* toward ~487 — not toward
583. So the residual is confirmed the **first-order-element convergence
rate** at the interface field peak, which `h` cannot fix. No failing
gate added (the ε_r=10.2 reconciliation stays a non-failing diagnostic);
FR-4 + uniform + homogeneous gates stand, all values bit-identical
(additive). **step-5.5 = p-refinement (second-order Nedelec)** is now
justified by data — but it is a large implementation for a marginal-case
(high-contrast) payoff, so it is **deprioritized**: the cross-section
eigensolver is production-validated at FR-4, and the recommended next
move is to rotate to a higher-breadth ROADMAP track rather than continue
this chain. CI note: the inline convergence study lifts
`eigensolver_inhomogeneous` to ~85 s (dense cutoff-pencil selection) — a
sparse selection / release-example move is a hygiene follow-on.

## References

* Jin, *FEM in EM* 3rd ed. §9 (h/p refinement, edge elements).
* Babuška & Suri (h-p convergence rates).
* ADR-0054 as-built (the discretization diagnosis + plateau evidence).
* Step 5.4 spec + plan.
* `crates/yee-mom/src/eigensolver/{mesh,solve,reference}.rs`.
