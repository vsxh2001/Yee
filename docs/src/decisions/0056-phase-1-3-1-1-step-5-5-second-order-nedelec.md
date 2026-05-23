# ADR-0056 — Phase 1.3.1.1 step 5.5: second-order Nedelec elements (reverses the ADR-0055 deprioritisation)

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.5 (high-contrast inhomogeneous accuracy)

## Context

ADR-0055 (step 5.4) closed the §4 inhomogeneous gate at FR-4, proved the
ε_r=10.2 residual is first-order-element discretisation (graded `h`
plateaus), and **deprioritised** p-refinement as "marginal-case,"
recommending a breadth-rotation instead.

On re-evaluation at the cross-section-chain milestone, the breadth
rotation was investigated and found unviable: every genuinely-different
high-value track is **grind-risky and poorly dispatchable** —
- **mom-002 numerical-microstrip-port:** a 10-track forensic quagmire
  entangled with the `MultilayerGreens` placeholder + a subtle
  2D-cross-section→planar-MoM coupling (the Numerical2D arm was validated
  only for homogeneous waveguide-TE10, not microstrip);
- **FDTD Q6/Q7:** an open energy-balance-closure problem (75-79% drift);
- **FEM real-waveguide-port:** entangled with the deprioritised
  fem-eig-006 modal-projection grind.
An agent dispatched on any of these would likely flail.

## Decision

**Reverse ADR-0055's deprioritisation: do step-5.5 (second-order
Nedelec/nodal elements) now.** It is the one **high-confidence,
cleanly-dispatchable, won't-grind** option remaining: the fix is textbook
(Jin §9.4, Webb hierarchal elements), the DoD is crisp (ε_r=10.2 ≤5% vs
the verified reference), step 5.4 *proved* it is the right lever, and the
interface-aligned mesh guarantees fast p-convergence (the slab interface
is a straight, element-aligned line — each element is within one smooth
material, no corner singularity).

## Rationale

(1) **Dispatchability is the deciding factor.** In an autonomous loop, a
marginal-but-certain, cleanly-executable win beats a high-variance
quagmire gamble that would burn many iterations flailing. step-5.5 has a
clear textbook target + a sharp correctness anchor (p=2 reproduces
analytic TE10 on the homogeneous guide); the breadth alternatives are
under-specified minefields.

(2) **Not a thrash.** This reverses a one-tick-old decision, but on the
basis of *new analysis* (the breadth alternatives were investigated this
tick and found grind-risky/un-dispatchable), not vacillation. The
ADR-0055 deprioritisation assumed clean breadth options existed; they do
not.

(3) **Genuine value beyond ε_r=10.2.** Higher-order elements are a
reusable FEM capability; the high-contrast close is the concrete
validation, but the element family serves future high-accuracy needs.

(4) **Bounded.** The plan front-loads the correctness anchors (J1
independent-quadrature element-matrix pin, J3 homogeneous-TE10 anchor)
before the high-contrast case, with a firm escape-hatch: a broken
element formulation stops early; a p=2-still-short result is a documented
finding (the p-convergence curve vs the p1 plateau is valuable evidence
either way), not a grind.

## Consequences

* New p=2 element matrices + a triangle Gauss rule + a higher-order DoF
  map; first-order stays the default (no regression to FR-4/uniform/
  homogeneous gates — p=2 is selected only for the high-contrast case).
* Reuses the step-5.3 sparse β-direct solve order-agnostically.
* If p=2 closes ε_r=10.2 ≤5%, the inhomogeneous §4 gate is closed across
  contrasts. If it falls short, step-5.6 (deeper modeling) is queued with
  the p-convergence evidence.
* The breadth-rotation tracks (mom-002 numerical-port, FDTD Q6, real-port)
  remain open but are documented as grind-risky — to be approached with
  tight bounded-experiment framing when chosen, not as open dives.

## References

* Jin, *FEM in EM* 3rd ed. §9.4. Webb, IEEE TAP 1999 (hierarchal vector
  bases). ADR-0054/0055. Step 5.5 spec + plan.
