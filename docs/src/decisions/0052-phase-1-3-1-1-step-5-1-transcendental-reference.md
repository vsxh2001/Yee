# ADR-0052 — Phase 1.3.1.1 step 5.1: published transcendental reference for the slab-loaded-guide gate

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5.1 (validation hardening)

## Context

Step 5 shipped the inhomogeneous cross-section eigensolver with its
inhomogeneous β gated only by a monotonic physics bracket
(`β_air < β_loaded < β_full`) + a self-referential regression value —
the DoD-V2′ escape-hatch, taken because a first attempt at the published
transcendental reference found no root corroborating the numerical β.
CLAUDE.md §4 requires a published-benchmark validation case; the bracket
is a necessary-condition sanity bound, not a benchmark. The numerical
machinery is otherwise well-validated (homogeneous canary to 4e-14,
coupling sign/scale pinned by an independent-quadrature unit test), so
the non-corroboration most likely reflects a *reference*-implementation
error (LSE/LSM mode family, bracket), not a solver error — but that must
be demonstrated, not assumed.

## Decision

Implement the **slab-loaded rectangular-waveguide transverse-resonance
transcendental dispersion** (Pozar §3.6 / Collin §6) as the published
reference, **verify it independently against a textbook-tabulated value**
before comparing, then reconcile with the numerical β:

* if they agree (≤5%), the transcendental comparison becomes the
  primary published-benchmark gate and the §4 gap closes;
* if they disagree *after the reference is independently verified*, the
  discrepancy is a root-caused **finding** (not a tolerance relaxation):
  the V2′ bracket stays as the floor, the reference comparison ships as
  a reported non-failing diagnostic, and step-5.2 is queued.

An unverified reference is **not** shipped — it would be worse than the
honest bracket gate.

## Rationale

(1) The bracket gate is genuinely weak (an 80%-wide window passes almost
any formulation error). Closing the §4 gap requires an independent
published value, and the slab-loaded guide has one.

(2) Verifying the reference against a textbook-tabulated number *before*
comparing to our solver is the only way to distinguish "reference wrong"
(the likely prior-attempt failure mode) from "solver wrong". Skipping
that step is what made the first attempt inconclusive.

(3) Treating a verified-reference disagreement as a finding rather than a
fallback preserves the §4 discipline: the validation question is
answered honestly either way, and a solver-accuracy issue (if any) is
surfaced rather than masked by a loose bracket.

## Consequences

* New reference dispersion solver (`reference.rs` or a test helper) +
  an independent textbook-value unit test.
* The horizontal-slab gate either upgrades to a published benchmark or
  carries a documented, root-caused solver-vs-reference diagnostic.
* If a solver-accuracy issue surfaces, it is a separate reviewed change
  (out of this step's lane).
* The vertical-slab LSM dual may need its own dispersion (step-5.2).

## References

* Pozar, *Microwave Engineering* 4th ed., §3.6.
* Collin, *Field Theory of Guided Waves* 2nd ed., §6.
* ADR-0051 (step-5 mixed solver + DoD-V2′ escape-hatch).
* Step 5.1 spec + plan
  `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-1-transcendental-reference-design.md`,
  `docs/superpowers/plans/2026-05-23-phase-1-3-1-1-step-5-1-transcendental-reference.md`.
