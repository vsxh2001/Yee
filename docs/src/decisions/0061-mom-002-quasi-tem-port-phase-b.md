# ADR-0061 — mom-002 quasi-TEM numerical wave-port: Phase B (now unblocked)

**Status:** Accepted (production wiring + bounded-experiment continuation)
**Date:** 2026-05-24
**Context Phase:** MoM beachhead follow-on (continues ADR-0059)

## Context

ADR-0059's bounded experiment stopped at Phase A with a precise finding:
the cross-section eigensolver could not find the microstrip quasi-TEM
mode, so the `Numerical2D` wave-port could not be excited for mom-002.
ADR-0060 (Phase 1.3.1.2) closed exactly that capability gap — a
quasi-TEM selection path (`solve_dense_mixed_quasi_tem`), HJ-validated to
1.2%. Phase B is therefore now unblocked: the numerical microstrip port
can finally be wired to the mom-002 line, which was the whole point of
the ADR-0059 experiment.

Two distinct pieces remain. (A) The quasi-TEM mode is currently
lib-internal — `NumericalCrossSection::solve` only runs the closed-guide
`solve_dense_mixed` (First/Second order), so the `mode_profile` the
`Numerical2D` arm + `e_tangential_at` consume is never populated from a
quasi-TEM solve. (B) The ADR-0059 Z_in comparison (delta-gap 674 Ω vs
HJ 51 Ω) was never produced because Phase A blocked.

## Decision

Run **Phase B** as a production wiring step (Part A) + a bounded
experiment (Part B):

- **Part A (reusable capability):** add a quasi-TEM solve path to
  `NumericalCrossSection` that runs `solve_dense_mixed_quasi_tem` and
  caches `mode_profile`, making the quasi-TEM mode usable for wave-port
  excitation (mirrors the existing First-path caching; the closed-guide
  `solve` stays the bit-identical default).
- **Part B (the experiment payoff):** feed the mom-002 microstrip
  cross-section's quasi-TEM modal field to the mom-002 line via the
  `Numerical2D` arm, extract `|Z_in|`, report vs 674 Ω (delta-gap) +
  51 Ω (HJ) as a NON-FAILING diagnostic. Hard ~30-min cap on the
  cross-section→RWG modal-RHS coupling for a microstrip; if it needs glue
  that doesn't exist, document the specific blocker + stop. Do NOT
  re-open the mom-002 kernel/Greens forensics; do NOT re-gate mom-002.

## Rationale

(1) **The payoff of the whole eigensolver chain.** Part A makes the
quasi-TEM capability actually usable (not lib-internal); Part B answers
the ADR-0059 question — does a principled numerical microstrip port beat
the delta-gap on the PRIMARY beachhead gate? Highest on-mission value.

(2) **Tractable now.** `solve_dense_mixed_quasi_tem` already returns the
same `MixedEigenSolution` the First path scatters into `mode_profile`, so
Part A is a small, well-targeted consumer change; the `Numerical2D` arm +
the mom-002 experiment infra already exist.

(3) **Grind-risk still bounded.** The microstrip cross-section↔planar-MoM
coupling has real subtlety (the arm was validated on waveguide TE10), and
mom-002 is a known quagmire. Feasibility-first + the hard cap + "either
branch is a deliverable" keep it from flailing — per the standing
bounded-experiment rule.

## Consequences

* `NumericalCrossSection` gains a quasi-TEM solve path (caches
  `mode_profile`); the closed-guide First/Second `solve` stays
  bit-identical (WR-90 / FR-4 / homogeneous + the HJ quasi-TEM gate
  guard it).
* A `|Z_in|` numerical-port-vs-delta-gap-vs-HJ diagnostic, OR a
  documented coupling-glue blocker finding. mom-002 gate + 674 Ω
  tripwire untouched.
* If the numerical port clearly moves `|Z_in|` toward `Z_0`, a follow-on
  adopts it as the mom-002 production excitation + re-gates.
* If it does not (or points back at the kernel/Greens), the experiment
  confirms the residual is not (only) the port — a useful narrowing; no
  forensics re-open.

## References

* ADR-0059 (the experiment + Phase-A finding), ADR-0060 (the quasi-TEM
  capability this consumes).
* `crates/yee-mom/src/ports.rs` (`NumericalCrossSection`, `Numerical2D`
  arm, `e_tangential_at`), `crates/yee-mom/src/eigensolver/solve.rs`
  (`solve_dense_mixed_quasi_tem`), `crates/yee-mom/tests/mom_002_numerical_waveport.rs`.
* Phase-B spec + plan (2026-05-24).

## Outcome (2026-05-24, merge `966d853`)

**Both parts delivered.** Part A: `NumericalCrossSection::with_quasi_tem()`
+ `solve` dispatch shipped; the closed-guide path is reviewer-verified
byte-identical when the flag is unset; quasi-TEM + Second-order returns
`Unimplemented`. Part B: the diagnostic runs end-to-end —
`|Z_in| = 378.1 Ω` (verified in `--release`, 266.95 s), **closer to
`Z_0 ≈ 51 Ω` than the 674 Ω delta-gap baseline** (7.41× vs 13.22×; the
quasi-TEM port roughly halves the port-excitation error).

**Recommendation: residual-not-(only)-the-port + frame-mapping-glue-needed.**
The cross-section `ε_eff = 2.105` (vs Hammerstad-Jensen 3.33) reveals
*why* the gap only halves: the diagnostic builds the cross-section in the
MoM `(x, y)` longitudinal-sampling frame, **not** the physical microstrip
transverse plane, so the mode it finds is not the line's true quasi-TEM.
The principled next increment is a **2-D-cross-section → 2.5-D-RWG frame
mapping** (project the physical transverse-plane mode onto the MoM port
edges in the correct coordinate frame), tracked as a follow-on. The
mom-002 Greens kernel is *not* the bottleneck — it was independently
re-confirmed shipped (multi-image DCIM + Sommerfeld surface-wave, within
1.83% of HJ); the remaining accuracy is port-side. The diagnostic is
`#[ignore]`'d from default CI (post-eigensolve Sommerfeld fill ≈ 267 s
release / hours debug); run it with `--release -- --ignored --nocapture`.
