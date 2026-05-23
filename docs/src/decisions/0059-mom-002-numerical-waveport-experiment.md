# ADR-0059 — mom-002 numerical-microstrip-wave-port bounded experiment

**Status:** Accepted (as a bounded experiment)
**Date:** 2026-05-24
**Context Phase:** MoM beachhead follow-on (post cross-section eigensolver)

## Context

The cross-section eigensolver chain (Phase 1.3.1.1 steps 4→5.8) is
complete: the quasi-TEM cross-section solver is FR-4-validated and the
`WavePort` `Numerical2D` arm can inject a cross-section modal field into
the MoM port RHS. Rotating off that subsystem, the genuinely-different
high-value tracks are all grind-risky or blocked (FDTD Q6/Q7
energy-balance; mom-002/003 forensic + Greens; Phase 4.fem.eig.1+
undrafted ↔ fem-eig-006). The project memory's standing guidance is to
approach such tracks as **tight bounded experiments** (negative result =
a deliverable), never open dives.

mom-002 (the PRIMARY MoM deliverable) passes only loosely (`|Z_in|≈674 Ω`
vs `Z_0≈51 Ω`); 10 forensic tracks exonerated the kernel and localised
the residual to delta-gap port-excitation modeling. A numerical
microstrip wave-port (the now-validated cross-section mode) is the
principled fix and a fresh angle those tracks never tried.

## Decision

Run a **tight bounded experiment**: feed a microstrip cross-section's
numerical modal field to the mom-002 line via the `Numerical2D` wave-port
and compare `|Z_in|` to the delta-gap baseline + the HJ target. The
deliverable is the *result* — either the comparison (a non-failing
diagnostic; do NOT re-gate mom-002) or, if the 2-D-cross-section↔planar-MoM
coupling cannot be wired cleanly within a hard 30-min feasibility cap, a
documented finding of what the `Numerical2D` arm lacks for microstrip
ports. Explicitly do NOT re-open the mom-002 kernel/Greens forensics or
edit the eigensolver.

## Rationale

(1) **Highest on-mission value among the rotation candidates** — mom-002
is the primary beachhead gate; the residual is the port, and this is the
principled, never-tried port fix.

(2) **Leverages the just-completed eigensolver** + is a genuinely
different subsystem (MoM port excitation, not eigensolver internals).

(3) **Grind-risk bounded by design** — the feasibility-first hard cap +
"negative result = deliverable" framing prevent a flail on the known
mom-002 quagmire / the subtle cross-formulation coupling. Per the memory
rule (value + dispatchability; tight bounded experiments for the
breadth tracks).

## Consequences

* A new non-failing diagnostic (`mom_002_numerical_waveport.rs`) OR a
  documented port-infra-glue finding; the mom-002 gate + tripwire band
  are untouched (this is an experiment, not a re-gate).
* If the numerical port clearly improves `|Z_in|` toward `Z_0`, a
  follow-on track adopts it as the mom-002 production excitation.
* If it does not (or points back at the kernel/Greens), the experiment
  documents that the residual is not (only) the port — itself a useful
  narrowing.
* The eigensolver is consumed read-only; no kernel/Greens re-analysis.

## References

* ADR-0036 / ADR-0037 (mom-002 validation reframe + R1 retraction).
* The completed cross-section eigensolver (ADRs 0050–0058).
* `crates/yee-validation/src/lib.rs`, `crates/yee-mom/src/ports.rs`,
  Track GGGGGGGG (the `Numerical2D` arm).
* Experiment spec + plan (2026-05-24).
