# ADR-0203: R.5c — the studio designs antennas too

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0198/0199 (the filter panel + verify loop whose pattern this
applies), ADR-0190..0193 (the A-track fixture and findings this packages).

## Decision

The studio gains the antenna panel (`antenna.rs` + `AntennaDesignPanel`),
completing the project's original goal statement — *design antennas AND
filters with the engine* — at the UI:

- **`design_antenna`** (instant): Balanis patch dims + inset-fed layout +
  byte-checked Gerber artifacts. The default inset is **0.25·L, the
  A.3-measured optimum** — not the closed-form seed (whose G₁-only slot
  model overestimates R_edge, the documented A.1 finding); the fraction is
  an exposed knob.
- **`verify_antenna`** (one solve): the A.1 single-run directional |S11|
  under the A.2 open-top per-face CPML boundary, streaming
  `verify://progress`, returning the curve + the dip (f, depth). Rendered
  by `SparamPlot`, which learned single-trace mode.

## Gates

`studio-antenna-e2e-001` (headless, in `studio-build` CI): design dims equal
the Balanis closed forms and the default inset is 0.25·L; Gerbers
byte-identical to `yee_export` for the same layout; unphysical inset
rejected with a designer-grade message; a reduced-fidelity verify-pipe run
streams progress and returns a finite, ≤ 0 dB curve (physics stays gated by
engine-antenna-001..004). Plus vitest DOM gates (16 total).

## Consequences

Spec → dims → Gerber → measured |S11| for patches, next to the filter flow.
Follow-ons: the far-field pattern cut in the verify response (the A.2 NTFF
machinery is protocol-ready, CPU-only), and an inset design-loop button
(the A.3 scan, ~5 solves).
