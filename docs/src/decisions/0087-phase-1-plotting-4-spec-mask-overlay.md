# ADR-0087: Phase 1.plotting.4 — S-parameter spec-mask overlay (`yee-plotters`)

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0063 (multi-trace S-param plotting), ADR-0069 (Smith CLI
multi-trace); `FILTER-DESIGN-ROADMAP.md` Stage 6 (full-wave verification)

---

## Context

The filter design flow's verification stage (FILTER-DESIGN-ROADMAP Stage 6)
compares an EM-simulated S-parameter response against a **spec mask** — passband
floors/return-loss and stopband rejection limits. `yee-plotters` can draw
multi-trace S-parameter magnitude plots (ADR-0063) but cannot overlay the
mask's forbidden regions, so a designer can't see at a glance whether the
response passes. This is a reusable plotting primitive, useful well beyond the
filter flow (any spec-compliance plot).

## Decision

Add a **spec-mask overlay** to `yee-plotters` (no new dependency — `plotters`
is already in tree):

- A plain, decoupled mask type:
  `MaskRegion { f_lo_hz, f_hi_hz, kind: MaskKind, limit_db }` where
  `MaskKind ∈ { Ceiling, Floor }` — `Ceiling` = the trace must stay **below**
  `limit_db` in `[f_lo, f_hi]` (e.g. stopband rejection); `Floor` = must stay
  **above** (e.g. passband min |S21|). Plain f64 data — NOT coupled to
  `yee-filter`'s `SpecMask`, so any caller can use it.
- `draw_sparam_with_mask(path, freqs_hz, traces, regions, opts) -> Result<()>`:
  draws each labeled dB trace (reusing the existing multi-trace style) with the
  forbidden side of each `MaskRegion` shaded (translucent red box from the limit
  to the plot edge across the region's frequency span), traces drawn on top.
- A `mask_violations(freqs_hz, trace_db, regions) -> Vec<usize>` helper returns
  the sample indices that violate any region (pure, unit-testable — also handy
  for a pass/fail badge).

## Consequences

**Ships:** the `MaskRegion`/`MaskKind` types, `draw_sparam_with_mask`, and
`mask_violations` in `yee-plotters`. Gates: a render smoke test (PNG produced,
non-empty — matching the ADR-0081 VSWR render-test precedent) and a
`mask_violations` unit test (a trace dipping into a Ceiling region returns the
expected violating indices; a compliant trace returns empty).

**Not in scope:** GUI (`yee-gui`) integration — deferred to keep this off the
wgpu build; CLI wiring; coupling to `yee-filter::SpecMask` (a thin adapter can
be added later where both crates meet).

**No new dependency.** Lane: `crates/yee-plotters/**`.

---

## References
- ADR-0063, ADR-0069, ADR-0081 (plotters render-test precedent).
- `FILTER-DESIGN-ROADMAP.md` Stage 6; spec/plan dated 2026-05-29.
