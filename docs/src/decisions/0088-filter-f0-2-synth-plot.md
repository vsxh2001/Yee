# ADR-0088: Filter Phase F0.2 — `yee filter synth --plot`

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0084 (F0 synthesis core), ADR-0087 (spec-mask plot overlay),
`FILTER-DESIGN-ROADMAP.md`

---

## Context

F0 (`yee filter synth`) computes a synthesized filter's closed-form ideal
response and grades it against the `FilterSpec`'s `SpecMask`, but the only
visual artifact is a Touchstone `.s2p`. 1.plotting.4 (ADR-0087) shipped
`yee-plotters::draw_sparam_with_mask` (shaded forbidden regions). Nothing yet
connects the synthesis `SpecMask` to that overlay, so the designer cannot *see*
the synthesized response against the spec. This is the "where the crates meet"
adapter flagged as pending after the parallel F0.1/F1.0/1.plotting.4 batch.

## Decision

Add a `--plot <path>` flag to `yee filter synth`. When given, after the response
sweep, render the **|S21|** magnitude (dB) with the spec mask overlaid via
`draw_sparam_with_mask`, then write the image (PNG/SVG by extension).

A small `SpecMask → Vec<MaskRegion>` adapter lives in `yee-cli` (the integration
point — keeps `yee-filter` and `yee-plotters` decoupled, no new inter-crate dep):
- **Passband** `[f1, f2] = f0·(1 ∓ FBW/2)` → one `MaskKind::Floor` at
  `limit_db = −passband_ripple_db` (|S21| must stay above −ripple in band).
- **Each stopband point** `(f_s, reject_db)` → a `MaskKind::Ceiling` at
  `−reject_db` over a ±2% band around `f_s` (|S21| must stay below −reject).

Only |S21| is plotted; return-loss / |S11| is an S11-plane constraint whose
mask would visually conflict with the S21 passband on one axis — deferred to a
future S11 plot. `yee filter synth` without `--plot` is unchanged.

## Consequences

**Ships:** `--plot` on `filter synth`; the `spec_mask_regions(&FilterSpec) ->
Vec<MaskRegion>` adapter in `yee-cli`; gates — an adapter unit test (a known
spec → the expected passband Floor + stopband Ceiling regions) and a CLI test
(`yee filter synth --plot <tmp.png>` exits 0 and writes a non-empty PNG).

**Not in scope:** S11/return-loss mask plotting; GUI; any EM; the F1.1 FDTD
coupling-extraction loop (next, separately).

**No new dependency** (`yee-cli` already deps `yee-plotters`/`yee-filter`/`yee-io`).
Lane: `crates/yee-cli/**`.

---

## References
- ADR-0084, ADR-0087; `FILTER-DESIGN-ROADMAP.md` Stage 6 / "SpecMask→MaskRegion adapter".
- `docs/superpowers/specs/2026-05-29-filter-f0-2-synth-plot-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f0-2-synth-plot.md`.
