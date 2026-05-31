# App.2.6 — Multi-technique response overlay — Design Spec

**ADR:** ADR-0143 · **Date:** 2026-05-31 · **Status:** Accepted
**Origin:** the visual companion to the App.2.5 compare table (the maintainer's
"deepen the flows — optimize/compare" direction). The table compares board size +
metrics; this overlays the swept `|S21|` responses on one chart so the user *sees* the
response differences for the current spec.

## Problem

The Compare panel (App.2.5) tabulates per-technique board size + verdict + metrics, but
shows no **response shape**. The per-flow Synthesis stages each plot one response; there
is no single chart overlaying the techniques' responses for the current spec.

## Honesty constraint (the curves that genuinely differ)

Edge-coupled and hairpin share the **same** coupled-resonator synthesis (identical
coupling matrix → identical ideal `|S21|`); they differ only *physically* (board layout
/ size — already shown in the compare table), **not** in the ideal response. The only
genuinely distinct swept responses for a band-pass spec are therefore:
- the **coupled-resonator ideal** (`Designed.sweep`, shared by edge-coupled + hairpin),
- the **lumped realized-ladder** (`LumpedDesigned.sweep` = `ladder_s21`, distinct).

For a low-pass spec, the single **stepped-impedance** ideal (`SteppedLowpassDesigned.sweep`).
The overlay must label these **truthfully** (not pretend three separate band-pass curves)
— the distributed techniques are one shared ideal curve, the lumped is a second.

## Method

### Engine (`engine.rs`)

```rust
/// One labelled response curve for the overlay.
pub struct OverlayCurve {
    pub label: String,                 // honest: "Coupled-resonator (edge-coupled / hairpin) — ideal", etc.
    pub sweep: Vec<SweepPoint>,
    pub realizable: bool,
}
/// The distinct swept responses to overlay for `spec` (deduplicating the
/// coupled-resonator ideal shared by edge-coupled + hairpin), plus the mask.
pub fn overlay_curves(spec: &FilterSpec) -> Vec<OverlayCurve>;
```

`overlay_curves` keys on `spec.response`:
- `Bandpass | Bandstop` → `[ coupled-resonator ideal (design_demo_from(.., EdgeCoupled).sweep,
  label naming both edge-coupled + hairpin), lumped realized (design_lumped_from(spec)
  → .sweep; `realizable=false`/empty on `Err`) ]`.
- `Lowpass` → `[ stepped ideal (design_stepped_from(spec).sweep) ]`.
- `Highpass` → `[]`.

Pure; sweeps are real engine output on the shared `sweep_freqs` grid.

### SVG (`svg.rs`)

`response_overlay(curves: &[(&str, &[SweepPoint], &str /*color*/)], bands: &[MaskBand])
-> String` — mirror the existing `response_plot`: the shaded mask bands, then one
`<polyline>` of `|S21|(f)` per curve in a distinct stroke colour, plus a small legend
(label + colour swatch). Reuse `response_plot`'s axis / scaling / band code.

### UI (`stages.rs`)

Render the overlay in the Compare panel (below the table): `response_overlay` over
`overlay_curves(&spec())` with `mask_bands(&spec())`. Each curve a distinct colour, the
legend naming the technique(s) honestly. Empty (high-pass) → reuse the panel's "no live
technique" note (no empty chart).

## Changes

- `crates/yee-studio-web/src/engine.rs` — `OverlayCurve`, `overlay_curves` (pure, doc) + a test.
- `crates/yee-studio-web/src/svg.rs` — `response_overlay`.
- `crates/yee-studio-web/src/stages.rs` — render the overlay in the Compare panel.

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): `overlay_curves` on a
   band-pass demo spec returns **2** curves (coupled-resonator ideal + lumped realized),
   the two sweeps are on the same frequency grid, and they are **NOT identical** (the
   lumped realized `|S21|` differs from the coupled-resonator ideal at ≥1 sweep point — a
   real difference, the whole point of the overlay). A low-pass spec → 1 curve
   (stepped); a high-pass spec → `[]`. Each curve's `sweep` equals the corresponding
   design's `.sweep` (real, not synthetic). A constant/empty `overlay_curves` fails.
2. `response_overlay` renders one `<polyline>` per curve + a legend (assert in a unit
   test or via the built bundle); `dx build --platform web --release` EXIT 0.
3. Existing tests pass; `cargo clippy ... -D warnings` + `cargo fmt --check` clean;
   `cargo check --workspace` green.

## Out of scope

A *realized* distributed response (needs EM — the deferred ADR-0133 wall); interactive
zoom/cursor; overlaying `|S11|`. No new physics.

## Why

Completes the compare view — table (board + metrics) **+** chart (response shape) — and
honestly shows that the distributed techniques share one ideal response (differing only
physically) while the lumped realized ladder is a distinct curve, all vs the spec mask.
