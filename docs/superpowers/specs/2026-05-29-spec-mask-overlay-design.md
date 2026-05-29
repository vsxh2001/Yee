# Phase 1.plotting.4 — S-parameter spec-mask overlay — Design Spec

**Phase:** 1.plotting.4 · **ADR:** ADR-0087 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
Overlay a spec mask (forbidden regions) on `yee-plotters` S-parameter magnitude
plots so spec compliance is visible at a glance, plus a pure `mask_violations`
helper for pass/fail. `yee-plotters` only; no new dependency.

## API (yee-plotters, new public items)
```rust
pub enum MaskKind { Ceiling, Floor }   // trace must stay below / above limit
pub struct MaskRegion { pub f_lo_hz: f64, pub f_hi_hz: f64,
    pub kind: MaskKind, pub limit_db: f64 }

/// Indices of `trace_db` samples (paired with `freqs_hz`) that violate ANY
/// region: a sample at freq in [f_lo,f_hi] is a violation if it is ABOVE a
/// Ceiling limit or BELOW a Floor limit. Pure; no I/O.
pub fn mask_violations(freqs_hz: &[f64], trace_db: &[f64], regions: &[MaskRegion]) -> Vec<usize>;

/// Render labeled dB traces with mask regions shaded on their forbidden side;
/// traces drawn on top. PNG or SVG by extension, matching existing draw_* fns.
pub fn draw_sparam_with_mask(
    path: &std::path::Path,
    freqs_hz: &[f64],
    traces: &[(&str, &[f64])],   // (label, magnitude in dB)
    regions: &[MaskRegion],
    title: &str,
) -> Result<(), <existing yee-plotters error type>>;
```
Match the existing `yee-plotters` draw-fn signatures / error type / multi-trace
colour cycle (read the crate before writing). `#![forbid(unsafe_code)]` /
`#![warn(missing_docs)]` already set on the crate; doc the new public items.

## Behaviour
- Frequency axis spans `freqs_hz`; dB y-axis auto-ranged over the traces with a
  small margin (reuse the existing magnitude-plot range logic if present).
- For each `MaskRegion`: shade a translucent red rectangle over the **forbidden**
  area — `Ceiling`: from `limit_db` up to the plot top, across `[f_lo,f_hi]`;
  `Floor`: from plot bottom up to `limit_db`. Draw regions BEFORE traces so the
  data sits on top.
- `mask_violations`: for each sample whose freq ∈ [f_lo,f_hi] of a region,
  flag if `kind==Ceiling && db>limit_db` or `kind==Floor && db<limit_db`. Return
  sorted unique indices.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-plotters --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-plotters` exit 0.
4. `mask_violations` unit test: a trace with one sample at −15 dB inside a
   `Ceiling{limit_db:−20}` region over its freq span returns that sample's index;
   an all-compliant trace returns `vec![]`; a `Floor` case symmetric.
5. Render smoke test (matching the ADR-0081 VSWR test): `draw_sparam_with_mask`
   to a temp PNG returns `Ok`, the file exists and is non-empty (> 0 bytes).

## Out of scope
`yee-gui` integration (avoid the wgpu build); CLI wiring; any `yee-filter`
coupling. Adapter from `yee-filter::SpecMask` → `Vec<MaskRegion>` is a later
increment where the crates meet.
