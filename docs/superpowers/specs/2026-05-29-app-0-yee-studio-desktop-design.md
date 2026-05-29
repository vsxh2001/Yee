# App.0 — `yee-studio` desktop skeleton — Design Spec

**Phase:** App.0 · **ADR:** ADR-0090 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
A native `eframe` desktop app wiring the shipped light flow (spec → synthesis →
spec-mask plot) into stage-gated panels. The first product increment of the
filter-design app (ADR-0089). No EM, no server, no new external dependency.

## Crate `yee-studio` (lib + bin)
`Cargo.toml`: deps `yee-synth`, `yee-filter`, `egui`, `eframe` (`["wgpu"]`),
`egui_plot` (all workspace), plus `num-complex`. `[lints.rust] unsafe_code =
"forbid"`, `missing_docs = "warn"`. `[[bin]] name = "yee-studio"`, path `src/main.rs`.

### `src/lib.rs` — testable state (NO egui types)
```rust
/// Editable design state + everything derived from the current spec.
pub struct StudioState {
    pub spec: yee_filter::FilterSpec,
    // derived (recomputed by `recompute`):
    pub project: yee_filter::FilterProject,
    pub freqs_hz: Vec<f64>,
    pub s21_db: Vec<f64>,
    pub mask_regions: Vec<MaskRegionView>,   // local plain struct (f_lo,f_hi,kind,limit)
    pub mask_pass: bool,
    pub mask_notes: Vec<String>,
}
impl StudioState {
    pub fn from_spec(spec: yee_filter::FilterSpec) -> Self;  // builds + recomputes
    pub fn recompute(&mut self);  // synthesize + sweep + s21_db + regions + check_mask
}
/// Plain mask-region view (kept egui-free + WASM-safe; mirrors the F0.2 mapping:
/// passband Floor at -ripple, per-stopband Ceiling at -reject over +-2%).
pub struct MaskRegionView { pub f_lo_hz: f64, pub f_hi_hz: f64, pub floor: bool, pub limit_db: f64 }
```
`recompute`: `synthesize(&spec)` → `freqs` (reuse a sweep like yee-cli's 401-pt
`f0·(1±6·fbw/2)`) → `ideal_response` → `s21_db` (20·log10, floored 1e-12) →
`mask_regions` (closed-form from `spec.mask`) → `check_mask` → `mask_pass`/notes.

### `src/app.rs` — `StudioApp` (`impl eframe::App`)
- Left `SidePanel`: editable spec fields (f0 GHz, FBW, order, ripple dB,
  return-loss dB, ±stopband points, approximation Butterworth/Chebyshev). On any
  change → `state.recompute()`.
- Central `CentralPanel`: g-values, the coupling matrix (grid), Qe_in/Qe_out, and
  a coloured **PASS/FAIL** mask verdict + notes.
- An `egui_plot::Plot`: `s21_db` vs frequency (GHz) line; shade each
  `MaskRegionView` as a translucent box on its forbidden side.

### `src/main.rs`
Thin: `eframe::run_native("Yee Filter Studio", opts, |_cc| Ok(Box::new(StudioApp::default())))`.
`StudioApp::default()` seeds a satisfiable Chebyshev BPF (reuse the
`cheb_bpf.toml`-equivalent values).

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-studio --all-targets -- -D warnings` exit 0.
3. `cargo build -p yee-studio` exit 0 (the windowed bin builds; not run in CI).
4. `cargo test -p yee-studio` exit 0. Test `studio_state_recompute_*`:
   - a satisfiable Chebyshev 0.5 dB N=5 BPF (f0=2 GHz, FBW=0.10, stopband
     (2.4 GHz, 40 dB)) → `state.mask_pass == true`; `project.prototype.g[1..=5]`
     within 1e-3 of `[1.7058,1.2296,2.5408,1.2296,1.7058]`; `s21_db.len() ==
     freqs_hz.len()`; `mask_regions` has one Floor (≈[1.9,2.1] GHz, −0.5 dB) and
     one Ceiling (≈[2.352,2.448] GHz, −40 dB).
   - a too-low order (N=2) on the same mask → `mask_pass == false` (negative control).
5. Workspace `Cargo.toml` gains `"crates/yee-studio"`.

## Out of scope
Layout preview (needs F1.2); EM/server (App.2); WASM build (App.1); save/load +
export (App.3). Keep `StudioState` egui-free + WASM-safe.
