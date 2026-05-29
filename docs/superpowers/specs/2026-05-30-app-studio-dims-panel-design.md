# App — yee-studio physical-dimensions panel — Design Spec

**ADR:** ADR-0099 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Show F1.2.0's physical dimensions live in `yee-studio` — a numeric panel, no
graphical canvas. Keep `StudioState` egui-free + WASM-safe (ADR-0089).

## Changes (`crates/yee-studio/**` ONLY)
- `Cargo.toml`: add `yee-layout = { workspace = true }` (serde-only, WASM-safe).
- `src/lib.rs` (`StudioState`):
  - Add fields `eps_r: f64`, `h_m: f64` (defaults `4.4`, `1.6e-3` — set in the
    constructor / `from_spec`).
  - Add a derived field `dims: Result<yee_filter::EdgeCoupledDimensions, String>`.
  - In `recompute`/`apply_derived` (wherever the synthesized `project` is
    available), compute `dims = dimension_edge_coupled(&project,
    &yee_layout::Substrate{ eps_r, h_m, .. }).map_err(|e| e.to_string())`. Read
    the `Substrate` field names/units from `yee-layout`. Keep `StudioState` free
    of any `egui`/`eframe` type (the `Result<_, String>` keeps it so).
- `src/app.rs` (behind `#[cfg(any(feature="desktop", feature="web"))]`): a
  "Physical dimensions" section — editable `εr` and `h` (mm) inputs that update
  `StudioState` + trigger `recompute`, and a read-out of `line_width`,
  `resonator_length`, and each gap (mm), or the error string on `Err`.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-studio --all-targets -- -D warnings` exit 0 (native
   default features).
3. `cargo test -p yee-studio` exit 0 (headless lib tests; no GUI, no FDTD).
4. **`studio_state_dims` test (lib):** default Chebyshev N=5 spec + FR-4 substrate
   → `state.dims` is `Ok` with `line_width_m > 0`, `resonator_length_m > 0`, every
   `gaps_m[i] > 0`.
5. **WASM-safety invariant holds:** `cargo check -p yee-studio
   --no-default-features --target wasm32-unknown-unknown` exit 0 and egui absent
   from that tree (`cargo tree -p yee-studio --no-default-features --target
   wasm32-unknown-unknown -i egui` errors "not found"). `StudioState` + `dims`
   compile egui-free.

## Out of scope
Graphical polygon/SVG canvas; Gerber/KiCad; `qe`→feed gap (F1.2.1); the wasm
`--features web` UI rendering of the panel (the panel code compiles under `web`
but is exercised only on native here). No `dimension_*` algorithm change.
