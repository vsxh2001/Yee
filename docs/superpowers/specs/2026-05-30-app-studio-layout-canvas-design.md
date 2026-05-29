# App — yee-studio layout preview canvas — Design Spec

**ADR:** ADR-0101 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
A top-view layout preview in `yee-studio`: draw the F1.2.0 `Layout` trace polygons
in an `egui_plot::Plot`. Keep `StudioState` egui-free + WASM-safe. Reuses the
existing `egui_plot::Polygon` pattern (the spec-mask boxes in `show_response_plot`).

## Changes (`crates/yee-studio/**` ONLY)
- `src/lib.rs` (`StudioState`):
  - Add derived field `pub layout: Result<yee_layout::Layout, String>`.
  - Initialise it in `from_spec` (like `dims: Err(String::new())`) and compute it
    in `apply_derived`, next to `dims`:
    `self.layout = dimension_edge_coupled_layout(&self.project, &substrate).map_err(|e| e.to_string());`
    (the `substrate` local already built there for `dims`). Import
    `dimension_edge_coupled_layout`. Keep `StudioState` free of egui/eframe types.
- `src/app.rs` (behind `#[cfg(any(feature="desktop", feature="web"))]`):
  - A `show_layout(ui, state)` fn: an `egui_plot::Plot::new("studio_layout_plot")`
    with `.data_aspect(1.0)` (equal aspect — geometry not distorted) + a legend
    optional. For each polygon in `state.layout`'s `traces`, build
    `PlotPoints::from(verts)` with vertices in **mm** (`p.x*1e3`, `p.y*1e3` — read
    `Point2` from yee-layout; coords are metres) and draw
    `plot_ui.polygon(Polygon::new(<name>, pts))`. On `Err`, show the error label.
  - Call `show_layout` from the central panel (same place `show_dimensions` /
    `show_response_plot` are invoked).

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-studio --all-targets -- -D warnings` exit 0 (native).
3. `cargo test -p yee-studio` exit 0 (headless lib tests; no GUI/FDTD).
4. **`studio_state_layout` test (lib):** default Chebyshev N=5 spec + FR-4 →
   `state.layout` is `Ok`, with `traces` non-empty and each polygon having ≥ 3
   vertices.
5. **WASM-safety invariant holds:** `cargo check -p yee-studio
   --no-default-features --target wasm32-unknown-unknown` exit 0 and egui absent
   (`cargo tree … -i egui` errors "not found"). `StudioState` + `layout` egui-free.

## Out of scope
Substrate/ground-plane fill, dimension rulers/labels, export buttons, 3-D, hover
tooltips. The `egui_plot` API specifics (0.35) are the implementer's to get right;
if `Polygon` fill is awkward, fall back to `Line` outlines per polygon. No
`dimension_*` change.
