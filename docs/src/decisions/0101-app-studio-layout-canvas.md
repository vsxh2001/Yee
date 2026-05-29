# ADR-0101: App — yee-studio layout preview canvas (egui_plot top view)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0099 (studio dims panel), ADR-0097 (F1.2.0 →
`dimension_edge_coupled_layout`), ADR-0086 (`yee-layout` `Layout`/`Polygon`),
ADR-0089 (WASM-safe light flow), `FILTER-DESIGN-ROADMAP.md`

---

## Context

The studio shows the synthesized `|S21|`/mask plot (ADR) and a numeric dimensions
panel (ADR-0099). The user's goal is the interactive app; the natural next visual
is to **see the actual filter geometry** — the edge-coupled microstrip layout
F1.2.0 produces. `yee-layout::Layout` already carries the trace polygons, and
`app.rs` already draws filled `egui_plot::Polygon` items (the spec-mask boxes), so
a top-view layout canvas reuses an in-place pattern with no new dependency beyond
the `yee-layout` dep the studio already has.

## Decision

Add a **layout preview canvas** to `yee-studio`:

- `StudioState` (egui-free) gains a derived
  `layout: Result<yee_layout::Layout, String>`, computed in `apply_derived`
  alongside `dims` via `yee_filter::dimension_edge_coupled_layout(&project,
  &substrate)` (`Err` mapped to a display string — keeps `StudioState`
  `egui`-free + WASM-safe).
- `app.rs` (behind `desktop`/`web`) adds a "Layout" panel: an `egui_plot::Plot`
  with equal data aspect (`data_aspect(1.0)`) drawing each `layout.traces`
  polygon as an `egui_plot::Polygon` (coordinates in mm, top view), or the error
  string on `Err`.

No new dependency (`yee-layout` already present); no `StudioState` egui leak; the
`Layout` type is WASM-safe.

## Consequences

**Ships:** a live top-view layout canvas. Editing the spec or substrate
re-derives and re-draws the geometry.

**Gate (`yee-studio` lib tests — headless, no GUI, no FDTD):** `studio_state_layout`
— the default Chebyshev N=5 spec + FR-4 substrate yields `layout` `Ok` with ≥ 1
trace polygon (and the polygons are non-degenerate). The WASM-safety invariant
(App.1.0/1.1) still holds: `--no-default-features --target wasm32` compiles with
egui absent and `StudioState`/`layout` egui-free.

**Not in scope:** substrate/ground rendering, dimension annotations/rulers,
export buttons, 3-D — later increments. The drawing itself is exercised only on
native (the panel compiles under `web` but the headless gate tests `StudioState`).

---

## References
- ADR-0097 (`dimension_edge_coupled_layout`); ADR-0099 (dims panel sibling).
- `docs/superpowers/specs/2026-05-30-app-studio-layout-canvas-design.md`;
  `docs/superpowers/plans/2026-05-30-app-studio-layout-canvas.md`.
