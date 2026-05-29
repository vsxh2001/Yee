# ADR-0099: App — yee-studio physical-dimensions panel (F1.2.0 in the studio)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0090/0092 (yee-studio App.0/1.0), ADR-0097 (F1.2.0 dimensional
synthesis), ADR-0089 (app architecture — WASM-safe light flow),
`FILTER-DESIGN-ROADMAP.md`

---

## Context

The studio (`yee-studio`) takes a spec → synthesis → `|S21|`/mask verdict, live.
The layout/dimensions preview was explicitly **deferred pending the F1.2 dims
mapping** — which F1.2.0 (ADR-0097) now provides. So the studio can finally show
the **physical microstrip dimensions** a synthesized filter implies, closing the
"spec → geometry" loop in the product surface.

## Decision

Wire F1.2.0 into `yee-studio` as a **numeric dimensions panel** (no graphical
canvas — that is a later increment):

- `StudioState` (the egui-free logic layer) gains substrate fields `eps_r: f64`
  and `h_m: f64` (FR-4 defaults `4.4` / `1.6e-3`), and a derived
  `dims: Result<EdgeCoupledDimensions, String>` computed in `recompute` by
  calling `yee_filter::dimension_edge_coupled(&project, &Substrate{..})`
  (`Err` mapped to a display string so `StudioState` stays `egui`-free and
  WASM-safe).
- `app.rs` (behind the `desktop`/`web` feature) renders a "Physical dimensions"
  section: editable `εr` / `h` inputs + the computed line width, resonator
  length, and per-section gaps (or the error string when the coupling is
  unrealizable on the chosen substrate).

`yee-studio` gains a `yee-layout` dependency — `yee-layout` is `serde`-only and
WASM-safe, so the WASM-safety invariant (ADR-0089) holds; the new derived state
lives in the egui-free `StudioState`, not the UI layer.

## Consequences

**Ships:** a live dimensions panel in the studio. Editing the spec or the
substrate re-derives the dimensions immediately.

**Gate (`yee-studio` lib tests — headless, no GUI, no FDTD):** a
`studio_state_dims_*` test that the default Chebyshev N=5 spec + FR-4 substrate
yields `Ok(dims)` with `line_width_m > 0`, `resonator_length_m > 0`, every gap
`> 0`; and that `StudioState` still compiles/derives with
`--no-default-features` (egui-free, WASM-safe — App.1.0/1.1 invariant intact).

**Not in scope:** a graphical polygon/SVG canvas; Gerber/KiCad export; the
`qe`→feed gap (F1.2.1). Numeric panel only.

---

## References
- ADR-0097 (`dimension_edge_coupled`); ADR-0089 (WASM-safe light flow).
- `docs/superpowers/specs/2026-05-30-app-studio-dims-panel-design.md`;
  `docs/superpowers/plans/2026-05-30-app-studio-dims-panel.md`.
