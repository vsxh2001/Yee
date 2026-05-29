# ADR-0068 — Smith Chart Constant-R/X Arc Overlays + GUI Multi-trace

**Status:** Accepted  
**Date:** 2026-05-25  
**Context Phase:** 1.gui.5 / 1.plotting.2

## Context

ADR-0063 (multi-trace S-parameter plotting) and ADR-0065 (GUI multi-trace
overlay) both explicitly deferred two Smith chart enhancements as follow-ons:

1. **Constant-R/X arc reference overlay** — the family of circles/arcs that
   transform the complex Γ-plane into the readable Smith chart grid used for
   impedance-matching design (Pozar §5.1).
2. **Multi-trace support in the GUI Smith tab** — multiple S-entries on one
   canvas with a legend, analogous to what the dB tab already does.

The current `show_smith_chart` (yee-gui) and `draw_smith` (yee-plotters) render
only the outer unit circle plus the S11 trajectory — enough to show the locus
but not enough to read off impedances visually.

## Decision

Land the two deferred enhancements in Phase 1.gui.5 / 1.plotting.2:

### Track A — yee-gui (`crates/yee-gui/src/**`)

* Add pure math helpers `smith_r_circle_points(r, n)` and
  `smith_x_arc_points(x, n)` (unit-disk clipped).
* Add `SmithSeries { label, points }` + `build_smith_series(file, selection)`.
* Update `show_smith_chart` to accept `&[SmithSeries]` and draw the R/X arc
  grid in light grey behind the data traces.
* Update `app.rs` Smith tab to build `SmithSeries` + "Show all entries"
  checkbox (reusing the `show_all_entries` flag already used by the dB tab).

### Track B — yee-plotters (`crates/yee-plotters/src/**`)

* Add private `smith_r_circle_points_plot` / `smith_x_arc_points_plot`.
* Add `SmithTrace { label, values }` (parallel to `SparamTrace`).
* Add `plot_smith_chart_multi(traces, out_path, config)`.
* Wrap existing `plot_smith_chart(s11, ...)` as a single-trace delegator to
  `plot_smith_chart_multi` — zero call-site changes needed.
* Update `draw_smith_multi` to draw R/X arcs + multi-trace with 8-colour
  palette.

### Reference overlay values

```
R circles: r ∈ {0.2, 0.5, 1, 2, 5}
X arcs:    x ∈ {±0.2, ±0.5, ±1, ±2, ±5}
```

All arc segments rendered in `RGBColor(210, 210, 210)` (light grey) so they
do not compete with the data traces.

## Consequences

* The Smith chart in the desktop GUI is now a functional impedance-matching
  tool, not just a locus viewer.
* `plot_smith_chart` remains API-stable; existing callers are unaffected.
* Both tracks carry full unit tests for the math helpers (centre/radius formula
  + unit-disk containment).
* **Deferred** (still out of scope): arc impedance labels, VSWR circles,
  admittance chart, `yee-cli --smith --all` multi-trace.

## References

* Spec: `docs/superpowers/specs/2026-05-25-smith-chart-arcs-design.md`
* Plan: `docs/superpowers/plans/2026-05-25-smith-chart-arcs.md`
* ADR-0063 `docs/src/decisions/0063-multi-trace-sparam-plotting.md`
* ADR-0065 `docs/src/decisions/0065-gui-multi-trace-overlay.md`
* Pozar, *Microwave Engineering*, 4th ed., §5.1.
