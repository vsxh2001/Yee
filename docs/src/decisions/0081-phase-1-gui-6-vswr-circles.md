# ADR-0081: Phase 1.gui.6 VSWR circles on the Smith chart

**Date:** 2026-05-29  
**Status:** Accepted  
**Deciders:** orchestrator  
**References:** ADR-0068 (Phase 1.gui.5 Smith chart arcs + multi-trace),
ADR-0069 (Phase 1.plotting.3 CLI Smith multi-trace)

---

## Context

ADR-0068 delivered constant-R circles and constant-X arcs on the Smith chart
in both `yee-gui` and `yee-plotters`, but explicitly deferred VSWR circles,
arc labels, and the admittance chart as follow-ons.

VSWR (Voltage Standing Wave Ratio) circles are among the most widely used
Smith chart overlays: every commercial RF tool (Keysight ADS, AWR MWO, QUCS,
scikit-rf) draws them by default. On the normalised Γ-plane a constant-VSWR
locus is a circle centred at the origin with radius
`ρ = (VSWR − 1) / (VSWR + 1)`. Standard values are VSWR ∈ {1.5, 2, 3, 5, 10}.

Reference: Pozar, *Microwave Engineering*, 4th ed., §2.5.

---

## Decision

Add VSWR circles to the Smith chart in **both** `yee-gui` (egui) and
`yee-plotters` (plotters).

**Implementation:**

- `pub fn smith_vswr_circle_points(vswr: f64, n: usize) -> Vec<[f64; 2]>`
  in `yee-gui/src/plots.rs` (unit-testable pure function).
- Private mirror `fn smith_vswr_circle_pts(vswr: f64, n: usize) -> Vec<(f64,f64)>`
  in `yee-plotters/src/lib.rs`.
- Both draw VSWR ∈ {1.5, 2, 3, 5, 10} in a muted light-blue-grey colour,
  before data traces so they don't obscure measurements.
- 6 new unit tests (4 in yee-gui, 2 in yee-plotters).

**Not in scope for this ADR:**
- Arc labels (text annotations on R/X/VSWR circles)
- Admittance chart overlay
- CLI `--vswr` flag (yee plot)

---

## Consequences

- The Yee Smith chart now matches the VSWR-circle convention used by every
  major RF EDA tool.
- No new dependencies; zero physics code touched; no validation gates affected.
- Arc labels, admittance overlay, and CLI multi-trace remain as follow-on items.
