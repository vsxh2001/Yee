# Phase 1.gui.6 — VSWR Circles on the Smith Chart: Design Spec

**Date:** 2026-05-29  
**Phase:** 1.gui.6  
**ADR:** 0081  
**Status:** Proposed  

---

## 1. Problem statement

The Smith chart in both `yee-gui` (egui `show_smith_chart`) and `yee-plotters`
(`plot_smith_chart_multi` / `draw_smith_multi`) currently renders:

1. Unit circle `|Γ| = 1`  
2. Constant-R circles for r ∈ {0.2, 0.5, 1.0, 2.0, 5.0}  
3. Constant-X arcs for x ∈ {±0.2, ±0.5, ±1.0, ±2.0, ±5.0}  
4. Data traces  

ADR-0068 (Phase 1.gui.5) deferred VSWR circles, arc labels, and the
admittance chart as follow-ons. This spec closes the VSWR-circles deferral.

VSWR (Voltage Standing Wave Ratio) circles are among the most commonly-used
Smith chart overlays in RF/microwave engineering.  Every commercial tool
(Keysight ADS, AWR Microwave Office, QUCS, scikit-rf) draws them by default.
Their absence makes the Yee Smith chart look incomplete relative to peers.

---

## 2. Background

### 2.1 VSWR on the Smith chart

For a load with reflection coefficient Γ:

```
VSWR = (1 + |Γ|) / (1 - |Γ|)
```

Equivalently, a constant-VSWR locus on the Γ-plane (Smith chart) is a
**circle centred at the origin** with radius:

```
ρ = |Γ| = (VSWR − 1) / (VSWR + 1)
```

Standard VSWR circles drawn on industry Smith charts:
VSWR ∈ {1.5, 2, 3, 5, 10}, giving ρ ∈ {1/4, 1/3, 1/2, 2/3, 9/11}.

Reference: Pozar, *Microwave Engineering*, 4th ed., §2.5 (Smith chart
fundamentals); Keysight Technologies, *Electronic Calibration Module (ECal)
Help*, Smith chart description (freely available).

### 2.2 Scope of this phase

- **In scope:** VSWR circle overlay in `yee-gui` + `yee-plotters`; unit
  tests for the circle-generation helper; ROADMAP and ADR updates.
- **Out of scope:** arc labels (R/X/VSWR text), admittance chart overlay,
  CLI `--vswr` flag — these remain in the ADR-0068 deferred list.

---

## 3. Design

### 3.1 Pure geometry helper (both crates)

```
vswr_circle_points(vswr: f64, n: usize) -> Vec<[f64; 2]>
```

- `vswr` must be > 1.0 (debug_assert).
- Returns `n + 1` points forming a closed circle centred at the Γ-plane
  origin with radius `ρ = (vswr − 1) / (vswr + 1)`.
- Convention: `points[0] == points[n]` (last point closes the loop).
- Points are in Γ-space: x = Re Γ, y = Im Γ.

```rust
pub fn smith_vswr_circle_points(vswr: f64, n: usize) -> Vec<[f64; 2]> {
    debug_assert!(vswr > 1.0, "vswr must be > 1.0 (got {vswr})");
    let rho = (vswr - 1.0) / (vswr + 1.0);
    (0..=n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            [rho * theta.cos(), rho * theta.sin()]
        })
        .collect()
}
```

In `yee-plotters`, a private mirror with `(f64, f64)` tuples:

```rust
fn smith_vswr_circle_pts(vswr: f64, n: usize) -> Vec<(f64, f64)>
```

### 3.2 GUI: `yee-gui/src/plots.rs`

In `show_smith_chart`, after step 3 (constant-X arcs), add step 3.5:

```rust
// 3.5 VSWR circles.
let vswr_style = egui::Color32::from_rgb(180, 180, 220);
for &vswr in &[1.5_f64, 2.0, 3.0, 5.0, 10.0] {
    let pts: PlotPoints =
        smith_vswr_circle_points(vswr, 128).into_iter().collect();
    plot_ui.line(
        Line::new(format!("VSWR={vswr}"), pts).color(vswr_style),
    );
}
```

(Note: `egui_plot::Line::color` exists in the 0.35 API. The label string
`"VSWR=1.5"` etc. appears in the Legend.)

### 3.3 Plotters: `yee-plotters/src/lib.rs`

In `draw_smith_multi`, after constant-X arcs, add:

```rust
// VSWR circles (very light blue-grey).
let vswr_style = RGBColor(190, 190, 220).stroke_width(1);
for &vswr in &[1.5_f64, 2.0, 3.0, 5.0, 10.0] {
    let pts = smith_vswr_circle_pts(vswr, 128);
    chart
        .draw_series(LineSeries::new(pts, vswr_style))
        .map_err(map_render_err)?;
}
```

---

## 4. Validation

Unit tests in `yee-gui/src/plots.rs` `#[cfg(test)]` block:

1. `smith_vswr_circle_points_has_n_plus_1_points` — length check.
2. `smith_vswr_circle_points_all_on_circle` — every point at radius ρ ± 1e-12.
3. `smith_vswr_circle_points_closes` — `first == last`.
4. `smith_vswr_circle_points_rho_for_vswr_2` — ρ = 1/3 for VSWR=2.

Unit tests in `yee-plotters/src/lib.rs` `#[cfg(test)]` block (within the
existing `tests` module):

1. `smith_vswr_circle_pts_radius_vswr_2` — ρ = 1/3 within 1e-12.
2. `smith_vswr_circle_pts_closed` — first ≈ last within 1e-12.

Smoke test (egui-free, plotters only):

3. `plot_smith_chart_multi_with_vswr_circles_renders` — call
   `plot_smith_chart_multi` with a synthetic `SmithTrace` and verify the
   output file is non-empty (mirrors `test_plot_smith_chart_writes_svg`).

**Verification command:**

```bash
cargo test -p yee-gui -p yee-plotters -- --nocapture 2>&1 | tail -5
# expected: test result: ok. N passed; 0 failed; 0 ignored
```

---

## 5. Risks

| Risk | Mitigation |
|------|-----------|
| `egui_plot::Line::color` API may differ from 0.35 | Check the egui_plot 0.35 docs; the method exists in 0.34+ |
| VSWR circles obscure data traces | Circles are drawn BEFORE data traces; use light muted colour |
| Out-of-lane edit (e.g. touching yee-cli) | Scope is strictly `crates/yee-gui/src/** + crates/yee-plotters/src/**` |
| VSWR < 1.0 or = 1.0 passed to helper | `debug_assert!(vswr > 1.0)` guards this |

---

## 6. Follow-ons (not this phase)

- Arc labels (R, X, VSWR values rendered as text on the chart)
- Admittance chart overlay (G and B circles)
- CLI `yee plot --vswr` flag to draw VSWR circles from the command line
- Marker annotation for best-match frequency on the VSWR circle nearest Γ=0

---

*Spec author: orchestrator run 2026-05-29*
