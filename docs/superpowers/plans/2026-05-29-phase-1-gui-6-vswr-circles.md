# Phase 1.gui.6 — VSWR Circles on the Smith Chart: Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-29-phase-1-gui-6-vswr-circles-design.md`  
**ADR:** `docs/src/decisions/0081-phase-1-gui-6-vswr-circles.md`  
**Base SHA:** `36dbf69681837f2522519501b9daab829bf46608`  
**Branch:** `feature/phase-1-gui-6-vswr-circles`  
**Lane:** `crates/yee-gui/src/**, crates/yee-plotters/src/**`  

---

## Step-by-step

### Step 1 — `yee-plotters`: add private helper + draw in `draw_smith_multi`

File: `crates/yee-plotters/src/lib.rs`

1a. After `fn smith_x_arc_pts`, add:

```rust
/// Compute sample points on the constant-VSWR circle centred at the
/// Γ-plane origin with radius `ρ = (vswr − 1) / (vswr + 1)`.
///
/// Returns `n + 1` points; `pts[0] == pts[n]` closes the loop.
/// `vswr` must be > 1.0.
fn smith_vswr_circle_pts(vswr: f64, n: usize) -> Vec<(f64, f64)> {
    debug_assert!(vswr > 1.0, "smith_vswr_circle_pts: vswr must be > 1.0");
    let rho = (vswr - 1.0) / (vswr + 1.0);
    (0..=n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            (rho * theta.cos(), rho * theta.sin())
        })
        .collect()
}
```

1b. In `draw_smith_multi`, after the constant-X arc loop (step 6) and
    BEFORE the data-trace legend (step 7), add:

```rust
// VSWR circles (light blue-grey, drawn before data traces).
let vswr_style = RGBColor(190, 190, 220).stroke_width(1);
for &vswr in &[1.5_f64, 2.0, 3.0, 5.0, 10.0] {
    let pts = smith_vswr_circle_pts(vswr, 128);
    chart
        .draw_series(LineSeries::new(pts, vswr_style))
        .map_err(map_render_err)?;
}
```

1c. Add unit tests at the bottom of the `#[cfg(test)]` block:

```rust
#[test]
fn smith_vswr_circle_pts_radius_vswr_2() {
    let pts = smith_vswr_circle_pts(2.0, 128);
    let rho_expected = 1.0_f64 / 3.0;
    for (x, y) in &pts {
        let r = (x * x + y * y).sqrt();
        assert!((r - rho_expected).abs() < 1e-12,
                "VSWR=2 point off-circle: r={r}, expected {rho_expected}");
    }
}

#[test]
fn smith_vswr_circle_pts_closed() {
    let pts = smith_vswr_circle_pts(3.0, 64);
    assert_eq!(pts.len(), 65);
    let (x0, y0) = pts[0];
    let (xn, yn) = pts[64];
    assert!((x0 - xn).abs() < 1e-12 && (y0 - yn).abs() < 1e-12,
            "VSWR circle not closed");
}
```

**Verify:** `cargo test -p yee-plotters -- --nocapture` → all pass, 0 failed.

---

### Step 2 — `yee-gui`: add public helper + draw in `show_smith_chart`

File: `crates/yee-gui/src/plots.rs`

2a. After `pub fn smith_x_arc_points`, add:

```rust
/// Compute sample points on the constant-VSWR circle centred at the
/// Γ-plane origin with radius `ρ = (vswr − 1) / (vswr + 1)`.
///
/// Returns `n + 1` points; `pts[0] == pts[n]` closes the loop.
/// Used by [`show_smith_chart`] and unit-testable without a live egui context.
///
/// # Panics
///
/// Panics (debug build) if `vswr ≤ 1.0`.
pub fn smith_vswr_circle_points(vswr: f64, n: usize) -> Vec<[f64; 2]> {
    debug_assert!(vswr > 1.0, "smith_vswr_circle_points: vswr must be > 1.0");
    let rho = (vswr - 1.0) / (vswr + 1.0);
    (0..=n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            [rho * theta.cos(), rho * theta.sin()]
        })
        .collect()
}
```

2b. In `show_smith_chart`, after step 3 (constant-X arcs) and BEFORE step 4
    (data traces), add:

```rust
// 3.5. VSWR circles (light blue-grey, before data traces).
let vswr_colour = egui::Color32::from_rgb(180, 180, 220);
for &vswr in &[1.5_f64, 2.0, 3.0, 5.0, 10.0] {
    let pts: PlotPoints =
        smith_vswr_circle_points(vswr, 128).into_iter().collect();
    plot_ui.line(
        Line::new(format!("VSWR={vswr}"), pts)
            .color(vswr_colour),
    );
}
```

2c. Add unit tests in the `#[cfg(test)]` block, after the existing
    `smith_x_arc_points` tests:

```rust
// -------------------------------------------------------------------------
// smith_vswr_circle_points
// -------------------------------------------------------------------------

/// `smith_vswr_circle_points(vswr, n)` returns exactly `n + 1` points.
#[test]
fn smith_vswr_circle_points_has_n_plus_1_points() {
    let pts = smith_vswr_circle_points(2.0, 64);
    assert_eq!(pts.len(), 65);
}

/// Every point returned by `smith_vswr_circle_points` must lie on the
/// circle of radius `ρ = (VSWR−1)/(VSWR+1)`.
#[test]
fn smith_vswr_circle_points_all_on_circle() {
    let vswr = 3.0_f64;
    let rho = (vswr - 1.0) / (vswr + 1.0);
    let pts = smith_vswr_circle_points(vswr, 128);
    for [x, y] in &pts {
        let r = (x * x + y * y).sqrt();
        assert!((r - rho).abs() < 1e-12,
                "VSWR={vswr} point off-circle: r={r}, expected {rho}");
    }
}

/// The first and last points must coincide (closed loop).
#[test]
fn smith_vswr_circle_points_closes() {
    let pts = smith_vswr_circle_points(2.0, 64);
    assert_eq!(pts.first(), pts.last(), "VSWR circle not closed");
}

/// For VSWR = 2, ρ = 1/3 exactly.
#[test]
fn smith_vswr_circle_points_rho_for_vswr_2() {
    let pts = smith_vswr_circle_points(2.0, 128);
    let rho_expected = 1.0_f64 / 3.0;
    // The first point lies at angle 0: [ρ, 0].
    let [x, y] = pts[0];
    assert!((y).abs() < 1e-12, "y[0] should be 0, got {y}");
    assert!((x - rho_expected).abs() < 1e-12,
            "x[0] = {x}, expected {rho_expected}");
}
```

**Verify:** `cargo test -p yee-gui -- --nocapture` → all pass, 0 failed.

---

### Step 3 — Lint check

```bash
cargo clippy -p yee-gui -p yee-plotters --all-targets -- -D warnings
cargo fmt --check -p yee-gui -p yee-plotters
```

Both must exit 0.

---

### Step 4 — ADR and docs update

4a. Write `docs/src/decisions/0081-phase-1-gui-6-vswr-circles.md`.
4b. Add ADR entry to `docs/src/SUMMARY.md` in the Decisions section.
4c. Update `ROADMAP.md`: add Phase 1.gui.6 bullet to the "Shipped" section.
4d. Update `CLAUDE.md` `*Last updated*` line.

These are in the **docs lane** — write them in the same commit or a follow-on
commit on the same branch. They should NOT be in a separate worktree since
the lane is strictly `crates/yee-gui/src/**, crates/yee-plotters/src/**`
plus docs.

---

## DoD checklist

- [ ] `smith_vswr_circle_pts` / `smith_vswr_circle_points` functions present
      in both crates
- [ ] VSWR circles (1.5, 2, 3, 5, 10) drawn in `show_smith_chart` (egui)
- [ ] VSWR circles (1.5, 2, 3, 5, 10) drawn in `draw_smith_multi` (plotters)
- [ ] 4 unit tests in `yee-gui` pass
- [ ] 2 unit tests in `yee-plotters` pass
- [ ] `cargo clippy -p yee-gui -p yee-plotters -- -D warnings` exits 0
- [ ] `cargo fmt --check -p yee-gui -p yee-plotters` exits 0
- [ ] ADR-0081 written and SUMMARY.md updated
- [ ] ROADMAP.md and CLAUDE.md updated

## Verification command

```bash
cargo test -p yee-gui -p yee-plotters 2>&1 | tail -5
# Expected: test result: ok. N passed; 0 failed; 0 ignored
```

Exit code 0.

## Escape hatch

Blocked > 15 min on any single step → surface the finding and stop.
Do NOT edit out-of-lane (i.e., do NOT touch `yee-cli`, `yee-fdtd`,
`yee-mom`, or any other crate except `yee-gui` and `yee-plotters`).
