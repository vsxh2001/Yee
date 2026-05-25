# Phase 1.gui.5 / 1.plotting.2 — Smith Chart Constant-R/X Arc Overlays + GUI Multi-trace

**Date:** 2026-05-25  
**Phase:** 1.gui.5 (GUI arc overlay) / 1.plotting.2 (plotters arc overlay)  
**ADR:** [0068](../../../docs/src/decisions/0068-smith-chart-constant-rx-arcs.md)  
**Status:** Design accepted — pending implementation

---

## 1. Background and motivation

ADR-0063 (multi-trace S-parameter plotting, merge `185335c`) and ADR-0065 (GUI
multi-trace overlay, merge `ab1425d`) both explicitly deferred the Smith chart
enhancements:

> *"Smith-chart multi-trace + constant-R/X arcs remain follow-ons."*  
> — ADR-0063 / ADR-0065

The current Smith chart renders exactly:

| Location | What it renders |
|----------|----------------|
| `yee-gui::plots::show_smith_chart` | unit circle + S11 trajectory |
| `yee-plotters::draw_smith` | unit circle + crosshair + S11 trajectory |

Neither renders the constant-resistance circles or constant-reactance arcs that
are the primary visual reference on a real Smith chart (Pozar §5.1, Figure
5.1). Without the arc grid, the chart is usable for visualising the S11 locus
but provides no impedance-matching guidance.

This design closes the gap by adding:

1. **Constant-R circles and constant-X arcs** in both `yee-gui` and
   `yee-plotters` (pure visual reference, drawn in light grey).
2. **Multi-trace support in the GUI Smith chart** tab (analogous to the dB plot
   multi-trace from ADR-0065), so a multi-port file can show S11, S22, etc. on
   the same canvas with a legend.

---

## 2. Smith chart mathematics

The Γ-plane (reflection-coefficient plane) is the natural domain of the Smith
chart. The normalised load impedance `z_L = Z_L / Z_0 = r + jx` maps to:

```
Γ = (z_L − 1) / (z_L + 1)
```

### 2.1 Constant-R circles

A constant normalised resistance `r` traces a circle in the Γ-plane:

```
centre:  (r / (r+1),  0)
radius:  1 / (r+1)
```

At `r = 0` the circle degenerates to the unit circle (the outer boundary).
At `r → ∞` the circle shrinks to the point `(1, 0)` (short circuit).

Reference R values to render: `r ∈ {0.2, 0.5, 1, 2, 5}` (plus `r = 0` —
the unit circle already drawn as the outer boundary).

### 2.2 Constant-X arcs

A constant normalised reactance `x` traces a circle in the Γ-plane:

```
centre:  (1,  1/x)
radius:  1/|x|
```

Only the arc inside the unit disk (`|Γ| ≤ 1`) is rendered. The upper half-plane
(`Im Γ > 0`) corresponds to inductive loads (`x > 0`); the lower to capacitive
(`x < 0`).

Reference X values to render: `x ∈ {±0.2, ±0.5, ±1, ±2, ±5}`.

### 2.3 Clipping to the unit disk

Each candidate point `(re, im)` on a constant-X circle is admitted only if
`re² + im² ≤ 1 + ε` (with a small epsilon for floating-point tolerance).
Points outside the disk are discarded; the resulting polyline may have fewer
samples than requested, which is fine — `egui_plot` and `plotters` both handle
sparse polylines.

---

## 3. Scope and lane boundaries

| Track | Lane | Allowed paths |
|-------|------|--------------|
| A | `yee-gui` | `crates/yee-gui/src/**` |
| B | `yee-plotters` | `crates/yee-plotters/src/**` |

These two tracks are **disjoint**: no shared source file is touched by both.
The math helpers are duplicated (a small, pure, dependency-free set of
functions) rather than extracted to a shared crate — this keeps the
`yee-core` API surface clean and avoids a cross-lane dep bump.

---

## 4. Track A — yee-gui

### 4.1 New pure helpers in `plots.rs`

```rust
/// Constant-R circle on the Smith chart.
///
/// Returns `n + 1` points forming a closed circle of centre
/// `(r/(r+1), 0)` and radius `1/(r+1)`.  All points lie on or
/// inside the unit disk (by construction: the r=0 circle *is* the
/// unit circle; every r>0 circle is strictly inside it).
pub fn smith_r_circle_points(r: f64, n: usize) -> Vec<[f64; 2]>

/// Constant-X arc on the Smith chart, clipped to the unit disk.
///
/// The full constant-X circle has centre `(1, 1/x)` and radius
/// `1/|x|`.  Only points inside `|Γ| ≤ 1` (within a floating-point
/// tolerance of 1e-9) are returned; the arc may be shorter than `n`
/// points when `|x|` is small (large arcs cross the unit-disk
/// boundary quickly).  Panics if `x == 0`.
pub fn smith_x_arc_points(x: f64, n: usize) -> Vec<[f64; 2]>
```

### 4.2 Multi-trace series for Smith

```rust
/// A single labelled S-parameter trace in the complex Γ-plane, ready
/// for [`show_smith_chart`].
///
/// `label` is the legend string (e.g. `"S11"`, `"S22"`).
/// Each element of `points` is `[Re(Γ), Im(Γ)]` at one frequency.
#[derive(Debug, Clone)]
pub struct SmithSeries {
    pub label: String,
    pub points: Vec<[f64; 2]>,
}

/// Build labelled Smith-chart series from a loaded Touchstone file.
///
/// Reuses [`Selection`] from the dB multi-trace path.  Each entry
/// selected becomes one [`SmithSeries`]; the raw complex S-parameter
/// value is used directly (no dB conversion).
pub fn build_smith_series(file: &TsFile, selection: &Selection) -> Vec<SmithSeries>
```

### 4.3 Updated `show_smith_chart`

Signature change:

```rust
// Before
pub fn show_smith_chart(ui: &mut egui::Ui, s11: &[Complex64])

// After
pub fn show_smith_chart(ui: &mut egui::Ui, series: &[SmithSeries])
```

Rendering order (back-to-front):

1. Outer unit circle (medium grey, width 1).
2. Constant-R circles for `r ∈ {0.2, 0.5, 1, 2, 5}` (light grey, width 1).
3. Constant-X arcs for `x ∈ {±0.2, ±0.5, ±1, ±2, ±5}` (light grey, width 1).
4. Origin crosshair (light grey, width 1, length = diameter).
5. Each `SmithSeries` as a coloured polyline (8-colour cycling palette from
   ADR-0065, `show_sparams_db_plot`), with `egui_plot::Legend` enabled.

### 4.4 Updated `app.rs`

The Smith tab changes analogously to the dB tab:

```rust
TabKind::Smith => {
    if f.n_ports > 1 {
        ui.checkbox(self.show_all_entries, "Show all entries");
    }
    let selection = if *self.show_all_entries && f.n_ports > 1 {
        Selection::All
    } else {
        Selection::Diagonal(0)
    };
    let series = build_smith_series(f, &selection);
    show_smith_chart(ui, &series);
}
```

Note: `show_all_entries` is the **same** boolean flag already in `AppState`
for the dB tab (ADR-0065). Both tabs read it so toggling affects both views.

---

## 5. Track B — yee-plotters

### 5.1 New pure math helpers (private to crate)

Same formulas as Track A; private (`fn`, not `pub`):

```rust
fn smith_r_circle_points(r: f64, n: usize) -> Vec<(f64, f64)>
fn smith_x_arc_points(x: f64, n: usize) -> Vec<(f64, f64)>
```

(Using `(f64, f64)` tuples to match the `plotters::LineSeries` expected input
rather than `[f64; 2]` arrays used by `egui_plot::PlotPoints`.)

### 5.2 New `SmithTrace` type and multi-trace function

```rust
/// A single labelled S-parameter trace for the Smith chart.
///
/// Parallel to [`SparamTrace`] used by [`plot_sparams_db`].
pub struct SmithTrace {
    /// Legend label, e.g. `"S11"`.
    pub label: String,
    /// Raw complex S-parameter values (one per frequency sample).
    pub values: Vec<Complex64>,
}

/// Plot multiple S-parameter traces on a Smith-chart canvas with R/X
/// arc reference overlays.
///
/// Each [`SmithTrace`] is rendered as a connected polyline in a distinct
/// colour (cycling through a fixed 8-colour palette matching
/// `plot_sparams_db`).  A legend is drawn if `traces.len() > 1`.
pub fn plot_smith_chart_multi(
    traces: &[SmithTrace],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error>
```

### 5.3 Back-compat for `plot_smith_chart`

The existing single-trace `plot_smith_chart(s11, out_path, config)` becomes a
thin wrapper:

```rust
pub fn plot_smith_chart(
    s11: &[Complex64],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    let trace = SmithTrace { label: "S11".to_string(), values: s11.to_vec() };
    plot_smith_chart_multi(&[trace], out_path, config)
}
```

No call site changes required (single-trace callers continue to work).

### 5.4 Updated `draw_smith` internal helper

`draw_smith` gains a `traces: &[SmithTrace]` parameter and draws:

1. Origin crosshair (light grey).
2. Unit circle (medium grey).
3. Constant-R circles for `r ∈ {0.2, 0.5, 1, 2, 5}` (light grey).
4. Constant-X arcs for `x ∈ {±0.2, ±0.5, ±1, ±2, ±5}` (light grey).
5. Each trace as a coloured polyline (8-colour palette: red, blue, green,
   orange, purple, brown, pink, grey — same palette as `draw_multi_trace`).

---

## 6. Verification commands

```bash
# Track A — yee-gui unit tests
cargo test -p yee-gui --lib -- 2>&1 | tail -5
# expected: test result: ok. N passed; 0 failed

# Track A — lint
cargo clippy -p yee-gui --lib -- -D warnings 2>&1 | tail -5
# expected: 0 warnings, exit 0

# Track B — yee-plotters unit tests
cargo test -p yee-plotters --lib -- 2>&1 | tail -5
# expected: test result: ok. N passed; 0 failed

# Track B — lint
cargo clippy -p yee-plotters --lib -- -D warnings 2>&1 | tail -5
# expected: 0 warnings, exit 0

# Full workspace sanity (run after both tracks merged)
cargo test --workspace --lib 2>&1 | tail -10
cargo clippy --workspace --all-targets -- -D warnings 2>&1 | tail -5
cargo fmt --check --all 2>&1 | tail -5
```

---

## 7. Definition of Done

### Track A (yee-gui)

- [ ] `smith_r_circle_points` — all returned points lie on the expected circle
      (centre/radius formula); test for `r ∈ {0, 1}`.
- [ ] `smith_r_circle_points` — all points satisfy `|Γ| ≤ 1 + 1e-9` (unit-disk
      containment); test for `r = 0.5` (strictest interior case).
- [ ] `smith_x_arc_points` — returned points lie on the expected circle
      (centre/radius formula); test for `x = 1`.
- [ ] `smith_x_arc_points` — all returned points satisfy `|Γ| ≤ 1 + 1e-9`; test
      for `x = 0.2` (largest arc, most points outside).
- [ ] `build_smith_series` — `Selection::All` on a 2-port file returns 4 series
      with correct labels and point counts.
- [ ] `show_smith_chart` compiles and passes the smoke test (call with empty
      series slice and with one-element series).
- [ ] `app.rs` uses `build_smith_series` for the Smith tab; "Show all entries"
      checkbox visible for multi-port files.
- [ ] `cargo test -p yee-gui --lib` → 0 failed.
- [ ] `cargo clippy -p yee-gui -- -D warnings` → 0 warnings.

### Track B (yee-plotters)

- [ ] `smith_r_circle_points` unit-tested (centre/radius + containment).
- [ ] `smith_x_arc_points` unit-tested (on-circle + containment).
- [ ] `plot_smith_chart_multi` with a 2-trace input writes a non-empty SVG file.
- [ ] `plot_smith_chart` (single-trace back-compat) still passes the existing
      `test_plot_smith_chart_writes_svg` test.
- [ ] `cargo test -p yee-plotters --lib` → 0 failed.
- [ ] `cargo clippy -p yee-plotters -- -D warnings` → 0 warnings.

---

## 8. Out of scope for this phase

- Annotating R/X arcs with impedance labels (text along the arc — needs
  `egui_plot` custom painter or `plotters` rotated text; deferred).
- VSWR circles (constant `|Γ|` circles — trivially added later as a separate
  overlay pass).
- `yee-cli` `--smith` flag multi-trace support (the CLI delegates to
  `yee-plotters::plot_smith_chart`; the back-compat wrapper keeps it working;
  a `--all` flag for Smith is a follow-on).
- Admittance chart (rotated Smith chart) — deferred.

---

## 9. References

- Pozar, *Microwave Engineering*, 4th ed., §5.1 (Smith Chart figure 5.1).
- ADR-0063 `docs/src/decisions/0063-multi-trace-sparam-plotting.md`
- ADR-0065 `docs/src/decisions/0065-gui-multi-trace-overlay.md`
- ADR-0068 `docs/src/decisions/0068-smith-chart-constant-rx-arcs.md`
