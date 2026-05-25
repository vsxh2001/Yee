# Implementation Plan — Phase 1.gui.5 / 1.plotting.2 Smith Chart Arcs

**Spec:** `docs/superpowers/specs/2026-05-25-smith-chart-arcs-design.md`  
**ADR:** `docs/src/decisions/0068-smith-chart-constant-rx-arcs.md`  
**Date:** 2026-05-25

---

## Tracks

Two disjoint implementation tracks run in parallel:

| Track | Worktree | Branch | Lane |
|-------|----------|--------|------|
| A | `worktrees/smith-gui/` | `feature/1.gui.5-smith-arcs-gui` | `crates/yee-gui/src/**` |
| B | `worktrees/smith-plotters/` | `feature/1.plotting.2-smith-arcs-plotters` | `crates/yee-plotters/src/**` |

---

## Track A — yee-gui (`crates/yee-gui/src/**`)

### A1 — Math helpers in `plots.rs`

Add two pure helper functions **before** the existing `unit_circle_points`:

```rust
pub fn smith_r_circle_points(r: f64, n: usize) -> Vec<[f64; 2]>
pub fn smith_x_arc_points(x: f64, n: usize) -> Vec<[f64; 2]>
```

- `smith_r_circle_points`: parametrise `θ ∈ [0, 2π]` in `n+1` steps (closing
  the loop); compute `Γ = centre + radius·(cos θ, sin θ)`.
- `smith_x_arc_points`: parametrise `θ ∈ [0, 2π]` in `n` samples; compute the
  full circle, then keep only points with `re²+im² ≤ 1.0 + 1e-9`.

### A2 — `SmithSeries` struct and `build_smith_series`

Add after `build_sparam_series`:

```rust
pub struct SmithSeries { pub label: String, pub points: Vec<[f64; 2]> }

pub fn build_smith_series(file: &TsFile, selection: &Selection) -> Vec<SmithSeries>
```

Implementation pattern: mirror `build_sparam_series`, but instead of
`[freq_ghz, db]` pairs, use `[z.re, z.im]` directly.

### A3 — `show_smith_chart` rewrite

Change signature from `(ui, s11: &[Complex64])` to `(ui, series: &[SmithSeries])`.

Drawing order:

1. Outer unit circle — `unit_circle_points(256)` — `RGBColor(160,160,160)`.
2. R circles — for each r in `[0.2, 0.5, 1.0, 2.0, 5.0]`:
   `smith_r_circle_points(r, 128)` — `RGBColor(210,210,210)`.
3. X arcs — for each x in `[0.2, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0]`:
   `smith_x_arc_points(x, 256)` — `RGBColor(210,210,210)`.
4. Each `SmithSeries` as a `Line` with a colour from the 8-colour palette.
5. `Legend::default()` if `series.len() > 1`.

Palette (matching `show_sparams_db_plot`):

```rust
const PALETTE: &[egui::Color32] = &[
    egui::Color32::from_rgb(228, 26, 28),   // red
    egui::Color32::from_rgb(55, 126, 184),  // blue
    egui::Color32::from_rgb(77, 175, 74),   // green
    egui::Color32::from_rgb(255, 127, 0),   // orange
    egui::Color32::from_rgb(152, 78, 163),  // purple
    egui::Color32::from_rgb(166, 86, 40),   // brown
    egui::Color32::from_rgb(247, 129, 191), // pink
    egui::Color32::from_rgb(153, 153, 153), // grey
];
```

Keep the `Plot::data_aspect(1.0)` lock.

### A4 — Update `app.rs`

Change the `TabKind::Smith` arm to:

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

### A5 — Unit tests

Add to the `#[cfg(test)] mod tests` block:

1. `smith_r_circle_returns_correct_n_plus_1_points` — for `r = 1`, `n = 64`:
   - length = 65
   - all points on the circle (centre `(0.5, 0)`, radius `0.5`) within 1e-9
   - first == last (closed loop)
2. `smith_r_circle_contained_in_unit_disk` — for `r = 0.5`, all `|Γ| ≤ 1 + 1e-9`
3. `smith_x_arc_points_on_circle` — for `x = 1.0`, `n = 256`: every returned
   point satisfies `|Γ - (1, 1)| ≈ 1` within 1e-9
4. `smith_x_arc_points_inside_unit_disk` — for `x = 0.2` (widest arc), all
   points satisfy `|Γ| ≤ 1 + 1e-9`
5. `build_smith_series_all_two_port` — Selection::All on a 2-port file:
   4 series, labels `["S11","S12","S21","S22"]`, each with 5 points
   (`[re, im]` match the raw `Complex64` values)

### A6 — Verification

```bash
cd worktrees/smith-gui
cargo test -p yee-gui --lib -- 2>&1 | tail -5
cargo clippy -p yee-gui --lib -- -D warnings 2>&1 | tail -5
cargo fmt --check -p yee-gui 2>&1 | tail -3
```

All must exit 0.

---

## Track B — yee-plotters (`crates/yee-plotters/src/**`)

### B1 — Private math helpers

Add private helper functions in `lib.rs` (near `draw_smith`):

```rust
fn smith_r_circle_points_plot(r: f64, n: usize) -> Vec<(f64, f64)>
fn smith_x_arc_points_plot(x: f64, n: usize) -> Vec<(f64, f64)>
```

Same formulas as Track A; return `(f64, f64)` tuples for `plotters::LineSeries`.

### B2 — `SmithTrace` type

Add after `SparamTrace`:

```rust
/// A labelled S-parameter trace for the Smith chart.
///
/// Parallel to [`SparamTrace`] for [`plot_sparams_db`].
pub struct SmithTrace {
    /// Legend label, e.g. `"S11"`.
    pub label: String,
    /// Raw complex values (one per frequency point).
    pub values: Vec<Complex64>,
}
```

### B3 — `plot_smith_chart_multi`

```rust
pub fn plot_smith_chart_multi(
    traces: &[SmithTrace],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error>
```

Delegates to updated `draw_smith_multi`.

### B4 — Back-compat `plot_smith_chart` wrapper

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

### B5 — Update `draw_smith` → `draw_smith_multi`

Rename (or add alongside) the internal draw function. New signature:

```rust
fn draw_smith_multi<DB>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    traces: &[SmithTrace],
) -> Result<(), Error>
```

Drawing order:

1. Origin crosshair — `RGBColor(200, 200, 200)`.
2. Unit circle — `RGBColor(160, 160, 160)`.
3. R circles for `r ∈ [0.2, 0.5, 1.0, 2.0, 5.0]` — `RGBColor(210, 210, 210)`.
4. X arcs for `x ∈ ±[0.2, 0.5, 1.0, 2.0, 5.0]` — `RGBColor(210, 210, 210)`.
5. Each trace: iterate `traces.iter().enumerate()`, use an 8-colour palette.

Palette (tuples of `(r, g, b)`):

```
(228, 26, 28), (55, 126, 184), (77, 175, 74), (255, 127, 0),
(152, 78, 163), (166, 86, 40), (247, 129, 191), (153, 153, 153)
```

### B6 — Unit tests

Add to the test block:

1. `smith_r_circle_plot_correct_geometry` — for `r = 1.0`, `n = 64`:
   all 65 points on circle with centre `(0.5, 0)`, radius `0.5`, within 1e-9.
2. `smith_x_arc_plot_inside_unit_disk` — for `x = 0.2`, all `re²+im² ≤ 1 + 1e-9`.
3. `plot_smith_chart_multi_two_traces_writes_svg` — 2 traces, 4 points each
   → a non-empty SVG file is written.
4. `plot_smith_chart_back_compat_still_works` — existing
   `test_plot_smith_chart_writes_svg` unchanged (zero modification to that test).

### B7 — Verification

```bash
cd worktrees/smith-plotters
cargo test -p yee-plotters --lib -- 2>&1 | tail -5
cargo clippy -p yee-plotters --lib -- -D warnings 2>&1 | tail -5
cargo fmt --check -p yee-plotters 2>&1 | tail -3
```

All must exit 0.

---

## Merge order

1. Review Track A and Track B independently (separate reviewer agents if
   possible).
2. Merge Track A → main (`git merge --no-ff feature/1.gui.5-smith-arcs-gui`).
3. Merge Track B → main (`git merge --no-ff feature/1.plotting.2-smith-arcs-plotters`).
4. Run full workspace check:
   ```
   cargo test --workspace --lib
   cargo clippy --workspace --all-targets -- -D warnings
   cargo fmt --check --all
   ```
5. Update ROADMAP.md and ADR-0068.
6. Push main.

---

## Escape hatch

**Blocked > 15 min → surface the blocker and stop.**

Known potential blockers:
- egui_plot `Line::new` API has a different signature in egui 0.34 from older
  versions. Pattern: match `show_sparams_db_plot` in the same file for the
  correct call shape.
- The `plotters` `LineSeries::new` API requires `(f64, f64)` tuples; Track B
  should match the existing `draw_smith` style exactly.
