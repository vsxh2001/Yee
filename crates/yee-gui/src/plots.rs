//! Plotting helpers for the Phase 1.gui.0 walking-skeleton studio shell.
//!
//! All functions here are pure data → screen-space transforms; nothing in this
//! module touches global state or owns a GUI runtime, which makes the math
//! easy to unit-test.
//!
//! ## Smith chart drawing strategy
//!
//! `egui_plot` does not ship a Smith chart widget, so we synthesise one by:
//! 1. Drawing the unit circle as a polyline with [`unit_circle_points`].
//! 2. Plotting `(Re S11, Im S11)` as a second polyline on the same axes.
//! 3. Locking the plot's data aspect ratio to 1:1 so the unit circle is round.
//!
//! ## dB plot strategy
//!
//! Standard `egui_plot::Plot` with frequency on X (in GHz for human-friendly
//! scaling) and `20·log10(|S11|)` on Y. We use [`db_from_complex`] to compute
//! the magnitude in dB so the math is testable independent of egui.
//!
//! ## Multi-trace overlay
//!
//! [`build_sparam_series`] is a **pure** (no egui) helper that converts a
//! loaded [`yee_io::touchstone::File`] plus a [`Selection`] into a
//! `Vec<`[`SparamSeries`]`>`.  The egui draw calls live in
//! [`show_sparams_db_plot`] and are only smoke-tested (they require a live
//! egui context). The series-building logic is unit-tested in the `tests`
//! block below.

use egui_plot::{Legend, Line, Plot, PlotPoints};
use num_complex::Complex64;
use yee_io::touchstone::File as TsFile;

/// Convert a complex S-parameter to its magnitude in decibels:
/// `20 · log10 |z|`.
///
/// Returns `f64::NEG_INFINITY` for `z = 0` rather than panicking, so callers
/// can hand the result straight to `egui_plot` (which clamps to the visible
/// Y range).
pub fn db_from_complex(z: Complex64) -> f64 {
    let mag = z.norm();
    if mag > 0.0 {
        20.0 * mag.log10()
    } else {
        f64::NEG_INFINITY
    }
}

// ---------------------------------------------------------------------------
// Multi-trace series building (pure — no egui types)
// ---------------------------------------------------------------------------

/// A single labelled S-parameter trace ready for [`show_sparams_db_plot`].
///
/// `label` appears in the egui_plot legend (e.g. `"S11"`, `"S21"`).
/// `points` is a parallel array of `[freq_ghz, db]` pairs.
#[derive(Debug, Clone)]
pub struct SparamSeries {
    /// Legend label, e.g. `"S11"` or `"S21"` (1-based, matching Touchstone).
    pub label: String,
    /// Plot points: each element is `[frequency_ghz, magnitude_db]`.
    pub points: Vec<[f64; 2]>,
}

/// Which S-matrix entries to include in a multi-trace overlay.
///
/// The default for a freshly-opened file is [`Selection::Diagonal`]`(0)`,
/// which reproduces the pre-existing single-trace (`S11`) behaviour.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Selection {
    /// A single diagonal entry `S[i][i]` (0-based).
    /// `Diagonal(0)` → `S11`, `Diagonal(1)` → `S22`, etc.
    Diagonal(usize),
    /// An explicit set of `(row, col)` entries (0-based).
    Entries(Vec<(usize, usize)>),
    /// Every entry of the `n × n` S-matrix in row-major order.
    All,
}

impl Default for Selection {
    fn default() -> Self {
        Selection::Diagonal(0)
    }
}

/// Build labelled (freq-GHz, dB) series from a loaded Touchstone file.
///
/// This is the **pure** (no egui) series-building helper.  It mirrors the
/// extraction idiom from `yee-cli/src/plot.rs`:
/// `flat_idx = row * n_ports + col`, with 1-based `S<ij>` labels.
///
/// - Out-of-range diagonal indices or explicit `(row, col)` pairs return an
///   empty `Vec` for that entry (they are silently skipped rather than
///   panicking, so the egui layer does not need error handling at the
///   drawing site).
/// - `data[k]` is row-major: `data[k][r * n + c] = S_{r,c}` at frequency
///   index `k`.
///
/// Returns one [`SparamSeries`] per selected entry, in the order they appear
/// in the selection (row-major for [`Selection::All`]).
pub fn build_sparam_series(file: &TsFile, selection: &Selection) -> Vec<SparamSeries> {
    let n = file.n_ports;

    // Expand the selection into (row, col) index pairs (0-based).
    let pairs: Vec<(usize, usize)> = match selection {
        Selection::Diagonal(i) => {
            if *i < n {
                vec![(*i, *i)]
            } else {
                vec![]
            }
        }
        Selection::Entries(pairs) => pairs
            .iter()
            .filter(|&&(r, c)| r < n && c < n)
            .copied()
            .collect(),
        Selection::All => (0..n).flat_map(|r| (0..n).map(move |c| (r, c))).collect(),
    };

    pairs
        .into_iter()
        .map(|(r, c)| {
            let flat_idx = r * n + c;
            let points: Vec<[f64; 2]> = file
                .freq_hz
                .iter()
                .zip(file.data.iter())
                .map(|(&f_hz, s_matrix)| {
                    let db = db_from_complex(s_matrix[flat_idx]);
                    [f_hz * 1.0e-9, db]
                })
                .collect();
            // 1-based label to match Touchstone / CLI convention.
            let label = format!("S{}{}", r + 1, c + 1);
            SparamSeries { label, points }
        })
        .collect()
}

/// A single labelled S-parameter trace in the complex Γ-plane.
///
/// Each element of `points` is `[Re(S), Im(S)]` at one frequency sample.
#[derive(Debug, Clone)]
pub struct SmithSeries {
    /// Legend label, e.g. `"S11"` or `"S22"`.
    pub label: String,
    /// Plot points: each is `[Re(S_ij), Im(S_ij)]`.
    pub points: Vec<[f64; 2]>,
}

/// Build labelled Smith-chart series from a loaded Touchstone file.
///
/// Mirrors [`build_sparam_series`], but stores raw `[Re, Im]` rather
/// than `[freq_ghz, dB]` pairs.
pub fn build_smith_series(file: &TsFile, selection: &Selection) -> Vec<SmithSeries> {
    let n = file.n_ports;

    // Expand the selection into (row, col) index pairs (0-based).
    let pairs: Vec<(usize, usize)> = match selection {
        Selection::Diagonal(i) => {
            if *i < n {
                vec![(*i, *i)]
            } else {
                vec![]
            }
        }
        Selection::Entries(pairs) => pairs
            .iter()
            .filter(|&&(r, c)| r < n && c < n)
            .copied()
            .collect(),
        Selection::All => (0..n).flat_map(|r| (0..n).map(move |c| (r, c))).collect(),
    };

    pairs
        .into_iter()
        .map(|(r, c)| {
            let flat_idx = r * n + c;
            let points: Vec<[f64; 2]> = file
                .data
                .iter()
                .map(|s_matrix| {
                    let z = s_matrix[flat_idx];
                    [z.re, z.im]
                })
                .collect();
            // 1-based label to match Touchstone / CLI convention.
            let label = format!("S{}{}", r + 1, c + 1);
            SmithSeries { label, points }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// egui drawing helpers
// ---------------------------------------------------------------------------

/// Sample `n` points evenly around the complex unit circle, returned as
/// `(x, y)` pairs suitable for [`egui_plot::PlotPoints`].
///
/// The circle is closed (the last point equals the first) so the polyline
/// renders without a visible gap.
pub fn unit_circle_points(n: usize) -> Vec<[f64; 2]> {
    assert!(n >= 3, "unit_circle_points needs at least 3 samples");
    let mut pts = Vec::with_capacity(n + 1);
    for i in 0..n {
        let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
        pts.push([theta.cos(), theta.sin()]);
    }
    // Close the loop.
    pts.push(pts[0]);
    pts
}

/// Constant-R circle on the Smith chart in the Γ-plane.
///
/// Returns `n + 1` points forming a closed circle with centre
/// `(r/(r+1), 0)` and radius `1/(r+1)`.  All points lie on or
/// inside the unit disk for any `r ≥ 0`.
pub fn smith_r_circle_points(r: f64, n: usize) -> Vec<[f64; 2]> {
    let centre_re = r / (r + 1.0);
    let radius = 1.0 / (r + 1.0);
    let mut pts = Vec::with_capacity(n + 1);
    for i in 0..n {
        let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
        pts.push([centre_re + radius * theta.cos(), radius * theta.sin()]);
    }
    // Close the loop.
    pts.push(pts[0]);
    pts
}

/// Constant-X arc on the Smith chart, clipped to the unit disk.
///
/// The full constant-X circle has centre `(1, 1/x)` and radius `1/|x|`.
/// Only points satisfying `re² + im² ≤ 1.0 + 1e-9` are returned.
/// Panics if `x == 0.0`.
pub fn smith_x_arc_points(x: f64, n: usize) -> Vec<[f64; 2]> {
    assert!(x != 0.0, "smith_x_arc_points: x must be non-zero");
    let centre_re = 1.0_f64;
    let centre_im = 1.0 / x;
    let radius = 1.0 / x.abs();
    (0..n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            [
                centre_re + radius * theta.cos(),
                centre_im + radius * theta.sin(),
            ]
        })
        .filter(|&[re, im]| re * re + im * im <= 1.0 + 1e-9)
        .collect()
}

/// Overlay multiple S-parameter traces (pre-built by [`build_sparam_series`])
/// as magnitude-dB lines with an automatic `egui_plot` legend.
///
/// Each [`SparamSeries`] is rendered as a distinct `Line`; `egui_plot` assigns
/// colours automatically based on the line name, giving consistent colouring
/// per label across repaints.
///
/// When `series` is empty the plot area is shown but contains no data.
pub fn show_sparams_db_plot(ui: &mut egui::Ui, series: &[SparamSeries]) {
    Plot::new("sparams_db_plot")
        .x_axis_label("Frequency (GHz)")
        .y_axis_label("|S| (dB)")
        .legend(Legend::default())
        .show(ui, |plot_ui| {
            for s in series {
                let pts: PlotPoints = s.points.iter().copied().collect();
                plot_ui.line(Line::new(s.label.clone(), pts));
            }
        });
}

/// Plot `20·log10|S11|` against frequency (in GHz on the X axis).
///
/// `freq_hz` and `s11` are parallel slices; they must have the same length.
/// The caller owns the [`egui::Ui`] this is rendered into.
///
/// This is a thin single-trace wrapper kept for backward compatibility with
/// call sites that already have a raw `s11` slice. For multi-trace use,
/// prefer [`build_sparam_series`] + [`show_sparams_db_plot`].
pub fn show_s11_db_plot(ui: &mut egui::Ui, freq_hz: &[f64], s11: &[Complex64]) {
    assert_eq!(
        freq_hz.len(),
        s11.len(),
        "freq_hz and s11 must have equal length"
    );
    let points: PlotPoints = freq_hz
        .iter()
        .zip(s11.iter())
        .map(|(f, z)| [f * 1.0e-9, db_from_complex(*z)])
        .collect();
    let line = Line::new("|S11| (dB)", points);
    Plot::new("s11_db_plot")
        .x_axis_label("Frequency (GHz)")
        .y_axis_label("|S11| (dB)")
        .show(ui, |plot_ui| {
            plot_ui.line(line);
        });
}

/// Plot S-parameter trajectories on a Smith-chart canvas.
///
/// Draws (in order):
/// 1. The unit circle `|Γ| = 1` as a reference boundary.
/// 2. Constant-R circles for `r ∈ [0.2, 0.5, 1.0, 2.0, 5.0]`.
/// 3. Constant-X arcs for `x ∈ [±0.2, ±0.5, ±1.0, ±2.0, ±5.0]`.
/// 4. One coloured line per [`SmithSeries`] in `series`.
///
/// The plot uses a locked 1:1 data aspect ratio so circles appear round, and
/// an `egui_plot` legend so trace labels are visible.
pub fn show_smith_chart(ui: &mut egui::Ui, series: &[SmithSeries]) {
    const R_VALUES: &[f64] = &[0.2, 0.5, 1.0, 2.0, 5.0];
    const X_VALUES: &[f64] = &[0.2, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0];

    Plot::new("smith_chart")
        .data_aspect(1.0)
        .x_axis_label("Re")
        .y_axis_label("Im")
        .legend(Legend::default())
        .show(ui, |plot_ui| {
            // 1. Unit circle.
            let unit: PlotPoints = unit_circle_points(256).into_iter().collect();
            plot_ui.line(Line::new("|Γ|=1", unit));

            // 2. Constant-R circles.
            for &r in R_VALUES {
                let pts: PlotPoints = smith_r_circle_points(r, 128).into_iter().collect();
                plot_ui.line(Line::new(format!("r={r}"), pts));
            }

            // 3. Constant-X arcs.
            for &x in X_VALUES {
                let pts: PlotPoints = smith_x_arc_points(x, 256).into_iter().collect();
                plot_ui.line(Line::new(format!("x={x}"), pts));
            }

            // 4. Data traces.
            for s in series {
                let pts: PlotPoints = s.points.iter().copied().collect();
                plot_ui.line(Line::new(s.label.clone(), pts));
            }
        });
}

// ----------------------------------------------------------------------------
// Unit tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;
    use yee_io::touchstone::{File as TsFile, Format, FreqUnit};

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /// Build a minimal synthetic 2-port `TsFile` for series-builder tests.
    ///
    /// Layout: 5 frequency points from 1–5 GHz. The S-matrix at each point
    /// has `S11 = 0.5+0j`, `S21 = 0.3+0j`, `S12 = 0.2+0j`, `S22 = 0.4+0j`.
    fn two_port_file() -> TsFile {
        let n = 5usize;
        let freq_hz: Vec<f64> = (1..=n).map(|i| i as f64 * 1.0e9).collect();
        let s11 = Complex64::new(0.5, 0.0);
        let s21 = Complex64::new(0.3, 0.0);
        let s12 = Complex64::new(0.2, 0.0);
        let s22 = Complex64::new(0.4, 0.0);
        // Row-major: [S11, S12, S21, S22] (indices [0,1,2,3])
        let data: Vec<Vec<Complex64>> = (0..n).map(|_| vec![s11, s12, s21, s22]).collect();
        TsFile {
            n_ports: 2,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz,
            data,
            comments: vec![],
        }
    }

    // -------------------------------------------------------------------------
    // db_from_complex
    // -------------------------------------------------------------------------

    #[test]
    fn test_dipole_db_from_complex() {
        // |0.5 + 0j| = 0.5 → 20·log10(0.5) ≈ -6.0206 dB
        let s11 = Complex64::new(0.5, 0.0);
        let db = db_from_complex(s11);
        assert!((db - (-6.020_599_913_279_624)).abs() < 1e-9, "db = {db}");
    }

    #[test]
    fn test_db_from_complex_zero_is_neg_inf() {
        // Edge case: |0| → -inf dB (callers must not panic on this).
        let db = db_from_complex(Complex64::new(0.0, 0.0));
        assert!(db.is_infinite() && db < 0.0, "db = {db}");
    }

    // -------------------------------------------------------------------------
    // unit_circle_points
    // -------------------------------------------------------------------------

    #[test]
    fn test_smith_unit_circle_has_n_points() {
        let n = 64;
        let pts = unit_circle_points(n);
        // n samples + 1 closing duplicate.
        assert_eq!(pts.len(), n + 1);
        // Every point should sit on the unit circle.
        for [x, y] in &pts {
            let r = (x * x + y * y).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "off-circle point: ({x}, {y})");
        }
        // Loop closes.
        assert_eq!(pts.first().unwrap(), pts.last().unwrap());
    }

    // -------------------------------------------------------------------------
    // build_sparam_series — pure series-builder unit tests
    // -------------------------------------------------------------------------

    /// `Selection::All` on a 2-port file → 4 series labelled S11/S12/S21/S22
    /// (row-major order) each with the correct point count.
    #[test]
    fn series_all_two_port_produces_four_series() {
        let file = two_port_file();
        let series = build_sparam_series(&file, &Selection::All);
        assert_eq!(series.len(), 4, "expected 4 series for 2-port All");

        let labels: Vec<&str> = series.iter().map(|s| s.label.as_str()).collect();
        assert_eq!(labels, vec!["S11", "S12", "S21", "S22"]);

        for s in &series {
            assert_eq!(
                s.points.len(),
                file.freq_hz.len(),
                "series '{}' has wrong point count",
                s.label
            );
        }
    }

    /// `Selection::Diagonal(0)` → a single `S11` series; default is unchanged.
    #[test]
    fn series_diagonal_zero_produces_s11_only() {
        let file = two_port_file();
        let series = build_sparam_series(&file, &Selection::Diagonal(0));
        assert_eq!(series.len(), 1);
        assert_eq!(series[0].label, "S11");
        assert_eq!(series[0].points.len(), file.freq_hz.len());
    }

    /// `Selection::Diagonal(1)` → a single `S22` series.
    #[test]
    fn series_diagonal_one_produces_s22() {
        let file = two_port_file();
        let series = build_sparam_series(&file, &Selection::Diagonal(1));
        assert_eq!(series.len(), 1);
        assert_eq!(series[0].label, "S22");
    }

    /// `Selection::Entries` with specific pairs → matching labels in request order.
    #[test]
    fn series_entries_specific_pairs() {
        let file = two_port_file();
        // Request S21 then S11 (non-alphabetical order to verify ordering is preserved).
        let sel = Selection::Entries(vec![(1, 0), (0, 0)]);
        let series = build_sparam_series(&file, &sel);
        assert_eq!(series.len(), 2);
        assert_eq!(series[0].label, "S21");
        assert_eq!(series[1].label, "S11");
    }

    /// Out-of-range diagonal index → empty result, no panic.
    #[test]
    fn series_out_of_range_diagonal_is_empty() {
        let file = two_port_file();
        // n_ports=2, so index 5 is out of range.
        let series = build_sparam_series(&file, &Selection::Diagonal(5));
        assert!(
            series.is_empty(),
            "out-of-range diagonal should yield empty Vec"
        );
    }

    /// Out-of-range entry in `Entries` is silently skipped; in-range entries
    /// are still returned.
    #[test]
    fn series_out_of_range_entry_is_skipped() {
        let file = two_port_file();
        // (99, 99) is out of range; (0, 1) = S12 is valid.
        let sel = Selection::Entries(vec![(99, 99), (0, 1)]);
        let series = build_sparam_series(&file, &sel);
        assert_eq!(series.len(), 1, "only in-range entry should be returned");
        assert_eq!(series[0].label, "S12");
    }

    /// X-axis values are in GHz (freq_hz * 1e-9) and Y-axis is dB.
    /// For S11 = 0.5+0j the dB value is ≈ -6.0206 dB.
    #[test]
    fn series_points_have_correct_freq_ghz_and_db_values() {
        let file = two_port_file();
        let series = build_sparam_series(&file, &Selection::Diagonal(0));
        assert_eq!(series.len(), 1);
        let pts = &series[0].points;
        for (k, &[x, y]) in pts.iter().enumerate() {
            // X in GHz.
            let expected_ghz = file.freq_hz[k] * 1.0e-9;
            assert!(
                (x - expected_ghz).abs() < 1e-12,
                "point {k}: x={x}, expected {expected_ghz}"
            );
            // Y = db_from_complex(S11) = 20*log10(0.5) ≈ -6.0206.
            let expected_db = db_from_complex(file.data[k][0]);
            assert!(
                (y - expected_db).abs() < 1e-9,
                "point {k}: y={y}, expected {expected_db}"
            );
        }
    }

    // -------------------------------------------------------------------------
    // smith_r_circle_points
    // -------------------------------------------------------------------------

    /// `smith_r_circle_points(1.0, 64)` returns exactly 65 points (n + 1),
    /// the loop is closed (first == last), and every point lies on the circle
    /// with centre (0.5, 0) and radius 0.5.
    #[test]
    fn smith_r_circle_returns_correct_n_plus_1_points() {
        let pts = smith_r_circle_points(1.0, 64);
        assert_eq!(pts.len(), 65, "expected n+1 = 65 points");
        assert_eq!(pts.first().unwrap(), pts.last().unwrap(), "loop not closed");
        let centre = [0.5_f64, 0.0_f64];
        let radius = 0.5_f64;
        for &[re, im] in &pts {
            let dist = ((re - centre[0]).powi(2) + (im - centre[1]).powi(2)).sqrt();
            assert!(
                (dist - radius).abs() < 1e-9,
                "point ({re},{im}) not on circle: dist={dist}"
            );
        }
    }

    /// Every point of a constant-R circle must lie inside the unit disk.
    #[test]
    fn smith_r_circle_contained_in_unit_disk() {
        let pts = smith_r_circle_points(0.5, 128);
        for &[re, im] in &pts {
            let mag2 = re * re + im * im;
            assert!(
                mag2 <= 1.0 + 1e-9,
                "point ({re},{im}) outside unit disk: |Γ|²={mag2}"
            );
        }
    }

    // -------------------------------------------------------------------------
    // smith_x_arc_points
    // -------------------------------------------------------------------------

    /// Every point returned by `smith_x_arc_points(1.0, 256)` must lie on the
    /// circle with centre (1, 1) and radius 1.
    #[test]
    fn smith_x_arc_points_on_circle() {
        let pts = smith_x_arc_points(1.0, 256);
        assert!(!pts.is_empty(), "expected non-empty arc for x=1");
        let centre = [1.0_f64, 1.0_f64];
        let radius = 1.0_f64;
        for &[re, im] in &pts {
            let dist = ((re - centre[0]).powi(2) + (im - centre[1]).powi(2)).sqrt();
            assert!(
                (dist - radius).abs() < 1e-9,
                "point ({re},{im}) not on circle: dist={dist}"
            );
        }
    }

    /// Every point returned by `smith_x_arc_points` must be inside the unit disk.
    #[test]
    fn smith_x_arc_points_inside_unit_disk() {
        let pts = smith_x_arc_points(0.2, 512);
        for &[re, im] in &pts {
            let mag2 = re * re + im * im;
            assert!(
                mag2 <= 1.0 + 1e-9,
                "point ({re},{im}) outside unit disk: |Γ|²={mag2}"
            );
        }
    }

    // -------------------------------------------------------------------------
    // build_smith_series
    // -------------------------------------------------------------------------

    /// `Selection::All` on a 2-port file → 4 series labelled S11/S12/S21/S22,
    /// each with 5 points; the first series first point is (Re=0.5, Im=0.0).
    #[test]
    fn build_smith_series_all_two_port() {
        let file = two_port_file();
        let series = build_smith_series(&file, &Selection::All);
        assert_eq!(series.len(), 4, "expected 4 series for 2-port All");

        let labels: Vec<&str> = series.iter().map(|s| s.label.as_str()).collect();
        assert_eq!(labels, vec!["S11", "S12", "S21", "S22"]);

        for s in &series {
            assert_eq!(
                s.points.len(),
                file.freq_hz.len(),
                "series '{}' has wrong point count",
                s.label
            );
        }

        // S11 = 0.5+0j → first point is [0.5, 0.0].
        let &[re, im] = &series[0].points[0];
        assert!((re - 0.5).abs() < 1e-12, "S11 Re: expected 0.5, got {re}");
        assert!(im.abs() < 1e-12, "S11 Im: expected 0.0, got {im}");
    }

    /// `build_smith_series` on a 1-port file with `Selection::All` → 1 series
    /// labelled `"S11"` with `[Re, Im]` points (not `[freq_ghz, dB]`).
    #[test]
    fn build_smith_series_one_port_produces_one_series_with_re_im_points() {
        let n_pts = 3usize;
        let file = TsFile {
            n_ports: 1,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz: (1..=n_pts).map(|i| i as f64 * 1.0e9).collect(),
            data: (0..n_pts).map(|_| vec![Complex64::new(0.5, 0.3)]).collect(),
            comments: vec![],
        };
        let series = build_smith_series(&file, &Selection::All);
        assert_eq!(series.len(), 1);
        assert_eq!(series[0].label, "S11");
        assert_eq!(series[0].points.len(), n_pts);
        // Each point must be [Re, Im], not [freq_ghz, dB].
        let &[re, im] = &series[0].points[0];
        assert!((re - 0.5).abs() < 1e-12, "Re: expected 0.5, got {re}");
        assert!((im - 0.3).abs() < 1e-12, "Im: expected 0.3, got {im}");
    }
}
