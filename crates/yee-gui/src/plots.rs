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

use egui_plot::{Line, Plot, PlotPoints};
use num_complex::Complex64;

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

/// Plot `20·log10|S11|` against frequency (in GHz on the X axis).
///
/// `freq_hz` and `s11` are parallel slices; they must have the same length.
/// The caller owns the [`egui::Ui`] this is rendered into.
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

/// Plot the `S11` trajectory on a Smith-chart-style canvas: the unit circle
/// is drawn for reference, and the data points are plotted in the complex
/// plane with a locked 1:1 aspect ratio.
pub fn show_smith_chart(ui: &mut egui::Ui, s11: &[Complex64]) {
    let unit: PlotPoints = unit_circle_points(256).into_iter().collect();
    let traj: PlotPoints = s11.iter().map(|z| [z.re, z.im]).collect();

    Plot::new("smith_chart")
        .data_aspect(1.0)
        .x_axis_label("Re")
        .y_axis_label("Im")
        .show(ui, |plot_ui| {
            plot_ui.line(Line::new("|Γ| = 1", unit));
            plot_ui.line(Line::new("S11(f)", traj));
        });
}

// ----------------------------------------------------------------------------
// Unit tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

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
}
