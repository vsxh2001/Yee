//! Static PNG/SVG export of S-parameter plots via the [`plotters`] crate.
//!
//! This crate is the "save the plot to a file" counterpart to the live
//! `egui_plot` views in `yee-gui`. It is consumed by:
//!
//! * Validation harnesses (per-case S-parameter sanity plots).
//! * Examples (golden-image generation).
//! * CI pipelines (artifact uploads — PNG for humans, SVG for diffing).
//!
//! ## Public surface
//!
//! Three entry points, one per plot type:
//!
//! * [`plot_s11_db`] — `|S₁₁|` in decibels vs. frequency.
//! * [`plot_s11_phase`] — `arg(S₁₁)` in degrees vs. frequency.
//! * [`plot_smith_chart`] — `S₁₁` on the complex unit disk (Smith-style).
//!
//! All three share a [`PlotConfig`] (size, title, output [`PlotFormat`]) and
//! dispatch to either a `BitMapBackend` (PNG) or an `SVGBackend` (SVG)
//! depending on `config.format`.

use std::path::Path;

use num_complex::Complex64;
use plotters::prelude::*;

/// Minimum dB value we will report. `|S₁₁| = 0` would otherwise become
/// `-∞ dB`, which both the renderer and downstream consumers (humans!)
/// dislike. We clamp instead.
pub const MIN_DB: f64 = -200.0;

/// Output image format for [`PlotConfig`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotFormat {
    /// Raster PNG via `plotters::backend::BitMapBackend`.
    Png,
    /// Vector SVG via `plotters::backend::SVGBackend`.
    Svg,
}

/// Common rendering knobs for all `plot_*` functions in this crate.
#[derive(Debug, Clone)]
pub struct PlotConfig {
    /// Image width in pixels (or SVG user units).
    pub width_px: u32,
    /// Image height in pixels (or SVG user units).
    pub height_px: u32,
    /// Plot title shown above the chart area.
    pub title: String,
    /// Output format — picks the `plotters` backend.
    pub format: PlotFormat,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            width_px: 800,
            height_px: 600,
            title: String::new(),
            format: PlotFormat::Png,
        }
    }
}

/// Error type returned by every `plot_*` function in this crate.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Filesystem or path-level error (e.g. parent directory missing).
    #[error("io error: {0}")]
    Io(String),
    /// Rendering error bubbled up from a `plotters` backend.
    #[error("render error: {0}")]
    Render(String),
}

/// Map a `plotters` `DrawingAreaErrorKind` into our [`Error`] type.
///
/// We collapse the backend-specific error variant into a string, because
/// `DrawingAreaErrorKind<E>` is generic over the backend and we want a
/// single [`Error`] enum the caller can match against without knowing
/// whether the underlying writer is a bitmap or an SVG.
fn map_render_err<E: std::error::Error + Send + Sync + 'static>(
    e: plotters::drawing::DrawingAreaErrorKind<E>,
) -> Error {
    Error::Render(format!("{e}"))
}

/// Convert a complex `S₁₁` sample to magnitude in dB, clamped at [`MIN_DB`].
///
/// `20 · log10 |z|` blows up to `-∞` for `|z| = 0`; instead of letting
/// that propagate through the plot we clamp to a sane floor.
#[inline]
pub fn db_clamped(z: Complex64) -> f64 {
    let mag = z.norm();
    if mag <= 0.0 {
        MIN_DB
    } else {
        let db = 20.0 * mag.log10();
        db.max(MIN_DB)
    }
}

/// Convert a complex `S₁₁` sample to phase in degrees on `[-180, 180]`.
#[inline]
pub fn phase_degrees(z: Complex64) -> f64 {
    z.arg().to_degrees()
}

/// Plot `|S₁₁|` in dB vs. frequency.
///
/// * X-axis: linear, labelled `"frequency (GHz)"` — values are
///   `freq_hz[i] / 1e9`.
/// * Y-axis: linear dB, labelled `"|S₁₁| (dB)"`. Zero magnitude is
///   clamped to [`MIN_DB`] to keep the y-range finite.
///
/// `freq_hz` and `s11` must have the same length.
pub fn plot_s11_db(
    freq_hz: &[f64],
    s11: &[Complex64],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    assert_eq!(
        freq_hz.len(),
        s11.len(),
        "freq_hz and s11 must have equal length"
    );

    let xs_ghz: Vec<f64> = freq_hz.iter().map(|f| f * 1.0e-9).collect();
    let ys_db: Vec<f64> = s11.iter().map(|z| db_clamped(*z)).collect();

    let (x_min, x_max) = finite_range(&xs_ghz, 0.0, 1.0);
    let (y_min_raw, y_max_raw) = finite_range(&ys_db, MIN_DB, 0.0);
    // Pad the y-range so the trace doesn't sit on the frame.
    let y_pad = ((y_max_raw - y_min_raw).abs() * 0.05).max(1.0);
    let y_min = y_min_raw - y_pad;
    let y_max = y_max_raw + y_pad;

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_xy_line(
                &root,
                &config.title,
                "frequency (GHz)",
                "|S₁₁| (dB)",
                &xs_ghz,
                &ys_db,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_xy_line(
                &root,
                &config.title,
                "frequency (GHz)",
                "|S₁₁| (dB)",
                &xs_ghz,
                &ys_db,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
    }
}

/// Plot the phase of `S₁₁` in degrees vs. frequency.
///
/// * X-axis: `"frequency (GHz)"`.
/// * Y-axis: `"phase (deg)"`, fixed to `[-180, 180]` because that's the
///   natural range of `Complex64::arg().to_degrees()`.
pub fn plot_s11_phase(
    freq_hz: &[f64],
    s11: &[Complex64],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    assert_eq!(
        freq_hz.len(),
        s11.len(),
        "freq_hz and s11 must have equal length"
    );

    let xs_ghz: Vec<f64> = freq_hz.iter().map(|f| f * 1.0e-9).collect();
    let ys_deg: Vec<f64> = s11.iter().map(|z| phase_degrees(*z)).collect();

    let (x_min, x_max) = finite_range(&xs_ghz, 0.0, 1.0);
    // Phase range is bounded by definition.
    let (y_min, y_max) = (-180.0_f64, 180.0_f64);

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_xy_line(
                &root,
                &config.title,
                "frequency (GHz)",
                "phase (deg)",
                &xs_ghz,
                &ys_deg,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_xy_line(
                &root,
                &config.title,
                "frequency (GHz)",
                "phase (deg)",
                &xs_ghz,
                &ys_deg,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
    }
}

/// Plot `S₁₁` on the complex unit disk (Smith-chart style).
///
/// We don't render the full Smith-chart constant-resistance / constant-
/// reactance arc family — for now the chart is just:
///
/// 1. A light-grey reference unit circle (200 samples).
/// 2. An origin crosshair through `(0, 0)`.
/// 3. The `S₁₁(f)` trajectory as a connected red line.
///
/// The aspect ratio is locked to 1:1 by giving both axes the range
/// `[-1.1, 1.1]`; downstream callers should pick `width_px ≈ height_px`
/// for a round circle.
pub fn plot_smith_chart(
    s11: &[Complex64],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_smith(&root, &config.title, s11)
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_smith(&root, &config.title, s11)
        }
    }
}

// ---------------------------------------------------------------------------
// Backend-agnostic drawing helpers
// ---------------------------------------------------------------------------

/// Render an `(xs, ys)` line plot into `root` with explicit axis ranges.
///
/// Factored out so the `PNG` and `SVG` arms of [`plot_s11_db`] /
/// [`plot_s11_phase`] don't repeat the chart-builder boilerplate.
#[allow(clippy::too_many_arguments)]
fn draw_xy_line<DB>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    xs: &[f64],
    ys: &[f64],
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
) -> Result<(), Error>
where
    DB: DrawingBackend,
    DB::ErrorType: 'static,
{
    root.fill(&WHITE).map_err(map_render_err)?;

    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(map_render_err)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()
        .map_err(map_render_err)?;

    let series: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
    chart
        .draw_series(LineSeries::new(series, RED.stroke_width(2)))
        .map_err(map_render_err)?;

    root.present().map_err(map_render_err)?;
    Ok(())
}

/// Render a Smith-style chart (reference unit circle + crosshair + S₁₁
/// trajectory) into `root`.
fn draw_smith<DB>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    s11: &[Complex64],
) -> Result<(), Error>
where
    DB: DrawingBackend,
    DB::ErrorType: 'static,
{
    root.fill(&WHITE).map_err(map_render_err)?;

    let lim = 1.1_f64;
    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-lim..lim, -lim..lim)
        .map_err(map_render_err)?;

    chart
        .configure_mesh()
        .x_desc("Re S₁₁")
        .y_desc("Im S₁₁")
        .draw()
        .map_err(map_render_err)?;

    // Origin crosshair (light grey).
    let crosshair_style = RGBColor(200, 200, 200).stroke_width(1);
    chart
        .draw_series(LineSeries::new(
            [(-1.0_f64, 0.0_f64), (1.0, 0.0)],
            crosshair_style,
        ))
        .map_err(map_render_err)?;
    chart
        .draw_series(LineSeries::new(
            [(0.0_f64, -1.0_f64), (0.0, 1.0)],
            crosshair_style,
        ))
        .map_err(map_render_err)?;

    // Reference unit circle (200 samples, closed).
    let n = 200usize;
    let unit: Vec<(f64, f64)> = (0..=n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            (theta.cos(), theta.sin())
        })
        .collect();
    chart
        .draw_series(LineSeries::new(
            unit,
            RGBColor(160, 160, 160).stroke_width(1),
        ))
        .map_err(map_render_err)?;

    // S₁₁ trajectory in red.
    let traj: Vec<(f64, f64)> = s11.iter().map(|z| (z.re, z.im)).collect();
    chart
        .draw_series(LineSeries::new(traj, RED.stroke_width(2)))
        .map_err(map_render_err)?;

    root.present().map_err(map_render_err)?;
    Ok(())
}

/// Compute `(min, max)` over a slice, ignoring non-finite values and
/// falling back to `(default_min, default_max)` if the slice is empty
/// or entirely non-finite.
fn finite_range(xs: &[f64], default_min: f64, default_max: f64) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &x in xs {
        if x.is_finite() {
            if x < min {
                min = x;
            }
            if x > max {
                max = x;
            }
        }
    }
    if !min.is_finite() || !max.is_finite() || min == max {
        if min.is_finite() && min == max {
            // Single-valued slice — give it a unit-wide window.
            return (min - 0.5, max + 0.5);
        }
        return (default_min, default_max);
    }
    (min, max)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::NamedTempFile;

    /// Synthetic 21-point S₁₁ sweep with a resonant dip near the middle.
    /// Used as input for the file-writing tests.
    fn synthetic_sweep() -> (Vec<f64>, Vec<Complex64>) {
        let n = 21;
        let f0 = 1.0e9;
        let f1 = 3.0e9;
        let fr = 2.0e9;
        let bw = 0.2e9;
        let freq: Vec<f64> = (0..n)
            .map(|i| f0 + (f1 - f0) * (i as f64) / ((n - 1) as f64))
            .collect();
        let s11: Vec<Complex64> = freq
            .iter()
            .map(|&f| {
                // Lorentzian-ish dip in magnitude with a phase swing.
                let x = (f - fr) / bw;
                let mag = (x * x / (1.0 + x * x)).sqrt(); // 0 at fr, → 1 far away
                let phase = (-2.0 * x).atan();
                Complex64::from_polar(mag, phase)
            })
            .collect();
        (freq, s11)
    }

    #[test]
    fn test_plot_s11_db_writes_png() {
        let (freq, s11) = synthetic_sweep();
        let tmp = NamedTempFile::with_suffix(".png").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 640,
            height_px: 480,
            title: "S11 dB test".to_string(),
            format: PlotFormat::Png,
        };
        plot_s11_db(&freq, &s11, tmp.path(), &cfg).expect("plot_s11_db");
        let len = fs::metadata(tmp.path()).expect("metadata").len();
        assert!(len > 1024, "PNG file is too small: {len} bytes");
    }

    #[test]
    fn test_plot_smith_chart_writes_svg() {
        let (_freq, s11) = synthetic_sweep();
        let tmp = NamedTempFile::with_suffix(".svg").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 600,
            height_px: 600,
            title: "Smith test".to_string(),
            format: PlotFormat::Svg,
        };
        plot_smith_chart(&s11, tmp.path(), &cfg).expect("plot_smith_chart");
        let body = fs::read_to_string(tmp.path()).expect("read svg");
        assert!(body.contains("<svg"), "SVG missing <svg tag: {body:.200}");
        assert!(body.len() > 256, "SVG body too short: {} bytes", body.len());
    }

    #[test]
    fn test_plot_s11_phase_within_range() {
        let (_freq, s11) = synthetic_sweep();
        for z in &s11 {
            let p = phase_degrees(*z);
            assert!((-180.0..=180.0).contains(&p), "phase out of range: {p}");
        }
    }

    #[test]
    fn test_db_clamp_at_zero_magnitude() {
        // |0+0j| -> MIN_DB, not -inf.
        let db = db_clamped(Complex64::new(0.0, 0.0));
        assert!(db.is_finite(), "clamped dB must be finite, got {db}");
        assert_eq!(db, MIN_DB, "expected exact clamp at MIN_DB, got {db}");

        // Tiny but non-zero magnitudes should also stay clamped.
        let db_tiny = db_clamped(Complex64::new(1.0e-300, 0.0));
        assert!(db_tiny.is_finite() && db_tiny >= MIN_DB - 1e-9);
        assert_eq!(db_tiny, MIN_DB);
    }

    #[test]
    fn test_db_clamp_normal_value() {
        // Sanity: |0.5+0j| → 20·log10(0.5) ≈ -6.0206 dB (well above MIN_DB).
        let db = db_clamped(Complex64::new(0.5, 0.0));
        assert!((db - (-6.020_599_913_279_624)).abs() < 1e-9, "db = {db}");
    }
}
