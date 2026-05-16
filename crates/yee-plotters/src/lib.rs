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
/// Stub — implementation lands in the next commit.
pub fn plot_smith_chart(
    _s11: &[Complex64],
    _out_path: &Path,
    _config: &PlotConfig,
) -> Result<(), Error> {
    Err(Error::Render("plot_smith_chart: not yet implemented".into()))
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
