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
//! Six entry points:
//!
//! * [`plot_s11_db`] — `|S₁₁|` in decibels vs. frequency (single trace).
//! * [`plot_s11_phase`] — `arg(S₁₁)` in degrees vs. frequency (single trace).
//! * [`plot_smith_chart`] — `S₁₁` on the complex unit disk (Smith-style),
//!   back-compat single-trace wrapper around [`plot_smith_chart_multi`].
//! * [`plot_smith_chart_multi`] — overlay multiple [`SmithTrace`] values on a
//!   Smith chart with constant-R/X arc reference overlays. [`SmithTrace`] is
//!   the Smith equivalent of [`SparamTrace`].
//! * [`plot_sparams_db`] — overlay multiple S-parameter traces in dB vs.
//!   frequency, each labelled and colour-coded with a legend.
//! * [`plot_sparams_phase`] — overlay multiple S-parameter traces (phase in
//!   degrees) with a legend.
//!
//! All share a [`PlotConfig`] (size, title, output [`PlotFormat`]) and dispatch
//! to either a `BitMapBackend` (PNG) or an `SVGBackend` (SVG) depending on
//! `config.format`.

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

/// A single named S-parameter trace for use with [`plot_sparams_db`] and
/// [`plot_sparams_phase`].
///
/// `label` appears in the legend (e.g. `"S11"`, `"S21"`).
/// `values` must have the same length as the `freq_hz` slice passed to the
/// plot function.
#[derive(Debug, Clone)]
pub struct SparamTrace {
    /// Human-readable label shown in the plot legend (e.g. `"S11"`, `"S21"`).
    pub label: String,
    /// Complex S-parameter samples, one per frequency point.
    pub values: Vec<Complex64>,
}

/// A single labelled S-parameter trace for Smith-chart export.
///
/// Parallel to [`SparamTrace`] for [`plot_sparams_db`].
/// Used with [`plot_smith_chart_multi`] to render multiple overlaid
/// trajectories on the Smith chart.
#[derive(Debug, Clone)]
pub struct SmithTrace {
    /// Legend label (e.g. `"S11"`, `"S22"`).
    pub label: String,
    /// Raw complex S-parameter values, one per frequency sample.
    pub values: Vec<Complex64>,
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
            let root =
                SVGBackend::new(out_path, (config.width_px, config.height_px)).into_drawing_area();
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
            let root =
                SVGBackend::new(out_path, (config.width_px, config.height_px)).into_drawing_area();
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
/// Back-compat single-trace wrapper around [`plot_smith_chart_multi`].
/// The chart includes:
///
/// 1. Constant-R circles and constant-X arc reference overlays (light grey).
/// 2. A light-grey reference unit circle (200 samples).
/// 3. An origin crosshair through `(0, 0)`.
/// 4. The `S₁₁(f)` trajectory as a connected coloured line.
///
/// The aspect ratio is locked to 1:1 by giving both axes the range
/// `[-1.1, 1.1]`; downstream callers should pick `width_px ≈ height_px`
/// for a round circle.
pub fn plot_smith_chart(
    s11: &[Complex64],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    let trace = SmithTrace {
        label: "S11".to_string(),
        values: s11.to_vec(),
    };
    plot_smith_chart_multi(&[trace], out_path, config)
}

/// Plot multiple S-parameter traces on a Smith-chart canvas with
/// constant-R/X arc reference overlays.
///
/// Each [`SmithTrace`] is rendered as a connected polyline. When only one
/// trace is provided the behaviour matches the original `plot_smith_chart`.
///
/// The reference arc grid is always drawn (light grey, `RGBColor(210,210,210)`):
/// - constant-R circles for r ∈ {0.2, 0.5, 1.0, 2.0, 5.0}
/// - constant-X arcs for x ∈ {±0.2, ±0.5, ±1.0, ±2.0, ±5.0}
///
/// The aspect ratio is locked to 1:1 by giving both axes the range
/// `[-1.1, 1.1]`; downstream callers should pick `width_px ≈ height_px`
/// for a round circle.
pub fn plot_smith_chart_multi(
    traces: &[SmithTrace],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_smith_multi(&root, &config.title, traces)
        }
        PlotFormat::Svg => {
            let root =
                SVGBackend::new(out_path, (config.width_px, config.height_px)).into_drawing_area();
            draw_smith_multi(&root, &config.title, traces)
        }
    }
}

/// Overlay multiple S-parameter traces as magnitude-dB lines with a legend.
///
/// Each [`SparamTrace`] in `traces` is drawn as a separate labelled line in a
/// distinct colour from a small fixed palette (up to 8 colours; traces beyond
/// the palette cycle back to the first colour). A legend entry is added for
/// each trace so the output is self-documenting.
///
/// * X-axis: `"frequency (GHz)"`.
/// * Y-axis: `"|S| (dB)"`. The y-range covers all traces with a 5% pad.
///
/// `freq_hz` must have the same length as every `trace.values` slice.
///
/// # Errors
///
/// Returns [`Error::Render`] if the `plotters` backend fails.
pub fn plot_sparams_db(
    freq_hz: &[f64],
    traces: &[SparamTrace],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    for t in traces {
        assert_eq!(
            freq_hz.len(),
            t.values.len(),
            "freq_hz and trace '{}' must have equal length",
            t.label
        );
    }

    let xs_ghz: Vec<f64> = freq_hz.iter().map(|f| f * 1.0e-9).collect();
    let all_db: Vec<Vec<f64>> = traces
        .iter()
        .map(|t| t.values.iter().map(|z| db_clamped(*z)).collect())
        .collect();

    let (x_min, x_max) = finite_range(&xs_ghz, 0.0, 1.0);
    let (y_min, y_max) = multi_trace_y_range(&all_db, MIN_DB, 0.0);

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_multi_trace(
                &root,
                &config.title,
                "frequency (GHz)",
                "|S| (dB)",
                &xs_ghz,
                traces,
                &all_db,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
        PlotFormat::Svg => {
            let root =
                SVGBackend::new(out_path, (config.width_px, config.height_px)).into_drawing_area();
            draw_multi_trace(
                &root,
                &config.title,
                "frequency (GHz)",
                "|S| (dB)",
                &xs_ghz,
                traces,
                &all_db,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
    }
}

/// Overlay multiple S-parameter traces as phase-in-degrees lines with a legend.
///
/// Each [`SparamTrace`] in `traces` is drawn as a separate labelled line in a
/// distinct colour. The y-axis is fixed to `[-180, 180]` degrees (the natural
/// range of `Complex64::arg().to_degrees()`). A legend entry is added for each
/// trace.
///
/// * X-axis: `"frequency (GHz)"`.
/// * Y-axis: `"phase (deg)"`, fixed to `[-180, 180]`.
///
/// `freq_hz` must have the same length as every `trace.values` slice.
///
/// # Errors
///
/// Returns [`Error::Render`] if the `plotters` backend fails.
pub fn plot_sparams_phase(
    freq_hz: &[f64],
    traces: &[SparamTrace],
    out_path: &Path,
    config: &PlotConfig,
) -> Result<(), Error> {
    for t in traces {
        assert_eq!(
            freq_hz.len(),
            t.values.len(),
            "freq_hz and trace '{}' must have equal length",
            t.label
        );
    }

    let xs_ghz: Vec<f64> = freq_hz.iter().map(|f| f * 1.0e-9).collect();
    let all_deg: Vec<Vec<f64>> = traces
        .iter()
        .map(|t| t.values.iter().map(|z| phase_degrees(*z)).collect())
        .collect();

    let (x_min, x_max) = finite_range(&xs_ghz, 0.0, 1.0);
    let (y_min, y_max) = (-180.0_f64, 180.0_f64);

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(out_path, (config.width_px, config.height_px))
                .into_drawing_area();
            draw_multi_trace(
                &root,
                &config.title,
                "frequency (GHz)",
                "phase (deg)",
                &xs_ghz,
                traces,
                &all_deg,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
        PlotFormat::Svg => {
            let root =
                SVGBackend::new(out_path, (config.width_px, config.height_px)).into_drawing_area();
            draw_multi_trace(
                &root,
                &config.title,
                "frequency (GHz)",
                "phase (deg)",
                &xs_ghz,
                traces,
                &all_deg,
                x_min,
                x_max,
                y_min,
                y_max,
            )
        }
    }
}

// ---------------------------------------------------------------------------
// Backend-agnostic drawing helpers
// ---------------------------------------------------------------------------

/// Fixed colour palette for multi-trace plots.
///
/// Up to 8 distinct colours; traces beyond index 7 wrap around. All colours
/// are chosen for reasonable contrast on a white background.
fn trace_colour(idx: usize) -> RGBColor {
    const PALETTE: &[RGBColor] = &[
        RGBColor(0xE6, 0x19, 0x4B), // red
        RGBColor(0x43, 0x63, 0xD8), // blue
        RGBColor(0x3C, 0xB4, 0x4B), // green
        RGBColor(0xFF, 0x7F, 0x00), // orange
        RGBColor(0x91, 0x1E, 0xB4), // purple
        RGBColor(0x42, 0xD4, 0xF4), // cyan
        RGBColor(0xF0, 0x32, 0xE6), // magenta
        RGBColor(0x80, 0x80, 0x00), // olive
    ];
    PALETTE[idx % PALETTE.len()]
}

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

/// Render multiple overlaid traces with a legend into `root`.
///
/// Each trace in `traces` is drawn as a coloured line using [`trace_colour`];
/// a legend entry (coloured swatch + label) is appended to each series.
#[allow(clippy::too_many_arguments)]
fn draw_multi_trace<DB>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    xs: &[f64],
    traces: &[SparamTrace],
    all_ys: &[Vec<f64>],
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
        // Reserve right margin for the legend so it does not overlap the chart.
        .right_y_label_area_size(80)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(map_render_err)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()
        .map_err(map_render_err)?;

    for (i, (trace, ys)) in traces.iter().zip(all_ys.iter()).enumerate() {
        let colour = trace_colour(i);
        let series_data: Vec<(f64, f64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
        let label = trace.label.clone();
        chart
            .draw_series(LineSeries::new(series_data, colour.stroke_width(2)))
            .map_err(map_render_err)?
            .label(label)
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], colour.stroke_width(2))
            });
    }

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .map_err(map_render_err)?;

    root.present().map_err(map_render_err)?;
    Ok(())
}

/// Compute sample points on the constant-resistance circle for normalised
/// resistance `r` on the Smith chart.
///
/// Centre = `(r/(r+1), 0)`, radius = `1/(r+1)`.
/// Returns `n + 1` points forming a closed loop (last point equals first).
fn smith_r_circle_pts(r: f64, n: usize) -> Vec<(f64, f64)> {
    let centre_re = r / (r + 1.0);
    let radius = 1.0 / (r + 1.0);
    (0..=n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            (centre_re + radius * theta.cos(), radius * theta.sin())
        })
        .collect()
}

/// Compute sample points on the constant-reactance arc for normalised
/// reactance `x` (`x ≠ 0`) on the Smith chart, clipped to the unit disk.
///
/// Full circle: centre = `(1, 1/x)`, radius = `1/|x|`.
/// Only points satisfying `re² + im² ≤ 1.0 + 1e-9` are returned.
fn smith_x_arc_pts(x: f64, n: usize) -> Vec<(f64, f64)> {
    debug_assert!(x != 0.0, "smith_x_arc_pts: x must be non-zero");
    let centre_im = 1.0 / x;
    let radius = 1.0 / x.abs();
    (0..n)
        .filter_map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            let re = 1.0 + radius * theta.cos();
            let im = centre_im + radius * theta.sin();
            if re * re + im * im <= 1.0 + 1e-9 {
                Some((re, im))
            } else {
                None
            }
        })
        .collect()
}

/// Compute sample points on the constant-VSWR circle in the Γ-plane.
///
/// The locus VSWR = (1+|Γ|)/(1−|Γ|) = const is a circle centred at the
/// origin with radius ρ = (VSWR−1)/(VSWR+1) (Pozar §2.5).
/// Returns `n + 1` points forming a closed loop (last point equals first).
fn smith_vswr_circle_pts(vswr: f64, n: usize) -> Vec<(f64, f64)> {
    debug_assert!(vswr > 1.0, "smith_vswr_circle_pts: vswr must be > 1.0");
    let rho = (vswr - 1.0) / (vswr + 1.0);
    let first = (rho, 0.0_f64);
    let mut pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n as f64);
            (rho * theta.cos(), rho * theta.sin())
        })
        .collect();
    pts.push(first);
    pts
}

/// Fixed colour palette for multi-trace Smith chart plots.
///
/// Eight distinct colours; traces beyond index 7 wrap around. Matches the
/// palette used by the `draw_multi_trace` family of functions.
const SMITH_PALETTE: &[RGBColor] = &[
    RGBColor(228, 26, 28),
    RGBColor(55, 126, 184),
    RGBColor(77, 175, 74),
    RGBColor(255, 127, 0),
    RGBColor(152, 78, 163),
    RGBColor(166, 86, 40),
    RGBColor(247, 129, 191),
    RGBColor(153, 153, 153),
];

/// Render multiple Smith-chart traces with constant-R/X arc overlays into
/// `root`.
///
/// Drawing order:
/// 1. White fill.
/// 2. Chart frame and mesh.
/// 3. Origin crosshair (light grey).
/// 4. Unit circle (light grey).
/// 5. Constant-R circles (very light grey, stroke 1).
/// 6. Constant-X arcs (very light grey, stroke 1).
/// 7. VSWR circles for VSWR ∈ {1.5, 2, 3, 5, 10} (light blue-grey, stroke 1).
/// 8. Data traces with legend.
fn draw_smith_multi<DB>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    title: &str,
    traces: &[SmithTrace],
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
    let n_unit = 200usize;
    let unit: Vec<(f64, f64)> = (0..=n_unit)
        .map(|i| {
            let theta = (i as f64) * std::f64::consts::TAU / (n_unit as f64);
            (theta.cos(), theta.sin())
        })
        .collect();
    chart
        .draw_series(LineSeries::new(
            unit,
            RGBColor(160, 160, 160).stroke_width(1),
        ))
        .map_err(map_render_err)?;

    // Constant-R circles.
    let arc_style = RGBColor(210, 210, 210).stroke_width(1);
    for &r in &[0.2_f64, 0.5, 1.0, 2.0, 5.0] {
        let pts = smith_r_circle_pts(r, 128);
        chart
            .draw_series(LineSeries::new(pts, arc_style))
            .map_err(map_render_err)?;
    }

    // Constant-X arcs.
    for &x in &[0.2_f64, 0.5, 1.0, 2.0, 5.0, -0.2, -0.5, -1.0, -2.0, -5.0] {
        let pts = smith_x_arc_pts(x, 256);
        chart
            .draw_series(LineSeries::new(pts, arc_style))
            .map_err(map_render_err)?;
    }

    // VSWR circles (light blue-grey, before data traces).
    let vswr_style = RGBColor(190, 190, 220).stroke_width(1);
    for &vswr in &[1.5_f64, 2.0, 3.0, 5.0, 10.0] {
        let pts = smith_vswr_circle_pts(vswr, 128);
        chart
            .draw_series(LineSeries::new(pts, vswr_style))
            .map_err(map_render_err)?;
    }

    // Data traces with legend entries.
    for (i, trace) in traces.iter().enumerate() {
        let colour = SMITH_PALETTE[i % SMITH_PALETTE.len()];
        let traj: Vec<(f64, f64)> = trace.values.iter().map(|z| (z.re, z.im)).collect();
        let label = trace.label.clone();
        chart
            .draw_series(LineSeries::new(traj, colour.stroke_width(2)))
            .map_err(map_render_err)?
            .label(label)
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], colour.stroke_width(2))
            });
    }

    if traces.len() > 1 {
        chart
            .configure_series_labels()
            .background_style(WHITE.mix(0.8))
            .border_style(BLACK)
            .draw()
            .map_err(map_render_err)?;
    }

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

/// Compute a padded `(min, max)` y-range that covers all traces.
///
/// Iterates over every value in every `ys` sub-slice, finds the global
/// finite min/max, and applies a 5% pad (minimum 1 dB) so traces don't
/// sit on the frame. Falls back to `(default_min, default_max)` when no
/// finite values are found.
fn multi_trace_y_range(all_ys: &[Vec<f64>], default_min: f64, default_max: f64) -> (f64, f64) {
    let flat: Vec<f64> = all_ys.iter().flatten().copied().collect();
    let (y_min_raw, y_max_raw) = finite_range(&flat, default_min, default_max);
    let y_pad = ((y_max_raw - y_min_raw).abs() * 0.05).max(1.0);
    (y_min_raw - y_pad, y_max_raw + y_pad)
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

    // --- Multi-trace tests ---------------------------------------------------

    /// Build two synthetic traces: S11 (same as `synthetic_sweep`) and a
    /// flat S21 at −10 dB.
    fn two_trace_input() -> (Vec<f64>, Vec<SparamTrace>) {
        let (freq, s11) = synthetic_sweep();
        let n = freq.len();
        let s21: Vec<Complex64> = (0..n)
            .map(|_| Complex64::from_polar(0.316_227_766, 0.0)) // −10 dB
            .collect();
        let traces = vec![
            SparamTrace {
                label: "S11".to_string(),
                values: s11,
            },
            SparamTrace {
                label: "S21".to_string(),
                values: s21,
            },
        ];
        (freq, traces)
    }

    /// `plot_sparams_db` with two traces → SVG exists, is non-trivial, and
    /// contains both legend labels.
    #[test]
    fn test_plot_sparams_db_writes_svg_with_legend() {
        let (freq, traces) = two_trace_input();
        let tmp = NamedTempFile::with_suffix(".svg").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 800,
            height_px: 600,
            title: "S11+S21 dB overlay test".to_string(),
            format: PlotFormat::Svg,
        };
        plot_sparams_db(&freq, &traces, tmp.path(), &cfg).expect("plot_sparams_db");

        let body = fs::read_to_string(tmp.path()).expect("read svg");
        assert!(body.contains("<svg"), "SVG missing <svg tag");
        assert!(body.len() > 512, "SVG body too short: {} bytes", body.len());
        // Both legend labels should appear in the SVG text nodes.
        assert!(
            body.contains("S11"),
            "SVG should contain legend label 'S11'"
        );
        assert!(
            body.contains("S21"),
            "SVG should contain legend label 'S21'"
        );
    }

    /// `plot_sparams_db` with two traces → PNG exists and is non-trivial in
    /// size (bitmap content check — no string assertions for PNG).
    #[test]
    fn test_plot_sparams_db_writes_png() {
        let (freq, traces) = two_trace_input();
        let tmp = NamedTempFile::with_suffix(".png").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 800,
            height_px: 600,
            title: "S11+S21 dB PNG test".to_string(),
            format: PlotFormat::Png,
        };
        plot_sparams_db(&freq, &traces, tmp.path(), &cfg).expect("plot_sparams_db PNG");
        let len = fs::metadata(tmp.path()).expect("metadata").len();
        assert!(len > 1024, "PNG file is too small: {len} bytes");
    }

    /// `plot_sparams_phase` with two traces → SVG exists, is non-trivial, and
    /// contains both legend labels.
    #[test]
    fn test_plot_sparams_phase_writes_svg_with_legend() {
        let (freq, traces) = two_trace_input();
        let tmp = NamedTempFile::with_suffix(".svg").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 800,
            height_px: 600,
            title: "S11+S21 phase overlay test".to_string(),
            format: PlotFormat::Svg,
        };
        plot_sparams_phase(&freq, &traces, tmp.path(), &cfg).expect("plot_sparams_phase");

        let body = fs::read_to_string(tmp.path()).expect("read svg");
        assert!(body.contains("<svg"), "SVG missing <svg tag");
        assert!(body.len() > 512, "SVG body too short: {} bytes", body.len());
        assert!(
            body.contains("S11"),
            "SVG should contain legend label 'S11'"
        );
        assert!(
            body.contains("S21"),
            "SVG should contain legend label 'S21'"
        );
    }

    // --- Smith chart arc / multi-trace tests --------------------------------

    /// `smith_r_circle_pts(1.0, 64)` returns 65 points forming a closed loop
    /// all lying on the circle centred at (0.5, 0.0) with radius 0.5.
    #[test]
    fn smith_r_circle_pts_correct_geometry() {
        let pts = smith_r_circle_pts(1.0, 64);
        assert_eq!(pts.len(), 65, "expected n+1 = 65 points");

        // First and last must be identical (closed loop).
        assert!(
            (pts[0].0 - pts[64].0).abs() < 1e-12 && (pts[0].1 - pts[64].1).abs() < 1e-12,
            "first and last point must be equal: {:?} vs {:?}",
            pts[0],
            pts[64]
        );

        // All points must lie on the circle: centre=(0.5,0), radius=0.5.
        for (re, im) in &pts {
            let dist = ((re - 0.5).powi(2) + im.powi(2)).sqrt();
            assert!(
                (dist - 0.5).abs() < 1e-9,
                "point ({re}, {im}) is not on r=1 circle (dist={dist})"
            );
        }
    }

    /// Every point of `smith_x_arc_pts(1.0, 256)` lies on the circle with
    /// centre (1.0, 1.0) and radius 1.0 — verifies the constant-X arc formula.
    #[test]
    fn smith_x_arc_pts_on_circle() {
        let pts = smith_x_arc_pts(1.0, 256);
        assert!(!pts.is_empty(), "expected non-empty arc for x=1.0");
        for (re, im) in &pts {
            let dist = ((re - 1.0).powi(2) + (im - 1.0).powi(2)).sqrt();
            assert!(
                (dist - 1.0).abs() < 1e-9,
                "point ({re}, {im}) not on circle: dist={dist}"
            );
        }
    }

    /// All points returned by `smith_x_arc_pts(0.2, 512)` must lie inside or
    /// on the unit disk.
    #[test]
    fn smith_x_arc_pts_inside_unit_disk() {
        let pts = smith_x_arc_pts(0.2, 512);
        assert!(!pts.is_empty(), "expected some points inside the unit disk");
        for (re, im) in &pts {
            let r2 = re * re + im * im;
            assert!(
                r2 <= 1.0 + 1e-9,
                "point ({re}, {im}) is outside the unit disk (r²={r2})"
            );
        }
    }

    // --- smith_vswr_circle_pts tests ----------------------------------------

    /// All points of `smith_vswr_circle_pts(2.0, 64)` must lie at radius
    /// ρ = (2−1)/(2+1) = 1/3.
    #[test]
    fn smith_vswr_circle_pts_radius_vswr_2() {
        let rho = 1.0_f64 / 3.0;
        let pts = smith_vswr_circle_pts(2.0, 64);
        for (re, im) in &pts {
            let r = (re * re + im * im).sqrt();
            assert!(
                (r - rho).abs() < 1e-12,
                "point ({re}, {im}) has radius {r}, expected {rho}"
            );
        }
    }

    /// `smith_vswr_circle_pts(3.0, 64)` returns `n+1 = 65` points and the
    /// first and last points are equal within 1e-12 (closed loop).
    #[test]
    fn smith_vswr_circle_pts_closed() {
        let n = 64usize;
        let pts = smith_vswr_circle_pts(3.0, n);
        assert_eq!(pts.len(), n + 1, "expected n+1 = {} points", n + 1);
        assert!(
            (pts[0].0 - pts[n].0).abs() < 1e-12 && (pts[0].1 - pts[n].1).abs() < 1e-12,
            "circle must be closed: first={:?}, last={:?}",
            pts[0],
            pts[n]
        );
    }

    /// `plot_smith_chart_multi` with two traces writes a non-empty SVG file.
    #[test]
    fn plot_smith_chart_multi_two_traces_writes_svg() {
        let traces = vec![
            SmithTrace {
                label: "S11".to_string(),
                values: vec![
                    Complex64::new(0.5, 0.3),
                    Complex64::new(0.2, -0.1),
                    Complex64::new(-0.1, 0.0),
                    Complex64::new(0.0, 0.5),
                ],
            },
            SmithTrace {
                label: "S22".to_string(),
                values: vec![
                    Complex64::new(-0.3, 0.2),
                    Complex64::new(0.1, 0.4),
                    Complex64::new(0.4, -0.2),
                    Complex64::new(0.6, 0.0),
                ],
            },
        ];
        let tmp = NamedTempFile::with_suffix(".svg").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 600,
            height_px: 600,
            title: "Smith multi-trace test".to_string(),
            format: PlotFormat::Svg,
        };
        plot_smith_chart_multi(&traces, tmp.path(), &cfg).expect("plot_smith_chart_multi");
        let body = fs::read_to_string(tmp.path()).expect("read svg");
        assert!(body.contains("<svg"), "SVG missing <svg tag");
        assert!(body.len() > 256, "SVG body too short: {} bytes", body.len());
    }

    /// `plot_smith_chart_multi` renders VSWR circles unconditionally; verify
    /// the output file is non-empty with a single trace.
    #[test]
    fn plot_smith_chart_multi_with_vswr_circles_renders() {
        let traces = vec![SmithTrace {
            label: "S11".to_string(),
            values: vec![
                Complex64::new(0.3, 0.1),
                Complex64::new(0.0, 0.0),
                Complex64::new(-0.2, 0.3),
            ],
        }];
        let tmp = NamedTempFile::with_suffix(".svg").expect("tempfile");
        let cfg = PlotConfig {
            width_px: 600,
            height_px: 600,
            title: "VSWR circles smoke test".to_string(),
            format: PlotFormat::Svg,
        };
        plot_smith_chart_multi(&traces, tmp.path(), &cfg)
            .expect("plot_smith_chart_multi with vswr circles");
        let body = fs::read_to_string(tmp.path()).expect("read svg");
        assert!(body.contains("<svg"), "SVG missing <svg tag");
        assert!(body.len() > 256, "SVG body too short: {} bytes", body.len());
    }
}
