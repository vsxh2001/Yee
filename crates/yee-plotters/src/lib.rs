//! Static PNG/SVG export of S-parameter plots via the [`plotters`] crate.
//!
//! Phase 1.plotting.0 walking-skeleton: this commit lands the public
//! surface (types + function signatures). Real implementations follow.

use std::path::Path;

use num_complex::Complex64;

/// Minimum dB value we will report. `|S₁₁| = 0` would otherwise become
/// `-∞ dB`; we clamp here instead.
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
    /// Filesystem or path-level error.
    #[error("io error: {0}")]
    Io(String),
    /// Rendering error from a `plotters` backend.
    #[error("render error: {0}")]
    Render(String),
}

/// Plot `|S₁₁|` in dB vs. frequency. Skeleton — implementation in a later commit.
pub fn plot_s11_db(
    _freq_hz: &[f64],
    _s11: &[Complex64],
    _out_path: &Path,
    _config: &PlotConfig,
) -> Result<(), Error> {
    Err(Error::Render("plot_s11_db: not yet implemented".into()))
}

/// Plot the phase of `S₁₁` in degrees vs. frequency. Skeleton.
pub fn plot_s11_phase(
    _freq_hz: &[f64],
    _s11: &[Complex64],
    _out_path: &Path,
    _config: &PlotConfig,
) -> Result<(), Error> {
    Err(Error::Render("plot_s11_phase: not yet implemented".into()))
}

/// Plot `S₁₁` on the complex unit disk (Smith-chart style). Skeleton.
pub fn plot_smith_chart(
    _s11: &[Complex64],
    _out_path: &Path,
    _config: &PlotConfig,
) -> Result<(), Error> {
    Err(Error::Render("plot_smith_chart: not yet implemented".into()))
}
