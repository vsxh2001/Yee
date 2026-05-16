//! `yee plot` handler.
//!
//! Reads a Touchstone file via `yee-io` and emits a PNG or SVG plot through
//! `yee-plotters`. The output image format is chosen from the file extension:
//!
//! | extension | backend            |
//! |-----------|--------------------|
//! | `.png`    | `BitMapBackend`    |
//! | `.svg`    | `SVGBackend`       |
//! | (none)    | PNG                |
//! | other     | error (exit 1)     |
//!
//! The plot kind is selected by `--format` (legacy alias `--kind`).
//! `db`, `smith`, and `phase` each emit a single file at `--output`; `both`
//! emits two files derived from `--output` by inserting `-db` / `-smith`
//! between the stem and the extension. For multi-port Touchstone files,
//! `--port` selects the diagonal entry `S[port][port]` to plot.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use num_complex::Complex64;

use crate::PlotKind;

/// Arguments for [`run_plot`]. Kept as its own struct so the dispatch site
/// in `main.rs` doesn't have to thread eight positional parameters.
pub(crate) struct PlotArgs {
    pub input: PathBuf,
    pub kind: PlotKind,
    pub output: PathBuf,
    pub width: u32,
    pub height: u32,
    pub title: Option<String>,
    pub port: usize,
}

/// Dispatch a `yee plot` invocation. Returns `Ok(ExitCode::SUCCESS)` on a
/// successful write, surfaces an `anyhow::Error` for IO/parse/render
/// failures (the binary's `main` turns that into exit code 1).
pub(crate) fn run_plot(args: PlotArgs) -> Result<ExitCode> {
    let file = yee_io::touchstone::read(&args.input)
        .with_context(|| format!("touchstone read: {}", args.input.display()))?;

    let n = file.n_ports;
    if args.port >= n {
        anyhow::bail!("port index {} out of range (file has {n} ports)", args.port);
    }

    // Diagonal entry S[port][port]. `file.data[k]` is row-major n×n.
    let idx = args.port * n + args.port;
    let s: Vec<Complex64> = file.data.iter().map(|row| row[idx]).collect();

    let format = pick_format(&args.output)?;

    let title = args.title.unwrap_or_else(|| {
        args.input
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("S-parameters")
            .to_string()
    });
    let config = yee_plotters::PlotConfig {
        width_px: args.width,
        height_px: args.height,
        title,
        format,
    };

    match args.kind {
        PlotKind::Db => {
            yee_plotters::plot_s11_db(&file.freq_hz, &s, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Smith => {
            yee_plotters::plot_smith_chart(&s, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Phase => {
            yee_plotters::plot_s11_phase(&file.freq_hz, &s, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Both => {
            let db_path = suffixed_path(&args.output, "-db");
            let smith_path = suffixed_path(&args.output, "-smith");
            yee_plotters::plot_s11_db(&file.freq_hz, &s, &db_path, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            yee_plotters::plot_smith_chart(&s, &smith_path, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!(
                "yee plot: wrote {} and {}",
                db_path.display(),
                smith_path.display()
            );
        }
    }

    Ok(ExitCode::SUCCESS)
}

/// Insert `suffix` between an output path's stem and its extension.
///
/// `foo/out.png` + `-db` becomes `foo/out-db.png`. Paths without an
/// extension append the suffix to the file name (`foo/out` → `foo/out-db`).
/// Paths with no file name (a bare `/`) fall through to a literal join,
/// which is consistent with how the rest of the handler treats malformed
/// outputs — the underlying plotter call will surface the IO error.
fn suffixed_path(output: &Path, suffix: &str) -> PathBuf {
    let parent = output.parent().unwrap_or_else(|| Path::new(""));
    let stem = output
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let new_name = match output.extension().and_then(|s| s.to_str()) {
        Some(ext) => format!("{stem}{suffix}.{ext}"),
        None => format!("{stem}{suffix}"),
    };
    parent.join(new_name)
}

/// Map an output path's extension to a [`yee_plotters::PlotFormat`].
///
/// `.png` and missing extensions both pick `Png`; `.svg` picks `Svg`;
/// anything else is rejected with a friendly error.
fn pick_format(output: &Path) -> Result<yee_plotters::PlotFormat> {
    match output.extension().and_then(|s| s.to_str()) {
        Some("svg") | Some("SVG") => Ok(yee_plotters::PlotFormat::Svg),
        Some("png") | Some("PNG") | None => Ok(yee_plotters::PlotFormat::Png),
        Some(other) => Err(anyhow::anyhow!(
            "unsupported output extension '.{other}' — use .png or .svg"
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pick_format_png_for_png_extension() {
        let p = PathBuf::from("/tmp/x.png");
        assert_eq!(pick_format(&p).unwrap(), yee_plotters::PlotFormat::Png);
    }

    #[test]
    fn pick_format_svg_for_svg_extension() {
        let p = PathBuf::from("/tmp/x.svg");
        assert_eq!(pick_format(&p).unwrap(), yee_plotters::PlotFormat::Svg);
    }

    #[test]
    fn pick_format_png_for_missing_extension() {
        let p = PathBuf::from("/tmp/noext");
        assert_eq!(pick_format(&p).unwrap(), yee_plotters::PlotFormat::Png);
    }

    #[test]
    fn pick_format_rejects_unknown_extension() {
        let p = PathBuf::from("/tmp/x.jpg");
        let err = pick_format(&p).unwrap_err().to_string();
        assert!(err.contains("jpg"), "error should mention extension: {err}");
    }

    #[test]
    fn suffixed_path_inserts_before_extension() {
        let p = PathBuf::from("/tmp/out.png");
        assert_eq!(suffixed_path(&p, "-db"), PathBuf::from("/tmp/out-db.png"));
        assert_eq!(
            suffixed_path(&p, "-smith"),
            PathBuf::from("/tmp/out-smith.png")
        );
    }

    #[test]
    fn suffixed_path_handles_missing_extension() {
        let p = PathBuf::from("/tmp/out");
        assert_eq!(suffixed_path(&p, "-db"), PathBuf::from("/tmp/out-db"));
    }

    #[test]
    fn suffixed_path_preserves_parent_directory() {
        let p = PathBuf::from("a/b/c.svg");
        assert_eq!(suffixed_path(&p, "-db"), PathBuf::from("a/b/c-db.svg"));
    }
}
