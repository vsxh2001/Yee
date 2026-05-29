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
//! between the stem and the extension.
//!
//! ## Port selection
//!
//! **Single-trace (default):** `--port <i>` selects the diagonal entry
//! `S[i][i]` (0-based). Default is 0 (`S₁₁`). The existing single-trace
//! entry points (`plot_s11_db`, `plot_s11_phase`, `plot_smith_chart`) are used
//! in this path — their signatures are unchanged.
//!
//! **Multi-trace overlay:** pass one or more `--entry <ij>` flags (e.g.
//! `--entry 11 --entry 21`) to select off-diagonal S-matrix entries and overlay
//! them in one plot. Alternatively pass `--all` to overlay every entry of a
//! multi-port file. When multi-trace mode is active, `plot_sparams_db` /
//! `plot_sparams_phase` / `plot_smith_chart_multi` are called as appropriate.
//! Out-of-range indices are rejected with a clean error.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use num_complex::Complex64;

use crate::PlotKind;

/// Arguments for [`run_plot`]. Kept as its own struct so the dispatch site
/// in `main.rs` doesn't have to thread ten positional parameters.
pub(crate) struct PlotArgs {
    /// Path to the input Touchstone file.
    pub input: PathBuf,
    /// What kind of plot to produce (`db`, `phase`, `smith`, `both`).
    pub kind: PlotKind,
    /// Output image path; extension selects PNG vs SVG.
    pub output: PathBuf,
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
    /// Optional plot title; defaults to the input file stem.
    pub title: Option<String>,
    /// Single diagonal port index (0-based) for the default single-trace path.
    /// Ignored when `entries` is non-empty or `all_entries` is true.
    pub port: usize,
    /// S-matrix entries to overlay, each expressed as a two-digit string like
    /// `"11"`, `"21"`, `"12"`, etc. (1-based, matching the Touchstone
    /// convention). When non-empty, the multi-trace overlay path is taken
    /// regardless of `port`.
    pub entries: Vec<String>,
    /// When true, overlay every entry of the S-matrix (all `n×n` entries).
    /// Takes precedence over individual `entries` items.
    pub all_entries: bool,
}

/// Dispatch a `yee plot` invocation. Returns `Ok(ExitCode::SUCCESS)` on a
/// successful write, surfaces an `anyhow::Error` for IO/parse/render
/// failures (the binary's `main` turns that into exit code 1).
pub(crate) fn run_plot(args: PlotArgs) -> Result<ExitCode> {
    let file = yee_io::touchstone::read(&args.input)
        .with_context(|| format!("touchstone read: {}", args.input.display()))?;

    let n = file.n_ports;

    // Decide whether we are in multi-trace mode.
    let multi_mode = args.all_entries || !args.entries.is_empty();

    if multi_mode {
        run_multi_trace(&args, &file, n)
    } else {
        run_single_trace(&args, &file, n)
    }
}

/// Execute the single-trace path (the pre-existing behaviour, unchanged).
///
/// Uses the `--port` index to select the diagonal entry `S[port][port]` and
/// calls the single-trace plotter functions (`plot_s11_db`, `plot_s11_phase`,
/// `plot_smith_chart`, or both).
fn run_single_trace(
    args: &PlotArgs,
    file: &yee_io::touchstone::File,
    n: usize,
) -> Result<ExitCode> {
    if args.port >= n {
        anyhow::bail!("port index {} out of range (file has {n} ports)", args.port);
    }

    // Diagonal entry S[port][port]. `file.data[k]` is row-major n×n.
    let idx = args.port * n + args.port;
    let s: Vec<Complex64> = file.data.iter().map(|row| row[idx]).collect();

    let format = pick_format(&args.output)?;
    let title = resolve_title(args);
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

/// Execute the multi-trace overlay path.
///
/// Resolves the requested S-matrix entries to (row, col) index pairs, extracts
/// each trace from the row-major `file.data`, labels them `S<row><col>` (e.g.
/// `S11`, `S21`), and calls the appropriate multi-trace plotter:
/// - `db`    → [`yee_plotters::plot_sparams_db`]
/// - `phase` → [`yee_plotters::plot_sparams_phase`]
/// - `smith` → [`yee_plotters::plot_smith_chart_multi`]
/// - `both`  → dB file (`out-db.<ext>`) + Smith file (`out-smith.<ext>`)
fn run_multi_trace(args: &PlotArgs, file: &yee_io::touchstone::File, n: usize) -> Result<ExitCode> {
    // Build the list of (row, col) pairs (0-based internally; 1-based on CLI).
    let pairs: Vec<(usize, usize)> = if args.all_entries {
        (0..n).flat_map(|r| (0..n).map(move |c| (r, c))).collect()
    } else {
        // Parse each entry string like "11", "21", "12" as (row-1, col-1).
        let mut parsed = Vec::with_capacity(args.entries.len());
        for entry in &args.entries {
            let (r1, c1) = parse_entry(entry)?;
            // Validate against n_ports.
            if r1 == 0 || c1 == 0 {
                anyhow::bail!(
                    "S-matrix entry '{}' uses 0-based indices; \
                     use 1-based indices matching the Touchstone convention (e.g. 11, 21)",
                    entry
                );
            }
            let r = r1 - 1;
            let c = c1 - 1;
            if r >= n || c >= n {
                anyhow::bail!(
                    "S-matrix entry '{}' out of range (file has {n} ports, \
                     so valid indices are 1..={n})",
                    entry
                );
            }
            parsed.push((r, c));
        }
        parsed
    };

    if pairs.is_empty() {
        anyhow::bail!("no S-matrix entries selected; pass --entry <ij> or --all");
    }

    // Extract S-parameter traces (label + raw complex values). Build as
    // SparamTrace; convert field-for-field to SmithTrace when needed.
    let sparam_traces: Vec<yee_plotters::SparamTrace> = pairs
        .iter()
        .map(|&(r, c)| {
            let flat_idx = r * n + c;
            let values: Vec<Complex64> = file.data.iter().map(|row| row[flat_idx]).collect();
            // Label uses 1-based, matches Touchstone notation.
            let label = format!("S{}{}", r + 1, c + 1);
            yee_plotters::SparamTrace { label, values }
        })
        .collect();

    let format = pick_format(&args.output)?;
    let title = resolve_title(args);
    let config = yee_plotters::PlotConfig {
        width_px: args.width,
        height_px: args.height,
        title,
        format,
    };

    match args.kind {
        PlotKind::Db => {
            yee_plotters::plot_sparams_db(&file.freq_hz, &sparam_traces, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Phase => {
            yee_plotters::plot_sparams_phase(&file.freq_hz, &sparam_traces, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Smith => {
            let smith_traces: Vec<yee_plotters::SmithTrace> = sparam_traces
                .into_iter()
                .map(|t| yee_plotters::SmithTrace {
                    label: t.label,
                    values: t.values,
                })
                .collect();
            yee_plotters::plot_smith_chart_multi(&smith_traces, &args.output, &config)
                .map_err(|e| anyhow::anyhow!("plot: {e}"))?;
            eprintln!("yee plot: wrote {}", args.output.display());
        }
        PlotKind::Both => {
            let db_path = suffixed_path(&args.output, "-db");
            let smith_path = suffixed_path(&args.output, "-smith");
            let smith_traces: Vec<yee_plotters::SmithTrace> = sparam_traces
                .iter()
                .map(|t| yee_plotters::SmithTrace {
                    label: t.label.clone(),
                    values: t.values.clone(),
                })
                .collect();
            yee_plotters::plot_sparams_db(&file.freq_hz, &sparam_traces, &db_path, &config)
                .map_err(|e| anyhow::anyhow!("plot (db): {e}"))?;
            yee_plotters::plot_smith_chart_multi(&smith_traces, &smith_path, &config)
                .map_err(|e| anyhow::anyhow!("plot (smith): {e}"))?;
            eprintln!(
                "yee plot: wrote {} and {}",
                db_path.display(),
                smith_path.display()
            );
        }
    }

    Ok(ExitCode::SUCCESS)
}

/// Parse an entry string like `"11"`, `"21"`, `"12"`, `"22"` into
/// `(row_1based, col_1based)`.
///
/// Accepts exactly two decimal digits (no spaces, no separator). Returns an
/// error for strings that don't match this shape.
fn parse_entry(entry: &str) -> Result<(usize, usize)> {
    let trimmed = entry.trim();
    if trimmed.len() != 2 {
        anyhow::bail!(
            "invalid entry '{}': expected exactly two digits, e.g. '11', '21'",
            entry
        );
    }
    let row_ch = trimmed.as_bytes()[0];
    let col_ch = trimmed.as_bytes()[1];
    if !row_ch.is_ascii_digit() || !col_ch.is_ascii_digit() {
        anyhow::bail!(
            "invalid entry '{}': both characters must be ASCII digits",
            entry
        );
    }
    let r = (row_ch - b'0') as usize;
    let c = (col_ch - b'0') as usize;
    Ok((r, c))
}

/// Resolve the plot title: use the explicit `--title` if provided, otherwise
/// fall back to the input file stem.
fn resolve_title(args: &PlotArgs) -> String {
    args.title.clone().unwrap_or_else(|| {
        args.input
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("S-parameters")
            .to_string()
    })
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

    #[test]
    fn parse_entry_ok_11() {
        assert_eq!(parse_entry("11").unwrap(), (1, 1));
    }

    #[test]
    fn parse_entry_ok_21() {
        assert_eq!(parse_entry("21").unwrap(), (2, 1));
    }

    #[test]
    fn parse_entry_ok_12() {
        assert_eq!(parse_entry("12").unwrap(), (1, 2));
    }

    #[test]
    fn parse_entry_rejects_too_long() {
        assert!(parse_entry("211").is_err());
    }

    #[test]
    fn parse_entry_rejects_non_digits() {
        assert!(parse_entry("ab").is_err());
    }

    #[test]
    fn parse_entry_rejects_one_char() {
        assert!(parse_entry("1").is_err());
    }
}
