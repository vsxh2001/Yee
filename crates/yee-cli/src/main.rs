//! # yee CLI
//!
//! Top-level command-line tool for the Yee electromagnetic-simulation studio.
//!
//! Phase 0 wires the following subcommands:
//!
//! - `yee validate <mom|fdtd|all> [--json]` — runs the `yee-validation`
//!   aggregator (real mom-001 NEC-4 gate; Phase-deferred FDTD cases report
//!   `Skipped`). Filters by target prefix (`mom-*`, FDTD-family) and
//!   exits 1 if any included case failed.
//! - `yee mesh <path>` — constructs a `yee_mesh::Session`. Without the
//!   `gmsh` feature this exits with code 2 and a guidance message.
//! - `yee export <input> --format <touchstone|hdf5> <output>` — reads/writes
//!   Touchstone via `yee-io`. `hdf5` is not yet enabled and exits with code 2.
//! - `yee plot <input> --kind <db|smith|phase> --output <out>` — reads a
//!   Touchstone file and emits a PNG/SVG plot via `yee-plotters` (the format
//!   is picked from the output file extension).
//! - `yee completions <shell>` — emits a shell completion script to stdout
//!   (`bash`, `zsh`, `fish`).
//! - `yee run <project>` — Phase 0 stub from the scaffold.

use std::io;
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{Shell, generate};

mod plot;

#[derive(Parser, Debug)]
#[command(
    name = "yee",
    version,
    about = "Yee — GPU-accelerated electromagnetic simulation studio"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Run validation cases against published benchmarks.
    ///
    /// Invokes the `yee-validation` aggregator and filters by `target`
    /// prefix. `mom` includes `mom-*` cases; `fdtd` includes
    /// `cpml-*` / `ntff-*` / `dispersive-*` / `fdtd-*`; `all` runs
    /// the full report. Default output is a human-readable table; pass
    /// `--json` to emit the serialized [`yee_validation::Report`].
    ///
    /// Exit code is non-zero iff any *included* case has status
    /// `Failed`. `Skipped` cases (Phase-deferred placeholders) do not
    /// fail the run.
    Validate {
        /// Which solver to validate.
        #[arg(value_enum, default_value_t = ValidateTarget::All)]
        target: ValidateTarget,
        /// Emit JSON report to stdout (default: human-readable table).
        #[arg(long)]
        json: bool,
    },
    /// Mesh a geometry file via Gmsh.
    Mesh {
        /// Input geometry path (.step / .iges / .kicad_pcb).
        input: PathBuf,
    },
    /// Run a simulation defined by a project file (Phase 0 stub).
    Run {
        /// Path to the project TOML.
        project: PathBuf,
    },
    /// Export results to Touchstone or HDF5.
    Export {
        /// Path to the input results file (e.g. a Touchstone `.s*p` file).
        input: PathBuf,
        /// Output format.
        #[arg(long, value_enum, default_value_t = ExportFormat::Touchstone)]
        format: ExportFormat,
        /// Output file path.
        output: PathBuf,
    },
    /// Plot S-parameters from a Touchstone file.
    ///
    /// The output format (PNG vs SVG) is chosen from the `--output` file
    /// extension; `.png` and `.svg` are accepted (no extension defaults to
    /// PNG).
    Plot {
        /// Input Touchstone path (.s1p, .s2p, etc.).
        input: PathBuf,
        /// What to plot.
        #[arg(long, value_enum, default_value_t = PlotKind::Db)]
        kind: PlotKind,
        /// Output file path; extension picks PNG vs SVG.
        #[arg(long, short)]
        output: PathBuf,
        /// Width in pixels.
        #[arg(long, default_value_t = 800)]
        width: u32,
        /// Height in pixels.
        #[arg(long, default_value_t = 600)]
        height: u32,
        /// Plot title; defaults to the input file stem.
        #[arg(long)]
        title: Option<String>,
        /// Which port (index into the S-matrix, 0-based). Default 0 (S₁₁).
        #[arg(long, default_value_t = 0)]
        port: usize,
    },
    /// Generate a shell completion script on stdout.
    ///
    /// Pre-generated scripts live in `crates/yee-cli/completions/`.
    Completions {
        /// Target shell (`bash`, `zsh`, `fish`, ...).
        shell: Shell,
    },
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum ExportFormat {
    /// Touchstone v1.1 (.s1p/.s2p/.s3p/.s4p).
    Touchstone,
    /// HDF5 (not yet enabled).
    Hdf5,
}

/// What `yee plot` should draw from the S-parameter sweep.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum PlotKind {
    /// `|S|` in dB vs frequency.
    Db,
    /// `S` on the Smith chart.
    Smith,
    /// `phase(S)` in degrees vs frequency.
    Phase,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum ValidateTarget {
    /// Method-of-moments planar solver.
    Mom,
    /// Finite-difference time-domain solver (Phase 2).
    Fdtd,
    /// Run every available solver's validation suite.
    All,
}

fn main() -> ExitCode {
    tracing_subscriber::fmt().with_target(false).init();

    let cli = Cli::parse();
    match run(cli) {
        Ok(code) => code,
        Err(err) => {
            eprintln!("error: {err:#}");
            ExitCode::from(1)
        }
    }
}

fn run(cli: Cli) -> Result<ExitCode> {
    match cli.command {
        Command::Validate { target, json } => run_validate(target, json),
        Command::Mesh { input } => Ok(run_mesh(&input)),
        Command::Run { project } => {
            println!("yee run {} — Phase 0 stub.", project.display());
            Ok(ExitCode::SUCCESS)
        }
        Command::Export {
            input,
            format,
            output,
        } => run_export(&input, format, &output),
        Command::Plot {
            input,
            kind,
            output,
            width,
            height,
            title,
            port,
        } => plot::run_plot(plot::PlotArgs {
            input,
            kind,
            output,
            width,
            height,
            title,
            port,
        }),
        Command::Completions { shell } => {
            let mut cmd = Cli::command();
            let bin_name = cmd.get_name().to_string();
            generate(shell, &mut cmd, bin_name, &mut io::stdout());
            Ok(ExitCode::SUCCESS)
        }
    }
}

fn run_export(
    input: &std::path::Path,
    format: ExportFormat,
    output: &std::path::Path,
) -> Result<ExitCode> {
    match format {
        ExportFormat::Touchstone => {
            let file = yee_io::touchstone::read(input)
                .map_err(|e| anyhow::anyhow!("failed to read Touchstone file: {e}"))?;
            yee_io::touchstone::write(output, &file)
                .map_err(|e| anyhow::anyhow!("failed to write Touchstone file: {e}"))?;
            Ok(ExitCode::SUCCESS)
        }
        ExportFormat::Hdf5 => {
            // Diagnostic messages go to stderr — keeps stdout clean for any
            // future success output (e.g. a confirmation path) and matches the
            // convention used by `run_mesh` for `NotEnabled`.
            eprintln!("hdf5 not yet enabled");
            Ok(ExitCode::from(2))
        }
    }
}

/// Construct a [`yee_mesh::Session`] for `input`. Without the `gmsh` feature
/// the underlying crate returns [`yee_mesh::Error::NotEnabled`]; we surface
/// this to the user with exit code 2.
fn run_mesh(_input: &std::path::Path) -> ExitCode {
    match yee_mesh::Session::new() {
        Ok(_session) => {
            // Phase 1 wires `import_step` and `mesh` against the real Gmsh
            // FFI; for Phase 0 simply constructing the session is the
            // smoke-test contract.
            println!("yee mesh: session opened (Phase 0 stub).");
            ExitCode::SUCCESS
        }
        Err(yee_mesh::Error::NotEnabled) => {
            eprintln!("mesh feature not enabled; rebuild with --features gmsh");
            ExitCode::from(2)
        }
        Err(err) => {
            eprintln!("mesh error: {err}");
            ExitCode::from(1)
        }
    }
}

/// Decide whether a case id belongs to the selected solver target.
///
/// Lives at the top of the validate-handler module so the prefix list is
/// in one place and matches the brief: `mom-*` for [`ValidateTarget::Mom`];
/// `cpml-*` / `ntff-*` / `dispersive-*` / `fdtd-*` for
/// [`ValidateTarget::Fdtd`]; everything for [`ValidateTarget::All`].
fn case_matches_target(case_id: &str, target: ValidateTarget) -> bool {
    match target {
        ValidateTarget::All => true,
        ValidateTarget::Mom => case_id.starts_with("mom-"),
        ValidateTarget::Fdtd => {
            case_id.starts_with("cpml-")
                || case_id.starts_with("ntff-")
                || case_id.starts_with("dispersive-")
                || case_id.starts_with("fdtd-")
        }
    }
}

/// Run the validation aggregator, filter by target, and print to stdout.
///
/// Returns [`ExitCode::FAILURE`] iff any *included* case has status
/// [`yee_validation::CaseStatus::Failed`]. `Skipped` cases (Phase-deferred
/// placeholders, see CLAUDE.md §10) never fail the run.
fn run_validate(target: ValidateTarget, json: bool) -> Result<ExitCode> {
    let mut report = yee_validation::Report::run_all();
    report.cases.retain(|c| case_matches_target(&c.id, target));

    if json {
        let s = report
            .to_json()
            .map_err(|e| anyhow::anyhow!("serializing report: {e}"))?;
        println!("{s}");
    } else {
        print_human_report(&report);
    }

    if report.has_failures() {
        Ok(ExitCode::FAILURE)
    } else {
        Ok(ExitCode::SUCCESS)
    }
}

/// Render a [`yee_validation::Report`] as a 4-column fixed-width table.
///
/// Columns: `CASE` (left), `STATUS` (left), `WALL_TIME (s)` (right),
/// `NOTES` (left). Widths are computed from the data so the right-aligned
/// wall-time column lines up regardless of how many cases passed vs
/// skipped. Wall-time for `Skipped` rows is rendered as an em-dash to
/// avoid implying meaningful timing on a no-op.
fn print_human_report(report: &yee_validation::Report) {
    use yee_validation::CaseStatus;

    const H_CASE: &str = "CASE";
    const H_STATUS: &str = "STATUS";
    const H_TIME: &str = "WALL_TIME (s)";
    const H_NOTES: &str = "NOTES";

    // Pre-format rows so column widths can be derived from the
    // final cell strings (status names, em-dashes, etc.).
    struct Row {
        id: String,
        status: String,
        time: String,
        notes: String,
    }
    let rows: Vec<Row> = report
        .cases
        .iter()
        .map(|c| Row {
            id: c.id.clone(),
            status: match c.status {
                CaseStatus::Passed => "Passed".to_string(),
                CaseStatus::Failed => "Failed".to_string(),
                CaseStatus::Skipped => "Skipped".to_string(),
            },
            time: match c.status {
                CaseStatus::Skipped => "—".to_string(),
                _ => format!("{:.1}", c.wall_time_seconds),
            },
            notes: c.notes.clone(),
        })
        .collect();

    let w_case = rows
        .iter()
        .map(|r| r.id.len())
        .max()
        .unwrap_or(0)
        .max(H_CASE.len());
    let w_status = rows
        .iter()
        .map(|r| r.status.len())
        .max()
        .unwrap_or(0)
        .max(H_STATUS.len());
    let w_time = rows
        .iter()
        .map(|r| r.time.chars().count())
        .max()
        .unwrap_or(0)
        .max(H_TIME.len());

    // Header
    println!("{H_CASE:<w_case$}  {H_STATUS:<w_status$}  {H_TIME:>w_time$}  {H_NOTES}");
    // The time column is right-aligned, which means simple padding via
    // {:>w_time$} on a multi-byte em-dash uses byte width — not what we
    // want. Pad with spaces explicitly using chars().count().
    for r in &rows {
        let time_pad = w_time.saturating_sub(r.time.chars().count());
        let pad = " ".repeat(time_pad);
        let id = &r.id;
        let status = &r.status;
        let time = &r.time;
        let notes = &r.notes;
        println!("{id:<w_case$}  {status:<w_status$}  {pad}{time}  {notes}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn target_filter_mom_only_matches_mom_prefix() {
        assert!(case_matches_target("mom-001", ValidateTarget::Mom));
        assert!(case_matches_target("mom-002", ValidateTarget::Mom));
        assert!(!case_matches_target("cpml-001", ValidateTarget::Mom));
        assert!(!case_matches_target("ntff-001", ValidateTarget::Mom));
    }

    #[test]
    fn target_filter_fdtd_matches_fdtd_families() {
        assert!(case_matches_target("cpml-001", ValidateTarget::Fdtd));
        assert!(case_matches_target("ntff-001", ValidateTarget::Fdtd));
        assert!(case_matches_target("dispersive-001", ValidateTarget::Fdtd));
        assert!(case_matches_target("fdtd-cavity", ValidateTarget::Fdtd));
        assert!(!case_matches_target("mom-001", ValidateTarget::Fdtd));
    }

    #[test]
    fn target_filter_all_matches_everything() {
        assert!(case_matches_target("mom-001", ValidateTarget::All));
        assert!(case_matches_target("cpml-001", ValidateTarget::All));
        assert!(case_matches_target("anything-else", ValidateTarget::All));
    }
}
