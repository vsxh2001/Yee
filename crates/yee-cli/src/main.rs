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
//! - `yee bench <target> [-- <criterion-args>]` — shells out to
//!   `cargo bench -p yee-bench` for the chosen benchmark binary (or all of
//!   them with `all`). Stdout/stderr are inherited so criterion's live
//!   progress output is preserved.
//! - `yee design <prompt> --output <p.toml> [--offline] [--model <id>]` —
//!   Phase 3.nl.0 natural-language → `yee.toml` surface. Default path is the
//!   deterministic offline parser (`yee_design::offline::parse`) when
//!   `--offline` is set or `ANTHROPIC_API_KEY` is unset; the live LLM path
//!   is provided by the `yee-py` Python sidecar (`yee.design.from_prompt_llm`)
//!   and PyO3 in-process embedding from this binary is deferred to
//!   Phase 3.nl.0.1 per the R5 escape hatch.
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
    /// The output image format (PNG vs SVG) is chosen from the `--output` file
    /// extension; `.png` and `.svg` are accepted (no extension defaults to
    /// PNG). The plot kind is selected with `--format` (or its legacy alias
    /// `--kind`) — `db`, `smith`, `phase`, or `both` (emits two files with
    /// `-db` / `-smith` suffixes inserted before the extension).
    Plot {
        /// Input Touchstone path (.s1p, .s2p, etc.).
        input: PathBuf,
        /// What to plot. `--kind` is accepted as a legacy alias.
        #[arg(long, visible_alias = "kind", value_enum, default_value_t = PlotKind::Db)]
        format: PlotKind,
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
    /// Run a yee-bench criterion benchmark.
    ///
    /// Shells out to `cargo bench -p yee-bench [--bench <name>]` and
    /// inherits stdout/stderr so criterion's progress and statistics
    /// display live. Extra arguments after `--` are forwarded as
    /// criterion CLI flags (e.g. `yee bench bo -- --warm-up-time 1`).
    Bench {
        /// Which bench to run.
        #[arg(value_enum)]
        target: BenchTarget,
        /// Pass through extra `--bench` args to cargo (e.g. `--baseline=foo`).
        #[arg(last = true)]
        extra: Vec<String>,
    },
    /// Synthesise a `yee.toml` project file from a natural-language prompt.
    ///
    /// Phase 3.nl.0 walking skeleton: turn `"2.4 GHz patch on FR4"` into a
    /// typed [`yee_design::DesignIntent`], run the Balanis Ch. 14 +
    /// Pozar §3.8 closed-form synthesis ([`yee_design::InitialEstimate::from_intent`]),
    /// then emit a deterministic project TOML via [`yee_design::emit()`]. Two
    /// artefacts are written: `<output>` (the TOML) and
    /// `<output>.intent.json` (the typed [`yee_design::DesignIntent`] for
    /// round-trip / provenance).
    ///
    /// Stage-1 (prompt → intent) selection:
    /// - If `--offline` is passed **or** `ANTHROPIC_API_KEY` is unset,
    ///   [`yee_design::offline::parse`] is invoked deterministically.
    /// - Otherwise this CLI surfaces a "live LLM path lives in the yee-py
    ///   sidecar" message and exits non-zero. Per the R5 escape hatch,
    ///   PyO3 in-process embedding from `yee-cli` is deferred to
    ///   Phase 3.nl.0.1; callers wanting the LLM path today must use
    ///   `python -c "import yee; ..."` against the `yee-py` bindings.
    Design {
        /// Free-form natural-language prompt, e.g.
        /// `"2.4 GHz inset-fed patch on RO4003C"`.
        prompt: String,
        /// Output project-TOML path. A sidecar `<output>.intent.json` is
        /// written next to it.
        #[arg(long, short)]
        output: PathBuf,
        /// Force the deterministic offline parser instead of the LLM path.
        ///
        /// The offline parser is always selected if `ANTHROPIC_API_KEY` is
        /// unset in the environment, regardless of this flag.
        #[arg(long)]
        offline: bool,
        /// Optional LLM model id to forward to the sidecar
        /// (e.g. `"claude-sonnet-4-5"`). Ignored in the offline path.
        #[arg(long)]
        model: Option<String>,
    },
    /// Run an FDTD simulation end-to-end and emit the radiation pattern as JSON.
    ///
    /// Composes [`yee_fdtd::FdtdDriver`] from a vacuum [`yee_fdtd::YeeGrid`]
    /// with the supplied grid / source / NTFF parameters, runs the time
    /// loop to completion, and writes the θ-cut of `|E_θ|` at `φ = 0` as
    /// JSON to `--output` (or stdout when unset). The JSON shape is
    /// `{"theta_deg": [...], "e_theta_phi0": [...]}` with both vectors of
    /// equal length; angles span `[0°, 180°]` in 5° steps.
    FdtdRun {
        /// Grid dimensions (Nx, Ny, Nz) — `--grid 60 60 60`.
        #[arg(long, num_args = 3, default_values_t = [60_usize, 60, 60])]
        grid: Vec<usize>,
        /// Cell size in meters.
        #[arg(long, default_value_t = 5.0e-3)]
        dx: f64,
        /// Number of timesteps.
        #[arg(long, default_value_t = 800)]
        steps: usize,
        /// Source center cell (i, j, k).
        #[arg(long, num_args = 3, default_values_t = [30_usize, 30, 30])]
        source: Vec<usize>,
        /// Dipole length in cells.
        #[arg(long, default_value_t = 5)]
        dipole_length: usize,
        /// Source frequency in Hz.
        #[arg(long, default_value_t = 1.0e9)]
        freq: f64,
        /// NTFF surface pad in cells.
        #[arg(long, default_value_t = 4)]
        ntff_pad: usize,
        /// CPML thickness in cells.
        #[arg(long, default_value_t = 10)]
        cpml: usize,
        /// Output JSON path. If unset, write to stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
}

/// Arguments to [`run_design`], mirroring the [`Command::Design`] variant.
///
/// Held in a struct so the dispatch in [`run`] stays readable and so the
/// handler can be unit-tested independently of clap parsing.
#[derive(Debug, Clone)]
struct DesignArgs {
    prompt: String,
    output: PathBuf,
    offline: bool,
    #[allow(dead_code)] // populated by clap; consumed by the LLM-path stub.
    model: Option<String>,
}

/// Arguments to [`run_fdtd`], mirroring the [`Command::FdtdRun`] variant.
///
/// Held in a struct so the handler signature stays manageable and so the
/// `clap`-parsed variant can be passed through one field at a time.
#[derive(Debug, Clone)]
struct FdtdArgs {
    grid: Vec<usize>,
    dx: f64,
    steps: usize,
    source: Vec<usize>,
    dipole_length: usize,
    freq: f64,
    ntff_pad: usize,
    cpml: usize,
    output: Option<PathBuf>,
}

/// Which `yee-bench` criterion target to invoke.
///
/// Each variant maps to a `[[bench]]` entry in `crates/yee-bench/Cargo.toml`;
/// see [`run_bench`] for the variant → `--bench <name>` translation.
/// `All` runs every bench by omitting the `--bench` flag entirely.
#[derive(ValueEnum, Clone, Debug)]
enum BenchTarget {
    /// MoM solve on dipole 8×8 single freq.
    Mom,
    /// FDTD step on 50³ vacuum grid.
    Fdtd,
    /// GMRES vs direct LU at 128×128.
    Gmres,
    /// Gaussian-process fit + fit_ml.
    Gp,
    /// Full Bayesian-optimization run.
    Bo,
    /// Run all yee-bench benches.
    All,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum ExportFormat {
    /// Touchstone v1.1 (.s1p/.s2p/.s3p/.s4p).
    Touchstone,
    /// HDF5 (not yet enabled).
    Hdf5,
}

/// What `yee plot` should draw from the S-parameter sweep.
///
/// `Both` is a convenience that emits a dB plot and a Smith chart in one
/// invocation; the two output paths are derived from `--output` by inserting
/// `-db` / `-smith` before the file extension.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum PlotKind {
    /// `|S|` in dB vs frequency.
    Db,
    /// `S` on the Smith chart.
    Smith,
    /// `phase(S)` in degrees vs frequency.
    Phase,
    /// Emit both a dB plot and a Smith chart (two output files).
    Both,
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
            format,
            output,
            width,
            height,
            title,
            port,
        } => plot::run_plot(plot::PlotArgs {
            input,
            kind: format,
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
        Command::Bench { target, extra } => run_bench(target, extra),
        Command::Design {
            prompt,
            output,
            offline,
            model,
        } => run_design(DesignArgs {
            prompt,
            output,
            offline,
            model,
        }),
        Command::FdtdRun {
            grid,
            dx,
            steps,
            source,
            dipole_length,
            freq,
            ntff_pad,
            cpml,
            output,
        } => run_fdtd(FdtdArgs {
            grid,
            dx,
            steps,
            source,
            dipole_length,
            freq,
            ntff_pad,
            cpml,
            output,
        }),
    }
}

/// Drive [`yee_fdtd::FdtdDriver`] end-to-end and emit the resulting
/// [`yee_fdtd::RadiationPattern`] as JSON.
///
/// The JSON shape is built by hand (not via `serde::Serialize` on the
/// pattern struct) so this handler is robust against future changes to the
/// `RadiationPattern` derives: only the two public `Vec<f64>` fields
/// (`theta_deg`, `e_theta_phi0`) are touched. When `args.output` is set
/// the JSON is written to that path and a confirmation line is printed;
/// otherwise the JSON is sent to stdout so callers can pipe it.
fn run_fdtd(args: FdtdArgs) -> Result<ExitCode> {
    use yee_fdtd::{FdtdDriver, FdtdDriverConfig, YeeGrid};

    let (nx, ny, nz) = (args.grid[0], args.grid[1], args.grid[2]);
    let (i, j, k) = (args.source[0], args.source[1], args.source[2]);
    let grid = YeeGrid::vacuum(nx, ny, nz, args.dx);
    let cfg = FdtdDriverConfig {
        n_steps: args.steps,
        dipole_center_cells: (i, j, k),
        dipole_length_cells: args.dipole_length,
        source_freq_hz: args.freq,
        ntff_surface_pad_cells: args.ntff_pad,
        cpml_thickness_cells: args.cpml,
    };
    let pattern = FdtdDriver::new(grid, cfg).run();

    let payload = serde_json::json!({
        "theta_deg": pattern.theta_deg,
        "e_theta_phi0": pattern.e_theta_phi0,
    });
    let text = serde_json::to_string_pretty(&payload)?;

    if let Some(path) = args.output {
        std::fs::write(&path, &text)?;
        println!("Wrote {}", path.display());
    } else {
        println!("{text}");
    }

    Ok(ExitCode::SUCCESS)
}

/// Phase 3.nl.0 walking-skeleton handler for `yee design`.
///
/// Pipeline:
///
/// 1. Stage-1: resolve the prompt → [`yee_design::DesignIntent`].
///    - `--offline` or no `ANTHROPIC_API_KEY` ⇒ deterministic
///      [`yee_design::offline::parse`].
///    - Otherwise: print a message pointing at the `yee-py` sidecar
///      (`yee.design.from_prompt_llm`) and exit non-zero. PyO3 in-process
///      embedding from `yee-cli` is deferred to Phase 3.nl.0.1 per the R5
///      escape hatch.
/// 2. Stage-3: [`yee_design::InitialEstimate::from_intent`] computes the
///    Balanis Ch. 14 + Pozar §3.8 closed-form dimensions.
/// 3. Stage-5: [`yee_design::emit()`] renders the deterministic project TOML +
///    `intent.json` sidecar. Both are written to disk.
/// 4. The resolved dimensions are echoed to stdout in millimetres for the
///    engineer to eyeball.
///
/// Returns `ExitCode::SUCCESS` iff every stage succeeds. Errors propagate as
/// `anyhow` errors and the binary's top-level handler renders them on stderr
/// with `ExitCode::FAILURE`. The "LLM path not wired" surface is also
/// non-zero exit so a script that forgets `--offline` does not silently no-op.
fn run_design(args: DesignArgs) -> Result<ExitCode> {
    // Stage 1 — prompt → DesignIntent.
    let use_offline = args.offline || std::env::var_os("ANTHROPIC_API_KEY").is_none();
    let intent = if use_offline {
        yee_design::offline::parse(&args.prompt)
            .map_err(|e| anyhow::anyhow!("offline parser: {e}"))?
    } else {
        // Per the R5 escape hatch, the live LLM path lives in `yee-py`'s
        // `yee.design.from_prompt_llm`. Embedding the Python interpreter
        // in-process from the `yee` binary is deferred to Phase 3.nl.0.1
        // (it requires opting `yee-cli` into PyO3 and a `python3` toolchain
        // at link time — a tech-stack change out of scope for R5).
        eprintln!(
            "yee design: live LLM path is provided by the yee-py sidecar \
             (yee.design.from_prompt_llm). Rerun with --offline, or invoke \
             the sidecar from Python yourself. PyO3 in-process embedding \
             from yee-cli is deferred to Phase 3.nl.0.1."
        );
        return Ok(ExitCode::from(2));
    };

    // Stage 3 — DesignIntent → InitialEstimate.
    let estimate = yee_design::InitialEstimate::from_intent(&intent)
        .map_err(|e| anyhow::anyhow!("initial-estimate synthesis: {e}"))?;

    // Stage 5 — emit project TOML + intent.json sidecar.
    let artefacts = yee_design::emit(&estimate, &intent);

    std::fs::write(&args.output, &artefacts.toml)
        .map_err(|e| anyhow::anyhow!("writing project TOML to {}: {e}", args.output.display()))?;

    let intent_path = intent_sidecar_path(&args.output);
    std::fs::write(&intent_path, &artefacts.intent_json)
        .map_err(|e| anyhow::anyhow!("writing intent sidecar to {}: {e}", intent_path.display()))?;

    // Stdout summary — millimetre-scaled, plus the resolved permittivity and
    // centre frequency so the engineer can eyeball the result without
    // having to re-open the TOML.
    println!("Wrote {}", args.output.display());
    println!("Wrote {}", intent_path.display());
    println!("Resolved design:");
    println!(
        "  center_frequency = {:.4} GHz",
        intent.target_frequency_hz / 1.0e9
    );
    println!("  substrate eps_r  = {:.4}", estimate.substrate.eps_r);
    println!("  width            = {:.4} mm", estimate.width_m * 1.0e3);
    println!("  length           = {:.4} mm", estimate.length_m * 1.0e3);
    println!(
        "  inset_offset     = {:.4} mm",
        estimate.inset_offset_m * 1.0e3
    );
    println!(
        "  feed_width       = {:.4} mm",
        estimate.feed_width_m * 1.0e3
    );

    Ok(ExitCode::SUCCESS)
}

/// Compute the `<output>.intent.json` sidecar path next to the project TOML.
///
/// Suffix-append — `foo/bar.toml` → `foo/bar.toml.intent.json`. We do not
/// replace the extension because a downstream consumer reading
/// `<project>.intent.json` from the project-file directory should not have
/// to know whether the project file ended in `.toml`, `.yee`, or nothing.
fn intent_sidecar_path(project: &std::path::Path) -> PathBuf {
    let mut s = project.as_os_str().to_owned();
    s.push(".intent.json");
    PathBuf::from(s)
}

/// Shell out to `cargo bench -p yee-bench [--bench <name>] [-- <extra>...]`.
///
/// Stdout/stderr are inherited (not captured) so users see criterion's
/// live progress, the warmup countdown, and the per-bench summary tables
/// as they're emitted. The function maps the [`BenchTarget`] variant to a
/// concrete `--bench` argument; `BenchTarget::All` deliberately omits the
/// flag so cargo runs every `[[bench]]` entry in `crates/yee-bench`.
///
/// Returns [`ExitCode::SUCCESS`] iff `cargo` exits zero. A non-zero exit
/// (failed bench, missing target, criterion arg error) is surfaced as
/// [`ExitCode::FAILURE`] without further interpretation.
fn run_bench(target: BenchTarget, extra: Vec<String>) -> Result<ExitCode> {
    let mut cmd = std::process::Command::new("cargo");
    cmd.args(["bench", "-p", "yee-bench"]);
    match target {
        BenchTarget::Mom => {
            cmd.args(["--bench", "mom_solve"]);
        }
        BenchTarget::Fdtd => {
            cmd.args(["--bench", "fdtd_step"]);
        }
        BenchTarget::Gmres => {
            cmd.args(["--bench", "gmres_vs_direct"]);
        }
        BenchTarget::Gp => {
            cmd.args(["--bench", "gp_fit"]);
        }
        BenchTarget::Bo => {
            cmd.args(["--bench", "bo_step"]);
        }
        BenchTarget::All => {}
    }
    if !extra.is_empty() {
        cmd.arg("--");
        cmd.args(&extra);
    }
    let status = cmd.status()?;
    Ok(if status.success() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
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
