//! # yee CLI
//!
//! Top-level command-line tool for the Yee electromagnetic-simulation studio.
//!
//! Phase 0 wires the following subcommands:
//!
//! - `yee validate <mom|fdtd|fem|synth|all> [--json]` — runs the
//!   `yee-validation` aggregator (real mom-001 NEC-4 gate; FEM-eig-001/002/004/005
//!   gates; filter-synthesis synth-001/synth-002/filt-001 gates;
//!   Phase-deferred FDTD cases report `Skipped`). Filters by target prefix
//!   (`mom-*`, FDTD-family, `fem-*`, `synth-*`/`filt-*`) and exits 1 if any
//!   included case failed.
//!   Pass `--list` to print the registered-case inventory (id, solver,
//!   policy, description) and exit 0 **without running any solver** — instant,
//!   unlike the default path which runs the ~7-8 min mom-001 solve. Add
//!   `--json` to `--list` to emit that inventory as a JSON array for CI/tooling.
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

mod filter;
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
    /// `cpml-*` / `ntff-*` / `dispersive-*` / `fdtd-*`; `fem` includes
    /// `fem-eig-*` cases (FEM eigenmode suite); `synth` includes the
    /// filter-synthesis gates (`synth-*` / `filt-*`, Filter Phase F0);
    /// `all` runs the full report. Default output is a human-readable
    /// table; pass `--json` to emit the serialized
    /// [`yee_validation::Report`].
    ///
    /// Exit code is non-zero iff any *included* case has status
    /// `Failed`. `Skipped` cases (Phase-deferred placeholders) do not
    /// fail the run.
    ///
    /// Pass `--list` to print the registered-case inventory (CASE /
    /// SOLVER / POLICY / DESCRIPTION) and exit 0 without running any
    /// solver. `--list` short-circuits before the aggregator, so it is
    /// instant; it still honours `target` (e.g. `yee validate fem
    /// --list` lists only `fem-*` cases). Combining `--list --json` emits
    /// the inventory as a JSON array (descriptors) instead of the table.
    Validate {
        /// Which solver to validate.
        #[arg(value_enum, default_value_t = ValidateTarget::All)]
        target: ValidateTarget,
        /// Emit JSON report to stdout (default: human-readable table).
        #[arg(long)]
        json: bool,
        /// List the registered cases (id, solver, policy, description)
        /// and exit 0 without running any solver. Filtered by `target`.
        #[arg(long)]
        list: bool,
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
    ///
    /// ## Multi-trace overlay
    ///
    /// Pass `--entry <ij>` one or more times (e.g. `--entry 11 --entry 21`) to
    /// overlay multiple S-matrix entries in one dB or phase plot. Use `--all`
    /// to overlay every entry of a multi-port file. Indices are 1-based and
    /// match the Touchstone convention. `--format smith` and `--format both`
    /// are not supported with `--entry`/`--all`.
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
        /// Ignored when --entry or --all is supplied.
        #[arg(long, default_value_t = 0)]
        port: usize,
        /// S-matrix entry to include in a multi-trace overlay (1-based, e.g.
        /// `11` for S₁₁, `21` for S₂₁). Repeat to overlay multiple entries.
        /// When present, activates the multi-trace path (`plot_sparams_db` /
        /// `plot_sparams_phase`). Incompatible with `--format smith` / `--format both`.
        #[arg(long = "entry", value_name = "IJ")]
        entries: Vec<String>,
        /// Overlay every S-matrix entry of a multi-port file. Activates the
        /// multi-trace path. Incompatible with `--format smith` / `--format both`.
        #[arg(long)]
        all: bool,
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
    /// Filter design (synthesis + ideal response + spec-mask check).
    ///
    /// Filter Phase F0 walking skeleton. The `synth` subcommand parses a
    /// [`yee_filter::FilterSpec`] TOML, synthesizes the lowpass prototype and
    /// all-pole coupling matrix ([`yee_filter::synthesize`]), sweeps the
    /// closed-form ideal response ([`yee_filter::ideal_response`]), writes the
    /// S-parameters as a Touchstone `.s2p` via `yee-io`, and grades the
    /// response against the spec mask ([`yee_filter::check_mask`]). Exit 0 on a
    /// PASS verdict, 1 on a mask FAIL.
    Filter {
        #[command(subcommand)]
        command: FilterCommand,
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

/// Subcommands under `yee filter`.
#[derive(Subcommand, Debug)]
enum FilterCommand {
    /// Synthesize a filter from a `FilterSpec` TOML and report the spec-mask
    /// verdict.
    ///
    /// Parses `<spec>`, synthesizes the prototype + coupling matrix, sweeps the
    /// closed-form ideal response, writes a Touchstone `.s2p`, and prints the
    /// mask PASS/FAIL verdict. Exit 0 on PASS, 1 on FAIL.
    Synth {
        /// Path to the `FilterSpec` TOML.
        spec: PathBuf,
        /// Output Touchstone path (`.s2p`). Defaults to the spec stem with a
        /// `.s2p` extension next to the spec file.
        #[arg(long, short)]
        output: Option<PathBuf>,
        /// Also render the synthesized |S21| response with the spec-mask
        /// forbidden regions overlaid, to this image path (`.png`/`.svg`).
        #[arg(long)]
        plot: Option<PathBuf>,
        /// Substrate relative permittivity `ε_r` for the F1.2.0 physical
        /// dimensioning (default FR-4 `4.4`).
        #[arg(long, default_value_t = 4.4)]
        eps_r: f64,
        /// Substrate dielectric height `h` in millimetres for the F1.2.0
        /// physical dimensioning (default FR-4 `1.6` mm).
        #[arg(long, default_value_t = 1.6)]
        h_mm: f64,
        /// Also write the synthesized edge-coupled layout as an SVG to this
        /// path (F1.2.0 `dimension_edge_coupled_layout`).
        #[arg(long)]
        layout_svg: Option<PathBuf>,
        /// Also write the synthesized edge-coupled layout as a single-copper
        /// Gerber to this path (F1.4.0 `yee_export::layout_to_gerber`).
        #[arg(long)]
        gerber: Option<PathBuf>,
        /// Also write the synthesized edge-coupled layout as a KiCad 7 board
        /// (`.kicad_pcb`) to this path (F1.4.1b `yee_export::layout_to_kicad_pcb`).
        #[arg(long)]
        kicad_pcb: Option<PathBuf>,
        /// Export the **lumped-LC** board (F2.2 `yee_filter::lumped_board`)
        /// instead of the planar edge-coupled layout. The `--layout-svg`,
        /// `--gerber`, and `--kicad-pcb` writers then emit the lumped board
        /// (signal line + ground rail + every L/C component pad) rather than
        /// the distributed half-wave resonators. The synthesized dimensions /
        /// Touchstone printout is unaffected.
        #[arg(long)]
        lumped: bool,
        /// SMD chip footprint for the `--lumped` board: `0402`, `0603`
        /// (default), or `0805`. Ignored unless `--lumped` is set.
        #[arg(long, value_enum, default_value_t = FootprintArg::Smd0603)]
        footprint: FootprintArg,
        /// Per-resonator unloaded quality factor `Q_u` (ADR-0161). When set,
        /// the exported Touchstone `.s2p` (and `--plot`) carries the realistic
        /// finite-Q lumped-LC response (`ladder_s21_lossy`) — with midband
        /// insertion loss and rounded passband corners — instead of the ideal
        /// lossless one. Independent of `--lumped`. Must be finite and `> 0`.
        #[arg(long)]
        q_unloaded: Option<f64>,
        /// Write the JLCPCB SMT-assembly upload set for the **lumped**
        /// realization into this directory (J4, ADR-0164): the lumped board's
        /// copper + outline Gerbers, `bom.csv` (autopicked LCSC parts,
        /// `Comment,Designator,Footprint,LCSC Part #`), and `cpl.csv`
        /// (`Designator,Mid X,Mid Y,Layer,Rotation`). Always uses the lumped LC
        /// board (the BOM/CPL describe the lumped parts), independent of
        /// `--lumped`; the directory is created if missing. Unrealizable values
        /// get a blank LCSC # (flagged, never dropped) and a one-line note.
        #[arg(long)]
        jlcpcb: Option<PathBuf>,
    },
}

/// SMD chip footprint selector for `yee filter synth --lumped`.
///
/// A thin clap-`ValueEnum` mirror of [`yee_filter::Footprint`] so the
/// `--footprint 0402|0603|0805` flag parses to a typed value (and rejects an
/// unknown size at the CLI boundary). Mapped to the library enum by
/// [`FootprintArg::to_footprint`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum FootprintArg {
    /// 0402 (1.0 × 0.5 mm body) chip.
    #[value(name = "0402")]
    Smd0402,
    /// 0603 (1.6 × 0.8 mm body) chip — the F2.2 default.
    #[value(name = "0603")]
    Smd0603,
    /// 0805 (2.0 × 1.25 mm body) chip.
    #[value(name = "0805")]
    Smd0805,
}

impl FootprintArg {
    /// Map the CLI selector to the [`yee_filter::Footprint`] library enum.
    fn to_footprint(self) -> yee_filter::Footprint {
        match self {
            FootprintArg::Smd0402 => yee_filter::Footprint::Smd0402,
            FootprintArg::Smd0603 => yee_filter::Footprint::Smd0603,
            FootprintArg::Smd0805 => yee_filter::Footprint::Smd0805,
        }
    }
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
    /// Finite-element method eigenmode suite (Phase 4).
    Fem,
    /// Filter-synthesis gates (`synth-*` / `filt-*`, Filter Phase F0).
    Synth,
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
        Command::Validate { target, json, list } => {
            if list {
                Ok(run_validate_list(target, json))
            } else {
                run_validate(target, json)
            }
        }
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
            entries,
            all,
        } => plot::run_plot(plot::PlotArgs {
            input,
            kind: format,
            output,
            width,
            height,
            title,
            port,
            entries,
            all_entries: all,
        }),
        Command::Completions { shell } => {
            let mut cmd = Cli::command();
            let bin_name = cmd.get_name().to_string();
            generate(shell, &mut cmd, bin_name, &mut io::stdout());
            Ok(ExitCode::SUCCESS)
        }
        Command::Bench { target, extra } => run_bench(target, extra),
        Command::Filter { command } => match command {
            FilterCommand::Synth {
                spec,
                output,
                plot,
                eps_r,
                h_mm,
                layout_svg,
                gerber,
                kicad_pcb,
                lumped,
                footprint,
                q_unloaded,
                jlcpcb,
            } => filter::run_synth(
                &spec,
                output.as_deref(),
                plot.as_deref(),
                eps_r,
                h_mm,
                layout_svg.as_deref(),
                gerber.as_deref(),
                kicad_pcb.as_deref(),
                lumped,
                footprint.to_footprint(),
                q_unloaded,
                jlcpcb.as_deref(),
            ),
        },
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
/// [`ValidateTarget::Fdtd`]; `fem-*` for [`ValidateTarget::Fem`];
/// `synth-*` / `filt-*` for [`ValidateTarget::Synth`]; everything for
/// [`ValidateTarget::All`].
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
        ValidateTarget::Fem => case_id.starts_with("fem-"),
        ValidateTarget::Synth => case_id.starts_with("synth-") || case_id.starts_with("filt-"),
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

/// List the registered validation cases (filtered by `target`) and
/// return [`ExitCode::SUCCESS`] — **no solver is executed**.
///
/// Backs `yee validate --list`. Pulls the case inventory from
/// [`yee_validation::list_cases`] (which maps the single-source
/// registry to descriptors without invoking any runner), filters with
/// the same [`case_matches_target`] used by the run path, and prints a
/// fixed-width `CASE | SOLVER | POLICY | DESCRIPTION` table. Listing
/// never "fails", so the exit code is always success.
fn run_validate_list(target: ValidateTarget, json: bool) -> ExitCode {
    use yee_validation::{ExecutionPolicy, Solver};

    let cases: Vec<yee_validation::CaseDescriptor> = yee_validation::list_cases()
        .into_iter()
        .filter(|d| case_matches_target(d.id, target))
        .collect();

    if json {
        // Emit the filtered inventory as a JSON array (ADR-0083). No solver
        // is run; each CaseDescriptor serializes as
        // `{"id","solver","description","policy"}`. CaseDescriptor + its enums
        // derive Serialize, so serializing our owned data cannot fail.
        let s = serde_json::to_string_pretty(&cases)
            .expect("CaseDescriptor serialization is infallible");
        println!("{s}");
        return ExitCode::SUCCESS;
    }

    const H_CASE: &str = "CASE";
    const H_SOLVER: &str = "SOLVER";
    const H_POLICY: &str = "POLICY";
    const H_DESC: &str = "DESCRIPTION";

    // Pre-format the categorical cells so column widths derive from the
    // final strings (mirrors the column-sizing idiom in
    // `print_human_report`).
    struct Row {
        id: &'static str,
        solver: &'static str,
        policy: &'static str,
        description: &'static str,
    }
    let rows: Vec<Row> = cases
        .iter()
        .map(|d| Row {
            id: d.id,
            solver: match d.solver {
                Solver::Mom => "MoM",
                Solver::Fdtd => "FDTD",
                Solver::Fem => "FEM",
                Solver::Synth => "Synth",
            },
            policy: match d.policy {
                ExecutionPolicy::Run => "Run",
                ExecutionPolicy::SkippedWallTime => "Skipped(wall-time)",
                ExecutionPolicy::SkippedGateOpen => "Skipped(gate-open)",
            },
            description: d.description,
        })
        .collect();

    let w_case = rows
        .iter()
        .map(|r| r.id.len())
        .max()
        .unwrap_or(0)
        .max(H_CASE.len());
    let w_solver = rows
        .iter()
        .map(|r| r.solver.len())
        .max()
        .unwrap_or(0)
        .max(H_SOLVER.len());
    let w_policy = rows
        .iter()
        .map(|r| r.policy.len())
        .max()
        .unwrap_or(0)
        .max(H_POLICY.len());

    println!("{H_CASE:<w_case$}  {H_SOLVER:<w_solver$}  {H_POLICY:<w_policy$}  {H_DESC}");
    for r in &rows {
        let id = r.id;
        let solver = r.solver;
        let policy = r.policy;
        let description = r.description;
        println!("{id:<w_case$}  {solver:<w_solver$}  {policy:<w_policy$}  {description}");
    }

    ExitCode::SUCCESS
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
    fn target_filter_fem_matches_fem_prefix() {
        assert!(case_matches_target("fem-eig-001", ValidateTarget::Fem));
        assert!(case_matches_target("fem-eig-006", ValidateTarget::Fem));
        assert!(!case_matches_target("mom-001", ValidateTarget::Fem));
        assert!(!case_matches_target("cpml-001", ValidateTarget::Fem));
    }

    #[test]
    fn target_filter_all_matches_everything() {
        assert!(case_matches_target("mom-001", ValidateTarget::All));
        assert!(case_matches_target("cpml-001", ValidateTarget::All));
        assert!(case_matches_target("fem-eig-001", ValidateTarget::All));
        assert!(case_matches_target("anything-else", ValidateTarget::All));
    }
}
