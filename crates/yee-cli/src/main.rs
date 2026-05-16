//! # yee CLI
//!
//! Top-level command-line tool for the Yee electromagnetic-simulation studio.
//!
//! Phase 0 wires the following subcommands:
//!
//! - `yee validate <mom|fdtd|all>` — prints planned validation cases.
//! - `yee mesh <path>` — constructs a `yee_mesh::Session`. Without the
//!   `gmsh` feature this exits with code 2 and a guidance message.
//! - `yee export <input> --format <touchstone|hdf5> <output>` — reads/writes
//!   Touchstone via `yee-io`. `hdf5` is not yet enabled and exits with code 2.
//! - `yee completions <shell>` — emits a shell completion script to stdout
//!   (`bash`, `zsh`, `fish`).
//! - `yee run <project>` — Phase 0 stub from the scaffold.

use std::io;
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::Result;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{Shell, generate};

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
    /// Run the validation suite for a given solver (Phase 0: prints planned cases).
    Validate {
        /// Which solver to validate.
        #[arg(value_enum, default_value_t = ValidateTarget::All)]
        target: ValidateTarget,
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
        Command::Validate { target } => {
            run_validate(target);
            Ok(ExitCode::SUCCESS)
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

fn run_validate(target: ValidateTarget) {
    match target {
        ValidateTarget::Mom => print_mom_report(),
        ValidateTarget::Fdtd => print_fdtd_report(),
        ValidateTarget::All => {
            print_mom_report();
            print_fdtd_report();
        }
    }
}

fn print_mom_report() {
    println!("yee validate mom (Phase 0)");
    println!("  planned cases:");
    println!("   - mom-001 half-wave dipole (Phase 1)");
    println!("   - mom-002 microstrip Z0  (Phase 1)");
    println!("   - mom-003 patch resonance (Phase 1)");
}

fn print_fdtd_report() {
    println!("yee validate fdtd (Phase 0)");
    println!("  Phase 2 deliverable — yee-fdtd not yet available");
}
