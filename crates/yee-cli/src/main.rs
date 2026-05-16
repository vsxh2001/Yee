//! # yee CLI
//!
//! Top-level command-line tool for the Yee electromagnetic-simulation studio.
//!
//! Phase 0 wires the following subcommands:
//!
//! - `yee validate <mom|fdtd|all>` — prints planned validation cases.
//! - `yee mesh <path>` — Phase 0 stub (wired in a later commit).
//! - `yee export <input> --format <touchstone|hdf5> <output>` — Phase 0
//!   stub (wired in a later commit).
//! - `yee run <project>` — Phase 0 stub from the scaffold.

use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};

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
        input: std::path::PathBuf,
    },
    /// Run a simulation defined by a project file (Phase 0 stub).
    Run {
        /// Path to the project TOML.
        project: std::path::PathBuf,
    },
    /// Export results to Touchstone or HDF5.
    Export {
        /// Path to the run results.
        results: std::path::PathBuf,
        /// Output format: `touchstone` or `hdf5`.
        #[arg(long, default_value = "touchstone")]
        format: String,
    },
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

fn main() -> Result<()> {
    tracing_subscriber::fmt().with_target(false).init();

    let cli = Cli::parse();
    match cli.command {
        Command::Validate { target } => run_validate(target),
        Command::Mesh { input } => {
            println!("yee mesh {} — Phase 0 stub.", input.display());
        }
        Command::Run { project } => {
            println!("yee run {} — Phase 0 stub.", project.display());
        }
        Command::Export { results, format } => {
            println!(
                "yee export {} --format {format} — Phase 0 stub.",
                results.display()
            );
        }
    }
    Ok(())
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
