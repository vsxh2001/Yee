//! # yee CLI
//!
//! `yee validate | mesh | run | export`. Phase 0: subcommands print their intended
//! behavior; only `validate` does meaningful work (calls into per-crate validation harnesses).

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "yee", version, about = "Yee — GPU-accelerated electromagnetic simulation studio")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Run the validation suite for a given solver (Phase 0: prints planned cases).
    Validate {
        /// Which solver to validate: `mom`, `fdtd`, or `all`.
        #[arg(default_value = "all")]
        solver: String,
    },
    /// Mesh a geometry file via Gmsh.
    Mesh {
        /// Input geometry path (.step / .iges / .kicad_pcb).
        input: std::path::PathBuf,
    },
    /// Run a simulation defined by a project file.
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

fn main() -> Result<()> {
    tracing_subscriber::fmt().with_target(false).init();

    let cli = Cli::parse();
    match cli.command {
        Command::Validate { solver } => {
            println!("yee validate {solver} — Phase 0 stub. See crates/*/validation/README.md.");
        }
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
