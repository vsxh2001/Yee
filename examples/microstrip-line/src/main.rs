//! Microstrip-line example: pure-mesh smoke test for the Gmsh import pipeline.
//!
//! This example exercises the `yee_mesh::Session` API surface:
//! `Session::new → import_step → mesh → tris`. The actual STEP file path can
//! be supplied via `--step <path>`; the default is a built-in placeholder so
//! the example still runs (and exits 0) without external assets.
//!
//! Without the `gmsh` cargo feature on `yee-mesh`, every session method
//! returns [`yee_mesh::Error::NotEnabled`]. Phase 0 ships this build path,
//! so this binary's *current* job is to demonstrate that the *typed*
//! pipeline compiles and to print a clear "feature not enabled" message,
//! exiting 0. Once Phase 1.mesh.0 wires in the bindgen FFI, the same driver
//! will report real triangle counts.

use anyhow::Result;
use std::path::PathBuf;

struct Args {
    step_path: Option<PathBuf>,
}

fn parse_args() -> Args {
    // Tiny hand-rolled parser to avoid pulling a CLI dependency into the
    // smoke example. We only support `--step <path>` (and `-h/--help`).
    let mut step_path: Option<PathBuf> = None;
    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--step" => {
                step_path = iter.next().map(PathBuf::from);
            }
            "-h" | "--help" => {
                println!("usage: microstrip-line [--step <path-to-step-file>]");
                std::process::exit(0);
            }
            other => {
                eprintln!("microstrip-line: ignoring unknown argument `{other}`");
            }
        }
    }
    Args { step_path }
}

fn main() -> Result<()> {
    // See note in half-wave-dipole: env-filter isn't enabled in the
    // workspace's `tracing-subscriber` feature set; use a plain fmt
    // subscriber. `try_init` is non-fatal if a subscriber is already set.
    let _ = tracing_subscriber::fmt::try_init();

    let args = parse_args();

    println!("microstrip-line: attempting to open a yee-mesh Session");
    match yee_mesh::Session::new() {
        Ok(mut session) => {
            // `gmsh` feature path. Currently `Session::new` returns Ok only
            // when the feature is on, and at Phase 0 the body is `todo!()`,
            // so we will not actually reach this branch in any shipped
            // build. We still wire the call site so that when Phase 1.mesh.0
            // lands the integration is exercised end-to-end.
            let step = args.step_path.unwrap_or_else(|| {
                PathBuf::from("examples/microstrip-line/reference/microstrip.step")
            });
            println!("microstrip-line: importing STEP file {}", step.display());
            session
                .import_step(&step)
                .map_err(|e| anyhow::anyhow!("import_step: {e}"))?;
            println!("microstrip-line: meshing surface (dim = 2)");
            session
                .mesh(2)
                .map_err(|e| anyhow::anyhow!("mesh: {e}"))?;
            let tris = session
                .tris()
                .map_err(|e| anyhow::anyhow!("tris: {e}"))?;
            println!(
                "microstrip-line: mesh has {} vertices, {} triangles",
                tris.vertices.len(),
                tris.n_tris()
            );
        }
        Err(yee_mesh::Error::NotEnabled) => {
            println!(
                "microstrip-line: yee-mesh built without the `gmsh` feature; \
                 mesh pipeline is gated until Phase 1.mesh.0."
            );
            if let Some(path) = args.step_path {
                println!(
                    "microstrip-line: requested STEP file `{}` will be processed once \
                     the gmsh feature is wired in.",
                    path.display()
                );
            }
            println!(
                "microstrip-line: rebuild with `--features yee-mesh/gmsh` (Phase 1.mesh.0) \
                 to exercise the real pipeline."
            );
        }
        Err(other) => {
            return Err(anyhow::anyhow!("unexpected mesh-session error: {other}"));
        }
    }

    println!("microstrip-line: done.");
    Ok(())
}
