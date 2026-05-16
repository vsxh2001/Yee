//! Half-wave dipole example: Phase 0 walking-skeleton driver for [`yee_mom::PlanarMoM`].
//!
//! Builds a minimal two-triangle planar mesh with tags `[1, 2]` (the port edge
//! is the shared edge between the differently-tagged triangles, per the basis
//! convention being established in Phase 1.0). The mesh is then fed to
//! [`yee_mom::PlanarMoM::run`] at a single resonance-band frequency.
//!
//! Phase 0 contract: `PlanarMoM::run` currently returns
//! [`yee_core::Error::Unimplemented`]. This example treats that as the
//! expected state — it prints a friendly message, writes a placeholder
//! Touchstone `.s1p` to `target/example-output/dipole.s1p`, and exits 0.
//! Once Track A merges Phase 1.0 MoM physics, the same driver will compute
//! the real free-space dipole impedance (Z ≈ 73 + j42 Ω).

use anyhow::{Context, Result};
use nalgebra::Vector3;
use num_complex::Complex64;
use yee_core::{FreqRange, Solver};
use yee_io::touchstone;
use yee_mesh::TriMesh;
use yee_mom::{PlanarMoM, SParameters};

/// Operating frequency: 300 MHz half-wave dipole (λ/2 ≈ 0.5 m).
const F0_HZ: f64 = 300.0e6;

/// Reference impedance for the Touchstone output.
const Z0_OHMS: f64 = 50.0;

fn build_dipole_mesh() -> Result<TriMesh> {
    // Two coplanar triangles sharing an edge along the y-axis. The shared
    // edge (between v1 and v2) is the port edge: its two adjacent triangles
    // carry different physical tags, which is the convention that the
    // upcoming basis machinery will use to identify the feed.
    //
    //        v3 ──── v2
    //         \  T0 /│
    //          \   / │
    //           \ /  │
    //            v1  │   (port edge along v1 → v2)
    //           / \  │
    //          /   \ │
    //         / T1  \│
    //        v0 ──── v4
    //
    // (Schematic only — the actual coordinates below place v1 at the origin
    // with v2 a quarter-wavelength up the y-axis, so the dipole spans λ/2
    // around the feed.)
    let quarter_wave = yee_core::units::C0 / (4.0 * F0_HZ);
    let vertices = vec![
        Vector3::new(0.0, -quarter_wave, 0.0), // v0: lower arm tip
        Vector3::new(0.0, 0.0, 0.0),           // v1: port edge — lower endpoint
        Vector3::new(0.0, quarter_wave, 0.0),  // v2: port edge — upper endpoint
        Vector3::new(0.001, quarter_wave, 0.0), // v3: thin width to give triangles area
        Vector3::new(0.001, -quarter_wave, 0.0), // v4: same, lower side
    ];
    let triangles = vec![
        [1, 2, 3], // upper arm triangle
        [0, 1, 4], // lower arm triangle
    ];
    // Different tags on the two arms identify the port edge between them.
    let tags = vec![1, 2];
    TriMesh::new(vertices, triangles, tags).context("building two-triangle dipole mesh")
}

fn placeholder_s1p() -> SParameters {
    // A single-port S-parameter container at the operating frequency. The
    // datum is `S11 = 0 + 0j` — explicitly a *placeholder*, not a physical
    // estimate. The intent is to exercise the SParameters → Touchstone
    // write path during Phase 0 so the I/O surface stays linked into the
    // example binary as the solver code lands.
    SParameters {
        freq_hz: vec![F0_HZ],
        data: vec![vec![Complex64::new(0.0, 0.0)]],
        n_ports: 1,
    }
}

fn write_placeholder(out_path: &std::path::Path) -> Result<()> {
    if let Some(parent) = out_path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output directory {}", parent.display()))?;
    }
    let mut file = placeholder_s1p().to_touchstone(Z0_OHMS);
    file.comments.push(
        " Phase 0 placeholder S-parameters from examples/half-wave-dipole.".to_string(),
    );
    file.comments
        .push(" Real values land once yee-mom Phase 1.0 ships.".to_string());
    touchstone::write(out_path, &file)
        .with_context(|| format!("writing Touchstone file {}", out_path.display()))?;
    Ok(())
}

fn main() -> Result<()> {
    // `tracing-subscriber` ships without the `env-filter` feature by default
    // in this workspace; use the plain fmt subscriber. Verbosity can still be
    // adjusted via `RUST_LOG` once env-filter is opted in upstream.
    let _ = tracing_subscriber::fmt::try_init();

    println!("half-wave-dipole: building two-triangle planar mesh");
    let mesh = build_dipole_mesh()?;
    println!(
        "half-wave-dipole: mesh has {} vertices, {} triangles, tags = {:?}",
        mesh.vertices.len(),
        mesh.n_tris(),
        mesh.tags
    );

    let band = FreqRange::new(F0_HZ, F0_HZ * 1.000_001, 1)
        .context("building single-point frequency range")?;
    println!(
        "half-wave-dipole: invoking PlanarMoM::run at {:.3} MHz",
        F0_HZ / 1.0e6
    );

    let solver = PlanarMoM::default();
    match solver.run(&mesh, band) {
        Ok(sparams) => {
            println!(
                "half-wave-dipole: PlanarMoM returned {} frequency points ({} ports)",
                sparams.freq_hz.len(),
                sparams.n_ports
            );
            let out =
                std::path::PathBuf::from("target/example-output/dipole.s1p");
            sparams
                .write_touchstone(&out, Z0_OHMS)
                .map_err(|e| anyhow::anyhow!("writing Touchstone: {e}"))?;
            println!("half-wave-dipole: wrote {}", out.display());
        }
        Err(yee_core::Error::Unimplemented(msg)) => {
            println!("half-wave-dipole: PlanarMoM::run is a Phase 0 stub ({msg}).");
            println!(
                "half-wave-dipole: writing placeholder S-parameters until Track A merges."
            );
            let out = std::path::PathBuf::from("target/example-output/dipole.s1p");
            write_placeholder(&out)?;
            println!("half-wave-dipole: wrote {}", out.display());
        }
        Err(other) => {
            return Err(anyhow::anyhow!("unexpected solver error: {other}"));
        }
    }

    println!("half-wave-dipole: done.");
    Ok(())
}
