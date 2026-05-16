//! 2.4 GHz rectangular patch on FR-4 — Phase 0→1 driver.
//!
//! Hand-builds a planar mesh of a ~50-triangle rectangular patch
//! (29.2 mm × 38.0 mm copper, tag = 1) over a ground plane (tag = 2) with
//! the port edge between them, then runs a 21-point sweep from 2.0 GHz to
//! 3.0 GHz through [`yee_mom::PlanarMoM`].
//!
//! The substrate (FR-4, εr ≈ 4.4, h ≈ 1.6 mm) is *not* modelled in the
//! mesh — accurate physics for this geometry depends on the multilayer
//! Green's function machinery that ships in Phase 1.1. Until then,
//! `PlanarMoM::run` returns [`yee_core::Error::Unimplemented`] and this
//! example exits 0 with a clear message.

use anyhow::{Context, Result};
use nalgebra::Vector3;
use yee_core::{FreqRange, Solver};
use yee_mesh::TriMesh;
use yee_mom::PlanarMoM;

/// Patch length along x (mm → m). Sized for 2.4 GHz resonance on a
/// 1.6 mm FR-4 substrate (Balanis, *Antenna Theory*, 4th ed., §14.2).
const PATCH_LX_M: f64 = 0.0292;
/// Patch width along y.
const PATCH_LY_M: f64 = 0.0380;
/// Ground-plane half-extent margin around the patch (each direction).
const GROUND_MARGIN_M: f64 = 0.0200;

/// Subdivision of the patch into a regular triangle grid. `NX × NY` cells,
/// 2 triangles per cell. The resulting count is `2 * NX * NY + 8` (the
/// ground-plane apron contributes 8 triangles, see below).
const NX: usize = 4;
const NY: usize = 5;

/// Build a hand-meshed rectangular patch + simple ground apron.
///
/// Layout (top-down, x → right, y → up):
///
/// ```text
///                ┌───────────────────┐
///                │   ground apron    │
///                │  ┌─────────────┐  │
///                │  │  copper     │  │
///                │  │  patch      │  │
///                │  │ (NX × NY    │  │
///                │  │  cells)     │  │
///                │  └─────────────┘  │
///                │                   │
///                └───────────────────┘
/// ```
///
/// All triangles carrying tag `1` are the radiating patch; tag `2` is the
/// ground apron. The basis-function convention being established for
/// Phase 1.0 treats any shared edge between a tag-1 and tag-2 triangle as
/// a port edge; here that is the bottom edge of the patch where it touches
/// the apron strip.
fn build_patch_mesh() -> Result<TriMesh> {
    let dx = PATCH_LX_M / (NX as f64);
    let dy = PATCH_LY_M / (NY as f64);

    let mut vertices: Vec<Vector3<f64>> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();

    // ---- Patch grid: (NX+1) x (NY+1) vertices, indexed row-major in y. ----
    // Patch vertex `(ix, iy)` lives at index `iy * (NX+1) + ix`; we inline
    // that arithmetic below rather than allocating an indexing closure.
    let patch_origin = Vector3::new(0.0, 0.0, 0.0);
    for iy in 0..=NY {
        for ix in 0..=NX {
            vertices.push(Vector3::new(
                patch_origin.x + (ix as f64) * dx,
                patch_origin.y + (iy as f64) * dy,
                0.0,
            ));
        }
    }
    let stride = (NX + 1) as u32;
    for iy in 0..NY as u32 {
        for ix in 0..NX as u32 {
            let v00 = iy * stride + ix;
            let v10 = v00 + 1;
            let v01 = v00 + stride;
            let v11 = v01 + 1;
            // Two triangles per cell, both tagged as patch copper.
            triangles.push([v00, v10, v11]);
            tags.push(1);
            triangles.push([v00, v11, v01]);
            tags.push(1);
        }
    }

    // ---- Ground apron: four trapezoidal strips around the patch. ----
    // The apron shares the patch's bottom-row vertices on its inner edge,
    // which is exactly the property we want for the port-edge convention
    // (shared edge ⇒ tag-1 / tag-2 boundary ⇒ port).
    //
    // We add four outer corner vertices to define the apron's outer
    // rectangle, then triangulate each side as two triangles.
    let outer_xmin = patch_origin.x - GROUND_MARGIN_M;
    let outer_xmax = patch_origin.x + PATCH_LX_M + GROUND_MARGIN_M;
    let outer_ymin = patch_origin.y - GROUND_MARGIN_M;
    let outer_ymax = patch_origin.y + PATCH_LY_M + GROUND_MARGIN_M;

    let outer_bl = vertices.len() as u32;
    vertices.push(Vector3::new(outer_xmin, outer_ymin, 0.0));
    let outer_br = vertices.len() as u32;
    vertices.push(Vector3::new(outer_xmax, outer_ymin, 0.0));
    let outer_tr = vertices.len() as u32;
    vertices.push(Vector3::new(outer_xmax, outer_ymax, 0.0));
    let outer_tl = vertices.len() as u32;
    vertices.push(Vector3::new(outer_xmin, outer_ymax, 0.0));

    // Patch corner indices.
    let patch_bl: u32 = 0;
    let patch_br: u32 = NX as u32;
    let patch_tr: u32 = (NY as u32) * stride + (NX as u32);
    let patch_tl: u32 = (NY as u32) * stride;

    // Bottom apron strip (two triangles): outer_bl, outer_br, patch_br, patch_bl.
    triangles.push([outer_bl, outer_br, patch_br]);
    tags.push(2);
    triangles.push([outer_bl, patch_br, patch_bl]);
    tags.push(2);
    // Right apron strip.
    triangles.push([patch_br, outer_br, outer_tr]);
    tags.push(2);
    triangles.push([patch_br, outer_tr, patch_tr]);
    tags.push(2);
    // Top apron strip.
    triangles.push([patch_tl, patch_tr, outer_tr]);
    tags.push(2);
    triangles.push([patch_tl, outer_tr, outer_tl]);
    tags.push(2);
    // Left apron strip.
    triangles.push([outer_bl, patch_bl, patch_tl]);
    tags.push(2);
    triangles.push([outer_bl, patch_tl, outer_tl]);
    tags.push(2);

    TriMesh::new(vertices, triangles, tags).context("building patch + apron mesh")
}

fn main() -> Result<()> {
    // See note in half-wave-dipole: env-filter isn't enabled in the
    // workspace's `tracing-subscriber` feature set; use a plain fmt
    // subscriber.
    let _ = tracing_subscriber::fmt::try_init();

    println!(
        "patch-2g4: building rectangular patch mesh ({:.1} mm × {:.1} mm, NX={NX}, NY={NY})",
        PATCH_LX_M * 1e3,
        PATCH_LY_M * 1e3,
    );
    let mesh = build_patch_mesh()?;
    println!(
        "patch-2g4: mesh has {} vertices, {} triangles ({} patch / {} ground)",
        mesh.vertices.len(),
        mesh.n_tris(),
        mesh.tags.iter().filter(|&&t| t == 1).count(),
        mesh.tags.iter().filter(|&&t| t == 2).count(),
    );

    let band = FreqRange::new(2.0e9, 3.0e9, 21).context("building 2–3 GHz sweep")?;
    println!(
        "patch-2g4: invoking PlanarMoM::run over [{:.2}, {:.2}] GHz, {} points",
        band.start_hz / 1e9,
        band.stop_hz / 1e9,
        band.n_points,
    );

    let solver = PlanarMoM::default();
    match solver.run(&mesh, band) {
        Ok(sparams) => {
            println!(
                "patch-2g4: PlanarMoM returned {} frequency points ({} ports)",
                sparams.freq_hz.len(),
                sparams.n_ports
            );
            let out = std::path::PathBuf::from("target/example-output/patch-2g4.s1p");
            if let Some(parent) = out.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating {}", parent.display()))?;
            }
            sparams
                .write_touchstone(&out, 50.0)
                .map_err(|e| anyhow::anyhow!("writing Touchstone: {e}"))?;
            println!("patch-2g4: wrote {}", out.display());
        }
        Err(yee_core::Error::Unimplemented(msg)) => {
            println!("patch-2g4: PlanarMoM::run is a Phase 0 stub ({msg}).");
            println!(
                "patch-2g4: accurate results require Phase 1.1 multilayer Green's functions \
                 (FR-4 εr ≈ 4.4, h ≈ 1.6 mm). Re-run once that lands."
            );
        }
        Err(other) => {
            return Err(anyhow::anyhow!("unexpected solver error: {other}"));
        }
    }

    println!("patch-2g4: done.");
    Ok(())
}
