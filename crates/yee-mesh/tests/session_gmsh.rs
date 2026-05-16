//! Integration tests that exercise the real Gmsh FFI bodies. These only
//! build when the `gmsh` feature is enabled and only execute when invoked
//! with `--ignored` (so default `cargo test --features gmsh` stays fast
//! and tolerant of an SDK that needs runtime resources).
//!
//! Phase 1.mesh.0 ships these `#[ignore]`-marked to keep CI green on hosts
//! that lack a runtime Gmsh library. A later phase will wire an
//! environment-gated runner that flips the `--ignored` flag when an SDK
//! is provisioned.

#![cfg(feature = "gmsh")]

use yee_mesh::Session;

/// Build a 1×1×1 cube via OCC, surface-mesh it, and assert the resulting
/// `TriMesh` has at least the analytical lower bound (12 triangles — two per
/// face) with vertex coordinates spanning [0, 1].
#[test]
#[ignore]
fn cube_mesh_via_gmsh_occ_box() {
    let mut sess = Session::new().expect("Session::new must succeed with the gmsh feature on");
    sess.add_occ_box(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        .expect("add_occ_box must succeed");
    sess.synchronize()
        .expect("synchronize must succeed after add_occ_box");
    sess.mesh(2).expect("mesh(2) must succeed");
    let mesh = sess.tris().expect("tris() must return a TriMesh");

    assert!(
        mesh.n_tris() >= 12,
        "a unit cube surface mesh should have at least 12 triangles (2 per face), got {}",
        mesh.n_tris()
    );

    let max_coord: f64 = mesh
        .vertices
        .iter()
        .map(|v| v.x.max(v.y).max(v.z))
        .fold(f64::NEG_INFINITY, f64::max);
    let min_coord: f64 = mesh
        .vertices
        .iter()
        .map(|v| v.x.min(v.y).min(v.z))
        .fold(f64::INFINITY, f64::min);
    assert!(
        (max_coord - 1.0).abs() < 1e-9,
        "max vertex coordinate should be 1.0, got {max_coord}"
    );
    assert!(
        min_coord.abs() < 1e-9,
        "min vertex coordinate should be 0.0, got {min_coord}"
    );
}
