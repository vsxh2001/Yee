//! Integration tests for [`yee_mesh::TetMesh3D::cavity_uniform`].
//!
//! Companion to `crates/yee-mesh/src/cavity.rs`. The Kuhn 6-tet
//! decomposition is asserted at the public-API level here so the
//! production-gate fem-eig-001 driver (Step T7) can consume it
//! without re-deriving its invariants.

use yee_mesh::{Error, TetMesh3D};

#[test]
fn cavity_uniform_cell_count() {
    // 2 × 2 × 2 bricks → 8 bricks × 6 Kuhn tets = 48 tets,
    // (2+1)^3 = 27 vertices.
    let m = TetMesh3D::cavity_uniform(1.0, 1.0, 1.0, 2, 2, 2).unwrap();
    assert_eq!(m.n_tets(), 48);
    assert_eq!(m.n_verts(), 27);
    assert_eq!(m.vertex_material.len(), 27);
    assert_eq!(m.tetrahedron_material.len(), 48);
}

#[test]
fn cavity_uniform_total_volume() {
    // Sum of signed_volume over all tets must equal a·b·d to floating
    // point precision. Uses a non-trivial set of extents so a coincidental
    // axis-product cancellation cannot mask a per-axis scaling bug.
    let a = 0.5;
    let b = 0.75;
    let d = 1.25;
    let m = TetMesh3D::cavity_uniform(a, b, d, 3, 4, 5).unwrap();
    let expected = a * b * d;
    let total: f64 = (0..m.n_tets()).map(|i| m.signed_volume(i)).sum();
    assert!(
        (total - expected).abs() < 1e-12,
        "total volume = {total}, expected {expected} within 1e-12"
    );
}

#[test]
fn cavity_uniform_all_tets_positive_volume() {
    // Construction must hand the consumer a mesh whose every tet has
    // strictly positive signed volume — the Kuhn table has three
    // sign-flipped entries, so this is the regression that catches
    // any missed re-orientation in `TetMesh3D::new`.
    let m = TetMesh3D::cavity_uniform(1.0, 2.0, 3.0, 2, 3, 4).unwrap();
    for i in 0..m.n_tets() {
        let v = m.signed_volume(i);
        assert!(
            v > 0.0,
            "tet {i} has non-positive signed volume {v}; Kuhn reorientation failed"
        );
    }
}

#[test]
fn cavity_uniform_wr90_dims_sanity() {
    // WR-90-based cavity dimensions a = 22.86 mm, b = 10.16 mm,
    // d = 30 mm with the v0 fem-eig-001 brick count 4 × 2 × 4
    // (selected so the production-gate solve fits the v0 sparse-eigen
    // budget; see Step T7 brief). Expected: 4 · 2 · 4 · 6 = 192 tets,
    // total volume 0.02286 · 0.01016 · 0.030 = 6.969768e-6 m³.
    let a = 0.02286;
    let b = 0.01016;
    let d = 0.030;
    let m = TetMesh3D::cavity_uniform(a, b, d, 4, 2, 4).unwrap();
    assert_eq!(m.n_tets(), 192);
    assert_eq!(m.n_verts(), 5 * 3 * 5);

    let expected = a * b * d;
    let total: f64 = (0..m.n_tets()).map(|i| m.signed_volume(i)).sum();
    assert!(
        (total - expected).abs() < 1e-12,
        "WR-90 cavity total volume = {total}, expected {expected} within 1e-12"
    );
    for i in 0..m.n_tets() {
        assert!(
            m.signed_volume(i) > 0.0,
            "WR-90 tet {i} has non-positive volume"
        );
    }
}

#[test]
fn cavity_uniform_invalid_zero_n_errors() {
    let err = TetMesh3D::cavity_uniform(1.0, 1.0, 1.0, 0, 2, 2).unwrap_err();
    match err {
        Error::Invalid(msg) => assert!(
            msg.contains("nx, ny, nz"),
            "expected message naming the offending arguments; got: {msg}"
        ),
        _ => panic!("expected Error::Invalid"),
    }
}
