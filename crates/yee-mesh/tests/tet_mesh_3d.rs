//! Integration tests for [`yee_mesh::TetMesh3D`].
//!
//! The bulk of fine-grained invariant checks (out-of-range indices,
//! material-tag mismatches, default tag fill) live as in-module
//! `#[cfg(test)]` tests next to the implementation in
//! `crates/yee-mesh/src/tetmesh.rs`. This file holds the high-level
//! integration cases that exercise the public API the way `yee-fem`
//! and its callers will use it.

use nalgebra::Vector3;
use yee_mesh::TetMesh3D;

/// Reference tet from Jin §9.4 used to back-out the local Nedelec
/// basis gradients (`∇λ_i = (face-normal opposite i) / (3·V)`). Volume
/// is `1/6` and the centroid is at `(1/4, 1/4, 1/4)`.
#[test]
fn reference_unit_tet_volume_and_centroid() {
    let m = TetMesh3D::new(
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ],
        vec![[0, 1, 2, 3]],
        None,
        None,
    )
    .unwrap();

    assert!((m.signed_volume(0) - 1.0 / 6.0).abs() < 1e-15);
    let c = m.centroid(0);
    assert!((c.x - 0.25).abs() < 1e-15);
    assert!((c.y - 0.25).abs() < 1e-15);
    assert!((c.z - 0.25).abs() < 1e-15);
}

/// A 5-tet decomposition of the unit cube `[0,1]^3` exercises:
///
/// * Construction over multiple tets, several of which arrive in
///   negative-orientation form and must be silently re-oriented.
/// * The total-volume sum matches the cube volume (1) to floating
///   point precision — this is the regression that catches any
///   sign-handling bug in `signed_volume` or in the reorientation
///   logic.
#[test]
fn unit_cube_five_tet_total_volume_is_one() {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0), // 0
        Vector3::new(1.0, 0.0, 0.0), // 1
        Vector3::new(1.0, 1.0, 0.0), // 2
        Vector3::new(0.0, 1.0, 0.0), // 3
        Vector3::new(0.0, 0.0, 1.0), // 4
        Vector3::new(1.0, 0.0, 1.0), // 5
        Vector3::new(1.0, 1.0, 1.0), // 6
        Vector3::new(0.0, 1.0, 1.0), // 7
    ];
    let tetrahedra = vec![
        [0, 1, 3, 4], // corner near v0
        [1, 2, 3, 6], // corner near v2
        [1, 5, 4, 6], // corner near v5
        [3, 4, 7, 6], // corner near v7
        [1, 3, 4, 6], // central
    ];
    let m = TetMesh3D::new(vertices, tetrahedra, None, None).unwrap();

    assert_eq!(m.n_tets(), 5);
    let total: f64 = (0..m.n_tets()).map(|i| m.signed_volume(i)).sum();
    assert!(
        (total - 1.0).abs() < 1e-12,
        "5-tet cube total volume = {total}, expected 1.0 within 1e-12"
    );
    for i in 0..m.n_tets() {
        assert!(
            m.signed_volume(i) > 0.0,
            "tet {i} has non-positive signed volume after construction"
        );
    }
}

/// Construction must auto-reorient a single negatively-oriented tet
/// without surfacing an error to the caller. This is the silent-fix
/// path documented on `TetMesh3D::new`.
#[test]
fn negatively_oriented_tet_is_silently_reoriented() {
    // Reference tet with v2 and v3 swapped → signed volume -1/6.
    let m = TetMesh3D::new(
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0), // formerly v3
            Vector3::new(0.0, 1.0, 0.0), // formerly v2
        ],
        vec![[0, 1, 2, 3]],
        None,
        None,
    )
    .unwrap();

    assert!(m.signed_volume(0) > 0.0);
    assert!((m.signed_volume(0) - 1.0 / 6.0).abs() < 1e-15);
}
