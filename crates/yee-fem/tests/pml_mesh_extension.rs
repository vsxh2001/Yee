//! Phase 4.fem.eig.3.5 P2 — unit tests for
//! [`yee_fem::extend_mesh_with_pml`].
//!
//! Verifies the PML mesh-extension contract per the spec §4.2 / plan
//! P2 brief:
//!
//! 1. `pml_zero_thickness_no_op` — `thickness_cells = 0` returns the
//!    original mesh tet count unchanged with every tet tagged
//!    [`yee_fem::PmlClass::Interior`]. This is the backward-compat
//!    canary the v3-equivalence assembly gate (P4) consumes.
//! 2. `pml_shell_tet_count_matches_brick_extension` — extending a
//!    `(nx, ny, nz)` cavity by `thickness_cells = t` on a single
//!    Cartesian face adds exactly `nx · ny · t · 6` PML tets (for an
//!    axial face) — the brick layer count times Kuhn-6 mesh density.
//! 3. `pml_inner_boundary_has_continuous_vertex_layer` — every vertex
//!    on the original cavity boundary face is preserved as a shared
//!    vertex in the extended mesh; no duplicates.
//! 4. `pml_class_depth_monotonic_outward` — for tets along the chosen
//!    PML axis, `d` is monotonically increasing as you move outward
//!    from the cavity inner boundary.
//!
//! Edge / corner wedges are explicitly out of scope for v3.5 per
//! ADR-0043 §4 — the multi-face stress test is queued for Phase
//! 4.fem.eig.3.5.1.

use yee_fem::{PmlAxis, PmlClass, extend_mesh_with_pml};
use yee_mesh::TetMesh3D;

/// 4 × 2 × 4 brick cavity matching the WR-90 stub aspect ratio.
fn wr90_stub_4_2_4() -> TetMesh3D {
    TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap()
}

#[test]
fn pml_zero_thickness_no_op() {
    let mesh = wr90_stub_4_2_4();
    let n_tets_in = mesh.tetrahedra.len();

    let (extended, classes, _faces) =
        extend_mesh_with_pml(&mesh, &[PmlAxis::ZMin], 0).expect("extend_mesh_with_pml");

    assert_eq!(
        extended.tetrahedra.len(),
        n_tets_in,
        "thickness_cells = 0 must preserve tet count"
    );
    assert_eq!(classes.len(), n_tets_in);
    assert!(
        classes.iter().all(|c| c.is_interior()),
        "thickness_cells = 0 must tag every tet Interior"
    );
}

#[test]
fn pml_shell_tet_count_matches_brick_extension() {
    let mesh = wr90_stub_4_2_4();
    let nx_in = 4;
    let ny_in = 2;
    let nz_in = 4;
    let n_tets_in = nx_in * ny_in * nz_in * 6;
    assert_eq!(mesh.tetrahedra.len(), n_tets_in);

    let t: usize = 3;
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&mesh, &[PmlAxis::ZMin], t).expect("extend_mesh_with_pml");

    // Single-axis ZMin extension adds nx · ny · t bricks → · 6 tets.
    let expected_total = (nx_in * ny_in * (nz_in + t)) * 6;
    assert_eq!(
        extended.tetrahedra.len(),
        expected_total,
        "extended mesh tet count must equal extended-brick count × 6"
    );
    assert_eq!(classes.len(), expected_total);

    let n_pml_tets = classes.iter().filter(|c| !c.is_interior()).count();
    let expected_pml = nx_in * ny_in * t * 6;
    assert_eq!(
        n_pml_tets, expected_pml,
        "PML-tagged tet count must equal added-brick-layer count × 6"
    );

    let n_interior_tets = classes.iter().filter(|c| c.is_interior()).count();
    assert_eq!(
        n_interior_tets, n_tets_in,
        "interior-tagged tet count must equal original cavity tet count"
    );
}

#[test]
fn pml_inner_boundary_has_continuous_vertex_layer() {
    let mesh = wr90_stub_4_2_4();
    let nvx = 5;
    let nvy = 3;
    let nvz = 5;
    assert_eq!(mesh.vertices.len(), nvx * nvy * nvz);

    let t: usize = 2;
    let (extended, _classes, _faces) =
        extend_mesh_with_pml(&mesh, &[PmlAxis::ZMin], t).expect("extend_mesh_with_pml");

    let expected_vert_count = nvx * nvy * (nvz + t);
    assert_eq!(
        extended.vertices.len(),
        expected_vert_count,
        "extended mesh vertex count must equal extended-brick vertex grid"
    );

    // Verify the cavity's z = 0 face vertices are shared in the
    // extended mesh — find every vertex with z ≈ 0 (relative to the
    // extended origin, which has shifted by `-t · dz`).
    let dz: f64 = 0.030 / 4.0;
    let expected_z0_extended = -(t as f64) * dz;
    let count_at_inner = extended
        .vertices
        .iter()
        .filter(|v| (v.z - 0.0).abs() < 1e-9)
        .count();
    let count_at_outer = extended
        .vertices
        .iter()
        .filter(|v| (v.z - expected_z0_extended).abs() < 1e-9)
        .count();
    // The original cavity face at z = 0 still has `nvx · nvy`
    // vertices in the extended mesh, and the new outer truncation
    // surface (at z = expected_z0_extended) has the same count.
    assert_eq!(
        count_at_inner,
        nvx * nvy,
        "inner cavity/PML interface must have nvx · nvy shared vertices"
    );
    assert_eq!(
        count_at_outer,
        nvx * nvy,
        "outer PML truncation surface must have nvx · nvy vertices"
    );
}

#[test]
fn pml_class_depth_monotonic_outward() {
    // Use a +z PML so depth increases with z.
    let mesh = wr90_stub_4_2_4();
    let t: usize = 4;
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&mesh, &[PmlAxis::ZMax], t).expect("extend_mesh_with_pml");

    // Compute centroid z for every PML-tagged tet; verify the
    // (centroid_z, depth) relationship is monotonic.
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    for (tet_idx, tet) in extended.tetrahedra.iter().enumerate() {
        let class = classes[tet_idx];
        let PmlClass::PmlZ { d } = class else {
            continue;
        };
        let v0 = extended.vertices[tet[0]];
        let v1 = extended.vertices[tet[1]];
        let v2 = extended.vertices[tet[2]];
        let v3 = extended.vertices[tet[3]];
        let centroid_z = (v0.z + v1.z + v2.z + v3.z) / 4.0;
        pairs.push((centroid_z, d));
    }
    assert!(
        !pairs.is_empty(),
        "ZMax PML must produce at least one PmlZ-tagged tet"
    );

    // Sort by centroid z and verify d is monotonically non-decreasing.
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut last_d = 0.0;
    for (_z, d) in &pairs {
        assert!(
            *d + 1e-9 >= last_d,
            "PML depth must be monotonically non-decreasing outward; \
             got d = {d} at z = {_z}, previous d = {last_d}"
        );
        last_d = *d;
    }
    // Last depth should be at most `(t - 1 + 0.5) * dz` — top brick
    // centroid in the t-deep shell.
    let dz: f64 = 0.030 / 4.0;
    let max_expected = ((t - 1) as f64 + 0.5) * dz + 1e-9;
    assert!(
        last_d <= max_expected,
        "max PML depth {last_d} must be ≤ ({t} - 0.5) · dz = {max_expected}"
    );
}
