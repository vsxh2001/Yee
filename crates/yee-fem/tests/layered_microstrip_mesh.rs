//! Gate test for [`yee_fem::layered_microstrip_mesh`] (FEM-EM brick 2,
//! ADR-0153).
//!
//! These are **geometric** assertions only — no `sweep`/LU solve, so the
//! whole test runs sub-second. It verifies that:
//!
//! 1. the tet count matches the `cavity_uniform` Kuhn-6 expectation
//!    exactly;
//! 2. the substrate/air repaint partitions the tets into the exact
//!    proportion implied by the z-split (substrate = bottom `n_sub / nz`
//!    fraction of the box height), and the tags are the right way round;
//! 3. the ground/trace predicates, fed through
//!    [`OpenBoundarySolver::interior_edges_matching`], each select a
//!    **non-empty** edge set (geometry + picker compose);
//! 4. the material database resolves FR-4 (`ε_r = 4.4`) and air
//!    (`ε_r = 1.0`).

use std::collections::HashMap;
use std::f64::consts::PI;

use yee_fem::{
    AIR_TAG, FR4_TAG, FaceKind, MaterialDatabase, OpenBoundarySolver, layered_microstrip_mesh,
};
use yee_mesh::TetMesh3D;

// ---------------------------------------------------------------------
// Fixture geometry. box_h = 4 mm, nz = 8 → dz = 0.5 mm, so the substrate
// (sub_h = 1 mm) is exactly n_sub = 2 cells: z = 1 mm lands on mesh plane
// 2 of 8 (a cell boundary, never mid-cell). Substrate therefore occupies
// the bottom 2/8 = 1/4 of the box height.
// ---------------------------------------------------------------------
const BOX_W: f64 = 4.0e-3;
const BOX_H: f64 = 4.0e-3;
const LINE_LEN: f64 = 10.0e-3;
const SUB_H: f64 = 1.0e-3;
const TRACE_W: f64 = 0.5e-3;
const NX: usize = 4;
const NY: usize = 10;
const NZ: usize = 8;

/// `n_sub`: number of z-cells inside the substrate (sub_h / dz).
const N_SUB: usize = 2;

/// Count exterior faces of a mesh (multiplicity-one face filter). Mirror
/// of the helper in `open_boundary_sweep_matrix.rs`; needed to size the
/// all-PEC `face_kinds` vector for the picker-only solver.
fn exterior_face_count(mesh: &TetMesh3D) -> usize {
    let mut face_map: HashMap<[usize; 3], usize> = HashMap::new();
    const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for tet in &mesh.tetrahedra {
        for &[a, b, c] in TET_FACES.iter() {
            let mut key = [tet[a], tet[b], tet[c]];
            key.sort_unstable();
            *face_map.entry(key).or_insert(0) += 1;
        }
    }
    face_map.values().filter(|&&c| c == 1).count()
}

#[test]
fn layered_microstrip_mesh_geometry_and_picker_compose() {
    let (mesh, db, ground_pred, trace_pred) =
        layered_microstrip_mesh(BOX_W, BOX_H, LINE_LEN, SUB_H, TRACE_W, NX, NY, NZ)
            .expect("fixture geometry must build (z = sub_h lands on a mesh plane)");

    // (1) Exact tet count: Kuhn-6 per brick.
    let expected_tets = NX * NY * NZ * 6;
    assert_eq!(
        mesh.n_tets(),
        expected_tets,
        "tet count must match cavity_uniform Kuhn-6 expectation"
    );
    assert_eq!(
        mesh.tetrahedron_material.len(),
        expected_tets,
        "per-tet material array must cover every tet"
    );

    // (2) Substrate/air partition. Because sub_h = n_sub * dz exactly, no
    // brick straddles the interface, so every tet in a z-layer k < n_sub
    // is FR-4 and every tet in k >= n_sub is air. The FR-4 count is
    // therefore EXACTLY nx*ny*n_sub*6, i.e. the bottom n_sub/nz fraction.
    let expected_fr4 = NX * NY * N_SUB * 6;
    let expected_air = NX * NY * (NZ - N_SUB) * 6;
    let n_fr4 = mesh
        .tetrahedron_material
        .iter()
        .filter(|&&t| t == FR4_TAG)
        .count();
    let n_air = mesh
        .tetrahedron_material
        .iter()
        .filter(|&&t| t == AIR_TAG)
        .count();
    assert_eq!(
        n_fr4, expected_fr4,
        "substrate tet count must be exactly nx*ny*n_sub*6 (bottom {N_SUB}/{NZ} of box height)"
    );
    assert_eq!(
        n_air, expected_air,
        "air tet count must be exactly nx*ny*(nz-n_sub)*6"
    );
    assert_eq!(
        n_fr4 + n_air,
        expected_tets,
        "every tet must be either FR-4 or air (tags partition the mesh)"
    );

    // Independent cross-check on the tag assignment itself: every FR-4 tet
    // really has its centroid below sub_h and every air tet at/above it
    // (proves the assertion is geometric, not a tautology on the counts).
    for (idx, tet) in mesh.tetrahedra.iter().enumerate() {
        let cz = 0.25
            * (mesh.vertices[tet[0]].z
                + mesh.vertices[tet[1]].z
                + mesh.vertices[tet[2]].z
                + mesh.vertices[tet[3]].z);
        if mesh.tetrahedron_material[idx] == FR4_TAG {
            assert!(
                cz < SUB_H,
                "tet {idx} tagged FR-4 but centroid.z = {cz} >= sub_h = {SUB_H}"
            );
        } else {
            assert!(
                cz >= SUB_H,
                "tet {idx} tagged air but centroid.z = {cz} < sub_h = {SUB_H}"
            );
        }
    }

    // (3) Predicates compose with the brick-1 interior-PEC picker and
    // select NON-EMPTY edge sets. Build a picker-only solver (all-PEC
    // boundary, no ports) — interior_edges_matching is a pure geometric
    // edge walk, so no solve happens.
    let n_exterior = exterior_face_count(&mesh);
    let solver = OpenBoundarySolver::new(
        &mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("picker-only solver must construct");

    let ground_edges = solver.interior_edges_matching(ground_pred);
    let trace_edges = solver.interior_edges_matching(trace_pred);
    assert!(
        !ground_edges.is_empty(),
        "ground_pred must select at least one edge on the z = 0 plane"
    );
    assert!(
        !trace_edges.is_empty(),
        "trace_pred must select at least one edge on the z = sub_h trace footprint"
    );

    // The ground plane is the full box footprint at z = 0; the trace is a
    // narrow strip at z = sub_h. The ground edge set must be strictly
    // larger than the trace edge set (sanity that the windows differ and
    // the trace is genuinely a sub-region, not the whole top plane).
    assert!(
        ground_edges.len() > trace_edges.len(),
        "ground footprint ({} edges) must exceed the narrower trace footprint ({} edges)",
        ground_edges.len(),
        trace_edges.len()
    );

    // (4) Material database resolves FR-4 and air. Evaluate at a
    // representative microwave ω (the value is irrelevant — FR-4 here has
    // no dispersive poles, so ε is flat).
    let omega = 2.0 * PI * 5.0e9; // 5 GHz
    assert_eq!(
        db.eps_at(FR4_TAG, omega).re,
        4.4,
        "FR-4 tag must resolve to ε_r = 4.4"
    );
    assert_eq!(
        db.eps_at(AIR_TAG, omega).re,
        1.0,
        "air tag must resolve to ε_r = 1.0"
    );
}

/// A `sub_h` that does not land on a mesh plane must be rejected — the
/// substrate/air interface has to be a cell boundary (ADR-0108), never
/// mid-cell. This guards the "clean tag proportion" contract the gate
/// above relies on.
#[test]
fn off_plane_substrate_is_rejected() {
    // dz = box_h / nz = 4 mm / 8 = 0.5 mm; sub_h = 0.7 mm is NOT a
    // multiple of 0.5 mm, so the interface falls mid-cell.
    let err = layered_microstrip_mesh(BOX_W, BOX_H, LINE_LEN, 0.7e-3, TRACE_W, NX, NY, NZ);
    assert!(
        err.is_err(),
        "an off-plane sub_h must be rejected so the interface stays a cell boundary"
    );
}
