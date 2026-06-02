//! FEM-EM brick 1 (ADR-0153) — interior-PEC wiring gate for
//! [`yee_fem::OpenBoundarySolver`].
//!
//! Brick 1 adds two surfaces:
//!
//! * [`OpenBoundarySolver::with_interior_pec_edges`] — unions
//!   caller-supplied **global edge IDs** into the private
//!   `pec_global_edges` (Dirichlet) set, so a conductor floating inside
//!   the mesh volume (e.g. a microstrip trace) is eliminated from the
//!   interior-DoF basis exactly like an exterior-PEC-face edge. The
//!   assembly primitive (`assemble_complex_with_pec_edges`) and the
//!   driven path (`assemble_driven_system` / `sweep` / `sweep_matrix`)
//!   already honour `pec_global_edges` verbatim, so interior PEC is
//!   inherited by every solve path with no further wiring.
//! * [`OpenBoundarySolver::interior_edges_matching`] — returns the
//!   global edge IDs whose endpoints satisfy a world-coordinate
//!   predicate, in the **same** canonical global-edge index space the
//!   assembly primitive filters against (it rebuilds the exact
//!   `EdgeKey` / `LOCAL_EDGES` dedup map).
//!
//! Gate (pure assembly, no LU sweep — runs in seconds):
//!
//! 1. [`interior_pec_removes_exactly_the_tagged_edges`] — on a tiny
//!    all-PEC-exterior `cavity_uniform(2, 2, 2)` air mesh, record
//!    `N0 = assemble_driven_system(ω).interior_edges.len()`. Pick a
//!    non-empty interior mid-plane edge set `E` via
//!    `interior_edges_matching`, assert `E` is disjoint from
//!    `pec_global_edges()` (genuinely interior), rebuild with
//!    `.with_interior_pec_edges(E)`, and assert the new
//!    `interior_edges.len() == N0 - |E|` (exact integer).
//! 2. [`retagging_exterior_pec_edge_is_idempotent`] — re-tagging an
//!    edge that is already exterior-PEC leaves `interior_edges.len()`
//!    unchanged (set union).
//!
//! References:
//! * Open-boundary face-kind assembly fixture style:
//!   `crates/yee-fem/tests/open_boundary_assembly.rs`.
//! * Canonical global-edge numbering:
//!   `crates/yee-fem/src/assembly.rs` `TetEdgeTable::build` /
//!   `assemble_complex_with_pec_edges`.

use std::collections::HashSet;
use std::f64::consts::PI;

use yee_fem::{FaceKind, MaterialDatabase, OpenBoundarySolver};
use yee_mesh::TetMesh3D;

/// Test angular frequency at 10 GHz in vacuum.
fn omega_10ghz() -> f64 {
    2.0 * PI * 10.0e9
}

/// Discover the exterior-face count of `mesh` from the authoritative
/// source — the constructor's own length-mismatch error — then build an
/// all-PEC-exterior [`OpenBoundarySolver`] over a unit-cube-style air
/// mesh (empty material DB → vacuum everywhere).
///
/// Passing a deliberately-wrong (length-0) `face_kinds` makes
/// [`OpenBoundarySolver::new`] report `"... exterior-face count N"`; we
/// parse the trailing integer rather than re-deriving the private
/// exterior-face enumeration in the test.
fn all_pec_solver(mesh: &TetMesh3D) -> OpenBoundarySolver<'_> {
    // `OpenBoundarySolver` is not `Debug`, so we can't use `expect_err`;
    // match the Result by hand.
    let msg = match OpenBoundarySolver::new(mesh, Vec::new(), Vec::new(), MaterialDatabase::new()) {
        Ok(_) => panic!("length-0 face_kinds must mismatch the exterior-face count"),
        Err(e) => e.to_string(),
    };
    let n_exterior: usize = msg
        .rsplit(|c: char| !c.is_ascii_digit())
        .find(|tok| !tok.is_empty())
        .and_then(|tok| tok.parse().ok())
        .unwrap_or_else(|| panic!("could not parse exterior-face count from error: {msg}"));
    assert!(
        n_exterior > 0,
        "parsed a zero exterior-face count from error: {msg}"
    );

    OpenBoundarySolver::new(
        mesh,
        vec![FaceKind::Pec; n_exterior],
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("all-PEC-exterior solver must build on a valid cavity mesh")
}

// ---------------------------------------------------------------------
// Gate 1 — interior-PEC removes exactly the tagged edges
// ---------------------------------------------------------------------

#[test]
fn interior_pec_removes_exactly_the_tagged_edges() {
    // 2 x 2 x 2 brick air cavity — small enough that the pure-assembly
    // path is sub-second, large enough to have genuine interior edges
    // (the z = d/2 mid-plane carries edges that touch no boundary face).
    let a = 1.0;
    let b = 1.0;
    let d = 1.0;
    let mesh = TetMesh3D::cavity_uniform(a, b, d, 2, 2, 2).unwrap();

    let solver = all_pec_solver(&mesh);
    let omega = omega_10ghz();

    // Baseline interior-DoF count with only exterior PEC.
    let n0 = solver
        .assemble_driven_system(omega)
        .unwrap()
        .interior_edges
        .len();
    assert!(
        n0 > 0,
        "all-PEC cavity must still have interior DoFs (edges off every boundary face)"
    );

    // Pick interior edges via the geometric picker: edges fully
    // contained in the z = d/2 mid-plane (both endpoints at z = d/2, so
    // they cannot lie on the z = 0 or z = d cavity walls). Some of those
    // mid-plane edges still graze a SIDE wall (x = 0, x = a, y = 0,
    // y = b), so per the picker's own contract we difference against the
    // exterior-PEC set to get the genuinely-interior subset. This is the
    // realistic caller pattern (pick by geometry, drop edges already PEC).
    let mid_z = d / 2.0;
    let pec_before: HashSet<usize> = solver.pec_global_edges().iter().copied().collect();
    let mid_plane: Vec<usize> = solver
        .interior_edges_matching(|p, q| (p.z - mid_z).abs() < 1e-9 && (q.z - mid_z).abs() < 1e-9);

    // The picker must return a strictly-ascending, de-duplicated list in
    // the canonical global-edge space.
    let unique: HashSet<usize> = mid_plane.iter().copied().collect();
    assert_eq!(
        unique.len(),
        mid_plane.len(),
        "interior_edges_matching returned duplicate IDs"
    );
    assert!(
        mid_plane.windows(2).all(|w| w[0] < w[1]),
        "interior_edges_matching must return a strictly-ascending, de-duplicated list"
    );
    assert!(
        !mid_plane.is_empty(),
        "mid-plane z = {mid_z} must contain at least one edge on a 2x2x2 cavity"
    );

    let e: Vec<usize> = mid_plane
        .iter()
        .copied()
        .filter(|gid| !pec_before.contains(gid))
        .collect();
    assert!(
        !e.is_empty(),
        "mid-plane z = {mid_z} must contain at least one GENUINELY interior edge \
         (mid-plane minus exterior-PEC) on a 2x2x2 cavity"
    );

    // Genuinely interior: the picked set must be disjoint from the
    // exterior-PEC Dirichlet set, otherwise removing them would not be a
    // clean -|E| delta.
    for &gid in &e {
        assert!(
            !pec_before.contains(&gid),
            "edge {gid} from the interior pick is already in the exterior-PEC set"
        );
    }

    // Rebuild with the interior edges folded into the PEC set. The
    // driven path inherits the interior PEC verbatim — no change to
    // assemble_driven_system was needed.
    let solver_pec = all_pec_solver(&mesh).with_interior_pec_edges(e.iter().copied());

    // The interior edges that were tagged must now be in the PEC set.
    let pec_after = solver_pec.pec_global_edges();
    for &gid in &e {
        assert!(
            pec_after.contains(&gid),
            "edge {gid} was passed to with_interior_pec_edges but is missing from pec_global_edges()"
        );
    }

    let n1 = solver_pec
        .assemble_driven_system(omega)
        .unwrap()
        .interior_edges
        .len();

    // Exact integer bookkeeping: each genuinely-interior edge tagged PEC
    // drops exactly one interior DoF.
    assert_eq!(
        n1,
        n0 - e.len(),
        "interior-DoF count after tagging |E| = {} interior edges must be N0 - |E| = {} - {} = {}, got {}",
        e.len(),
        n0,
        e.len(),
        n0 - e.len(),
        n1
    );

    // And the lift map must contain none of the newly-PEC edges.
    let system = solver_pec.assemble_driven_system(omega).unwrap();
    for &gid in &e {
        assert!(
            !system.interior_edges.contains(&gid),
            "edge {gid} is interior-PEC but still appears in the interior-edge lift map"
        );
    }
}

// ---------------------------------------------------------------------
// Gate 2 — re-tagging an existing exterior-PEC edge is idempotent
// ---------------------------------------------------------------------

#[test]
fn retagging_exterior_pec_edge_is_idempotent() {
    let mesh = TetMesh3D::cavity_uniform(1.0, 1.0, 1.0, 2, 2, 2).unwrap();

    let solver = all_pec_solver(&mesh);
    let omega = omega_10ghz();

    let n0 = solver
        .assemble_driven_system(omega)
        .unwrap()
        .interior_edges
        .len();

    // Grab an edge that is ALREADY in the exterior-PEC set.
    let existing: usize = *solver
        .pec_global_edges()
        .iter()
        .next()
        .expect("all-PEC cavity must have a non-empty exterior-PEC set");

    // Re-tagging it via with_interior_pec_edges must be a no-op on the
    // interior-DoF count (set union, not multiset).
    let n1 = all_pec_solver(&mesh)
        .with_interior_pec_edges([existing])
        .assemble_driven_system(omega)
        .unwrap()
        .interior_edges
        .len();

    assert_eq!(
        n1, n0,
        "re-tagging an already-PEC edge changed the interior-DoF count ({n0} -> {n1}); \
         with_interior_pec_edges must be an idempotent set union"
    );
}
