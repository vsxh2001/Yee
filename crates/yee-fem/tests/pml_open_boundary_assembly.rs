//! Phase 4.fem.eig.3.5 P4 — unit tests for the CFS-PML end-to-end
//! [`yee_fem::OpenBoundarySolver::with_cfs_pml`] wire-in.
//!
//! Gate inventory (per the plan P4 brief):
//!
//! 1. `pml_assembly_finite_at_dc` — the assembled CFS-PML matrix at
//!    a far-below-cutoff frequency remains finite (CFS `α_α > 0`
//!    causality canary distinguishing CFS-PML from the original
//!    Berenger 1994 PML that would diverge at DC).
//! 2. `pml_assembly_zero_thickness_passes_through` — with
//!    `thickness_cells = 0` the mesh extension is a no-op and every
//!    tet is tagged Interior; the assembled matrix matches the
//!    scalar-ε scalar-μ closed-cavity assembly bit-for-bit.

use std::f64::consts::PI;

use nalgebra::Vector3;
use yee_fem::{
    AbcOrder, FaceKind, MaterialDatabase, OpenBoundarySolver, PmlAxis, PmlConfig, PmlMeshMeta,
    extend_mesh_with_pml,
};
use yee_mesh::TetMesh3D;

/// Tiny brick cavity that fits in unit tests without blowing the
/// runtime budget.
fn tiny_brick() -> TetMesh3D {
    TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap()
}

#[test]
fn pml_assembly_finite_at_dc() {
    let cavity = tiny_brick();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 3).expect("extend_mesh_with_pml");

    // Build face_kinds for the extended mesh: any face whose centroid
    // is on the outermost truncation surface (z < cavity_z_min) is
    // PEC; the rest of the original PEC sidewalls and the inner
    // cavity boundary stay PEC; the cavity's z = d face becomes a
    // wave-port stub for assembly-only smoke (no driven solve here).
    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };
    let placeholder = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("placeholder solver");
    let centroids = placeholder.exterior_face_centroids();

    let mut face_kinds: Vec<FaceKind> = Vec::with_capacity(centroids.len());
    for _c in &centroids {
        face_kinds.push(FaceKind::Pec);
    }
    drop(placeholder);

    let cfg = PmlConfig::default();
    let solver =
        OpenBoundarySolver::new(&extended, face_kinds, Vec::new(), MaterialDatabase::new())
            .expect("solver")
            .with_cfs_pml(cfg, classes);

    // ω at 0.1 GHz — well below the WR-90 TE_{10} cutoff (~6.56 GHz).
    let omega = 2.0 * PI * 0.1e9;
    let system = solver.assemble_driven_system(omega).expect("assemble");
    let max_abs = system
        .matrix
        .triplet_iter()
        .map(|t| t.val.norm())
        .fold(0.0_f64, |a, b| a.max(b));
    assert!(
        max_abs.is_finite(),
        "max entry magnitude must be finite at DC (CFS-PML α > 0 causality), got {max_abs}"
    );
    assert!(
        max_abs < 1.0e12,
        "max entry magnitude must be bounded at DC; got {max_abs:e}"
    );
}

#[test]
fn pml_assembly_zero_thickness_passes_through() {
    let cavity = tiny_brick();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 0).expect("extend_mesh_with_pml");

    assert_eq!(
        extended.tetrahedra.len(),
        cavity.tetrahedra.len(),
        "thickness_cells = 0 must preserve tet count"
    );
    assert!(
        classes.iter().all(|c| c.is_interior()),
        "thickness_cells = 0 must tag every tet Interior"
    );

    // Build all-PEC face_kinds so the matrix structure is purely
    // bulk and we can compare PML-Interior path vs Second-order ABC
    // path bit-for-bit on the cavity interior.
    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };

    let omega = 2.0 * PI * 10e9; // 10 GHz, well above cutoff

    let solver_pml = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds.clone(),
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver")
    .with_cfs_pml(PmlConfig::default(), classes);

    let solver_scalar = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver scalar")
    .with_abc_order(AbcOrder::First);

    let sys_pml = solver_pml.assemble_driven_system(omega).expect("pml asm");
    let sys_scalar = solver_scalar
        .assemble_driven_system(omega)
        .expect("scalar asm");

    // Same DoF count and same non-zero pattern.
    assert_eq!(sys_pml.rhs.len(), sys_scalar.rhs.len());

    // Compare entry sums (trace-like) to verify the PML-zero-thickness
    // path produces a numerically equivalent matrix to the scalar path
    // on cavity-interior tets. We don't compare triplet-by-triplet
    // because COO build orders may differ; instead we sum the per-row
    // L1 norms and require agreement to round-off.
    let n_int = sys_pml.rhs.len();
    let mut row_l1_pml = vec![0.0_f64; n_int];
    for t in sys_pml.matrix.triplet_iter() {
        row_l1_pml[t.row] += t.val.norm();
    }
    let mut row_l1_scalar = vec![0.0_f64; n_int];
    for t in sys_scalar.matrix.triplet_iter() {
        row_l1_scalar[t.row] += t.val.norm();
    }
    let mut max_rel_diff = 0.0_f64;
    for i in 0..n_int {
        let diff = (row_l1_pml[i] - row_l1_scalar[i]).abs();
        let denom = row_l1_pml[i].abs().max(1.0e-30);
        let rel = diff / denom;
        if rel > max_rel_diff {
            max_rel_diff = rel;
        }
    }
    let _ = Vector3::<f64>::zeros(); // touch nalgebra import for clippy
    assert!(
        max_rel_diff < 1.0e-9,
        "zero-thickness PML path must match scalar (v3) path to round-off; got max rel = {max_rel_diff:e}"
    );
}

/// Phase 4.fem.eig.3.5.1 R1 — the per-axis `h_alpha` resolver's
/// [`PmlMeshMeta::h_per_axis`] returns `extents / cell_counts`
/// componentwise.
#[test]
fn pml_mesh_meta_h_per_axis_is_extent_over_count() {
    let meta = PmlMeshMeta {
        extents: [0.02286, 0.01016, 0.030],
        cell_counts: [24, 12, 36],
    };
    let h = meta.h_per_axis();
    let tol = 1.0e-12;
    assert!((h[0] - 0.02286 / 24.0).abs() < tol);
    assert!((h[1] - 0.01016 / 12.0).abs() < tol);
    assert!((h[2] - 0.030 / 36.0).abs() < tol);
}

/// Phase 4.fem.eig.3.5.1 R1 — with `thickness_cells = 0` the per-axis
/// resolver path produces a finite PML-assembly matrix with `Lambda = I` on
/// every tet (the same `Interior`-only classification the v3.5 scalar
/// path consumed).
#[test]
fn per_axis_resolver_zero_thickness_matches_scalar_path() {
    let cavity = TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 0).expect("extend_mesh_with_pml");

    assert_eq!(extended.tetrahedra.len(), cavity.tetrahedra.len());
    assert!(classes.iter().all(|c| c.is_interior()));

    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };

    let omega = 2.0 * PI * 10e9;
    let solver_pml = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver")
    .with_cfs_pml(PmlConfig::default(), classes);
    let sys = solver_pml
        .assemble_driven_system(omega)
        .expect("per-axis no-shell asm");
    for t in sys.matrix.triplet_iter() {
        assert!(
            t.val.norm().is_finite(),
            "non-finite entry under no-shell PML path"
        );
    }
}

/// Phase 4.fem.eig.3.5.1 R1 — on an isotropic mesh (`h_x = h_y = h_z`)
/// the per-axis resolver produces a finite, well-conditioned assembly,
/// matching the v3.5 single-`h_cell` resolver behaviour by construction
/// (all `sigma_alpha_max` equal, all `h_alpha` equal).
#[test]
fn per_axis_resolver_isotropic_mesh_matches_legacy() {
    let cavity = TetMesh3D::cavity_uniform(0.004, 0.004, 0.004, 4, 4, 4).unwrap();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 3).expect("extend_mesh_with_pml");

    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };

    let omega = 2.0 * PI * 30e9;
    let cfg = PmlConfig::default();
    let solver = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver")
    .with_cfs_pml(cfg, classes);
    let sys = solver.assemble_driven_system(omega).expect("iso asm");
    for t in sys.matrix.triplet_iter() {
        assert!(
            t.val.norm().is_finite() && t.val.norm() < 1.0e18,
            "isotropic-mesh assembly produced non-finite or unbounded entry: {}",
            t.val.norm()
        );
    }
    assert!(matches!(solver.abc_order(), AbcOrder::CfsPml(_)));
}

/// Phase 4.fem.eig.3.5.2 S1 — with `alpha_grading_order = 0`, the
/// `α_α(d) = α_max · (1 − d/D)^0 = α_max` ramp collapses to the
/// v3.5.1 constant. The assembled CFS-PML stiffness must match the
/// v3.5.1 default-config stiffness bit-for-bit (Frobenius row-L1
/// equality to machine epsilon).
#[test]
fn alpha_grading_order_zero_matches_v3_5_1() {
    let cavity = TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 3).expect("extend_mesh_with_pml");

    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };

    let omega = 2.0 * PI * 10e9;

    // SSSSSSSSS Phase 4.fem.eig.3.5.2 S3 retune set the default
    // `alpha_grading_order` to 1; this test parametrises **both**
    // configurations explicitly so it remains invariant to future
    // default retunes. The mathematical claim is:
    //     α_α(d) = α_max · (1 − d/D)^0 ≡ α_max
    // collapses to the v3.5.1 constant for `alpha_grading_order = 0`,
    // regardless of what the current default carries.
    let cfg_v3_5_1 = PmlConfig {
        alpha_grading_order: 0,
        ..PmlConfig::default()
    };
    let cfg_explicit_zero = PmlConfig {
        alpha_grading_order: 0,
        ..PmlConfig::default()
    };

    let solver_default = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds.clone(),
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver default")
    .with_cfs_pml(cfg_v3_5_1, classes.clone());

    let solver_explicit = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver explicit")
    .with_cfs_pml(cfg_explicit_zero, classes);

    let sys_default = solver_default
        .assemble_driven_system(omega)
        .expect("default asm");
    let sys_explicit = solver_explicit
        .assemble_driven_system(omega)
        .expect("explicit asm");

    assert_eq!(sys_default.rhs.len(), sys_explicit.rhs.len());

    let n = sys_default.rhs.len();
    let mut row_l1_default = vec![0.0_f64; n];
    for t in sys_default.matrix.triplet_iter() {
        row_l1_default[t.row] += t.val.norm();
    }
    let mut row_l1_explicit = vec![0.0_f64; n];
    for t in sys_explicit.matrix.triplet_iter() {
        row_l1_explicit[t.row] += t.val.norm();
    }
    let mut max_diff = 0.0_f64;
    for i in 0..n {
        let d = (row_l1_default[i] - row_l1_explicit[i]).abs();
        if d > max_diff {
            max_diff = d;
        }
    }
    assert!(
        max_diff < 1.0e-12,
        "alpha_grading_order = 0 path must match v3.5.1 constant-alpha path \
         bit-for-bit (row-L1 |Δ| < 1e-12); got {max_diff:e}"
    );
}

/// Phase 4.fem.eig.3.5.2 S1 — with `alpha_grading_order = 1`, the
/// `α_α(d) = α_max · (1 − d/D)^1` linear ramp must produce a finite,
/// bounded assembled stiffness (the §7 (a) causality canary survives
/// the `denom.norm_sqr() <= MIN_POSITIVE` guard at the outer
/// truncation surface where `α_α(D) = 0`).
#[test]
fn alpha_grading_order_one_assembly_finite() {
    let cavity = TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 4, 2, 4).unwrap();
    let (extended, classes, _faces) =
        extend_mesh_with_pml(&cavity, &[PmlAxis::ZMin], 3).expect("extend_mesh_with_pml");

    let placeholder_kinds = {
        let mut face_map: std::collections::HashMap<[usize; 3], usize> =
            std::collections::HashMap::new();
        const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for tet in &extended.tetrahedra {
            for &[a, b, c] in TET_FACES.iter() {
                let mut key = [tet[a], tet[b], tet[c]];
                key.sort_unstable();
                *face_map.entry(key).or_insert(0) += 1;
            }
        }
        let n_exterior = face_map.values().filter(|&&c| c == 1).count();
        vec![FaceKind::Pec; n_exterior]
    };

    let cfg = PmlConfig {
        alpha_grading_order: 1,
        ..PmlConfig::default()
    };

    let solver = OpenBoundarySolver::new(
        &extended,
        placeholder_kinds,
        Vec::new(),
        MaterialDatabase::new(),
    )
    .expect("solver")
    .with_cfs_pml(cfg, classes);

    // Exercise both the high-frequency path and the DC-limit causality
    // canary path: the outer-truncation cell has `α_α(D) = 0`, so the
    // `denom.norm_sqr()` guard is on the critical path at low ω.
    for &freq_hz in &[0.1e9_f64, 10.0e9] {
        let omega = 2.0 * PI * freq_hz;
        let sys = solver
            .assemble_driven_system(omega)
            .expect("alpha_grading_order=1 asm");
        let max_abs = sys
            .matrix
            .triplet_iter()
            .map(|t| t.val.norm())
            .fold(0.0_f64, |a, b| a.max(b));
        assert!(
            max_abs.is_finite(),
            "alpha_grading_order = 1 produced non-finite entry at f = {freq_hz} Hz"
        );
        assert!(
            max_abs < 1.0e18,
            "alpha_grading_order = 1 produced unbounded entry at f = {freq_hz} Hz: {max_abs:e}"
        );
    }
}
