//! Integration tests for the Phase 4 T4 global FEM assembly:
//! dimensions, symmetry, positive diagonal, and the free-space default
//! constructor.
//!
//! The geometry fixtures are kept tiny and hand-rolled — a single free
//! tet, a pair of tets sharing one face, and a four-tet "umbrella"
//! whose four ring tets share one interior edge — so the boundary-edge
//! classifier and orientation-aware scatter can be validated by direct
//! counting before any production-scale mesh is involved.

use nalgebra::Vector3;
use nalgebra_sparse::csr::CsrMatrix;
use yee_fem::FemEigenAssembly;
use yee_mesh::TetMesh3D;

/// Reference unit tet `[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]`.
fn single_tet() -> TetMesh3D {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let tetrahedra = vec![[0, 1, 2, 3]];
    TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
}

/// Two tets sharing a triangular face (face (0, 1, 2) on the z = 0 plane;
/// v3 is above the plane, v4 is below).
fn two_tets_shared_face() -> TetMesh3D {
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, -1.0),
    ];
    let tetrahedra = vec![[0, 1, 2, 3], [0, 1, 2, 4]];
    TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
}

/// Four-tet "umbrella" fixture sharing a single interior edge along
/// the z-axis between `v0 = (0, 0, 0)` and `v1 = (0, 0, 1)`. The four
/// outer ring vertices `v2..v5` are placed at the cardinal compass
/// points on the `z = 1/2` plane so every tet has positive signed
/// volume. The shared edge `(v0, v1)` is the only edge of the mesh
/// whose every incident face is shared by two tets — so it is the only
/// **interior** edge.
fn umbrella_four_tet() -> TetMesh3D {
    // Ring vertices on z = 1/2 at the four cardinal points.
    let vertices = vec![
        Vector3::new(0.0, 0.0, 0.0),  // v0 — bottom of shared edge
        Vector3::new(0.0, 0.0, 1.0),  // v1 — top of shared edge
        Vector3::new(1.0, 0.0, 0.5),  // v2 — east ring vertex
        Vector3::new(0.0, 1.0, 0.5),  // v3 — north
        Vector3::new(-1.0, 0.0, 0.5), // v4 — west
        Vector3::new(0.0, -1.0, 0.5), // v5 — south
    ];
    // Four tets walking CCW around the shared edge (v0 -> v1). Each
    // tet shares face `(0, 1, v_k)` with the next; the four ring tets
    // therefore close into a full revolution. Orientation is fixed up
    // silently by `TetMesh3D::new` if any of these turn out CW.
    let tetrahedra = vec![[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 5], [0, 1, 5, 2]];
    TetMesh3D::new(vertices, tetrahedra, None, None).unwrap()
}

/// Dense `(K - K^T)` Frobenius-style norm: `sqrt(Σ_{i,j} (K[i,j] - K[j,i])²)`.
/// Implemented dense-only because the fixtures here are at most a handful
/// of DoFs, so a dense scan is the simplest path.
fn anti_symmetry_norm(k: &CsrMatrix<f64>) -> f64 {
    let n = k.nrows();
    assert_eq!(n, k.ncols(), "expected square matrix");
    let mut dense = vec![vec![0.0f64; n]; n];
    for (i, j, v) in k.triplet_iter() {
        dense[i][j] += v;
    }
    let mut sum_sq = 0.0;
    #[allow(clippy::needless_range_loop)]
    for i in 0..n {
        for j in 0..n {
            let d = dense[i][j] - dense[j][i];
            sum_sq += d * d;
        }
    }
    sum_sq.sqrt()
}

/// Extract the diagonal of a CSR matrix as a dense vector.
fn diagonal(k: &CsrMatrix<f64>) -> Vec<f64> {
    let n = k.nrows();
    let mut diag = vec![0.0f64; n];
    for (i, j, v) in k.triplet_iter() {
        if i == j {
            diag[i] += v;
        }
    }
    diag
}

#[test]
fn assembled_k_dimensions_match_interior_edge_count_single_tet() {
    // A single free tet has 6 edges, all of them boundary (every face
    // belongs to exactly one tet). After PEC Dirichlet elimination the
    // interior-DoF basis is empty, so K and M are both 0×0.
    let mesh = single_tet();
    let asm = FemEigenAssembly::new_free_space(&mesh).assemble().unwrap();
    assert_eq!(asm.interior_edges.len(), 0);
    assert_eq!(asm.k.nrows(), 0);
    assert_eq!(asm.k.ncols(), 0);
    assert_eq!(asm.m.nrows(), 0);
    assert_eq!(asm.m.ncols(), 0);
}

#[test]
fn two_tet_shared_face_interior_edge_count_matches_classifier() {
    // Two tets sharing the face (0, 1, 2) on z = 0:
    //   * 9 total edges (3 shared-face edges + 6 tip edges to v3 / v4).
    //   * Every shared-face edge still lies on a *non-shared* face of one
    //     of the tets (the three other tet faces touch only one tet), so
    //     every edge of the mesh is boundary by the
    //     `face-touches-< 2-tets` classifier.
    //   * Consequently the interior DoF count is 0, and both K and M
    //     come out 0×0. The brief explicitly notes this match-the-
    //     classifier path; document the assertion accordingly.
    let mesh = two_tets_shared_face();
    let asm = FemEigenAssembly::new_free_space(&mesh).assemble().unwrap();
    assert_eq!(
        asm.interior_edges.len(),
        0,
        "every edge of the two-tet shared-face fixture lies on at least one boundary face"
    );
    assert_eq!(asm.k.nrows(), 0);
    assert_eq!(asm.m.nrows(), 0);
}

#[test]
fn assembled_k_is_symmetric_on_umbrella_fixture() {
    // The four-tet umbrella has exactly one interior edge — the
    // shared (v0, v1) z-axis edge — because every face incident on it
    // is shared by exactly two tets. K must be exactly symmetric
    // (real entries; the orientation-sign scatter cannot introduce
    // anti-symmetry).
    let mesh = umbrella_four_tet();
    let asm = FemEigenAssembly::new_free_space(&mesh).assemble().unwrap();
    // The umbrella fixture is designed so the shared z-axis edge is
    // the unique interior DoF after PEC elimination. Assert non-empty
    // first so the test fails loudly if the fixture geometry changes.
    assert!(
        !asm.interior_edges.is_empty(),
        "umbrella fixture should yield at least one interior edge after PEC elimination"
    );
    let asym = anti_symmetry_norm(&asm.k);
    assert!(
        asym < 1e-10,
        "K is not symmetric on the umbrella fixture: ||K - K^T||_F = {asym}"
    );
}

#[test]
fn assembled_m_is_symmetric_and_positive_diagonal() {
    let mesh = umbrella_four_tet();
    let asm = FemEigenAssembly::new_free_space(&mesh).assemble().unwrap();
    assert!(!asm.interior_edges.is_empty());
    let asym = anti_symmetry_norm(&asm.m);
    assert!(
        asym < 1e-10,
        "M is not symmetric on the umbrella fixture: ||M - M^T||_F = {asym}"
    );
    // The vector mass is a Gram matrix of the Nedelec basis weighted by
    // ε_r > 0; every diagonal entry must be strictly positive.
    let diag = diagonal(&asm.m);
    for (i, &d) in diag.iter().enumerate() {
        assert!(
            d > 0.0,
            "M[{i},{i}] = {d} should be strictly positive (mass diagonals are ε_r-weighted ||N_i||² ≥ 0)"
        );
    }
}

#[test]
fn free_space_default_constructor_produces_nonempty_matrices() {
    // The convenience `new_free_space` constructor should set ε_r = μ_r
    // = 1 on every tet and produce the same K / M as an explicit
    // construction with uniform-1 vectors.
    let mesh = umbrella_four_tet();
    let asm_default = FemEigenAssembly::new_free_space(&mesh).assemble().unwrap();
    assert!(!asm_default.interior_edges.is_empty());
    assert!(asm_default.k.nnz() > 0);
    assert!(asm_default.m.nnz() > 0);

    let n_tets = mesh.tetrahedra.len();
    let asm_explicit = FemEigenAssembly::new(&mesh, vec![1.0; n_tets], vec![1.0; n_tets])
        .unwrap()
        .assemble()
        .unwrap();
    assert_eq!(asm_default.k.nrows(), asm_explicit.k.nrows());
    assert_eq!(asm_default.k.nnz(), asm_explicit.k.nnz());
    // Diagonals match exactly (both built from the same per-tet material
    // and the same scatter order).
    let d_default = diagonal(&asm_default.k);
    let d_explicit = diagonal(&asm_explicit.k);
    assert_eq!(d_default.len(), d_explicit.len());
    for (a, b) in d_default.iter().zip(d_explicit.iter()) {
        assert!(
            (a - b).abs() < 1e-15,
            "default and explicit free-space K diagonals must agree: {a} vs {b}"
        );
    }
}
