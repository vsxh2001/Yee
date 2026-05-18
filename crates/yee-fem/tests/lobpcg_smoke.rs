//! Smoke tests for [`yee_fem::solve::SparseEigen`] —
//! shift-invert sparse eigensolve regression tests.
//!
//! These tests exercise the public [`yee_fem::InverseIterEigen`]
//! implementation (Phase 4 T5; see crate docs for the LOBPCG → inverse-
//! iteration escape-hatch rationale) against analytically-known
//! eigenpairs. They are intentionally cheap (≤ 64 DoF) so they run in
//! milliseconds on the default `cargo test` flow — the full fem-eig-001
//! production gate lives in `crates/yee-validation/tests/`.
//!
//! Test inventory:
//!
//! 1. `recovers_smallest_eigenvalue_on_known_dense_pencil` — 4×4
//!    diagonal pencil with eigenvalues `{0.5, 1.2, 3.4, 7.8}`; shift
//!    `σ = 0.1`, num_eigs = 3; returned eigenvalues within `1e-6`
//!    of `{0.5, 1.2, 3.4}`, sorted ascending.
//! 2. `scaled_identity_pencil` — `K = αI`, `M = βI` → every
//!    eigenvalue is `α/β`; shift-invert converges in one iteration.
//! 3. `eigenvectors_M_orthogonal` — `e_i^T M e_j ≈ δ_{ij}` within
//!    `1e-8` for the returned eigenvector basis.
//! 4. `converges_within_max_iter_for_3d_laplacian` — 7-point
//!    Dirichlet Laplacian on a 4³ grid (64 DoFs); smallest
//!    eigenvalue recovered within `max_iter = 100`, `tol = 1e-8`.

use nalgebra_sparse::coo::CooMatrix;
use nalgebra_sparse::csr::CsrMatrix;
use yee_fem::{InverseIterEigen, SparseEigen};

/// Build a diagonal CSR matrix from a slice.
fn diag_csr(diag: &[f64]) -> CsrMatrix<f64> {
    let n = diag.len();
    let mut coo = CooMatrix::new(n, n);
    for (i, &d) in diag.iter().enumerate() {
        if d != 0.0 {
            coo.push(i, i, d);
        }
    }
    CsrMatrix::from(&coo)
}

/// Build a CSR from a dense row-major slice.
fn csr_from_dense(rows: usize, cols: usize, data: &[f64]) -> CsrMatrix<f64> {
    assert_eq!(data.len(), rows * cols);
    let mut coo = CooMatrix::new(rows, cols);
    for r in 0..rows {
        for c in 0..cols {
            let v = data[r * cols + c];
            if v != 0.0 {
                coo.push(r, c, v);
            }
        }
    }
    CsrMatrix::from(&coo)
}

/// Helper: y = A x for CSR A.
fn csr_matvec(a: &CsrMatrix<f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let mut y = vec![0.0f64; n];
    let row_offsets = a.row_offsets();
    let col_indices = a.col_indices();
    let values = a.values();
    for row in 0..n {
        let start = row_offsets[row];
        let end = row_offsets[row + 1];
        let mut sum = 0.0f64;
        for k in start..end {
            sum += values[k] * x[col_indices[k]];
        }
        y[row] = sum;
    }
    y
}

#[test]
fn recovers_smallest_eigenvalue_on_known_dense_pencil() {
    let lambdas = [0.5, 1.2, 3.4, 7.8];
    let k = diag_csr(&lambdas);
    let m = diag_csr(&[1.0; 4]);

    let solver = InverseIterEigen::new(1000, 1e-10);
    let pairs = solver.solve(&k, &m, 3, 0.1).expect("solve must succeed");

    assert_eq!(pairs.k.len(), 3);
    assert!(
        (pairs.k[0] - 0.5).abs() < 1e-6,
        "expected 0.5, got {}",
        pairs.k[0]
    );
    assert!(
        (pairs.k[1] - 1.2).abs() < 1e-6,
        "expected 1.2, got {}",
        pairs.k[1]
    );
    assert!(
        (pairs.k[2] - 3.4).abs() < 1e-6,
        "expected 3.4, got {}",
        pairs.k[2]
    );
    for w in pairs.k.windows(2) {
        assert!(w[0] <= w[1], "eigenvalues must be sorted ascending");
    }
}

#[test]
fn scaled_identity_pencil() {
    let alpha = 3.7;
    let beta = 1.5;
    let k = diag_csr(&[alpha; 5]);
    let m = diag_csr(&[beta; 5]);

    let solver = InverseIterEigen::new(50, 1e-10);
    let pairs = solver.solve(&k, &m, 3, 0.5).expect("solve must succeed");

    let expected = alpha / beta;
    for &k_sq in &pairs.k {
        assert!(
            (k_sq - expected).abs() < 1e-8,
            "expected {expected}, got {k_sq}"
        );
    }
}

#[test]
#[allow(non_snake_case)]
fn eigenvectors_M_orthogonal() {
    let lambdas = [0.5, 1.2, 3.4, 7.8];
    let k = diag_csr(&lambdas);
    let m = diag_csr(&[1.0; 4]);

    let solver = InverseIterEigen::new(1000, 1e-10);
    let pairs = solver.solve(&k, &m, 3, 0.1).expect("solve must succeed");

    let n = pairs.e.nrows();
    let ncols = pairs.e.ncols();
    for i in 0..ncols {
        for j in 0..ncols {
            let col_j: Vec<f64> = (0..n).map(|r| pairs.e[(r, j)]).collect();
            let mxj = csr_matvec(&m, &col_j);
            let acc: f64 = (0..n).map(|r| pairs.e[(r, i)] * mxj[r]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (acc - expected).abs() < 1e-8,
                "M-orthogonality failed at ({i},{j}): got {acc}, expected {expected}"
            );
        }
    }
}

#[test]
fn converges_within_max_iter_for_3d_laplacian() {
    // 7-point Dirichlet Laplacian on a 4×4×4 interior grid (64 DoFs).
    let nx = 4;
    let n = nx * nx * nx;
    let idx = |i: usize, j: usize, k: usize| i + nx * (j + nx * k);

    let mut k_dense = vec![0.0; n * n];
    for i in 0..nx {
        for j in 0..nx {
            for kk in 0..nx {
                let p = idx(i, j, kk);
                k_dense[p * n + p] = 6.0;
                if i + 1 < nx {
                    k_dense[p * n + idx(i + 1, j, kk)] = -1.0;
                }
                if i >= 1 {
                    k_dense[p * n + idx(i - 1, j, kk)] = -1.0;
                }
                if j + 1 < nx {
                    k_dense[p * n + idx(i, j + 1, kk)] = -1.0;
                }
                if j >= 1 {
                    k_dense[p * n + idx(i, j - 1, kk)] = -1.0;
                }
                if kk + 1 < nx {
                    k_dense[p * n + idx(i, j, kk + 1)] = -1.0;
                }
                if kk >= 1 {
                    k_dense[p * n + idx(i, j, kk - 1)] = -1.0;
                }
            }
        }
    }
    let k_mat = csr_from_dense(n, n, &k_dense);
    let m_mat = diag_csr(&vec![1.0; n]);

    // Smallest 5-point-Dirichlet Laplacian eigenvalue on `nx` interior
    // points (per dimension, 3-D summed):
    //   λ_111 = 6 − 6 cos(π/(nx+1))
    let pi = std::f64::consts::PI;
    let expected = 6.0 - 6.0 * (pi / (nx as f64 + 1.0)).cos();

    let solver = InverseIterEigen::new(100, 1e-8);
    let pairs = solver
        .solve(&k_mat, &m_mat, 1, 0.1)
        .expect("3-D Laplacian solve must succeed");

    assert!(
        (pairs.k[0] - expected).abs() < 1e-5,
        "smallest 3-D Laplacian eigenvalue: expected {expected}, got {}",
        pairs.k[0]
    );
}
