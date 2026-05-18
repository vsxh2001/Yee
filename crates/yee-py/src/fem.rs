//! `yee.fem` submodule — 3-D FEM eigenmode binding.
//!
//! Python wrapper for the Phase 4.fem.eig.0 walking-skeleton pipeline:
//! a one-shot `solve_cavity(a, b, d, nx, ny, nz, num_eigs)` that builds a
//! Kuhn 6-tet uniform mesh of a closed PEC rectangular cavity via
//! [`yee_mesh::TetMesh3D::cavity_uniform`], assembles the Nedelec
//! curl-curl pencil via [`yee_fem::FemEigenAssembly::new_free_space`],
//! and solves `K e = k² M e` for the lowest `num_eigs` physical modes via
//! shift-invert deflated inverse-power iteration
//! ([`yee_fem::InverseIterEigen`]).
//!
//! Mirrors the Rust-side `yee_validation::run_fem_eig_001_rectangular_cavity`
//! driver (Phase 4 T7) end-to-end; the shift heuristic is the same
//! `σ = 2.5 · k₀_TE101²` value vetted on the (8, 6, 10) mesh in T7's
//! production-gate run. See the documented limitation in the validation
//! driver's source for why this shift sits above the lowest physical
//! mode; the [`yee_fem::SparseEigen`] trait keeps the choice behind an
//! abstraction so a future LOBPCG / ARPACK swap fixes the dependency in
//! one PR.
//!
//! ## Submodule registration
//!
//! Mirrors the `yee.touchstone` / `yee.eigensolver` pattern: we insert
//! the module into `sys.modules` from `lib.rs` so that
//! `from yee.fem import solve_cavity` works in addition to the
//! attribute-access form `yee.fem.solve_cavity`.

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use yee_core::units::C0;
use yee_fem::{FemEigenAssembly, InverseIterEigen, SparseEigen};
use yee_mesh::TetMesh3D;

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve_cavity, m)?)?;
    Ok(())
}

/// Solve the lowest-`num_eigs` resonant modes of a closed PEC
/// rectangular cavity.
///
/// Parameters
/// ----------
/// a, b, d : float
///     Cavity extents along the `x`, `y`, `z` axes, in metres. Must be
///     strictly positive and finite.
/// nx, ny, nz : int
///     Brick subdivisions along the `x`, `y`, `z` axes. Each axis-aligned
///     brick is decomposed into 6 tetrahedra (Kuhn decomposition), so the
///     total tet count is `nx * ny * nz * 6`. Must be ≥ 1.
/// num_eigs : int, optional
///     Number of lowest physical eigenvalues to return. Default 10.
///
/// Returns
/// -------
/// frequencies : list[float]
///     Resonant frequencies in Hz, length `num_eigs`, sorted ascending.
/// mode_coefficients : numpy.ndarray
///     `(num_interior_edges, num_eigs)` `float64` array of eigenvectors
///     on the interior-edge basis (PEC-eliminated). Column `n` is the
///     mode for `frequencies[n]`.
///
/// Notes
/// -----
/// The shift `σ = 2.5 · k₀_TE101²` is computed from the analytic
/// `TE_{101}` wavenumber `k₀_TE101² = (π/a)² + (π/d)²`, matching the
/// Phase 4 T7 production-gate driver. The current
/// `yee_fem::InverseIterEigen` is the Phase 4 T5 escape-hatch
/// implementation; a future LOBPCG / ARPACK swap on the
/// `yee_fem::SparseEigen` trait will remove the shift-sensitivity caveat
/// documented in `yee_validation::run_fem_eig_001_rectangular_cavity`.
///
/// Raises
/// ------
/// ValueError
///     If any extent is non-positive / non-finite, any subdivision is
///     zero, or `num_eigs == 0`.
/// RuntimeError
///     If the inner sparse LU or any eigenmode fails to converge in the
///     solver's iteration budget.
#[pyfunction]
#[pyo3(signature = (a, b, d, nx, ny, nz, num_eigs = 10))]
#[allow(clippy::too_many_arguments)]
fn solve_cavity(
    py: Python<'_>,
    a: f64,
    b: f64,
    d: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    num_eigs: usize,
) -> PyResult<(Vec<f64>, Py<PyArray2<f64>>)> {
    // ---- 1. Build the cavity mesh -----------------------------------
    let mesh =
        TetMesh3D::cavity_uniform(a, b, d, nx, ny, nz).map_err(crate::errors::yee_mesh_to_py)?;

    // ---- 2. Assemble free-space K, M with PEC Dirichlet -------------
    let assembled = FemEigenAssembly::new_free_space(&mesh)
        .assemble()
        .map_err(crate::errors::yee_to_py)?;

    // ---- 3. Shift heuristic — same as the Phase 4 T7 driver ---------
    //
    // `solve_cavity` cannot know which returned eigenvalue is TE_{101}
    // ahead of time, so the shift uses the analytic TE_{101} wavenumber
    // for the generic rectangular cavity: `k₀² = (π/a)² + (π/d)²` (the
    // `(m, n, p) = (1, 0, 1)` mode is insensitive to `b`). The factor
    // `2.5` lifts `σ` above the gradient-kernel cluster at `k² = 0`
    // and above the lowest physical mode, per the documented limitation
    // of the inverse-power escape-hatch solver — see
    // `yee_validation::run_fem_eig_001_rectangular_cavity` source.
    let pi = std::f64::consts::PI;
    let k0_te101_sq = (pi / a).powi(2) + (pi / d).powi(2);
    let sigma = 2.5 * k0_te101_sq;

    // ---- 4. Solve K e = k² M e --------------------------------------
    let pairs = InverseIterEigen::default()
        .solve(&assembled.k, &assembled.m, num_eigs, sigma)
        .map_err(crate::errors::yee_to_py)?;

    // ---- 5. Convert k² → frequency, sort with paired eigenvectors ---
    //
    // `pairs.k` is sorted ascending by `k²` already (see
    // `InverseIterEigen::solve` post-processing), but we sort defensively
    // here so the returned `frequencies` list contract is
    // implementation-agnostic. The eigenvectors are reordered in
    // lockstep so column `n` of `mode_coefficients` always matches
    // `frequencies[n]`.
    let n_modes = pairs.k.len();
    let n_dofs = pairs.e.nrows();
    let mut order: Vec<usize> = (0..n_modes).collect();
    order.sort_by(|&i, &j| pairs.k[i].total_cmp(&pairs.k[j]));

    let two_pi = 2.0 * pi;
    let frequencies: Vec<f64> = order
        .iter()
        .map(|&i| {
            let k_sq = pairs.k[i];
            let k_abs = if k_sq > 0.0 { k_sq.sqrt() } else { 0.0 };
            C0 * k_abs / two_pi
        })
        .collect();

    // Build a `(n_dofs, n_modes)` row-major buffer in the permuted
    // column order.
    let mut buf: Vec<f64> = Vec::with_capacity(n_dofs * n_modes);
    for row in 0..n_dofs {
        for &col in &order {
            buf.push(pairs.e[(row, col)]);
        }
    }
    let modes = numpy::ndarray::Array2::from_shape_vec((n_dofs, n_modes), buf)
        .expect("buffer length matches (n_dofs, n_modes) by construction")
        .into_pyarray(py)
        .unbind();

    Ok((frequencies, modes))
}
