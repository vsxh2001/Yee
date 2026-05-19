//! Dense generalized-eigensolve fallback for the 2-D cross-section
//! eigenproblem.
//!
//! Solves `S x = k_c² T x` (real-symmetric for the lossless case),
//! where `S` is the curl-curl stiffness and `T` is the ε_r-weighted
//! Nedelec mass — both produced by [`super::assembly::assemble_transverse`].
//! The propagation constant `β` for a given operating frequency
//! follows from the dispersion relation `β² = k_0² − k_c²`.
//!
//! **Spurious-mode handling.** First-order Nedelec edge elements admit
//! a large gradient null-space (every nodal-Lagrange gradient is in
//! their span and lies in the kernel of `curl`); these modes show up
//! with `k_c² ≈ 0` and are filtered out by a threshold relative to the
//! largest eigenvalue. The physical dominant mode is the **smallest
//! strictly-positive** `k_c²` eigenvalue.
//!
//! **Numerical method.** Reduce `S x = k_c² T x` to a standard
//! symmetric problem `M y = k_c² y` via the Cholesky factor `T = L Lᵀ`,
//! then `M = L⁻¹ S L⁻ᵀ`. Solve `M` with [`nalgebra::SymmetricEigen`]
//! (tridiagonal QR). `O(n³)` flops; viable only for coarse cross
//! sections (≤ a few hundred DoF). Sparse shift-and-invert with
//! `arpack-rs` is Phase 1.3.1.1 step 4 and is escape-hatched away from
//! the lossless TE10 validation gate.
//!
//! **Real-arithmetic restriction.** Lossless inputs only for Phase
//! 1.3.1.1 step 3. The [`super::assembly`] module produces
//! `DMatrix<Complex64>` to keep the API future-proof for lossy fills;
//! this routine takes the real parts and surfaces an
//! `Error::Unimplemented` if the imaginary parts exceed a small
//! threshold relative to the real norm. Complex extension lands in
//! Phase 1.3.1.2.

use nalgebra::{DMatrix, SymmetricEigen};
use num_complex::Complex64;

use super::assembly::AssembledTransverse;

/// Solved-eigensolution payload returned by [`solve_dense`].
pub(crate) struct EigenSolution {
    /// `β² = k_0² − k_c²` for the dominant propagating mode at the
    /// supplied frequency. Stored as `Complex64` to keep the API
    /// future-proof for the lossy / complex-symmetric path; the
    /// Phase 1.3.1.1 step 3 path always returns a real value.
    pub beta_sq: Complex64,
    /// Eigenvector for the dominant mode in the **interior-edge DoF**
    /// ordering of [`AssembledTransverse::interior_to_global`]. Length
    /// equals the interior-edge DoF count `n`. Real-valued on the
    /// lossless path (stored as `Complex64` to mirror the API contract
    /// of [`Self::beta_sq`]).
    ///
    /// **Sign convention.** The eigenvector returned by
    /// [`nalgebra::SymmetricEigen`] has an arbitrary global sign. To
    /// fix it deterministically, [`solve_dense`] post-processes so the
    /// largest-magnitude component is positive (`v[argmax |v|] > 0`).
    /// Callers that need a physically-meaningful sign (e.g. picking
    /// the positive-going wave for a wave-port RHS) should additionally
    /// renormalize against a known reference point — see
    /// [`crate::ports::NumericalCrossSection::e_tangential_at`] which
    /// fixes the sign by `E_y > 0` at the cross-section centroid.
    pub eigenvector: Vec<Complex64>,
}

/// Solve `S x = k_c² T x` densely on the lossless / real-symmetric path
/// and return `β² = k_0² − k_c²` for the dominant propagating mode at
/// `freq_hz`.
///
/// **Size limit.** The implementation builds dense `n×n` matrices and
/// runs an `O(n³)` Cholesky + symmetric-tridiagonal QR. Callers must
/// keep `n` below a few hundred. The eigensolver's only validation
/// case (WR-90 TE10 on a ~72-triangle mesh) lands at `n ≈ 84` interior
/// edges, well inside this envelope.
pub(crate) fn solve_dense(
    asm: &AssembledTransverse,
    freq_hz: f64,
) -> Result<EigenSolution, yee_core::Error> {
    let n = asm.s.nrows();
    if n == 0 {
        return Err(yee_core::Error::Numerical(
            "eigensolver: empty interior-edge DoF set (mesh has no interior edges?)".into(),
        ));
    }

    // Reject lossy-input case: this solver is lossless-only for Phase
    // 1.3.1.1 step 3. The `Complex64` storage in `assembly` is kept so
    // the API extends cleanly when the lossy / complex-symmetric path
    // lands.
    let im_norm = asm
        .s
        .iter()
        .chain(asm.t.iter())
        .map(|z| z.im.abs())
        .fold(0.0_f64, f64::max);
    let re_norm = asm
        .s
        .iter()
        .chain(asm.t.iter())
        .map(|z| z.re.abs())
        .fold(0.0_f64, f64::max);
    if im_norm > 1e-9 * re_norm.max(1.0) {
        return Err(yee_core::Error::Unimplemented(
            "eigensolver: complex (lossy) ε_r / μ_r is Phase 1.3.1.2; current path is lossless only",
        ));
    }

    let s_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.s[(i, j)].re);
    let t_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.t[(i, j)].re);

    // Cholesky factor T = L Lᵀ. For the lossless Nedelec mass matrix
    // with ε_r > 0 on the interior-edge DoF set this is SPD by
    // construction.
    let chol_t = t_re.clone().cholesky().ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: Nedelec mass matrix T is not SPD on the interior-edge DoF set".into(),
        )
    })?;
    let l = chol_t.l();
    // Compute M = L⁻¹ S L⁻ᵀ via two triangular solves preserving
    // symmetry. Step 1: solve L Y = S. Step 2: solve L Z = Yᵀ → M = Zᵀ.
    let y = l.clone().solve_lower_triangular(&s_re).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: L Y = S solve failed (L from Cholesky should be non-singular)".into(),
        )
    })?;
    let y_t = y.transpose();
    let z = l.solve_lower_triangular(&y_t).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: L Z = Yᵀ solve failed (L from Cholesky should be non-singular)".into(),
        )
    })?;
    let m = z.transpose();

    // Symmetrize explicitly to suppress floating-point drift and feed
    // the symmetric QR path, which is bulletproof for symmetric pencils.
    let m_sym = 0.5 * (&m + m.transpose());

    let eig = SymmetricEigen::new(m_sym);
    // SymmetricEigen returns real eigenvalues (T::RealField = f64 here).
    // Filter:
    //   - reject `k_c² <= 0` (gradient null-space lives here),
    //   - reject the spurious-mode floor `k_c² < spurious_floor` where
    //     `spurious_floor = max_eig · 1e-6` catches the
    //     numerically-non-zero gradient eigenvalues,
    // then take the **smallest** strictly-positive eigenvalue. Physical
    // dominant mode = lowest cutoff.
    let max_eig = eig.eigenvalues.iter().cloned().fold(0.0_f64, f64::max);
    let spurious_floor = max_eig * 1e-6;

    // Track the (eigenvalue, eigenvector-column-index) of the dominant
    // physical mode so we can recover its eigenvector after the
    // standard-form back-transform `x = L⁻ᵀ y`.
    let mut dominant: Option<(f64, usize)> = None;
    for (col, &lam) in eig.eigenvalues.iter().enumerate() {
        if !lam.is_finite() || lam <= spurious_floor {
            continue;
        }
        dominant = Some(match dominant {
            None => (lam, col),
            Some((curr, curr_col)) => {
                if lam < curr {
                    (lam, col)
                } else {
                    (curr, curr_col)
                }
            }
        });
    }
    let (k_c_sq, dom_col) = dominant.ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: no strictly-positive k_c² eigenvalue above the spurious-mode floor"
                .into(),
        )
    })?;

    // β² = k_0² − k_c². ε_r is folded into T (and therefore into k_c²);
    // the relation as written is correct for real, uniform ε_r and
    // exact for the WR-90 air-filled validation case. Heterogeneous /
    // lossy ε_r is Phase 1.3.1.2.
    let omega = std::f64::consts::TAU * freq_hz;
    let k0 = omega / yee_core::units::C0;
    let beta_sq_re = k0 * k0 - k_c_sq;
    if beta_sq_re <= 0.0 {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver: mode is evanescent at {freq_hz} Hz (k_c² = {k_c_sq}, k_0² = {})",
            k0 * k0
        )));
    }

    // Recover the eigenvector in the **original** (T-weighted) basis.
    // SymmetricEigen solves `M y = λ y` with `M = L⁻¹ S L⁻ᵀ`; the
    // generalized-problem eigenvector satisfies `S x = λ T x` with
    // `x = L⁻ᵀ y`. So back-transform with one upper-triangular solve.
    let y_dom = eig.eigenvectors.column(dom_col).clone_owned();
    let l_t = l.transpose();
    let x_dom = l_t.solve_upper_triangular(&y_dom).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: Lᵀ x = y back-transform failed (L from Cholesky should be non-singular)"
                .into(),
        )
    })?;

    // Fix the global sign deterministically: largest-magnitude
    // component positive. Downstream consumers (e.g.
    // `NumericalCrossSection::e_tangential_at`) renormalize against a
    // physical reference point.
    let argmax = x_dom
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.abs()
                .partial_cmp(&b.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    let sign = if x_dom[argmax] < 0.0 { -1.0 } else { 1.0 };
    let eigenvector: Vec<Complex64> = x_dom
        .iter()
        .map(|&v| Complex64::new(sign * v, 0.0))
        .collect();

    Ok(EigenSolution {
        beta_sq: Complex64::new(beta_sq_re, 0.0),
        eigenvector,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigensolver::{assembly::assemble_transverse, mesh::EdgeTable};
    use std::collections::HashMap;
    use yee_mesh::TriMesh2D;

    /// Build a structured `nx × ny` quad-grid of CCW triangles spanning
    /// `[0, a] × [0, b]`. Each quad cell is split along the
    /// `(low-x, low-y) → (high-x, high-y)` diagonal into two triangles.
    fn rectangular_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
        let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut triangles = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let v00 = idx(i, j);
                let v10 = idx(i + 1, j);
                let v11 = idx(i + 1, j + 1);
                let v01 = idx(i, j + 1);
                triangles.push([v00, v10, v11]);
                triangles.push([v00, v11, v01]);
            }
        }
        TriMesh2D::new(vertices, triangles, None, None).unwrap()
    }

    #[test]
    fn coarse_wr90_te10_sanity() {
        // Sanity smoke test: on a coarse WR-90 mesh, the dominant β² at
        // 10 GHz should be positive and within an order of magnitude of
        // the analytic value. The tight 1 % tolerance is the integration
        // gate in `tests/eigensolver_wr90.rs`.
        let a = 22.86e-3;
        let b = 10.16e-3;
        let mesh = rectangular_mesh(a, b, 4, 4);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_transverse(&mesh, &eps, &mu, &table);
        let sol = solve_dense(&asm, 10e9).unwrap();
        let beta_sq = sol.beta_sq.re;
        // Analytic TE10 β² = k_0² − (π/a)²
        let omega = std::f64::consts::TAU * 10e9;
        let k0 = omega / yee_core::units::C0;
        let kc = std::f64::consts::PI / a;
        let beta_sq_analytic = k0 * k0 - kc * kc;
        assert!(
            beta_sq > 0.0 && beta_sq.is_finite(),
            "β² = {beta_sq} not positive-finite"
        );
        // Within an order of magnitude — assembly correctness check.
        let ratio = beta_sq / beta_sq_analytic;
        assert!(
            ratio.is_finite() && ratio > 0.1 && ratio < 10.0,
            "β² = {beta_sq}, analytic = {beta_sq_analytic}: ratio {ratio} out of band"
        );
    }
}
