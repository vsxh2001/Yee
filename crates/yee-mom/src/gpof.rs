//! Generalized Pencil of Function (GPOF) fitting — Phase 1.1.1.0.
//!
//! Given complex samples `y_m = y(m · Δt)` of a function that is known
//! (or expected) to be expressible as a sum of `N` complex exponentials
//!
//! ```text
//!     y(t) = Σ_{n=1}^N  α_n · exp(β_n · t)
//! ```
//!
//! GPOF recovers the `(α_n, β_n)` pairs from `M ≥ 2N` uniformly-spaced
//! samples via a closed-form linear-algebra pipeline — no iteration,
//! no nonlinear solve. The algorithm is the matrix-pencil method of
//! Hua & Sarkar (1989); the present implementation follows the
//! Aksun-DCIM specialisation (1996) where the input samples come from
//! a spectral Green's function on a deformed integration contour and
//! the recovered `(α_n, β_n)` are turned into discrete complex image
//! coefficients.
//!
//! ## Numerical pipeline
//!
//! 1. Build two `(M − L) × L` Hankel matrices
//!    `Y_1[i, j] = y_{i + j}`, `Y_2[i, j] = y_{i + j + 1}` with
//!    `L = M − N` (so each block is `N × (M − N)` — the minimum that
//!    still resolves `N` poles from `M ≥ 2N` samples).
//! 2. SVD `Y_1 = U Σ V^H`; truncate to the top `N` components giving
//!    `U_N ∈ ℂ^{N × N}`, `Σ_N ∈ ℝ^{N × N}`, `V_N ∈ ℂ^{L × N}`.
//! 3. Form the `N × N` pencil matrix
//!    `Z = Σ_N^{-1} · U_N^H · Y_2 · V_N`. Its complex eigenvalues
//!    `z_n` are `exp(β_n · Δt)`.
//! 4. Solve the rectangular Vandermonde least-squares system
//!    `V_{m, n} = z_n^m`, `V · α = y` for the amplitudes `α_n` via
//!    an SVD-backed `solve_lstsq`.
//!
//! All four stages are pure linear algebra; the routine returns in
//! deterministic time with no convergence loop.
//!
//! ## References
//!
//! * Y. Hua and T. K. Sarkar, "Generalized pencil-of-function method
//!   for extracting poles of an EM system from its transient
//!   response," *IEEE Trans. Antennas Propag.*, vol. 37, no. 2,
//!   pp. 229–234, Feb 1989.
//! * M. I. Aksun, "A robust approach for the derivation of closed-form
//!   Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44,
//!   no. 5, pp. 651–658, May 1996.

#![allow(dead_code)]

use faer::Mat;
use faer::linalg::solvers::{SolveLstsq, Svd};
use num_complex::Complex64;

/// Errors that the GPOF fit can surface to its caller.
#[derive(Debug, thiserror::Error)]
pub enum GpofError {
    /// Not enough samples to resolve `N` poles. Requires `M ≥ 2 N`.
    #[error("GPOF requires M >= 2 N samples; got M={m}, N={n}")]
    TooFewSamples {
        /// Number of samples provided.
        m: usize,
        /// Requested number of poles.
        n: usize,
    },
    /// Requested zero poles. The trivial fit `y(t) = 0` is degenerate
    /// and the caller almost certainly does not want it.
    #[error("GPOF: n_poles must be >= 1; got 0")]
    ZeroPoles,
    /// faer's SVD or eigendecomposition step failed (effectively never
    /// in practice — the input matrices are small and dense).
    #[error("GPOF linear-algebra failure: {0}")]
    LinAlg(String),
}

/// Result of a successful GPOF fit. Each tuple is `(α_n, β_n)` such
/// that `y(t) ≈ Σ_n α_n · exp(β_n · t)`. The returned vector has
/// length `n_poles` requested by the caller.
pub type GpofPoles = Vec<(Complex64, Complex64)>;

/// Fit `y(m · dt) ≈ Σ_n α_n · exp(β_n · m · dt)` for `m = 0..M-1`.
///
/// `dt` is the sample spacing (the same units as the recovered `β_n`
/// inverses; if `t` carries no physical unit the recovery is in
/// "per-sample" units and the caller must rescale).
///
/// Returns `n_poles` pairs `(α_n, β_n)` of complex amplitudes and
/// exponential rates, or a [`GpofError`] if the input is malformed or
/// the linear-algebra steps failed.
///
/// Determinism: the routine is a fixed sequence of dense SVD,
/// eigendecomposition, and least-squares solves — no iteration, no
/// random initialisation, identical results across runs given
/// identical inputs.
pub fn gpof(samples: &[Complex64], dt: f64, n_poles: usize) -> Result<GpofPoles, GpofError> {
    if n_poles == 0 {
        return Err(GpofError::ZeroPoles);
    }
    let m = samples.len();
    if m < 2 * n_poles {
        return Err(GpofError::TooFewSamples { m, n: n_poles });
    }

    // Pencil parameter L. With L = M - N each Hankel block is
    // (M-L) × L = N × (M-N). Hua & Sarkar's noise-robustness analysis
    // favours L ≈ M/2; for the DCIM use case M = 2N so L = M/2 = N
    // — both bands collapse to the same value, matching the
    // minimum-data setup expected by the spec.
    let l = m - n_poles;
    let rows = m - l; // = n_poles
    let cols = l;

    // Build Y1 and Y2 (Hankel structure). faer is column-major so the
    // index ordering below does not matter functionally — what matters
    // is that Y1[i,j] = y_{i+j} and Y2[i,j] = y_{i+j+1}.
    let mut y1: Mat<Complex64> = Mat::zeros(rows, cols);
    let mut y2: Mat<Complex64> = Mat::zeros(rows, cols);
    for i in 0..rows {
        for j in 0..cols {
            y1[(i, j)] = samples[i + j];
            y2[(i, j)] = samples[i + j + 1];
        }
    }

    // Step 2: SVD Y1 (thin form). With rows=N and cols=M-N=N this is
    // a square N×N decomposition; we take the full thing.
    let svd = Svd::new(y1.as_ref()).map_err(|e| GpofError::LinAlg(format!("SVD(Y1): {e:?}")))?;
    let u = svd.U();
    let s = svd.S();
    let v = svd.V();

    // Truncate to top n_poles. With rows=cols=n_poles in the DCIM
    // baseline, the truncation is a no-op, but doing it explicitly
    // keeps the code correct for M > 2N cases too.
    let n = n_poles.min(rows).min(cols);

    // Form Z = Σ_N^{-1} · U_N^H · Y_2 · V_N ∈ ℂ^{N × N}.
    //
    // 1) tmp1 = U_N^H · Y_2   (n × cols)
    // 2) tmp2 = tmp1 · V_N    (n × n)
    // 3) Z[i,j] = tmp2[i,j] / σ_i
    let mut tmp1: Mat<Complex64> = Mat::zeros(n, cols);
    for i in 0..n {
        for j in 0..cols {
            let mut acc = Complex64::new(0.0, 0.0);
            for k in 0..rows {
                acc += u[(k, i)].conj() * y2[(k, j)];
            }
            tmp1[(i, j)] = acc;
        }
    }
    let mut z_mat: Mat<Complex64> = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let mut acc = Complex64::new(0.0, 0.0);
            for k in 0..cols {
                acc += tmp1[(i, k)] * v[(k, j)];
            }
            // Divide row i by σ_i (real positive). Guard against zero
            // singular values: a rank-deficient Y1 means the model has
            // fewer than n_poles independent modes — surface as a
            // numerical failure rather than propagate NaN.
            let sigma = s[i].re;
            if !sigma.is_finite() || sigma <= 0.0 {
                return Err(GpofError::LinAlg(format!(
                    "GPOF: degenerate singular value σ[{i}] = {sigma} — Y1 rank-deficient"
                )));
            }
            z_mat[(i, j)] = acc / Complex64::new(sigma, 0.0);
        }
    }

    // Step 3: eigenvalues of Z are exp(β_n · dt).
    let eig = faer::linalg::solvers::Eigen::new(z_mat.as_ref())
        .map_err(|e| GpofError::LinAlg(format!("Eigen(Z): {e:?}")))?;
    let s_diag = eig.S();

    let mut zs: Vec<Complex64> = (0..n).map(|i| s_diag[i]).collect();

    // β_n = ln(z_n) / dt. The principal branch of `ln` is used; this
    // is the conventional choice for GPOF and matches the way Aksun
    // 1996 unwraps the exponents into image locations. The branch
    // ambiguity is harmless for our DCIM input: the spectral function
    // decays along the integration contour, so the recovered β_n
    // necessarily have Re β_n < 0 and the principal-branch argument
    // is the physically-meaningful one.
    let betas: Vec<Complex64> = zs
        .iter()
        .map(|z| z.ln() / Complex64::new(dt, 0.0))
        .collect();

    // Step 4: solve Vandermonde least squares  V · α = y  with
    // V[m, n] = z_n^m for m = 0..M-1, n = 0..N-1. Vandermonde is
    // ill-conditioned for clustered `z_n`; an SVD-backed
    // `solve_lstsq` is the standard mitigation.
    let mut vmat: Mat<Complex64> = Mat::zeros(m, n);
    for row in 0..m {
        // z_n^row computed by repeated multiplication to avoid the
        // branch-cut ambiguity in `Complex64::powi`/`powc` for
        // |z| ≪ 1 (which happens for the higher-rate images).
        for col in 0..n {
            let mut acc = Complex64::new(1.0, 0.0);
            for _ in 0..row {
                acc *= zs[col];
            }
            vmat[(row, col)] = acc;
        }
    }
    let mut rhs: Mat<Complex64> = Mat::zeros(m, 1);
    for row in 0..m {
        rhs[(row, 0)] = samples[row];
    }
    let vsvd = Svd::new(vmat.as_ref()).map_err(|e| GpofError::LinAlg(format!("SVD(V): {e:?}")))?;
    vsvd.solve_lstsq_in_place(rhs.as_mut());
    let alphas: Vec<Complex64> = (0..n).map(|i| rhs[(i, 0)]).collect();

    // Suppress the unused-warning on `zs` by re-asserting it after
    // the borrow scope above; we use it only for the Vandermonde
    // construction, and the assignment is logically complete.
    zs.truncate(n);

    Ok(alphas.into_iter().zip(betas).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_pair_match(
        recovered: &[(Complex64, Complex64)],
        expected: &[(Complex64, Complex64)],
        tol: f64,
    ) -> bool {
        // Permutation-invariant comparison: for each expected pole,
        // find the recovered pole closest in (α, β) space and check
        // the joint relative error.
        if recovered.len() != expected.len() {
            return false;
        }
        let mut used = vec![false; recovered.len()];
        for (a_exp, b_exp) in expected {
            let mut best_err = f64::INFINITY;
            let mut best_idx = usize::MAX;
            for (i, (a_rec, b_rec)) in recovered.iter().enumerate() {
                if used[i] {
                    continue;
                }
                let err_a = (a_rec - a_exp).norm() / a_exp.norm().max(1e-30);
                let err_b = (b_rec - b_exp).norm() / b_exp.norm().max(1e-30);
                let err = err_a.max(err_b);
                if err < best_err {
                    best_err = err;
                    best_idx = i;
                }
            }
            if best_err > tol || best_idx == usize::MAX {
                eprintln!(
                    "GPOF mismatch: expected (α={:?}, β={:?}); best_err = {:.3e}",
                    a_exp, b_exp, best_err
                );
                return false;
            }
            used[best_idx] = true;
        }
        true
    }

    /// Synthetic three-image recovery: generate uniform samples from a
    /// known sum of three exponentials, fit GPOF, and check that the
    /// recovered `(α_n, β_n)` match the originals to within `1e-6`
    /// (relative). This is the spec DoD #5 gate.
    #[test]
    fn gpof_recovers_three_synthetic_poles() {
        let true_poles: Vec<(Complex64, Complex64)> = vec![
            (Complex64::new(1.0, 0.5), Complex64::new(-0.3, 1.2)),
            (Complex64::new(0.4, -0.2), Complex64::new(-1.1, 0.7)),
            (Complex64::new(2.0, 1.0), Complex64::new(-0.7, -1.5)),
        ];
        let n = true_poles.len();
        let m = 4 * n; // oversampled to keep the Vandermonde well-conditioned
        let dt = 0.1;

        let samples: Vec<Complex64> = (0..m)
            .map(|k| {
                let t = (k as f64) * dt;
                true_poles
                    .iter()
                    .map(|(a, b)| a * (b * Complex64::new(t, 0.0)).exp())
                    .sum()
            })
            .collect();

        let recovered = gpof(&samples, dt, n).expect("GPOF fit");
        assert!(
            approx_pair_match(&recovered, &true_poles, 1.0e-6),
            "recovered = {:?}\n expected = {:?}",
            recovered,
            true_poles
        );
    }

    /// Single-pole recovery to machine precision. The N=1 case
    /// exercises the same pipeline but degenerates to a 1×1 SVD plus a
    /// scalar eigenvalue — keep the path covered so the linear-algebra
    /// shortcut paths cannot regress unnoticed.
    #[test]
    fn gpof_recovers_one_pole() {
        let alpha = Complex64::new(1.5, -0.3);
        let beta = Complex64::new(-0.8, 0.9);
        let dt = 0.2;
        let m = 6;
        let samples: Vec<Complex64> = (0..m)
            .map(|k| alpha * (beta * Complex64::new((k as f64) * dt, 0.0)).exp())
            .collect();
        let recovered = gpof(&samples, dt, 1).expect("GPOF fit");
        assert_eq!(recovered.len(), 1);
        let (a, b) = recovered[0];
        assert!(
            (a - alpha).norm() / alpha.norm() < 1e-10,
            "α: got {a:?}, expected {alpha:?}"
        );
        assert!(
            (b - beta).norm() / beta.norm() < 1e-10,
            "β: got {b:?}, expected {beta:?}"
        );
    }

    /// Sample-count guard: GPOF must reject `M < 2 N` cleanly.
    #[test]
    fn gpof_rejects_under_sampled_input() {
        let samples = vec![Complex64::new(1.0, 0.0); 3];
        let err = gpof(&samples, 1.0, 2).unwrap_err();
        matches!(err, GpofError::TooFewSamples { .. });
    }

    /// Zero-pole guard: `n_poles = 0` is degenerate and rejected.
    #[test]
    fn gpof_rejects_zero_poles() {
        let samples = vec![Complex64::new(1.0, 0.0); 4];
        let err = gpof(&samples, 1.0, 0).unwrap_err();
        matches!(err, GpofError::ZeroPoles);
    }
}
