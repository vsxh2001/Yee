//! Iterative Krylov-subspace solvers for the MoM impedance matrix.
//!
//! Phase 1.6 walking-skeleton scope:
//! - Restarted GMRES(m) with default m = 30, max_restarts = 50.
//! - Diagonal Jacobi preconditioner (cheapest non-trivial preconditioner).
//!   Block-diagonal coming in Phase 1.6.1 when basis-block structure is
//!   exposed by `RwgBasis`.
//! - Convergence target: residual / b_norm <= 1e-10.
//!
//! Intended use case: n ≥ 50k, where the dense partial-pivot LU path used
//! since Phase 1.0 overflows GPU memory or wallclock budgets. For smaller
//! n the direct path remains preferable — GMRES has a higher constant
//! factor and only wins asymptotically.

use faer::Mat;
use num_complex::Complex64;

/// Tuning knobs for [`gmres_jacobi`].
#[derive(Debug, Clone, Copy)]
pub struct GmresParams {
    /// Restart parameter `m`: dimension of the Krylov subspace built between
    /// restarts. Larger `m` converges in fewer outer restarts but costs
    /// `O(m·n)` extra memory and `O(m²·n)` extra work per cycle.
    pub restart: usize,
    /// Maximum number of outer restart cycles before giving up.
    pub max_restarts: usize,
    /// Convergence threshold on the relative preconditioned residual
    /// `||M^{-1}(b - A x)||₂ / ||M^{-1} b||₂`.
    pub tolerance: f64,
}

impl Default for GmresParams {
    fn default() -> Self {
        Self {
            restart: 30,
            max_restarts: 50,
            tolerance: 1.0e-10,
        }
    }
}

/// Outcome of a [`gmres_jacobi`] call.
#[derive(Debug, Clone)]
pub struct GmresResult {
    /// Approximate solution `x ≈ A⁻¹ b`.
    pub x: Mat<Complex64>,
    /// Total number of Arnoldi iterations (summed across restart cycles).
    pub iterations: usize,
    /// Final relative preconditioned residual norm at termination.
    pub final_residual: f64,
    /// `true` iff `final_residual <= params.tolerance`.
    pub converged: bool,
}

/// Solve `A x = b` with restarted GMRES(m).
///
/// This commit lands the unpreconditioned restarted GMRES(m) core
/// (Modified Gram–Schmidt Arnoldi + accumulated Givens rotations). The
/// Jacobi preconditioner is folded in by the next commit; for now `M = I`,
/// so the function name is forward-compatible but the math is plain
/// GMRES(m).
///
/// Algorithm reference: Saad, *Iterative Methods for Sparse Linear
/// Systems* (2nd ed., 2003), §6.5 Algorithm 6.9 / 6.11.
pub fn gmres_jacobi(
    a: faer::MatRef<Complex64>,
    b: faer::MatRef<Complex64>,
    params: GmresParams,
) -> GmresResult {
    let n = a.nrows();

    // M = I (no preconditioning in this commit; Jacobi diagonal lands
    // in the follow-up commit).
    let m_inv: Vec<Complex64> = vec![Complex64::new(1.0, 0.0); n];

    // Initial guess x = 0. Then r = M⁻¹ (b - A·0) = b.
    let mut x = Mat::<Complex64>::zeros(n, 1);

    let mut mb_norm_sq = 0.0_f64;
    for i in 0..n {
        let v = m_inv[i] * b[(i, 0)];
        mb_norm_sq += v.norm_sqr();
    }
    let mb_norm = mb_norm_sq.sqrt().max(f64::MIN_POSITIVE);

    let restart = params.restart.max(1);

    let mut total_iters = 0_usize;
    let mut final_res = f64::INFINITY;
    let mut converged = false;

    // Working scratch for Arnoldi and Givens. Reused across restarts.
    // V is stored as (restart+1) column vectors of length n.
    let mut v: Vec<Vec<Complex64>> = vec![vec![Complex64::new(0.0, 0.0); n]; restart + 1];
    // Upper-Hessenberg H stored column-major as a flat Vec; H[i, j] lives
    // at h[j * (restart + 1) + i].
    let mut h: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); (restart + 1) * restart];
    // Givens rotation cosine (real, ≥0) and sine (complex) for each column.
    let mut cs: Vec<f64> = vec![0.0; restart];
    let mut sn: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); restart];
    // RHS of the projected least-squares problem (β e_1 rotated).
    let mut g: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); restart + 1];

    'outer: for _restart_idx in 0..params.max_restarts {
        // r = M⁻¹ (b - A x)
        let r = preconditioned_residual(a, b, x.as_ref(), &m_inv);
        let beta = vec_norm_slice(&r);
        final_res = beta / mb_norm;
        if final_res <= params.tolerance {
            converged = true;
            break;
        }

        // v_0 = r / beta
        let inv_beta = 1.0 / beta;
        let scale_beta = Complex64::new(inv_beta, 0.0);
        for (dst, &src) in v[0].iter_mut().zip(r.iter()) {
            *dst = src * scale_beta;
        }
        for slot in g.iter_mut() {
            *slot = Complex64::new(0.0, 0.0);
        }
        g[0] = Complex64::new(beta, 0.0);

        // Inner Arnoldi loop.
        let mut k_final = 0_usize;
        let mut inner_converged = false;
        for k in 0..restart {
            total_iters += 1;
            k_final = k;

            // w = M⁻¹ A v_k (left preconditioning).
            let mut w = matvec(a, &v[k]);
            for (w_i, &m_i) in w.iter_mut().zip(m_inv.iter()) {
                *w_i *= m_i;
            }

            // Modified Gram–Schmidt: orthogonalise w against v_0..v_k.
            for j in 0..=k {
                let hjk = cdot(&v[j], &w); // <v_j, w>
                h[k * (restart + 1) + j] = hjk;
                let neg = -hjk;
                for (w_i, &v_ji) in w.iter_mut().zip(v[j].iter()) {
                    *w_i += neg * v_ji;
                }
            }
            let h_kp1_k = vec_norm_slice(&w);
            h[k * (restart + 1) + (k + 1)] = Complex64::new(h_kp1_k, 0.0);

            // Normalise; if breakdown (h_{k+1,k} ≈ 0) we still apply the
            // Givens rotations to the existing Hessenberg column and stop
            // the inner loop — this is the "happy breakdown" case where
            // the Krylov subspace becomes A-invariant and the projected
            // problem is solvable exactly.
            if h_kp1_k > 0.0 {
                let inv = Complex64::new(1.0 / h_kp1_k, 0.0);
                for (dst, &src) in v[k + 1].iter_mut().zip(w.iter()) {
                    *dst = src * inv;
                }
            } else {
                for slot in v[k + 1].iter_mut() {
                    *slot = Complex64::new(0.0, 0.0);
                }
            }

            // Apply previous Givens rotations to column k of H.
            for j in 0..k {
                let a_top = h[k * (restart + 1) + j];
                let a_bot = h[k * (restart + 1) + (j + 1)];
                let c = Complex64::new(cs[j], 0.0);
                let s = sn[j];
                // Standard complex Givens: see Saad §6.5.3.
                //   top' =  c·top + s·bot
                //   bot' = -conj(s)·top + c·bot
                h[k * (restart + 1) + j] = c * a_top + s * a_bot;
                h[k * (restart + 1) + (j + 1)] = -s.conj() * a_top + c * a_bot;
            }

            // Form new Givens rotation that zeros h[k+1, k].
            let a_top = h[k * (restart + 1) + k];
            let a_bot = h[k * (restart + 1) + (k + 1)];
            let (c, s) = givens(a_top, a_bot);
            cs[k] = c;
            sn[k] = s;

            // Apply it to the kth column of H and to g.
            h[k * (restart + 1) + k] = Complex64::new(c, 0.0) * a_top + s * a_bot;
            h[k * (restart + 1) + (k + 1)] = Complex64::new(0.0, 0.0);

            let g_top = g[k];
            let g_bot = g[k + 1];
            g[k] = Complex64::new(c, 0.0) * g_top + s * g_bot;
            g[k + 1] = -s.conj() * g_top + Complex64::new(c, 0.0) * g_bot;

            // The current residual norm is |g[k+1]|.
            let cur_res = g[k + 1].norm() / mb_norm;
            final_res = cur_res;
            if cur_res <= params.tolerance {
                inner_converged = true;
                break;
            }
            if h_kp1_k == 0.0 {
                // Lucky breakdown — Krylov subspace is invariant. The
                // projected problem still has a solution; let the back-
                // substitution step recover it and exit.
                inner_converged = true;
                break;
            }
        }

        // Back-substitution: solve the upper-triangular m×m system
        // R y = g[0..m] where m = k_final + 1.
        let m = k_final + 1;
        let mut y = vec![Complex64::new(0.0, 0.0); m];
        for i in (0..m).rev() {
            let mut sum = g[i];
            for j in (i + 1)..m {
                sum -= h[j * (restart + 1) + i] * y[j];
            }
            let diag = h[i * (restart + 1) + i];
            if diag.norm() == 0.0 {
                // Should not happen after Givens unless m == 0; bail out.
                y[i] = Complex64::new(0.0, 0.0);
            } else {
                y[i] = sum / diag;
            }
        }

        // x ← x + V_m y
        for i in 0..n {
            let mut acc = Complex64::new(0.0, 0.0);
            for j in 0..m {
                acc += v[j][i] * y[j];
            }
            x[(i, 0)] += acc;
        }

        if inner_converged {
            // Recompute the true preconditioned residual to confirm.
            let r = preconditioned_residual(a, b, x.as_ref(), &m_inv);
            final_res = vec_norm_slice(&r) / mb_norm;
            if final_res <= params.tolerance {
                converged = true;
                break 'outer;
            }
        }
    }

    GmresResult {
        x,
        iterations: total_iters,
        final_residual: final_res,
        converged,
    }
}

/// `r = M⁻¹ (b - A x)` as a plain `Vec` — avoids allocating a `Mat`
/// every restart cycle.
fn preconditioned_residual(
    a: faer::MatRef<Complex64>,
    b: faer::MatRef<Complex64>,
    x: faer::MatRef<Complex64>,
    m_inv: &[Complex64],
) -> Vec<Complex64> {
    let n = a.nrows();
    let mut r = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n {
        let mut ax_i = Complex64::new(0.0, 0.0);
        for j in 0..n {
            ax_i += a[(i, j)] * x[(j, 0)];
        }
        r[i] = m_inv[i] * (b[(i, 0)] - ax_i);
    }
    r
}

/// `y = A x` for a dense `MatRef<Complex64>` and a `&[Complex64]` x.
fn matvec(a: faer::MatRef<Complex64>, x: &[Complex64]) -> Vec<Complex64> {
    let n = a.nrows();
    let mut y = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n {
        let mut acc = Complex64::new(0.0, 0.0);
        for j in 0..n {
            acc += a[(i, j)] * x[j];
        }
        y[i] = acc;
    }
    y
}

/// Complex Euclidean inner product `⟨u, v⟩ = Σ ū_i v_i`. Note the
/// conjugation on the left argument — this is the convention that makes
/// Modified Gram–Schmidt produce a true orthonormal basis under the
/// complex inner product.
fn cdot(u: &[Complex64], v: &[Complex64]) -> Complex64 {
    let mut s = Complex64::new(0.0, 0.0);
    for i in 0..u.len() {
        s += u[i].conj() * v[i];
    }
    s
}

fn vec_norm_slice(v: &[Complex64]) -> f64 {
    let mut s = 0.0_f64;
    for c in v {
        s += c.norm_sqr();
    }
    s.sqrt()
}

/// Complex Givens rotation that annihilates `b` against `a`:
///   [  c    s ] [a]   [r]
///   [-s̄    c ] [b] = [0]
/// with real, nonnegative `c` and complex `s` satisfying `c² + |s|² = 1`.
/// Reference: Saad §6.5.3 Eq. (6.80)–(6.81).
fn givens(a: Complex64, b: Complex64) -> (f64, Complex64) {
    let abs_a = a.norm();
    let abs_b = b.norm();
    if abs_b == 0.0 {
        return (1.0, Complex64::new(0.0, 0.0));
    }
    if abs_a == 0.0 {
        return (0.0, Complex64::new(1.0, 0.0));
    }
    let r = (abs_a * abs_a + abs_b * abs_b).sqrt();
    let c = abs_a / r;
    // s = (a / |a|) * conj(b) / r  — chosen so that c·a + s·b = (a/|a|)·r
    // is real-positive and the standard "c real, c² + |s|² = 1" identity
    // holds.
    let s = (a / Complex64::new(abs_a, 0.0)) * b.conj() / Complex64::new(r, 0.0);
    (c, s)
}
