//! Sparse generalised eigensolver — `K e = k² M e` via shift-invert
//! deflated inverse-power iteration on `faer` sparse LU.
//!
//! ## Public surface
//!
//! * [`SparseEigen`] — the load-bearing trait abstracting the
//!   `K e = k² M e` solve. Per Phase 4 spec §8 the trait is the
//!   load-bearing decision; the concrete implementation is the swap
//!   point.
//! * [`InverseIterEigen`] — the concrete implementation shipped in
//!   Phase 4 T5. **Pre-flight finding:** the `lobpcg` crate referenced
//!   in the T5 brief is not published on crates.io at base SHA
//!   `817955a` (2026-05-18). Per the documented escape hatch
//!   (`docs/superpowers/plans/2026-05-18-phase-4-fem-eigenmode.md`
//!   step T5) we ship a hand-rolled deflated inverse-power iteration
//!   on `faer` sparse LU one mode at a time. The downstream gate
//!   (`fem-eig-001`, T7) is unaffected — it consumes the trait.
//! * [`LobpcgEigen`] — the Phase 1.3.1.1 step 4 **block** LOBPCG
//!   (Knyazev 2001) implementation of the same [`SparseEigen`] trait.
//!   It computes the `num_eigs` smallest `k²` eigenpairs
//!   *simultaneously*, carrying an `n × b` block (`b = num_eigs +
//!   guard`) through a single Rayleigh-Ritz step per outer iteration,
//!   and reuses the very same shift-invert operator `(K − σM)⁻¹M`
//!   (factored once via the shared `build_shifted` + faer sparse LU)
//!   as its preconditioner. This resolves the clustered / degenerate
//!   spectra where the sequential `InverseIterEigen` deflation is
//!   weak (TE/TM degeneracies, `TE_{mn}`/`TE_{nm}` pairs). It adds
//!   **no** new dependency (ADR-0050): the small dense `3b × 3b`
//!   Rayleigh-Ritz subproblem reduces via a Cholesky of `SᵀMS` to a
//!   standard symmetric eigenproblem solved by `nalgebra`, which is
//!   already a workspace dep. `InverseIterEigen` remains the default
//!   for existing consumers; `LobpcgEigen` is selected by the caller.
//!   The complex arm (`ComplexLobpcgEigen`) is a step-4.1 follow-on:
//!   lossy dispersive cavities (`fem-eig-002`) keep
//!   [`ComplexInverseIterEigen`]. An optional `arpack` feature behind
//!   the same trait remains available if a >10⁵-DoF cross-section ever
//!   demands Krylov–Schur (ADR-0050 §rationale (4)).
//! * [`EigenpairList`] — the result type: eigenvalues `k²` (sorted
//!   ascending) plus the column-stacked eigenvectors on the
//!   interior-DoF basis the caller supplied.
//!
//! ## Algorithm
//!
//! Solving `K e = k² M e` for the smallest physical eigenvalues is
//! equivalent, after a *shift-invert* transformation at shift `σ`, to
//! finding the *largest* eigenvalues `θ` of the operator
//!
//! ```text
//!     T = (K − σM)^{-1} M.
//! ```
//!
//! The relationship `K e = k² M e` ⇔ `T e = θ e` with
//! `θ = 1 / (k² − σ)` gives back physical eigenvalues
//! `k² = σ + 1 / θ`.
//!
//! `T` is dense in principle but we never form it: the matvec
//! `y ← T x` is one sparse mass-matvec `r ← M x` followed by a
//! sparse triangular back-substitution `y ← (K − σM)^{-1} r` against
//! the LU factorisation produced once up-front.
//!
//! For each requested mode `n ∈ 0..num_eigs`:
//!
//! 1. Seed `x` random, normalize.
//! 2. Loop: `x ← T x`, deflate against converged eigenvectors with
//!    explicit `M`-orthogonalisation (Gram-Schmidt: subtract
//!    `(e_j^T M x) e_j` for each prior `e_j`), normalize in the
//!    `M`-inner product so `x^T M x = 1`. The Rayleigh quotient
//!    `θ = x^T T x / (x^T x)` is updated each iteration; convergence
//!    is `|θ - θ_prev| / |θ| < tol`.
//! 3. Convert `θ → k² = σ + 1 / θ`.
//!
//! Modes are extracted **one at a time** in decreasing-θ order, which
//! is the natural order for the inverse-power iteration. Because the
//! shift-invert spectrum maps small `k²` to large `θ`, this also is
//! the natural increasing-k² order: the most strongly-amplified
//! eigenvector of `T` is the one nearest `σ` in `k²`. Each newly-
//! converged eigenvector is added to the deflation set used for the
//! next mode, preserving `M`-orthogonality of the returned basis to
//! `≤ tol` per the spec §6 invariant.
//!
//! ## Why deflated inverse iteration rather than LOBPCG
//!
//! Inverse-power iteration is the textbook fallback (Golub & Van Loan
//! §8.2.2; Saad *Numerical Methods for Large Eigenvalue Problems*
//! §4.2). It converges linearly with rate proportional to
//! `|θ_{n+1} / θ_n|`; for shift-invert spectra of well-separated
//! interior modes this is rapid in practice. Sparsity is preserved
//! end-to-end: the LU factorisation of `(K − σM)` is the only sparse
//! linear-algebra dependency, and `faer` ships a battle-tested
//! sparse LU. The escape-hatch trade-off vs LOBPCG is one extra
//! Gram-Schmidt sweep per mode (cost `O(n · num_eigs)`), which is
//! cheap at the fem-eig-001 scale (~2000 interior DoFs × 10 modes).
//!
//! ## Convergence diagnostics
//!
//! `InverseIterEigen` returns [`yee_core::Error::Numerical`] when a
//! mode fails to converge within `max_iter` iterations. The error
//! message contains the mode index, the final residual, and the
//! requested tolerance, so callers can either raise `max_iter`,
//! relax `tol`, or move the shift `σ` to skip a known cluster.

use faer::linalg::solvers::SolveCore;
use faer::sparse::{SparseColMat, Triplet, linalg::solvers::Lu};
use nalgebra::DMatrix;
use nalgebra_sparse::csr::CsrMatrix;
use num_complex::Complex64;

/// Solved eigenpairs returned by [`SparseEigen::solve`].
///
/// Eigenvalues are stored as `k²` (the physical eigenvalues of the
/// generalised problem `K e = k² M e`, real-positive for the lossless
/// case), sorted **ascending**. Eigenvectors are column-stacked on the
/// interior-DoF basis the caller supplied; callers needing the full-
/// edge representation lift via the assembly's `interior_edges` map.
///
/// ## Invariants
///
/// * `k.len() == e.ncols()`.
/// * `e.nrows()` equals the dimension of the interior-DoF basis.
/// * Eigenvectors are `M`-orthonormalized to within the solver's
///   working tolerance: `e[:, i]^T M e[:, j] ≈ δ_{ij}` (≤ `tol`).
/// * `k.iter()` is monotonically non-decreasing.
#[derive(Debug, Clone)]
pub struct EigenpairList {
    /// Eigenvalues `k²`, sorted ascending. v0 is real (lossless); the
    /// complex extension lands with `fem-eig-002` per spec §13.
    pub k: Vec<f64>,
    /// Mode-coefficient vectors stacked column-wise on the
    /// interior-DoF basis: `e[:, n]` is the eigenvector for `k[n]`.
    pub e: DMatrix<f64>,
}

/// Trait abstracting the sparse generalized eigensolve
/// `K e = k² M e`.
///
/// Per Phase 4 spec §8 the trait is the load-bearing decision; the
/// concrete implementation is the swap point. v0 ships
/// [`InverseIterEigen`]; a future ARPACK / LOBPCG / SLEPc binding
/// would implement the same trait without touching downstream
/// consumers (`yee-fem`'s assembly, `yee-validation`'s `fem-eig-001`
/// gate, the optional `yee-py` binding).
pub trait SparseEigen {
    /// Solve `K e = k² M e` for the `num_eigs` smallest physical
    /// eigenvalues *near* the shift `σ` (`sigma`). Shift-invert
    /// converts the problem to `(K − σM)^{-1} M e = θ e` and
    /// recovers `k² = σ + 1 / θ`.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] if `k.nrows() != m.nrows()`
    /// or `num_eigs == 0` or `num_eigs > k.nrows()`. Returns
    /// [`yee_core::Error::Numerical`] if the inner sparse LU of
    /// `(K − σM)` fails or any mode fails to converge within the
    /// implementation's configured iteration budget.
    fn solve(
        &self,
        k: &CsrMatrix<f64>,
        m: &CsrMatrix<f64>,
        num_eigs: usize,
        sigma: f64,
    ) -> Result<EigenpairList, yee_core::Error>;
}

/// Deflated shift-invert **inverse-power** iteration on a `faer`
/// sparse LU factorisation of `(K − σM)`.
///
/// **Phase 4 T5 escape-hatch implementation.** The published `lobpcg`
/// crate is not on crates.io at base SHA `817955a` (2026-05-18); per
/// the plan's documented escape hatch we ship a hand-rolled deflated
/// inverse-power iteration instead. The [`SparseEigen`] trait keeps
/// the solver behind an abstraction so the eventual LOBPCG / ARPACK
/// swap is one PR.
///
/// ## Tuning
///
/// * `max_iter` — per-mode iteration budget. The fem-eig-001 gate
///   (Phase 4 T7) hits convergence in `< 50` iterations on the
///   `(8, 6, 10)` cavity mesh; `max_iter = 1000` is safe headroom.
/// * `tol` — relative Rayleigh-quotient convergence target. `1e-8`
///   reaches the fem-eig-001 ±0.3% bound with margin.
///
/// ## Example
///
/// ```ignore
/// use yee_fem::solve::{InverseIterEigen, SparseEigen};
/// let solver = InverseIterEigen::new(1000, 1e-8);
/// let pairs = solver.solve(&k, &m, 10, sigma_k2)?;
/// // pairs.k is sorted ascending; pairs.k[0] is the dominant mode.
/// ```
#[derive(Debug, Clone, Copy)]
pub struct InverseIterEigen {
    /// Per-mode iteration cap. Eigenpairs failing to converge in
    /// `max_iter` iterations cause [`SparseEigen::solve`] to return
    /// [`yee_core::Error::Numerical`] with the mode index and the
    /// last-seen residual.
    pub max_iter: usize,
    /// Relative Rayleigh-quotient convergence tolerance. Iteration
    /// stops when `|θ − θ_prev| / |θ| < tol`.
    pub tol: f64,
}

impl InverseIterEigen {
    /// Construct a configured solver. See type docs for tuning notes.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }
}

impl Default for InverseIterEigen {
    /// `max_iter = 1000`, `tol = 1e-8` — the defaults used by
    /// `fem-eig-001`. See type docs.
    fn default() -> Self {
        Self::new(1000, 1e-8)
    }
}

impl SparseEigen for InverseIterEigen {
    fn solve(
        &self,
        k: &CsrMatrix<f64>,
        m: &CsrMatrix<f64>,
        num_eigs: usize,
        sigma: f64,
    ) -> Result<EigenpairList, yee_core::Error> {
        // ---- Validate shapes --------------------------------------
        if k.nrows() != k.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "InverseIterEigen: K must be square, got {}×{}",
                k.nrows(),
                k.ncols()
            )));
        }
        if m.nrows() != m.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "InverseIterEigen: M must be square, got {}×{}",
                m.nrows(),
                m.ncols()
            )));
        }
        if k.nrows() != m.nrows() {
            return Err(yee_core::Error::Invalid(format!(
                "InverseIterEigen: K and M must have matching dimensions, got K = {}×{} M = {}×{}",
                k.nrows(),
                k.ncols(),
                m.nrows(),
                m.ncols()
            )));
        }
        let n = k.nrows();
        if num_eigs == 0 {
            return Err(yee_core::Error::Invalid(
                "InverseIterEigen: num_eigs must be >= 1".to_string(),
            ));
        }
        if num_eigs > n {
            return Err(yee_core::Error::Invalid(format!(
                "InverseIterEigen: num_eigs = {num_eigs} exceeds dimension {n}"
            )));
        }

        // ---- Build (K − σM) as a faer SparseColMat ----------------
        let shifted = build_shifted(k, m, sigma)?;

        // ---- Factor once via faer sparse LU -----------------------
        let lu: Lu<usize, f64> = shifted.sp_lu().map_err(|e| {
            yee_core::Error::Numerical(format!(
                "InverseIterEigen: sparse LU of (K − σM) failed: {e:?}"
            ))
        })?;

        // ---- Inverse-power iteration with deflation ---------------
        let mut eig_vals: Vec<f64> = Vec::with_capacity(num_eigs);
        let mut eig_vecs: Vec<Vec<f64>> = Vec::with_capacity(num_eigs);

        for mode_idx in 0..num_eigs {
            let mut x = seed_vector(n, mode_idx);
            // M-orthogonalize the seed against converged eigenvectors
            // so we land in the deflated subspace from the start.
            m_orthogonalize(&mut x, &eig_vecs, m);
            m_normalize(&mut x, m)?;

            let mut theta_prev = f64::NAN;
            let mut converged = false;
            let mut last_residual = f64::INFINITY;

            for _iter in 0..self.max_iter {
                // y = M x
                let mx = csr_matvec(m, &x);
                // x_new = (K − σM)^{-1} (M x)
                let mut x_new = lu_solve(&lu, &mx);
                // Deflate against converged eigenvectors:
                // x_new ← x_new − Σ_j (e_j^T M x_new) e_j.
                m_orthogonalize(&mut x_new, &eig_vecs, m);
                // M-normalize: x_new ← x_new / sqrt(x_new^T M x_new).
                m_normalize(&mut x_new, m)?;

                // Rayleigh quotient θ for the operator T = (K−σM)^{-1} M
                // on the M-normalized vector x_new:
                //   T x_new = x_new'  (one step of the same iteration)
                // For the *just-computed* x_new we have
                //   T (previous_x) = x_new (before normalization),
                // so θ ≈ x_new^T M (previous T x_new) / norm.
                // Cleanest path: θ = x_new^T M x_new_T_step.
                // Reconstruct one extra T step to get a current θ
                // estimate that doesn't drift with normalization.
                let mx_new = csr_matvec(m, &x_new);
                let t_x_new = lu_solve(&lu, &mx_new);
                let mut t_x_def = t_x_new.clone();
                m_orthogonalize(&mut t_x_def, &eig_vecs, m);
                // θ = x_new^T M (T x_new)
                let theta = dot(&x_new, &csr_matvec(m, &t_x_def));

                last_residual = if theta_prev.is_finite() && theta.abs() > 0.0 {
                    (theta - theta_prev).abs() / theta.abs()
                } else {
                    f64::INFINITY
                };

                x = x_new;
                if last_residual < self.tol && theta.is_finite() {
                    converged = true;
                    theta_prev = theta;
                    break;
                }
                theta_prev = theta;
            }

            if !converged {
                return Err(yee_core::Error::Numerical(format!(
                    "InverseIterEigen: mode {mode_idx} failed to converge in {} iterations \
                     (last relative residual = {last_residual:e}, tol = {:e})",
                    self.max_iter, self.tol
                )));
            }

            // theta_prev is the final converged θ. Convert to k² and
            // store. Guard against θ → 0 (would correspond to k² = ∞,
            // not physical for finite-energy resonances).
            if theta_prev.abs() < f64::EPSILON {
                return Err(yee_core::Error::Numerical(format!(
                    "InverseIterEigen: mode {mode_idx} converged to θ ≈ 0 \
                     (k² = ∞ — shift σ = {sigma} is far from the spectrum)"
                )));
            }
            let k_sq = sigma + 1.0 / theta_prev;
            eig_vals.push(k_sq);
            eig_vecs.push(x);
        }

        // ---- Sort ascending by k² --------------------------------
        let mut order: Vec<usize> = (0..num_eigs).collect();
        order.sort_by(|&a, &b| eig_vals[a].total_cmp(&eig_vals[b]));

        let sorted_k: Vec<f64> = order.iter().map(|&i| eig_vals[i]).collect();
        let mut e = DMatrix::<f64>::zeros(n, num_eigs);
        for (col, &i) in order.iter().enumerate() {
            for row in 0..n {
                e[(row, col)] = eig_vecs[i][row];
            }
        }

        Ok(EigenpairList { k: sorted_k, e })
    }
}

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

/// Build `(K − σM)` as a `faer::SparseColMat<usize, f64>` ready for
/// sparse LU factorisation.
///
/// The pencil is the standard *shifted operator* of inverse iteration:
/// at non-eigenvalue shifts `σ` it is non-singular and `faer`'s sparse
/// LU factors it once-per-call. Triplets from both `K` and `M` are
/// concatenated; `try_new_from_triplets` accumulates duplicates,
/// matching the entry-summing semantics of the assembly's COO build.
fn build_shifted(
    k: &CsrMatrix<f64>,
    m: &CsrMatrix<f64>,
    sigma: f64,
) -> Result<SparseColMat<usize, f64>, yee_core::Error> {
    let n = k.nrows();
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::with_capacity(k.nnz() + m.nnz());
    // K contributes +k_ij
    for (row, col, &val) in k.triplet_iter() {
        triplets.push(Triplet::new(row, col, val));
    }
    // -σM contributes -σ m_ij
    for (row, col, &val) in m.triplet_iter() {
        triplets.push(Triplet::new(row, col, -sigma * val));
    }
    SparseColMat::try_new_from_triplets(n, n, &triplets).map_err(|e| {
        yee_core::Error::Numerical(format!(
            "InverseIterEigen: failed to build sparse (K − σM): {e:?}"
        ))
    })
}

/// Sparse `y = A x` for a CSR matrix `A` and dense vector `x`. We do
/// this by hand rather than reaching into `nalgebra-sparse`'s
/// `*`-operator surface because we want a `Vec<f64>` output without
/// pulling in `nalgebra::DVector` for trivial workhorse calls.
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

/// Solve `(K − σM) y = b` in place via the pre-computed sparse LU.
///
/// `faer`'s `Lu::solve_in_place_with_conj` operates on a `MatMut<f64>`;
/// we wrap the rhs as a one-column dense matrix, dispatch, and read the
/// column out. Allocation cost is `O(n)` per call which is fine: the
/// inner loop is `max_iter ≤ 1000` invocations per mode, each
/// dominated by the actual triangular back-substitutions inside
/// `faer`.
fn lu_solve(lu: &Lu<usize, f64>, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut rhs = faer::Mat::<f64>::zeros(n, 1);
    for (i, &bi) in b.iter().enumerate() {
        rhs[(i, 0)] = bi;
    }
    lu.solve_in_place_with_conj(faer::Conj::No, rhs.as_mut());
    let mut out = vec![0.0f64; n];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = rhs[(i, 0)];
    }
    out
}

/// `x · y` (Euclidean dot product on dense vectors).
fn dot(x: &[f64], y: &[f64]) -> f64 {
    debug_assert_eq!(x.len(), y.len());
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum()
}

/// In-place `M`-orthogonalize `x` against the converged eigenvectors:
/// `x ← x − Σ_j (e_j^T M x) e_j`. Each `e_j` is already `M`-normalized,
/// so the projection coefficient is the bilinear form `e_j^T M x`.
fn m_orthogonalize(x: &mut [f64], eig_vecs: &[Vec<f64>], m: &CsrMatrix<f64>) {
    for ej in eig_vecs {
        let mx = csr_matvec(m, x);
        let coeff = dot(ej, &mx);
        for (xi, eji) in x.iter_mut().zip(ej.iter()) {
            *xi -= coeff * eji;
        }
    }
}

/// In-place `M`-normalize `x`: `x ← x / sqrt(x^T M x)`. Returns an
/// error if the `M`-norm of `x` underflows to zero (would indicate the
/// deflation killed the seed entirely; bump the seed or relax `tol`).
fn m_normalize(x: &mut [f64], m: &CsrMatrix<f64>) -> Result<(), yee_core::Error> {
    let mx = csr_matvec(m, x);
    let norm_sq = dot(x, &mx);
    if norm_sq <= 0.0 || !norm_sq.is_finite() {
        return Err(yee_core::Error::Numerical(format!(
            "InverseIterEigen: M-norm collapsed to {norm_sq} during deflation"
        )));
    }
    let inv_norm = 1.0 / norm_sq.sqrt();
    for xi in x.iter_mut() {
        *xi *= inv_norm;
    }
    Ok(())
}

/// Seed vector for the `mode_idx`-th inverse-power iteration.
///
/// We use a deterministic seed (a sawtooth shifted by mode index) so
/// the eigensolve is bit-reproducible across runs — critical for the
/// CI gate's pass/fail boundary. Each seed lives outside the span of
/// the previous converged eigenvectors with high probability, which is
/// all the deflation needs.
fn seed_vector(n: usize, mode_idx: usize) -> Vec<f64> {
    let mut x = vec![0.0f64; n];
    let phase = mode_idx as f64;
    for (i, xi) in x.iter_mut().enumerate() {
        // Mix two oscillating components so seeds for different modes
        // are visibly different in every coordinate without being
        // identical-up-to-sign for trivial mode_idx ≠ 0 cases.
        let t = (i as f64 + 1.0) / (n as f64 + 1.0);
        *xi = (1.0 + phase) * t + (1.0 + phase * 0.37).sin() * (t * 7.0).cos();
    }
    x
}

// =====================================================================
// Phase 1.3.1.1 step 4 — block LOBPCG (Knyazev 2001)
// =====================================================================

/// Block **LOBPCG** (Locally Optimal Block Preconditioned Conjugate
/// Gradient, Knyazev 2001) shift-invert eigensolver implementing the
/// [`SparseEigen`] trait.
///
/// Where [`InverseIterEigen`] iterates one mode at a time and deflates
/// sequentially, `LobpcgEigen` carries an `n × b` block
/// (`b = num_eigs + guard`) and resolves it *simultaneously* via a
/// single dense Rayleigh-Ritz step per outer iteration over the
/// search space `S = [X | W | P]` — the current block `X`, the
/// preconditioned residual `W = T·R`, and the previous block `P`. This
/// is exactly the structure that resolves **clustered / degenerate
/// spectra** (TE/TM degeneracies, `TE_{mn}`/`TE_{nm}` pairs on a
/// symmetric cross-section): the block subspace spans the degenerate
/// eigenspace directly, instead of accumulating Gram-Schmidt
/// orthogonality error across a cluster.
///
/// ## Preconditioner
///
/// The preconditioner is the *same* shift-invert operator inverse
/// iteration uses: `T = (K − σM)⁻¹M`, with `(K − σM)` factored exactly
/// once via the shared [`build_shifted`] + faer sparse LU. There is no
/// second factorisation; the only cost delta vs inverse iteration is
/// the small dense `3b × 3b` Rayleigh-Ritz eigensolve per outer
/// iteration, negligible for `b ≤ 20`.
///
/// ## Dense Rayleigh-Ritz (no new dependency — ADR-0050)
///
/// The reduced generalized symmetric problem `(SᵀKS) c = θ (SᵀMS) c`
/// is solved by the **Cholesky-reduction** path (spec §7 risk (b)):
/// `SᵀMS = L Lᵀ`, transform to the standard symmetric problem
/// `(L⁻¹ SᵀKS L⁻ᵀ) y = θ y`, run `nalgebra`'s symmetric eigensolver,
/// then back-transform `c = L⁻ᵀ y`. `nalgebra` is already a workspace
/// dependency, so this adds **zero** new `Cargo.toml` lines, matching
/// the project's pure-Rust LA ethos.
///
/// ## Soft-locking
///
/// Near convergence the previous block `P` collapses into `span(X)`
/// and the Gram matrix `SᵀMS` goes near-singular. Following Knyazev §4
/// ("soft locking") the `[X|W|P]` columns are M-orthonormalised by
/// modified Gram-Schmidt and any column whose post-orthogonalisation
/// M-norm falls below `√ε` is dropped, shrinking the search block that
/// iteration. The first outer iteration carries no `P`; it is
/// introduced from iteration two onward.
///
/// ## Determinism
///
/// The initial block `X₀` is seeded from the same deterministic
/// [`seed_vector`] generator the inverse-power path uses (one column
/// per block index), so the eigensolve is bit-reproducible across runs
/// — critical for the CI gate's pass/fail boundary. No thread RNG.
///
/// ## Tuning
///
/// * `max_iter` — outer-iteration budget for the whole block.
///   `LobpcgEigen` typically converges the leading `num_eigs` columns
///   in far fewer outer iterations than `InverseIterEigen` consumes
///   per mode, because the block step is locally optimal.
/// * `tol` — per-column relative residual target
///   `‖K xᵢ − k²ᵢ M xᵢ‖₂ / (k²ᵢ ‖M xᵢ‖₂) < tol` on the leading
///   `num_eigs` columns.
/// * `guard` — extra columns beyond `num_eigs` (block width
///   `b = num_eigs + guard`, capped at `n`). Guard columns accelerate
///   cluster resolution by giving the block room to separate nearly-
///   degenerate roots; `guard = 2` is the default.
///
/// ## Example
///
/// ```ignore
/// use yee_fem::solve::{LobpcgEigen, SparseEigen};
/// let solver = LobpcgEigen::new(1000, 1e-8, 2);
/// let pairs = solver.solve(&k, &m, 10, sigma_k2)?;
/// // pairs.k is sorted ascending and M-orthonormal — same
/// // postcondition contract as InverseIterEigen.
/// ```
#[derive(Debug, Clone, Copy)]
pub struct LobpcgEigen {
    /// Outer-iteration budget for the block. Failure to converge the
    /// leading `num_eigs` columns within `max_iter` outer iterations
    /// causes [`SparseEigen::solve`] to return
    /// [`yee_core::Error::Numerical`] with the worst-column residual.
    pub max_iter: usize,
    /// Per-column relative residual convergence tolerance on the
    /// leading `num_eigs` columns.
    pub tol: f64,
    /// Guard columns beyond `num_eigs`: block width
    /// `b = (num_eigs + guard).min(n)`. Improves cluster robustness.
    pub guard: usize,
}

impl LobpcgEigen {
    /// Construct a configured solver. See type docs for tuning notes.
    pub fn new(max_iter: usize, tol: f64, guard: usize) -> Self {
        Self {
            max_iter,
            tol,
            guard,
        }
    }
}

impl Default for LobpcgEigen {
    /// `max_iter = 1000`, `tol = 1e-8`, `guard = 2` — the defaults
    /// mirroring [`InverseIterEigen`] plus a two-column cluster guard.
    fn default() -> Self {
        Self::new(1000, 1e-8, 2)
    }
}

impl SparseEigen for LobpcgEigen {
    fn solve(
        &self,
        k: &CsrMatrix<f64>,
        m: &CsrMatrix<f64>,
        num_eigs: usize,
        sigma: f64,
    ) -> Result<EigenpairList, yee_core::Error> {
        // ---- Validate shapes (identical guards to InverseIterEigen) -
        if k.nrows() != k.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "LobpcgEigen: K must be square, got {}×{}",
                k.nrows(),
                k.ncols()
            )));
        }
        if m.nrows() != m.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "LobpcgEigen: M must be square, got {}×{}",
                m.nrows(),
                m.ncols()
            )));
        }
        if k.nrows() != m.nrows() {
            return Err(yee_core::Error::Invalid(format!(
                "LobpcgEigen: K and M must have matching dimensions, got K = {}×{} M = {}×{}",
                k.nrows(),
                k.ncols(),
                m.nrows(),
                m.ncols()
            )));
        }
        let n = k.nrows();
        if num_eigs == 0 {
            return Err(yee_core::Error::Invalid(
                "LobpcgEigen: num_eigs must be >= 1".to_string(),
            ));
        }
        if num_eigs > n {
            return Err(yee_core::Error::Invalid(format!(
                "LobpcgEigen: num_eigs = {num_eigs} exceeds dimension {n}"
            )));
        }

        // ---- Build (K − σM) and factor once via faer sparse LU ------
        let shifted = build_shifted(k, m, sigma)?;
        let lu: Lu<usize, f64> = shifted.sp_lu().map_err(|e| {
            yee_core::Error::Numerical(format!("LobpcgEigen: sparse LU of (K − σM) failed: {e:?}"))
        })?;
        // T x = (K − σM)^{-1} (M x): the shared shift-invert operator.
        let apply_t = |cols: &[Vec<f64>]| -> Vec<Vec<f64>> {
            cols.iter()
                .map(|c| {
                    let mc = csr_matvec(m, c);
                    lu_solve(&lu, &mc)
                })
                .collect()
        };

        // ---- Block width and deterministic, M-orthonormal seed ------
        let b = (num_eigs + self.guard).min(n);
        let mut x: Vec<Vec<f64>> = (0..b).map(|j| seed_vector(n, j)).collect();
        block_m_orthonormalize(&mut x, m)?;
        if x.len() < num_eigs {
            // Seeds were rank-deficient against M; this only happens for
            // pathological pencils. Surface rather than silently return
            // fewer modes than requested.
            return Err(yee_core::Error::Numerical(format!(
                "LobpcgEigen: initial block rank {} < num_eigs {num_eigs} after \
                 M-orthonormalisation (degenerate seed against M)",
                x.len()
            )));
        }

        // Previous block P, M-orthonormalised against X each iteration.
        let mut p: Vec<Vec<f64>> = Vec::new();
        // Current Ritz values (k²) for the block columns.
        let mut theta_k2: Vec<f64> = vec![sigma; x.len()];

        let mut last_max_res = f64::INFINITY;
        let mut converged = false;

        for _outer in 0..self.max_iter {
            // ---- Block residual R = K X − M X Λ  (Λ = current k²) ----
            let mut residual: Vec<Vec<f64>> = Vec::with_capacity(x.len());
            let mut max_res = 0.0f64;
            for (col, xi) in x.iter().enumerate() {
                let kx = csr_matvec(k, xi);
                let mx = csr_matvec(m, xi);
                let lam = theta_k2[col];
                let ri: Vec<f64> = kx
                    .iter()
                    .zip(mx.iter())
                    .map(|(&kxi, &mxi)| kxi - lam * mxi)
                    .collect();
                // Per-column relative residual on the leading num_eigs.
                if col < num_eigs {
                    let rnorm = dot(&ri, &ri).sqrt();
                    let mxnorm = dot(&mx, &mx).sqrt();
                    let denom = lam.abs() * mxnorm;
                    let rel = if denom > 0.0 { rnorm / denom } else { rnorm };
                    if rel > max_res {
                        max_res = rel;
                    }
                }
                residual.push(ri);
            }
            last_max_res = max_res;
            if max_res < self.tol {
                converged = true;
                break;
            }

            // ---- Preconditioned residual W = T R --------------------
            let w = apply_t(&residual);

            // ---- Search space S = [X | W | P], M-orthonormalised ----
            // Soft-locking: M-orthonormalise the whole stack by modified
            // Gram-Schmidt and drop near-null columns (Knyazev §4).
            let nx = x.len();
            let nw = w.len();
            let mut s: Vec<Vec<f64>> = Vec::with_capacity(nx + nw + p.len());
            s.extend(x.iter().cloned());
            s.extend(w.into_iter());
            s.extend(p.iter().cloned());
            block_m_orthonormalize(&mut s, m)?;
            let bb = s.len();
            if bb < num_eigs {
                return Err(yee_core::Error::Numerical(format!(
                    "LobpcgEigen: search space collapsed to rank {bb} < num_eigs \
                     {num_eigs} (basis ill-conditioning — raise guard or move σ)"
                )));
            }

            // ---- Dense Rayleigh-Ritz on S ---------------------------
            // Sᵀ K S and Sᵀ M S (bb × bb, symmetric). S is M-orthonormal
            // so SᵀMS ≈ I, but we form it explicitly for the Cholesky
            // reduction so rounding does not bias the Ritz values.
            let st_k_s = block_gram(&s, k);
            let st_m_s = block_gram(&s, m);
            let (ritz_vals, ritz_vecs) = dense_gen_sym_eigen(&st_k_s, &st_m_s)?;

            // The bb Ritz pairs are sorted ascending; the leading `nx`
            // (= block width) define the next X. Columns are the
            // coefficient vectors in the S basis.
            let take = nx.min(bb);

            // New X = S · C[:, 0..take].
            let new_x = combine(&s, &ritz_vecs, 0, take);
            // New P built from the W,P portion of the Ritz combination
            // (columns nx.. of S) so the conjugate direction is carried
            // forward (Knyazev's "local optimality").
            let new_p = if bb > nx {
                combine_rows(&s, &ritz_vecs, nx, bb, 0, take)
            } else {
                Vec::new()
            };

            theta_k2 = ritz_vals[0..take].to_vec();
            x = new_x;
            p = new_p;
            block_m_orthonormalize(&mut x, m)?;
            if !p.is_empty() {
                block_m_orthonormalize(&mut p, m)?;
            }
        }

        if !converged {
            return Err(yee_core::Error::Numerical(format!(
                "LobpcgEigen: block failed to converge in {} outer iterations \
                 (worst leading-column relative residual = {last_max_res:e}, tol = {:e})",
                self.max_iter, self.tol
            )));
        }

        // ---- Assemble the leading num_eigs pairs, sorted ascending --
        // theta_k2 is already ascending (Ritz order); take the leading
        // num_eigs and the matching block columns. Re-M-orthonormalise
        // the returned set so the EigenpairList postcondition holds to
        // working tolerance.
        let mut take_vecs: Vec<Vec<f64>> = x[0..num_eigs].to_vec();
        block_m_orthonormalize(&mut take_vecs, m)?;
        if take_vecs.len() < num_eigs {
            return Err(yee_core::Error::Numerical(format!(
                "LobpcgEigen: returned block rank {} < num_eigs {num_eigs} after final \
                 M-orthonormalisation",
                take_vecs.len()
            )));
        }

        // Recompute Ritz values on the final orthonormal vectors so the
        // reported k² is the Rayleigh quotient of the returned vector,
        // not a stale block value: k²ᵢ = xᵢᵀ K xᵢ / xᵢᵀ M xᵢ.
        let mut k_vals: Vec<f64> = Vec::with_capacity(num_eigs);
        for xi in &take_vecs {
            let kx = csr_matvec(k, xi);
            let mx = csr_matvec(m, xi);
            let num = dot(xi, &kx);
            let den = dot(xi, &mx);
            k_vals.push(num / den);
        }

        let mut order: Vec<usize> = (0..num_eigs).collect();
        order.sort_by(|&a, &c| k_vals[a].total_cmp(&k_vals[c]));
        let sorted_k: Vec<f64> = order.iter().map(|&i| k_vals[i]).collect();
        let mut e = DMatrix::<f64>::zeros(n, num_eigs);
        for (col, &i) in order.iter().enumerate() {
            for row in 0..n {
                e[(row, col)] = take_vecs[i][row];
            }
        }

        Ok(EigenpairList { k: sorted_k, e })
    }
}

// ---------------------------------------------------------------------
// Block helpers — peers of the single-vector helpers above, reusing
// `csr_matvec` / `dot` / `seed_vector` (no duplication).
// ---------------------------------------------------------------------

/// In-place block M-orthonormalisation by **modified Gram-Schmidt** in
/// the `M`-inner product: each column is M-orthogonalised against the
/// already-accepted columns, then M-normalised. Columns whose
/// post-orthogonalisation M-norm falls below `√ε` are **dropped** (the
/// returned block shrinks) — this is Knyazev's "soft locking" guard
/// against the near-singular `[X|W|P]` Gram matrix at convergence.
///
/// Postcondition: `block[i]ᵀ M block[j] ≈ δ_{ij}` for the surviving
/// columns. Reuses [`csr_matvec`] and [`dot`].
fn block_m_orthonormalize(
    block: &mut Vec<Vec<f64>>,
    m: &CsrMatrix<f64>,
) -> Result<(), yee_core::Error> {
    // Drop-tolerance on the M-norm after orthogonalisation. sqrt(eps)
    // is the standard soft-locking threshold (Knyazev §4).
    let drop_tol = f64::EPSILON.sqrt();
    let mut accepted: Vec<Vec<f64>> = Vec::with_capacity(block.len());
    for col in block.drain(..) {
        let mut v = col;
        // Modified Gram-Schmidt against accepted columns (each already
        // M-normalised, so the projection coefficient is e_jᵀ M v).
        for ej in &accepted {
            let mv = csr_matvec(m, &v);
            let coeff = dot(ej, &mv);
            for (vi, eji) in v.iter_mut().zip(ej.iter()) {
                *vi -= coeff * eji;
            }
        }
        let mv = csr_matvec(m, &v);
        let norm_sq = dot(&v, &mv);
        if !norm_sq.is_finite() {
            return Err(yee_core::Error::Numerical(format!(
                "LobpcgEigen: block M-norm went non-finite ({norm_sq}) during \
                 orthonormalisation"
            )));
        }
        if norm_sq <= drop_tol * drop_tol {
            // Soft-lock: this column collapsed into the accepted span;
            // drop it and continue (search block shrinks this iter).
            continue;
        }
        let inv = 1.0 / norm_sq.sqrt();
        for vi in v.iter_mut() {
            *vi *= inv;
        }
        accepted.push(v);
    }
    *block = accepted;
    Ok(())
}

/// Dense `b' × b'` Gram matrix `Sᵀ A S` for a block `S` (`Vec` of `n`-
/// vectors) and a sparse CSR `A`. Symmetric by construction for
/// symmetric `A`; we symmetrise the result to kill rounding asymmetry
/// before it feeds the symmetric dense eigensolver. Reuses
/// [`csr_matvec`] and [`dot`].
fn block_gram(s: &[Vec<f64>], a: &CsrMatrix<f64>) -> DMatrix<f64> {
    let bb = s.len();
    let a_s: Vec<Vec<f64>> = s.iter().map(|c| csr_matvec(a, c)).collect();
    let mut g = DMatrix::<f64>::zeros(bb, bb);
    for i in 0..bb {
        for j in i..bb {
            let v = dot(&s[i], &a_s[j]);
            g[(i, j)] = v;
            g[(j, i)] = v;
        }
    }
    g
}

/// Solve the small dense **generalized symmetric** eigenproblem
/// `A c = θ B c` (`A = SᵀKS`, `B = SᵀMS`) via the **Cholesky-reduction**
/// path (spec §7 risk (b)): `B = L Lᵀ`, transform to the standard
/// symmetric problem `(L⁻¹ A L⁻ᵀ) y = θ y`, solve with `nalgebra`'s
/// symmetric eigensolver, then back-transform the generalized
/// eigenvectors `c = L⁻ᵀ y`. Returns `(eigenvalues, eigenvectors)`
/// sorted **ascending** by eigenvalue, eigenvectors column-stacked in
/// the `S` basis.
///
/// `nalgebra` is already a workspace dependency, so this Rayleigh-Ritz
/// dense solve adds **no** new `Cargo.toml` line (ADR-0050).
fn dense_gen_sym_eigen(
    a: &DMatrix<f64>,
    bmat: &DMatrix<f64>,
) -> Result<(Vec<f64>, DMatrix<f64>), yee_core::Error> {
    let bb = a.nrows();
    // B = L Lᵀ. B is SᵀMS with M SPD on the surviving block, so the
    // Cholesky should succeed; a failure means the block is rank-
    // deficient despite the soft-lock drop — surface it.
    let chol = bmat.clone().cholesky().ok_or_else(|| {
        yee_core::Error::Numerical(
            "LobpcgEigen: Rayleigh-Ritz Cholesky of SᵀMS failed (near-singular \
             block Gram matrix; raise guard or move σ)"
                .to_string(),
        )
    })?;
    let l = chol.l();
    let l_inv = l.clone().try_inverse().ok_or_else(|| {
        yee_core::Error::Numerical(
            "LobpcgEigen: Rayleigh-Ritz triangular inverse of L failed".to_string(),
        )
    })?;
    // Standard symmetric problem matrix Ã = L⁻¹ A L⁻ᵀ. Symmetrise to
    // remove rounding asymmetry before the symmetric eigensolve.
    let mut a_tilde = &l_inv * a * l_inv.transpose();
    let a_sym = (&a_tilde + a_tilde.transpose()) * 0.5;
    a_tilde = a_sym;
    let eig = a_tilde.symmetric_eigen();
    // nalgebra does not guarantee sorted eigenvalues; sort ascending and
    // permute the eigenvectors to match.
    let mut idx: Vec<usize> = (0..bb).collect();
    idx.sort_by(|&i, &j| eig.eigenvalues[i].total_cmp(&eig.eigenvalues[j]));
    let vals: Vec<f64> = idx.iter().map(|&i| eig.eigenvalues[i]).collect();
    // Back-transform generalized eigenvectors c = L⁻ᵀ y, then permute.
    let l_inv_t = l_inv.transpose();
    let mut c = DMatrix::<f64>::zeros(bb, bb);
    for (new_col, &old_col) in idx.iter().enumerate() {
        let y = eig.eigenvectors.column(old_col);
        let cy = &l_inv_t * y;
        for row in 0..bb {
            c[(row, new_col)] = cy[row];
        }
    }
    Ok((vals, c))
}

/// Linear combination of block columns: `out[j] = Σ_r S[r] · C[r, c0+j]`
/// for `j ∈ 0..(take)`. Forms the `take` new physical-space vectors
/// from the Ritz coefficient matrix `C` (S-basis) columns `c0..c0+take`.
fn combine(s: &[Vec<f64>], c: &DMatrix<f64>, c0: usize, take: usize) -> Vec<Vec<f64>> {
    let n = s[0].len();
    let bb = s.len();
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(take);
    for col in c0..c0 + take {
        let mut v = vec![0.0f64; n];
        for (r, sr) in s.iter().enumerate().take(bb) {
            let coeff = c[(r, col)];
            if coeff == 0.0 {
                continue;
            }
            for (vi, &sri) in v.iter_mut().zip(sr.iter()) {
                *vi += coeff * sri;
            }
        }
        out.push(v);
    }
    out
}

/// Like [`combine`] but uses only the `S` rows `r0..r1` of the Ritz
/// coefficients (the `W`,`P` portion), forming the next conjugate
/// block `P` from the non-`X` part of the Ritz combination — Knyazev's
/// local-optimality direction. Columns `c0..c0+take` of `C`.
fn combine_rows(
    s: &[Vec<f64>],
    c: &DMatrix<f64>,
    r0: usize,
    r1: usize,
    c0: usize,
    take: usize,
) -> Vec<Vec<f64>> {
    let n = s[0].len();
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(take);
    for col in c0..c0 + take {
        let mut v = vec![0.0f64; n];
        for r in r0..r1 {
            let coeff = c[(r, col)];
            if coeff == 0.0 {
                continue;
            }
            for (vi, &sri) in v.iter_mut().zip(s[r].iter()) {
                *vi += coeff * sri;
            }
        }
        out.push(v);
    }
    out
}

// =====================================================================
// Phase 4.fem.eig.1 — complex-coefficient peer
// =====================================================================
//
// Per ADR-0039 / spec §8, Phase 4.fem.eig.1 introduces a complex peer
// `ComplexInverseIterEigen` of the v0 `InverseIterEigen` behind a new
// `SparseEigenComplex` trait sibling to `SparseEigen<f64>`. Two traits
// is cleaner than one parametric trait given the complex symmetric (not
// Hermitian) inner-product conventions for lossy materials and the
// pivoting differences in complex LU.
//
// Algorithm mirrors the v0 path: shift-invert
// `T = (K − σM)^{-1} M`, deflated inverse-power iteration mode-by-mode,
// Gram-Schmidt M-orthogonalisation against converged eigenvectors. The
// substitutions are:
//
// * Scalar field: `f64 → Complex64`.
// * Inner product: **transposed** `e^T M e` (not Hermitian `e^H M e`).
//   The complex symmetric pencil arising from FEM dispersive materials
//   is symmetric under transposition only — the Hellmann–Feynman
//   identity in plan step D5 lands in the natural transposed form
//   exactly because of this convention.
// * Normalisation: `x ← x / sqrt(x^T M x)` where the square root is
//   the principal complex square root.
// * Convergence: `|θ_new − θ_prev| / |θ_new|` with `Complex64::norm`.
// * Sort order: ascending by `Re(k²)` (the natural one-mode-at-a-time
//   ordering from inverse-power iteration on a complex-symmetric
//   pencil with mostly-real eigenvalues; see spec §11 risk register
//   for the mode-crossing caveat).
//
// `faer 0.23` exposes `Lu<usize, T>` for `T: ComplexField` and
// `num_complex::Complex<T: RealField>` is a `ComplexField`, so the
// sparse LU path is the same `try_new_from_triplets → sp_lu →
// solve_in_place_with_conj` chain as the v0 real path. No dense
// fallback is necessary at the base SHA `a08f0db`.

/// Complex-coefficient eigenpairs returned by [`SparseEigenComplex::solve`].
///
/// Eigenvalues are stored as **complex** `k²` (the physical
/// eigenvalues of the generalised problem `K e = k² M e` for the
/// dispersive case; `Im(k²) ≤ 0` for lossy media in the engineering
/// `exp(+jωt)` convention used throughout Yee). Sorted ascending by
/// `Re(k²)`. Eigenvectors are column-stacked on the interior-DoF basis
/// the caller supplied.
///
/// ## Invariants
///
/// * `k.len() == e.ncols()`.
/// * `e.nrows()` equals the dimension of the interior-DoF basis.
/// * Eigenvectors are `M`-orthonormalised in the *transposed* inner
///   product to within the solver's working tolerance:
///   `e[:, i]^T M e[:, j] ≈ δ_{ij}` (≤ `tol`). Note: **transposed**,
///   not Hermitian — see module-level comment for why.
/// * `k.iter().map(|z| z.re)` is monotonically non-decreasing.
#[derive(Debug, Clone)]
pub struct EigenpairListComplex {
    /// Complex eigenvalues `k²`, sorted ascending by `Re(k²)`.
    pub k: Vec<Complex64>,
    /// Mode-coefficient vectors stacked column-wise on the
    /// interior-DoF basis: `e[:, n]` is the eigenvector for `k[n]`.
    pub e: DMatrix<Complex64>,
}

/// Trait abstracting the **complex-coefficient** sparse generalised
/// eigensolve `K e = k² M e` for `K`, `M` ∈ `CsrMatrix<Complex64>`.
///
/// Sibling of [`SparseEigen`] (the real-coefficient trait); per
/// ADR-0039 §8 the two traits live alongside each other rather than
/// being unified parametrically because the complex symmetric (not
/// Hermitian) inner-product conventions and the pivoting differences
/// in complex LU make a single trait clumsy.
pub trait SparseEigenComplex {
    /// Solve `K e = k² M e` for the `num_eigs` eigenvalues *near* the
    /// complex shift `σ`. Shift-invert converts the problem to
    /// `(K − σM)^{-1} M e = θ e` and recovers `k² = σ + 1 / θ`.
    ///
    /// # Errors
    ///
    /// Returns [`yee_core::Error::Invalid`] for shape mismatches or
    /// `num_eigs == 0` / `num_eigs > k.nrows()`. Returns
    /// [`yee_core::Error::Numerical`] when the sparse LU of `(K − σM)`
    /// fails or a mode does not converge within the implementation's
    /// configured iteration budget.
    fn solve(
        &self,
        k: &CsrMatrix<Complex64>,
        m: &CsrMatrix<Complex64>,
        num_eigs: usize,
        sigma: Complex64,
    ) -> Result<EigenpairListComplex, yee_core::Error>;
}

/// Complex peer of [`InverseIterEigen`] — deflated shift-invert
/// inverse-power iteration over `faer::sparse` `Lu<usize, Complex64>`.
///
/// Mirrors the v0 algorithm with `Complex64` arithmetic throughout
/// and a **transposed** (not Hermitian) M-inner product per ADR-0039.
/// For real-valued `K`, `M` and real `σ` the results are bit-for-bit
/// equal to the v0 [`InverseIterEigen`] on the same input — the
/// `complex_lobpcg_smoke` integration test pins this invariant.
///
/// ## Tuning
///
/// * `max_iter` — per-mode iteration budget. Same defaults as the
///   real path (`1000` is safe headroom for fem-eig-002 scale).
/// * `tol` — relative Rayleigh-quotient convergence target in the
///   complex norm. `1e-8` is the v0 default.
#[derive(Debug, Clone, Copy)]
pub struct ComplexInverseIterEigen {
    /// Per-mode iteration cap.
    pub max_iter: usize,
    /// Relative complex Rayleigh-quotient convergence tolerance.
    pub tol: f64,
}

impl ComplexInverseIterEigen {
    /// Construct a configured solver. See type docs for tuning notes.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }
}

impl Default for ComplexInverseIterEigen {
    /// `max_iter = 1000`, `tol = 1e-8` — matches the v0
    /// [`InverseIterEigen`] defaults.
    fn default() -> Self {
        Self::new(1000, 1e-8)
    }
}

impl SparseEigenComplex for ComplexInverseIterEigen {
    fn solve(
        &self,
        k: &CsrMatrix<Complex64>,
        m: &CsrMatrix<Complex64>,
        num_eigs: usize,
        sigma: Complex64,
    ) -> Result<EigenpairListComplex, yee_core::Error> {
        // ---- Validate shapes --------------------------------------
        if k.nrows() != k.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "ComplexInverseIterEigen: K must be square, got {}×{}",
                k.nrows(),
                k.ncols()
            )));
        }
        if m.nrows() != m.ncols() {
            return Err(yee_core::Error::Invalid(format!(
                "ComplexInverseIterEigen: M must be square, got {}×{}",
                m.nrows(),
                m.ncols()
            )));
        }
        if k.nrows() != m.nrows() {
            return Err(yee_core::Error::Invalid(format!(
                "ComplexInverseIterEigen: K and M must have matching dimensions, \
                 got K = {}×{} M = {}×{}",
                k.nrows(),
                k.ncols(),
                m.nrows(),
                m.ncols()
            )));
        }
        let n = k.nrows();
        if num_eigs == 0 {
            return Err(yee_core::Error::Invalid(
                "ComplexInverseIterEigen: num_eigs must be >= 1".to_string(),
            ));
        }
        if num_eigs > n {
            return Err(yee_core::Error::Invalid(format!(
                "ComplexInverseIterEigen: num_eigs = {num_eigs} exceeds dimension {n}"
            )));
        }

        // ---- Build (K − σM) as a faer SparseColMat<usize, Complex64>
        let shifted = build_shifted_complex(k, m, sigma)?;

        // ---- Factor once via faer sparse LU -----------------------
        let lu: Lu<usize, Complex64> = shifted.sp_lu().map_err(|e| {
            yee_core::Error::Numerical(format!(
                "ComplexInverseIterEigen: sparse LU of (K − σM) failed: {e:?}"
            ))
        })?;

        // ---- Inverse-power iteration with deflation ---------------
        let mut eig_vals: Vec<Complex64> = Vec::with_capacity(num_eigs);
        let mut eig_vecs: Vec<Vec<Complex64>> = Vec::with_capacity(num_eigs);

        for mode_idx in 0..num_eigs {
            let mut x = seed_vector_complex(n, mode_idx);
            m_orthogonalize_complex(&mut x, &eig_vecs, m);
            m_normalize_complex(&mut x, m)?;

            let mut theta_prev = Complex64::new(f64::NAN, f64::NAN);
            let mut converged = false;
            let mut last_residual = f64::INFINITY;

            for _iter in 0..self.max_iter {
                let mx = csr_matvec_complex(m, &x);
                let mut x_new = lu_solve_complex(&lu, &mx);
                m_orthogonalize_complex(&mut x_new, &eig_vecs, m);
                m_normalize_complex(&mut x_new, m)?;

                // Reconstruct one extra T step so the Rayleigh quotient
                // is a fresh estimate at the normalised x_new (matches
                // the v0 algorithm bit-for-bit).
                let mx_new = csr_matvec_complex(m, &x_new);
                let t_x_new = lu_solve_complex(&lu, &mx_new);
                let mut t_x_def = t_x_new.clone();
                m_orthogonalize_complex(&mut t_x_def, &eig_vecs, m);
                // θ = x_new^T M (T x_new) — transposed, not Hermitian.
                let theta = dot_complex(&x_new, &csr_matvec_complex(m, &t_x_def));

                last_residual = if theta_prev.is_finite() && theta.norm() > 0.0 {
                    (theta - theta_prev).norm() / theta.norm()
                } else {
                    f64::INFINITY
                };

                x = x_new;
                if last_residual < self.tol && theta.is_finite() {
                    converged = true;
                    theta_prev = theta;
                    break;
                }
                theta_prev = theta;
            }

            if !converged {
                return Err(yee_core::Error::Numerical(format!(
                    "ComplexInverseIterEigen: mode {mode_idx} failed to converge in {} \
                     iterations (last relative residual = {last_residual:e}, tol = {:e})",
                    self.max_iter, self.tol
                )));
            }

            if theta_prev.norm() < f64::EPSILON {
                return Err(yee_core::Error::Numerical(format!(
                    "ComplexInverseIterEigen: mode {mode_idx} converged to θ ≈ 0 \
                     (k² = ∞ — shift σ = {sigma} is far from the spectrum)"
                )));
            }
            let k_sq = sigma + Complex64::new(1.0, 0.0) / theta_prev;
            eig_vals.push(k_sq);
            eig_vecs.push(x);
        }

        // ---- Sort ascending by Re(k²) -----------------------------
        let mut order: Vec<usize> = (0..num_eigs).collect();
        order.sort_by(|&a, &b| eig_vals[a].re.total_cmp(&eig_vals[b].re));

        let sorted_k: Vec<Complex64> = order.iter().map(|&i| eig_vals[i]).collect();
        let mut e = DMatrix::<Complex64>::zeros(n, num_eigs);
        for (col, &i) in order.iter().enumerate() {
            for row in 0..n {
                e[(row, col)] = eig_vecs[i][row];
            }
        }

        Ok(EigenpairListComplex { k: sorted_k, e })
    }
}

// ---------------------------------------------------------------------
// Complex helpers — peers of the real helpers above
// ---------------------------------------------------------------------

/// Build `(K − σM)` as a `faer::SparseColMat<usize, Complex64>` ready
/// for sparse LU factorisation. Peer of [`build_shifted`].
fn build_shifted_complex(
    k: &CsrMatrix<Complex64>,
    m: &CsrMatrix<Complex64>,
    sigma: Complex64,
) -> Result<SparseColMat<usize, Complex64>, yee_core::Error> {
    let n = k.nrows();
    let mut triplets: Vec<Triplet<usize, usize, Complex64>> = Vec::with_capacity(k.nnz() + m.nnz());
    for (row, col, &val) in k.triplet_iter() {
        triplets.push(Triplet::new(row, col, val));
    }
    for (row, col, &val) in m.triplet_iter() {
        triplets.push(Triplet::new(row, col, -sigma * val));
    }
    SparseColMat::try_new_from_triplets(n, n, &triplets).map_err(|e| {
        yee_core::Error::Numerical(format!(
            "ComplexInverseIterEigen: failed to build sparse (K − σM): {e:?}"
        ))
    })
}

/// Sparse `y = A x` for a complex CSR matrix. Peer of [`csr_matvec`].
fn csr_matvec_complex(a: &CsrMatrix<Complex64>, x: &[Complex64]) -> Vec<Complex64> {
    let n = a.nrows();
    let mut y = vec![Complex64::new(0.0, 0.0); n];
    let row_offsets = a.row_offsets();
    let col_indices = a.col_indices();
    let values = a.values();
    for row in 0..n {
        let start = row_offsets[row];
        let end = row_offsets[row + 1];
        let mut sum = Complex64::new(0.0, 0.0);
        for k in start..end {
            sum += values[k] * x[col_indices[k]];
        }
        y[row] = sum;
    }
    y
}

/// Solve `(K − σM) y = b` in place via the pre-computed complex sparse
/// LU. Peer of [`lu_solve`]. `Conj::No` selects the unconjugated solve
/// (matching the natural transposed `e^T M e` convention used by the
/// outer inverse-power iteration).
fn lu_solve_complex(lu: &Lu<usize, Complex64>, b: &[Complex64]) -> Vec<Complex64> {
    let n = b.len();
    let mut rhs = faer::Mat::<Complex64>::zeros(n, 1);
    for (i, &bi) in b.iter().enumerate() {
        rhs[(i, 0)] = bi;
    }
    lu.solve_in_place_with_conj(faer::Conj::No, rhs.as_mut());
    let mut out = vec![Complex64::new(0.0, 0.0); n];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = rhs[(i, 0)];
    }
    out
}

/// Transposed (not Hermitian) bilinear product `x^T y` on complex
/// dense vectors. Peer of [`dot`].
///
/// **Not conjugated** — this is the load-bearing choice for the
/// complex symmetric pencil of dispersive FEM eigenproblems (ADR-0039
/// / spec §11). For real-valued inputs it agrees with the real
/// Euclidean dot bit-for-bit (imaginary parts are identically zero).
fn dot_complex(x: &[Complex64], y: &[Complex64]) -> Complex64 {
    debug_assert_eq!(x.len(), y.len());
    let mut acc = Complex64::new(0.0, 0.0);
    for (a, b) in x.iter().zip(y.iter()) {
        acc += a * b;
    }
    acc
}

/// In-place complex M-orthogonalisation in the *transposed* inner
/// product: `x ← x − Σ_j (e_j^T M x) e_j`. Peer of [`m_orthogonalize`].
fn m_orthogonalize_complex(
    x: &mut [Complex64],
    eig_vecs: &[Vec<Complex64>],
    m: &CsrMatrix<Complex64>,
) {
    for ej in eig_vecs {
        let mx = csr_matvec_complex(m, x);
        let coeff = dot_complex(ej, &mx);
        for (xi, eji) in x.iter_mut().zip(ej.iter()) {
            *xi -= coeff * eji;
        }
    }
}

/// In-place complex M-normalisation in the *transposed* inner product:
/// `x ← x / sqrt(x^T M x)`. The square root is the principal complex
/// square root (`num_complex::Complex::sqrt`). Peer of [`m_normalize`].
///
/// Returns an error if `x^T M x` is zero or non-finite (would indicate
/// deflation killed the seed; bump the seed or relax `tol`).
fn m_normalize_complex(
    x: &mut [Complex64],
    m: &CsrMatrix<Complex64>,
) -> Result<(), yee_core::Error> {
    let mx = csr_matvec_complex(m, x);
    let norm_sq = dot_complex(x, &mx);
    if norm_sq.norm() == 0.0 || !norm_sq.is_finite() {
        return Err(yee_core::Error::Numerical(format!(
            "ComplexInverseIterEigen: M-norm (transposed) collapsed to {norm_sq} during \
             deflation"
        )));
    }
    let inv_norm = Complex64::new(1.0, 0.0) / norm_sq.sqrt();
    for xi in x.iter_mut() {
        *xi *= inv_norm;
    }
    Ok(())
}

/// Seed vector for the `mode_idx`-th complex inverse-power iteration.
///
/// Matches the v0 real `seed_vector` exactly in the real-component
/// generation (so the real-and-complex paths agree bit-for-bit on
/// real input) and adds a small linearly-varying imaginary tail so
/// the seeds for different `mode_idx` are visibly different in every
/// coordinate even when the K, M pencil is purely real (in which
/// case the imaginary tail is killed by the first M-orthogonalise
/// step against the prior converged real eigenvectors). Peer of
/// [`seed_vector`].
fn seed_vector_complex(n: usize, mode_idx: usize) -> Vec<Complex64> {
    let mut x = vec![Complex64::new(0.0, 0.0); n];
    let phase = mode_idx as f64;
    for (i, xi) in x.iter_mut().enumerate() {
        let t = (i as f64 + 1.0) / (n as f64 + 1.0);
        let re = (1.0 + phase) * t + (1.0 + phase * 0.37).sin() * (t * 7.0).cos();
        *xi = Complex64::new(re, 0.0);
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a CSR matrix from a dense slice (row-major). Convenience
    /// for the unit tests; not exported.
    fn csr_from_dense(rows: usize, cols: usize, data: &[f64]) -> CsrMatrix<f64> {
        use nalgebra_sparse::coo::CooMatrix;
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

    /// Diagonal CSR matrix.
    fn diag_csr(diag: &[f64]) -> CsrMatrix<f64> {
        use nalgebra_sparse::coo::CooMatrix;
        let n = diag.len();
        let mut coo = CooMatrix::new(n, n);
        for (i, &d) in diag.iter().enumerate() {
            if d != 0.0 {
                coo.push(i, i, d);
            }
        }
        CsrMatrix::from(&coo)
    }

    /// Test 1 — known 4×4 dense pencil with eigenvalues {0.5, 1.2, 3.4, 7.8}.
    ///
    /// Construction: `K = diag(λ_i)`, `M = I`. The generalized
    /// eigenproblem `K e = λ M e` then has eigenvalues exactly
    /// `λ_i` with eigenvectors `e_i`. The diagonal form is the
    /// simplest published-reference pencil with prescribed spectrum.
    #[test]
    fn recovers_smallest_eigenvalue_on_known_dense_pencil() {
        let lambdas = [0.5, 1.2, 3.4, 7.8];
        let k = diag_csr(&lambdas);
        let m = diag_csr(&[1.0; 4]);

        let solver = InverseIterEigen::new(1000, 1e-10);
        let pairs = solver.solve(&k, &m, 3, 0.1).expect("solve");

        assert_eq!(pairs.k.len(), 3);
        // Eigenvalues sorted ascending; match the three smallest of the spectrum.
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
        // Monotone ascending.
        for w in pairs.k.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }

    /// Test 2 — scaled-identity pencil. `K = αI`, `M = βI` →
    /// every eigenvalue is exactly `α/β`. Shift-invert converges
    /// in one iteration (the operator is a scalar multiple of `I`).
    #[test]
    fn scaled_identity_pencil() {
        let alpha = 3.7;
        let beta = 1.5;
        let k = diag_csr(&[alpha; 5]);
        let m = diag_csr(&[beta; 5]);

        let solver = InverseIterEigen::new(50, 1e-10);
        let pairs = solver.solve(&k, &m, 3, 0.5).expect("solve");

        let expected = alpha / beta;
        for &k_sq in &pairs.k {
            assert!(
                (k_sq - expected).abs() < 1e-8,
                "expected {expected}, got {k_sq}"
            );
        }
    }

    /// Test 3 — `M`-orthogonality of returned eigenvectors.
    #[test]
    fn eigenvectors_m_orthogonal() {
        // Same 4×4 diagonal pencil as Test 1.
        let lambdas = [0.5, 1.2, 3.4, 7.8];
        let k = diag_csr(&lambdas);
        let m = diag_csr(&[1.0; 4]);

        let solver = InverseIterEigen::new(1000, 1e-10);
        let pairs = solver.solve(&k, &m, 3, 0.1).expect("solve");

        // e_i^T M e_j ≈ δ_{ij}.
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

    /// Test 4 — small 3-D scalar Laplacian on a 4³ grid (5-point
    /// stencil generalized to 7-point in 3-D). Smallest eigenvalue
    /// of `-Δ` with Dirichlet BCs on `[0,1]³` is `3π² ≈ 29.608`;
    /// the discrete 5-point stencil at `h = 1/(n+1)` underestimates
    /// it, but the solver must recover the smallest *discrete*
    /// eigenvalue within max_iter=100, tol=1e-8.
    #[test]
    fn converges_within_max_iter_for_3d_laplacian() {
        // Build the 7-point Dirichlet Laplacian on a 4×4×4 interior
        // grid. n = 64 unknowns. Stencil at interior point (i,j,k):
        //   6 x_{ijk} − x_{i±1,j,k} − x_{i,j±1,k} − x_{i,j,k±1} = h² λ x.
        // We omit the h² scaling — the solver should converge to the
        // discrete eigenvalue of the integer-Laplacian operator.
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

        // The smallest 5-point-Dirichlet Laplacian eigenvalue on a
        // 4-interior-point grid in 3-D is
        //   λ_111 = 2 (3 − cos(π/(nx+1)) · 3)
        //         = 6 − 6 cos(π/5)
        //         = 6 (1 − 0.809…) ≈ 1.1459…
        let pi = std::f64::consts::PI;
        let expected = 6.0 - 6.0 * (pi / (nx as f64 + 1.0)).cos();

        let solver = InverseIterEigen::new(100, 1e-8);
        let pairs = solver.solve(&k_mat, &m_mat, 1, 0.1).expect("solve");

        assert!(
            (pairs.k[0] - expected).abs() < 1e-5,
            "smallest 3-D Laplacian eigenvalue: expected {expected}, got {}",
            pairs.k[0]
        );
    }
}
