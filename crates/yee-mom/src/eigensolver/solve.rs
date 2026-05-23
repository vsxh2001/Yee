//! Dense generalized-eigensolve fallback for the 2-D cross-section
//! eigenproblem.
//!
//! **β-direct formulation (Phase 1.3.1.1 step 5.2).** Solves
//! `(k_0² T_ε − S) x = β² T_1 x` (real-symmetric for the lossless case),
//! where `S` is the curl-curl stiffness, `T_ε = ∫ε_r N·N` the ε_r-weighted
//! Nedelec mass, and `T_1 = ∫N·N` the **unweighted** mass — all produced
//! by [`super::assembly::assemble_transverse`]. The eigenvalue is `β²`
//! **directly**: this is the discrete form of the transverse Helmholtz
//! equation `∇×(1/μ_r ∇×E_t) = (k_0² ε_r − β²) E_t`, where ε_r appears only
//! on the `k_0²` side. The earlier `S x = k_c² T_ε x` / `β² = k_0² − k_c²`
//! extraction (vacuum `k_0`, ε_r-weighted RHS) was correct only for
//! `ε_r ≡ 1` — for any `ε_r ≠ 1` it under-counted the dielectric; the
//! β-direct form fixes that.
//!
//! **Spurious-mode handling.** First-order Nedelec edge elements admit
//! a large gradient null-space (every nodal-Lagrange gradient is in
//! their span and lies in the kernel of `curl`). These curl-free modes
//! satisfy `S x ≈ 0`, so in the β-direct pencil `(k_0² T_ε − S) x = β² T_1 x`
//! they land at `β² ≈ k_0² ⟨ε_r⟩` (the **top** of the spectrum, not the
//! bottom). They are filtered by their vanishing **cutoff Rayleigh
//! quotient** `k_c² := (xᵀS x)/(xᵀT_ε x) ≈ 0` (the same physical quantity
//! the old pencil used as its eigenvalue), and the physical dominant mode
//! is then the **largest β²** (equivalently the lowest cutoff `k_c²`)
//! among the survivors.
//!
//! **Numerical method.** `T_1` is SPD (a Gram matrix with ε_r ≡ 1 > 0), so
//! reduce `(k_0² T_ε − S) x = β² T_1 x` to a standard symmetric problem
//! `M y = β² y` via the Cholesky factor `T_1 = L Lᵀ`, then
//! `M = L⁻¹ (k_0² T_ε − S) L⁻ᵀ` (the operator `k_0² T_ε − S` is symmetric
//! **indefinite**, which the symmetric QR handles). Solve `M` with
//! [`nalgebra::SymmetricEigen`] (tridiagonal QR). `O(n³)` flops; viable
//! only for coarse cross sections (≤ a few hundred DoF). Sparse
//! shift-and-invert with `arpack-rs` is Phase 1.3.1.1 step 4 and is
//! escape-hatched away from the lossless TE10 validation gate.
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

use super::assembly::{AssembledMixed, AssembledTransverse};

/// Solved-eigensolution payload returned by [`solve_dense`].
///
/// Retained as the transverse-only reference path (and exercised by the
/// homogeneous-guide regression tests / DoD-V1 canary) after
/// [`crate::ports::NumericalCrossSection::solve`] switched to the mixed
/// [`solve_dense_mixed`] path in Phase 1.3.1.1 step 5.
#[allow(dead_code)]
pub(crate) struct EigenSolution {
    /// `β²` for the dominant propagating mode at the supplied frequency,
    /// the **direct eigenvalue** of `(k_0² T_ε − S) x = β² T_1 x` (Phase
    /// 1.3.1.1 step 5.2). Stored as `Complex64` to keep the API
    /// future-proof for the lossy / complex-symmetric path; the lossless
    /// path always returns a real value.
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

/// Solve the β-direct generalized eigenproblem `(k_0² T_ε − S) x = β² T_1 x`
/// densely on the lossless / real-symmetric path and return `β²` for the
/// dominant propagating mode at `freq_hz` (Phase 1.3.1.1 step 5.2).
///
/// **Size limit.** The implementation builds dense `n×n` matrices and
/// runs an `O(n³)` Cholesky + symmetric-tridiagonal QR. Callers must
/// keep `n` below a few hundred. The eigensolver's only validation
/// case (WR-90 TE10 on a ~72-triangle mesh) lands at `n ≈ 84` interior
/// edges, well inside this envelope.
///
/// Retained as the transverse-only reference path after
/// [`crate::ports::NumericalCrossSection::solve`] switched to the mixed
/// [`solve_dense_mixed`] in Phase 1.3.1.1 step 5; exercised by the
/// homogeneous-guide regression tests (the DoD-V1 canary asserts the
/// mixed solve reproduces this path's β to machine precision).
#[allow(dead_code)]
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
    let t_eps_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.t[(i, j)].re);
    let t1_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.t1[(i, j)].re);

    // β-direct LHS operator: K := k_0² T_ε − S (symmetric **indefinite**).
    let omega = std::f64::consts::TAU * freq_hz;
    let k0 = omega / yee_core::units::C0;
    let k0_sq = k0 * k0;
    let k_op = k0_sq * &t_eps_re - &s_re;

    // Cholesky factor T_1 = L Lᵀ. The unweighted Nedelec mass is SPD (a
    // Gram matrix with ε_r ≡ 1 > 0) on the interior-edge DoF set.
    let chol_t1 = t1_re.clone().cholesky().ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: unweighted Nedelec mass T_1 is not SPD on the interior-edge DoF set"
                .into(),
        )
    })?;
    let l = chol_t1.l();
    // Compute M = L⁻¹ K L⁻ᵀ via two triangular solves preserving symmetry.
    // Step 1: solve L Y = K. Step 2: solve L Z = Yᵀ → M = Zᵀ.
    let y = l.clone().solve_lower_triangular(&k_op).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: L Y = K solve failed (L from Cholesky should be non-singular)".into(),
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
    // the symmetric QR path, which is bulletproof for symmetric pencils
    // (here `K = k_0² T_ε − S` is symmetric indefinite).
    let m_sym = 0.5 * (&m + m.transpose());

    let eig = SymmetricEigen::new(m_sym);
    // SymmetricEigen returns real eigenvalues `β²` (T::RealField = f64).
    // The curl-free gradient null-space satisfies `S x ≈ 0`, so those
    // spurious modes land at `β² ≈ k_0² ⟨ε_r⟩` (the TOP of the spectrum).
    // Filter them by their vanishing cutoff Rayleigh quotient
    // `k_c² := (xᵀ S x)/(xᵀ T_ε x)` — the same physical quantity the old
    // pencil used as its eigenvalue — then take the **largest** surviving
    // `β²` (the physical dominant mode = lowest cutoff). The eigenvector
    // is recovered per-candidate by the standard-form back-transform
    // `x = L⁻ᵀ y`.
    let l_t = l.transpose();
    // Cutoff-Rayleigh floor relative to the largest cutoff seen, mirroring
    // the old pencil's `k_c² ≤ max_eig · 1e-6` gradient filter.
    let mut k_c_sq_max = 0.0_f64;
    let mut cands: Vec<(f64, f64, usize)> = Vec::new(); // (β², k_c², col)
    for (col, &beta_sq) in eig.eigenvalues.iter().enumerate() {
        if !beta_sq.is_finite() {
            continue;
        }
        let y_col = eig.eigenvectors.column(col).clone_owned();
        let Some(x) = l_t.solve_upper_triangular(&y_col) else {
            continue;
        };
        let s_energy = (x.transpose() * &s_re * &x)[(0, 0)];
        let t_eps_energy = (x.transpose() * &t_eps_re * &x)[(0, 0)];
        if t_eps_energy <= 0.0 {
            continue;
        }
        let k_c_sq = s_energy / t_eps_energy;
        if k_c_sq > k_c_sq_max {
            k_c_sq_max = k_c_sq;
        }
        cands.push((beta_sq, k_c_sq, col));
    }
    let spurious_floor = k_c_sq_max * 1e-6;

    // Among modes with a genuine (above-floor) cutoff, pick the largest
    // β² = lowest cutoff = dominant propagating mode.
    let mut dominant: Option<(f64, usize)> = None;
    for &(beta_sq, k_c_sq, col) in &cands {
        if k_c_sq <= spurious_floor {
            continue; // curl-free gradient null-space mode
        }
        dominant = Some(match dominant {
            None => (beta_sq, col),
            Some((curr, curr_col)) => {
                if beta_sq > curr {
                    (beta_sq, col)
                } else {
                    (curr, curr_col)
                }
            }
        });
    }
    let (beta_sq_re, dom_col) = dominant.ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver: no mode with a strictly-positive cutoff above the spurious-mode floor"
                .into(),
        )
    })?;

    if beta_sq_re <= 0.0 {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver: dominant mode is evanescent at {freq_hz} Hz (β² = {beta_sq_re})"
        )));
    }

    // Recover the eigenvector in the **original** (T_1-weighted) basis.
    // SymmetricEigen solves `M y = β² y` with `M = L⁻¹ K L⁻ᵀ`; the
    // generalized-problem eigenvector satisfies `K x = β² T_1 x` with
    // `x = L⁻ᵀ y`. So back-transform with one upper-triangular solve
    // (`l_t` already formed above for the candidate eigenvector recovery).
    let y_dom = eig.eigenvectors.column(dom_col).clone_owned();
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

/// Solved-eigensolution payload returned by [`solve_dense_mixed`].
pub(crate) struct MixedEigenSolution {
    /// `β²` for the dominant quasi-TEM mode at the supplied frequency, the
    /// **direct eigenvalue** of `(k_0² B − A) x = β² B_1 x` (Phase 1.3.1.1
    /// step 5.2). Real-valued on the lossless path.
    pub beta_sq: Complex64,
    /// Transverse `E_t` eigenvector components in the **interior-edge
    /// DoF** ordering of [`AssembledMixed::interior_to_global_edges`]
    /// (length `n_t`).
    pub e_t: Vec<Complex64>,
    /// Longitudinal `E_z` eigenvector components in the **interior-vertex
    /// DoF** ordering of [`AssembledMixed::interior_to_global_verts`]
    /// (length `n_z`).
    pub e_z: Vec<Complex64>,
}

/// Solve the mixed `(E_t, E_z)` β-direct block pencil
/// `(k_0² B − A) x = β² B_1 x` densely on the lossless / real-symmetric
/// path and return `β²` for the dominant quasi-TEM mode at `freq_hz`
/// (Phase 1.3.1.1 step 5.2). `A` is the block-stiffness, `B` the
/// ε_r-weighted block-mass + coupling, and `B_1` the **unweighted**
/// block-mass + coupling (see [`AssembledMixed`] for why the `−β²` term
/// carries the unweighted mass).
///
/// **Method (Phase 1.3.1.1 step 5.2) — cutoff-pencil select + β-direct
/// Rayleigh quotient.** The physical mode is the dominant guided mode of
/// the cutoff pencil `A x = k_c² B x` (the step-5 pencil, unchanged +
/// validated), and its propagation constant is the **β-direct Rayleigh
/// quotient** of that eigenvector:
///
/// ```text
///   β² = R(x) = (xᵀ (k_0² B − A) x) / (xᵀ B_1 x)
/// ```
///
/// 1. **Select on `A x = k_c² B x`.** Reduce to `B⁻¹A` (one LU; `B`
///    nonsingular though indefinite) and take `complex_eigenvalues`. The
///    curl-free gradient null-space sits cleanly at `k_c² ≈ 0` (rejected
///    by the `k_c² ≤ max|k_c²| · 1e-6` floor) and the **transverse-energy
///    filter** (`‖e_t‖²/‖x‖² ≥` [`TRANSVERSE_ENERGY_FLOOR`]) removes
///    E_z-dominated contamination. The dominant guided mode is the
///    **smallest** valid `k_c²`; its eigenvector (with the genuine E_z
///    content of an inhomogeneous hybrid mode) is recovered by inverse
///    iteration on `(A − σ B)` — the step-5 recovery, *not* a
///    smallest-singular-vector null-space (which grabbed a spurious
///    `E_t`-only gradient direction).
/// 2. **Extract β² via the β-direct Rayleigh quotient `R(x)`.** Since
///    `A x = k_c² B x`, `R(x) = (k_0² − k_c²)·⟨ε_r⟩` with the
///    mode-resolved `⟨ε_r⟩ = (xᵀ B x)/(xᵀ B_1 x)`. This is **exact** on a
///    uniformly-filled guide (`B = ε_r B_1` ⇒ `⟨ε_r⟩ = ε_r` and
///    `R = ε_r k_0² − k_c0²` with `k_c0² = (xᵀA x)/(xᵀB_1 x)` the
///    unweighted cutoff — the analytic
///    `β = √(ε_r k_0² − (π/a)²)`, DoD-1), and reduces to `k_0² − k_c²` on
///    a homogeneous guide (`B_1 ≡ B`, the DoD-V1 canary). The step-5
///    `β² = k_0² − k_c²` with vacuum `k_0` dropped the `⟨ε_r⟩` factor and
///    so under-counted any `ε_r ≠ 1` fill.
///
/// **Why the Rayleigh quotient rather than the β-direct pencil's direct
/// eigenvalue (spec §3 option A).** Solving `K x = β² B_1 x` directly was
/// tried and *drifts off the physical mode*: its dominant eigenvalue near
/// the physical `β²` belongs to a spurious `E_z ≈ 0` branch (the curl-free
/// gradient null-space lands at `β² ≈ k_0² ⟨ε_r⟩`, interleaved with the
/// physical mode, and `(K − σ B_1) ≈ −A` is singular there). The cutoff
/// pencil cleanly isolates the gradient cluster at `k_c² ≈ 0`, so it
/// reliably picks the physical hybrid mode — verified by `‖E_z‖/‖E_t‖`
/// matching the published LSM-to-y reference for the horizontal slab. The
/// Rayleigh quotient on that correctly-selected eigenvector is the right
/// β² for the physical mode and is exact on the DoD-1 uniform anchor;
/// option A's mode-dependence concern is moot because β² is evaluated on
/// the physically-selected mode, not a generic vector.
///
/// `O(n³)` with `n = n_t + n_z`; returns in milliseconds at the
/// validation `n ≈ 121`. Revisit for a sparse / large-DoF
/// symmetric-indefinite solver.
pub(crate) fn solve_dense_mixed(
    asm: &AssembledMixed,
    freq_hz: f64,
) -> Result<MixedEigenSolution, yee_core::Error> {
    let n_t = asm.n_t;
    let n = n_t + asm.n_z;
    debug_assert_eq!(n, asm.a.nrows(), "mixed pencil size must be n_t + n_z");
    if n == 0 || n_t == 0 {
        return Err(yee_core::Error::Numerical(
            "eigensolver(mixed): empty DoF set (no interior edges?)".into(),
        ));
    }

    // Lossless-only, mirroring `solve_dense`.
    let im_norm = asm
        .a
        .iter()
        .chain(asm.b.iter())
        .chain(asm.b1.iter())
        .map(|z| z.im.abs())
        .fold(0.0_f64, f64::max);
    let re_norm = asm
        .a
        .iter()
        .chain(asm.b.iter())
        .chain(asm.b1.iter())
        .map(|z| z.re.abs())
        .fold(0.0_f64, f64::max);
    if im_norm > 1e-9 * re_norm.max(1.0) {
        return Err(yee_core::Error::Unimplemented(
            "eigensolver(mixed): complex (lossy) ε_r / μ_r is Phase 1.3.1.2; current path is lossless only",
        ));
    }

    let a_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.a[(i, j)].re);
    let b_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.b[(i, j)].re);
    let b1_re = DMatrix::<f64>::from_fn(n, n, |i, j| asm.b1[(i, j)].re);

    let omega = std::f64::consts::TAU * freq_hz;
    let k0 = omega / yee_core::units::C0;
    let k0_sq = k0 * k0;
    // β-direct LHS operator: K := k_0² B − A. The β-direct generalized
    // eigenproblem is `K x = β² B_1 x`.
    let k_op = k0_sq * &b_re - &a_re;

    // ── Stage 1: select the dominant guided mode on the cutoff pencil ──
    // `A x = k_c² B x` (the step-5 pencil). This is unchanged from step 5
    // and is fast + validated: its `k_c²` eigenvalue is the physically
    // meaningful cutoff (curl energy / ε-weighted field energy), so the
    // gradient null-space sits at `k_c² ≈ 0` (easy to filter) and the
    // transverse-energy filter removes E_z-dominated contamination.
    //
    // Why use the cutoff pencil for SELECTION rather than the β-direct
    // pencil directly: in the β-direct pencil the curl-free gradient
    // modes land at `β² ≈ k_0² ⟨ε_r⟩` (the TOP of the spectrum, mixed in
    // with the physical mode) and their shifted matrix `(K − σ B_1) ≈ −A`
    // is singular (A's gradient null-space), so a blind inverse-iteration
    // sweep over every β-direct eigenvalue thrashes on the gradient
    // cluster. Selecting on the cutoff pencil keeps the gradient cluster
    // cleanly at `k_c² ≈ 0` and runs inverse iteration only on genuine
    // candidates.
    let b_lu = b_re.clone().lu();
    let binv_a = b_lu.solve(&a_re).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver(mixed): block mass matrix B is singular on the interior DoF set".into(),
        )
    })?;
    let kc_evals = binv_a.complex_eigenvalues();
    let max_abs = kc_evals
        .iter()
        .map(|z| z.norm())
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let spurious_floor = max_abs * 1e-6;

    // Among real, above-floor cutoffs whose eigenvector is
    // transverse-energy-dominated, pick the SMALLEST k_c² (dominant guided
    // mode = lowest cutoff). Eigenvector recovered by inverse iteration on
    // `(A − σ B)` (σ off the cutoff eigenvalue), which is well away from
    // the gradient null-space and converges fast.
    let mut best: Option<(f64, Vec<f64>)> = None;
    for ev in kc_evals.iter() {
        if ev.im.abs() > 1e-6 * ev.re.abs().max(1.0) {
            continue;
        }
        let k_c_sq = ev.re;
        if !k_c_sq.is_finite() || k_c_sq <= spurious_floor {
            continue; // curl-free gradient null-space
        }
        let Some(x) = inverse_iterate(&a_re, &b_re, k_c_sq) else {
            continue;
        };
        let total: f64 = x.iter().map(|&v| v * v).sum();
        if total <= 0.0 {
            continue;
        }
        let trans: f64 = (0..n_t).map(|i| x[i] * x[i]).sum();
        if trans / total < TRANSVERSE_ENERGY_FLOOR {
            continue;
        }
        let take = match &best {
            None => true,
            Some((curr, _)) => k_c_sq < *curr,
        };
        if take {
            best = Some((k_c_sq, x));
        }
    }
    let (_k_c_sq, x_sel) = best.ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver(mixed): no transverse-energy-dominated real k_c² above the spurious floor"
                .into(),
        )
    })?;

    // ── Stage 2: β² = the β-direct Rayleigh quotient on the selected mode ──
    // β² = R(x) = (xᵀ(k_0² B − A)x)/(xᵀ B_1 x) (Phase 1.3.1.1 step 5.2).
    // Since `A x_sel = k_c² B x_sel`, this is `(k_0² − k_c²)·⟨ε_r⟩` with the
    // mode-resolved `⟨ε_r⟩ = (xᵀB x)/(xᵀB_1 x)`: exact on a uniform fill
    // (`⟨ε_r⟩ = ε_r`, the DoD-1 analytic anchor) and `= k_0² − k_c²` on a
    // homogeneous guide (`B_1 ≡ B`, the DoD-V1 canary). The step-5
    // `β² = k_0² − k_c²` with vacuum `k_0` dropped the `⟨ε_r⟩` factor.
    // Evaluating R on the cutoff pencil's correctly-selected hybrid mode
    // avoids spec §3 option A's drift onto the spurious `E_z≈0` β-direct
    // branch (see the function docstring).
    let rayleigh_beta_sq = |x: &[f64]| -> Option<f64> {
        let xv = DMatrix::<f64>::from_column_slice(n, 1, x);
        let num = (xv.transpose() * &k_op * &xv)[(0, 0)];
        let den = (xv.transpose() * &b1_re * &xv)[(0, 0)];
        (den.abs() > 0.0).then_some(num / den)
    };
    let x_dom = x_sel;
    let beta_sq_re = rayleigh_beta_sq(&x_dom).ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver(mixed): β-direct Rayleigh quotient has a zero B_1-norm denominator".into(),
        )
    })?;

    if beta_sq_re <= 0.0 {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver(mixed): dominant mode is evanescent at {freq_hz} Hz (β² = {beta_sq_re})"
        )));
    }

    // Fix the global sign deterministically off the transverse block:
    // largest-magnitude E_t component positive (matches `solve_dense`).
    let argmax = (0..n_t)
        .max_by(|&a, &b| {
            x_dom[a]
                .abs()
                .partial_cmp(&x_dom[b].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(0);
    let sign = if x_dom[argmax] < 0.0 { -1.0 } else { 1.0 };

    let e_t: Vec<Complex64> = (0..n_t)
        .map(|i| Complex64::new(sign * x_dom[i], 0.0))
        .collect();
    let e_z: Vec<Complex64> = (n_t..n)
        .map(|i| Complex64::new(sign * x_dom[i], 0.0))
        .collect();

    Ok(MixedEigenSolution {
        beta_sq: Complex64::new(beta_sq_re, 0.0),
        e_t,
        e_z,
    })
}

/// Recover the generalized eigenvector of `A x = λ B x` for the
/// eigenvalue nearest `lambda` by **inverse iteration** on the shifted
/// pencil `(A − σ B)` with `σ = lambda · (1 + δ)` (a small relative
/// shift `δ` off the eigenvalue so the shifted matrix is non-singular
/// and the iteration converges to *this* eigenvector rather than
/// stalling on an exact null space).
///
/// Returns the converged eigenvector (length `A.nrows()`), or `None` if
/// the shifted-pencil LU is singular or the iterate collapses to zero.
/// A handful of iterations is sufficient because the shift sits right on
/// top of the target eigenvalue (the dominant amplification factor
/// `1/(λ_i − σ)` is enormous for the nearest eigenvalue and small for
/// all others).
///
/// This replaces a smallest-singular-vector null-space recovery that
/// failed on the mixed pencil: `A_tt` carries the Nedelec curl
/// gradient null-space, so `(A − λ B)` has many spurious near-null
/// directions in the `E_t`-only subspace, and the global smallest
/// singular vector picked one of those (`E_z ≡ 0`) instead of the
/// physical mode.
///
/// **Assumes a simple (non-degenerate) eigenvalue.** For a degenerate
/// `λ` the iteration converges to *some* vector in the eigenspace (an
/// arbitrary combination of the degenerate modes), not a canonical
/// basis. The dominant guided mode is well-separated for the validation
/// cross sections; resolving an explicit degenerate cluster M-orthonormal
/// is a later step (cf. the FEM-side `LobpcgEigen` degenerate handling).
fn inverse_iterate(a: &DMatrix<f64>, b: &DMatrix<f64>, lambda: f64) -> Option<Vec<f64>> {
    let n = a.nrows();
    // Relative shift off the eigenvalue. Large enough to keep (A − σB)
    // well away from exact singularity, small enough that the target
    // eigenvalue still dominates the inverse-iteration amplification.
    let sigma = lambda * (1.0 + 1e-6) + 1e-6;
    let shifted = a - sigma * b;
    let lu = shifted.lu();

    // Seed with a deterministic non-symmetric vector (avoids accidental
    // orthogonality to the target eigenvector).
    let mut z = DMatrix::<f64>::from_fn(n, 1, |i, _| 1.0 + (i as f64) * 0.001);
    let mut last_norm = 0.0;
    for _ in 0..50 {
        // z_{k+1} = (A − σB)⁻¹ (B z_k)
        let rhs = b * &z;
        let y = lu.solve(&rhs)?;
        let norm = y.norm();
        if !norm.is_finite() || norm == 0.0 {
            return None;
        }
        z = y / norm;
        // Converged when the normalized iterate stops moving (Rayleigh
        // amplification has saturated on the dominant eigenvector).
        if (norm - last_norm).abs() <= 1e-10 * norm {
            break;
        }
        last_norm = norm;
    }
    Some(z.column(0).iter().copied().collect())
}

/// Minimum transverse-block **Euclidean** energy fraction
/// `‖e_t‖² / ‖x‖²` for a candidate to count as a physical quasi-TEM /
/// quasi-TE mode in [`solve_dense_mixed`]. Modes below this carry their
/// energy in the longitudinal `E_z` / nodal-gradient block and are
/// rejected as spurious. The fraction is the plain `ℓ²` (Euclidean)
/// ratio of the eigenvector's components, **not** a `B`-norm: `B` is
/// indefinite here (it carries the edge-node coupling), so a true
/// `B`-norm could be negative and is not a meaningful energy. `0.5` is a
/// wide margin: the homogeneous-guide dominant mode sits at fraction
/// `1.0` (E_z ≡ 0) and inhomogeneous quasi-TEM modes remain strongly
/// transverse-dominated for the validation cross sections.
const TRANSVERSE_ENERGY_FLOOR: f64 = 0.5;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigensolver::{
        assembly::{assemble_mixed, assemble_transverse},
        mesh::EdgeTable,
    };
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

    #[test]
    fn mixed_solve_reproduces_transverse_beta_on_homogeneous_guide() {
        // DoD-V1 canary (unit-level): on a homogeneous air-filled WR-90
        // cross-section the longitudinal E_z block contributes zero to
        // the dominant mode (E_z ≡ 0), so the mixed pencil must
        // reproduce the transverse-only β to high precision. A gross
        // block sign/placement error breaks this immediately.
        let a = 22.86e-3;
        let b = 10.16e-3;
        let freq_hz = 10e9;
        let mesh = rectangular_mesh(a, b, 6, 6);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);

        let asm_t = assemble_transverse(&mesh, &eps, &mu, &table);
        let sol_t = solve_dense(&asm_t, freq_hz).unwrap();
        let beta_t = sol_t.beta_sq.re.sqrt();

        let asm_m = assemble_mixed(&mesh, &eps, &mu, &table);
        let sol_m = solve_dense_mixed(&asm_m, freq_hz).unwrap();
        let beta_m = sol_m.beta_sq.re.sqrt();

        let rel = (beta_m - beta_t).abs() / beta_t;
        eprintln!(
            "homogeneous WR-90 β: transverse {beta_t:.6}, mixed {beta_m:.6}, rel err {rel:.3e}"
        );
        assert!(
            rel < 1e-3,
            "mixed β {beta_m} must reproduce transverse β {beta_t} within 0.1% (rel {rel:.3e})"
        );

        // On the homogeneous guide the dominant mode is purely
        // transverse: the recovered E_z block must be ~zero.
        let ez_norm: f64 = sol_m.e_z.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        let et_norm: f64 = sol_m.e_t.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        eprintln!("homogeneous WR-90: ‖E_z‖ = {ez_norm:.3e}, ‖E_t‖ = {et_norm:.3e}");
        assert!(
            ez_norm < 1e-6 * et_norm.max(1e-30),
            "homogeneous-guide E_z block should be ~zero: ‖E_z‖={ez_norm}, ‖E_t‖={et_norm}"
        );
    }

    /// Horizontal-slab WR-90 mesh: lower-y half tagged 1, rest tagged 0.
    fn horizontal_slab_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
        let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                vertices.push([a * (i as f64) / (nx as f64), b * (j as f64) / (ny as f64)]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut triangles = Vec::with_capacity(2 * nx * ny);
        let mut tags = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let v00 = idx(i, j);
                let v10 = idx(i + 1, j);
                let v11 = idx(i + 1, j + 1);
                let v01 = idx(i, j + 1);
                let yc = b * ((j as f64) + 0.5) / (ny as f64);
                let tag = if yc < b / 2.0 { 1u32 } else { 0u32 };
                triangles.push([v00, v10, v11]);
                tags.push(tag);
                triangles.push([v00, v11, v01]);
                tags.push(tag);
            }
        }
        TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
    }

    #[test]
    fn zeroing_coupling_changes_hybrid_mode() {
        // Step-5-review P1-1 load-bearing guard: on a horizontal-slab
        // (hybrid-mode) guide, zeroing ONLY the off-diagonal coupling
        // block B_tz (= B_ztᵀ) must measurably change the dominant
        // eigenpair — proving the coupling participates and pinning that
        // it is non-trivially placed. (On a homogeneous or vertical-slab
        // guide this delta is zero because the dominant mode has E_z = 0;
        // the horizontal slab is exactly where the coupling bites.)
        let a = 22.86e-3;
        let b = 10.16e-3;
        let freq_hz = 10e9;
        let mesh = horizontal_slab_mesh(a, b, 8, 8);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        eps.insert(1u32, Complex64::new(10.2, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        mu.insert(1u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);

        let asm = assemble_mixed(&mesh, &eps, &mu, &table);
        let beta_full = solve_dense_mixed(&asm, freq_hz).unwrap().beta_sq.re.sqrt();

        // Build a coupling-zeroed copy: drop the edge↔vertex blocks of B.
        let mut asm_nc = assemble_mixed(&mesh, &eps, &mu, &table);
        let n_t = asm_nc.n_t;
        let n = n_t + asm_nc.n_z;
        for i in 0..n_t {
            for j in n_t..n {
                asm_nc.b[(i, j)] = Complex64::new(0.0, 0.0);
                asm_nc.b[(j, i)] = Complex64::new(0.0, 0.0);
            }
        }
        let beta_nc = solve_dense_mixed(&asm_nc, freq_hz)
            .unwrap()
            .beta_sq
            .re
            .sqrt();

        let rel = (beta_full - beta_nc).abs() / beta_full;
        eprintln!(
            "coupling delta (horizontal slab): β with coupling {beta_full:.4}, \
             β without {beta_nc:.4}, rel Δ {rel:.3e}"
        );
        // The coupling shifts β by ~1% here; require a clearly non-zero,
        // not-floating-point-noise delta.
        assert!(
            rel > 1e-4,
            "zeroing B_tz must change the hybrid-mode β (rel Δ {rel:.3e}); \
             a zero delta means the coupling block is inert/misplaced"
        );
    }
}
