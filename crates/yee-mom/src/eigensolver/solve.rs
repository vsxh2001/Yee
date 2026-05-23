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

use faer::linalg::solvers::SolveCore;
use faer::sparse::linalg::solvers::Lu;
use faer::sparse::{SparseColMat, Triplet};
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
/// `(k_0² B − A) x = β² B_1 x` for `β²` and the eigenvector of the
/// dominant quasi-TEM mode at `freq_hz` (Phase 1.3.1.1 step 5.3). `A` is
/// the block-stiffness, `B` the ε_r-weighted block-mass + coupling, and
/// `B_1` the **unweighted** block-mass + coupling (see [`AssembledMixed`]
/// for why the `−β²` term carries the unweighted mass).
///
/// **Method (Phase 1.3.1.1 step 5.3) — direct β-direct sparse
/// shift-and-invert.** This is the production path consumed by
/// [`crate::ports::NumericalCrossSection::solve`]. It supersedes the
/// step-5.2 *hybrid* (cutoff-pencil select + β-direct Rayleigh quotient
/// on the *cutoff*-pencil eigenvector), which carried a mesh-stable
/// eigenvector-mismatch bias for inhomogeneous fills because its β² was a
/// Rayleigh quotient on the wrong (cutoff-pencil) eigenvector. The direct
/// solve recovers the **true β-direct eigenvector**, so its β² is exact
/// for that mode.
///
/// 1. **Physics-informed shift `σ₀`.** Run the dominant-mode selection
///    ([`select_dominant_mode`] — the cutoff-pencil candidate screen, then
///    the Phase 1.3.1.1 step 5.6 highest-β-direct-Rayleigh-quotient pick)
///    to obtain the dominant guided mode's eigenvector `x_c` and its
///    β-direct Rayleigh quotient
///    `σ₀ = R(x_c) = (x_cᵀ(k_0²B−A)x_c)/(x_cᵀB_1 x_c)`. Selecting on the
///    highest `R(x)` (= highest ε_eff) lands the physical dominant mode
///    even when the enlarged p=2 curl-free cluster contributes spurious
///    low-`k_c²` transverse-dominated candidates; the shift then targets
///    that eigenpair, which sits well below the spurious curl-free
///    gradient cluster at `β² ≈ k_0²⟨ε_r⟩`.
/// 2. **Sparse shift-and-invert.** Build `(K − σ₀ B_1)` with
///    `K = k_0² B − A` as a faer [`SparseColMat`] (mirroring the yee-fem
///    `build_shifted` pattern), factor once via `sp_lu`, and inverse-
///    iterate `z ← (K − σ₀B_1)⁻¹ B_1 z` to the eigenpair *nearest* `σ₀`
///    — the true β-direct eigenvector `x_true`. Then `β² = R(x_true)` is
///    exact (the Rayleigh quotient is stationary *at* its own eigenvector,
///    so there is no cutoff-vs-β-direct mismatch). Targeting `σ₀` rather
///    than sweeping every eigenvalue avoids the spurious-cluster thrash
///    that sank spec §3 "option A" (a blind global β-direct solve grabs
///    the `E_z ≈ 0` gradient branch, against which `(K − σB_1) ≈ −A` is
///    singular).
/// 3. **Spurious-mode screen.** The converged mode is rejected if it is
///    not transverse-energy-dominated (`‖e_t‖²/‖x‖² ≥`
///    [`TRANSVERSE_ENERGY_FLOOR`]); on the validation cross-sections the
///    direct solve lands the physical mode cleanly (`‖e_t‖²/‖x‖² ≈ 1`).
///
/// **Validation disposition (Phase 1.3.1.1 step 5.3, ADR-0054).** At the
/// FR-4 contrast (ε_r=4.4) the direct β lands within ≤5 % of the verified
/// `reference::slab_loaded_beta` (the §4 published-benchmark closure). At
/// the high-contrast stretch (ε_r=10.2) the direct β improves on the
/// hybrid only ~1 % (483 → 486 rad/m vs the reference 583) and **plateaus
/// under mesh refinement** (8×8 → 16×16 within ~0.6 %): this is decisive
/// evidence that the residual there is **discretization-dominated**
/// (first-order Nedelec/nodal elements under-resolving the field peak at
/// the high-contrast interface), *not* the eigenvector mismatch the
/// direct solve removes. Closing the ε_r=10.2 gap needs higher-order
/// elements — queued to step-5.4. See
/// `tests/eigensolver_inhomogeneous.rs` for both gates.
///
/// On a **homogeneous** guide (`B_1 ≡ B`) the β-direct pencil reduces to
/// the cutoff form (`β² = k_0² − k_c²`), so the WR-90 TE10 canary and the
/// uniform-fill analytic anchor are preserved exactly.
///
/// The sparse LU is `O(n^{3/2})`-ish on these 2-D meshes and the inverse
/// iteration is a handful of triangular back-substitutions, so the solve
/// scales to far finer meshes than the step-5.2 dense `O(n³)` path
/// (which becomes minutes-scale by 24×24, `n ≈ 2100`).
///
/// The function name is retained from step 5.2 to keep the public
/// boundary in [`crate::ports`] stable; the implementation is now the
/// sparse direct solve above (the dense Rayleigh-quotient is kept only as
/// the small-`n` reference in [`solve_dense_mixed_rq`]).
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

    let (a_re, b_re, b1_re) = mixed_real_blocks(asm)?;

    let omega = std::f64::consts::TAU * freq_hz;
    let k0 = omega / yee_core::units::C0;
    let k0_sq = k0 * k0;
    // β-direct LHS operator: K := k_0² B − A. The β-direct generalized
    // eigenproblem is `K x = β² B_1 x`.
    let k_op = k0_sq * &b_re - &a_re;

    // ── Step 1: gather propagating cutoff-pencil candidate shifts ──
    // The cutoff pencil isolates the curl-free gradient null-space at
    // `k_c² ≈ 0` (floored out), and each surviving candidate's β-direct
    // Rayleigh quotient `R(x_c)` is a physics-informed shift. At p=2 the
    // physical dominant mode's cutoff-pencil eigenvector is NOT transverse-
    // dominated (the enlarged curl-free cluster mixes it with the E_z /
    // gradient block), so we do NOT pre-filter on the cutoff-pencil
    // transverse fraction — every propagating candidate is kept as a shift
    // and the transverse screen is applied to the converged β-direct
    // eigenvector below. Candidates are ordered by `R(x_c)` descending
    // (the dominant quasi-TEM mode = slowest wave = highest ε_eff first).
    let candidates = cutoff_candidates(&a_re, &b_re, &k_op, &b1_re, n, k0_sq)?;
    if candidates.is_empty() {
        return Err(yee_core::Error::Numerical(
            "eigensolver(mixed): no propagating cutoff candidate above the spurious floor".into(),
        ));
    }

    // ── Step 2: shift-and-invert from each candidate; keep the highest-β²
    // transverse-dominated converged β-direct eigenpair (Phase 1.3.1.1 step
    // 5.6). The β-direct shift-and-invert recovers the TRUE eigenvector for
    // the mode nearest the shift and screens it for transverse-dominance on
    // that converged eigenvector; among the modes that pass, the physical
    // dominant mode is the largest β² (highest ε_eff). This rejects the
    // p=2 spurious gradient-cluster captures (low ε_eff) that defeated the
    // smallest-cutoff selection. ──
    let mut best: Option<(f64, Vec<f64>)> = None; // (β²_true, x_true)
    let n_candidates = candidates.len();
    let mut n_converged = 0usize; // shift-inverts that yielded a valid transverse mode
    for (sigma0, x_c) in &candidates {
        match beta_direct_shift_invert(&k_op, &b1_re, n, n_t, *sigma0, Some(x_c)) {
            Ok((beta_sq, x_true)) => {
                if beta_sq <= 0.0 || !beta_sq.is_finite() {
                    continue;
                }
                n_converged += 1;
                let take = match &best {
                    None => true,
                    Some((curr, _)) => beta_sq > *curr,
                };
                if take {
                    best = Some((beta_sq, x_true));
                }
            }
            // A candidate whose shift-invert converges to a non-transverse
            // (spurious E_z / gradient) mode, fails to converge, or is
            // evanescent is skipped — another candidate targets the
            // physical mode.
            Err(_) => continue,
        }
    }

    // Failure here is rare (no propagating transverse mode found at all); the
    // candidate-tried / converged counts make the message grep-able for CI
    // triage rather than requiring a debugger re-run (step-5.6 review P1-1).
    let (beta_sq_re, x_true) = best.ok_or_else(|| {
        yee_core::Error::Numerical(format!(
            "eigensolver(mixed): no candidate shift converged to a transverse-dominated \
             propagating mode ({n_converged}/{n_candidates} candidates yielded a valid mode)"
        ))
    })?;

    Ok(pack_mixed_solution(beta_sq_re, n_t, n, &x_true))
}

/// Real-valued mixed-pencil blocks `(A, B, B_1)` extracted from the
/// (future-proof complex) [`AssembledMixed`] on the lossless path.
type RealBlocks = (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>);

/// Extract the real parts of the mixed-pencil blocks `(A, B, B_1)`,
/// rejecting a lossy (complex) input (Phase 1.3.1.2). Shared by the
/// sparse production path [`solve_dense_mixed`] and the dense reference
/// [`solve_dense_mixed_rq`].
fn mixed_real_blocks(asm: &AssembledMixed) -> Result<RealBlocks, yee_core::Error> {
    let n = asm.a.nrows();
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
    Ok((a_re, b_re, b1_re))
}

/// β-direct Rayleigh quotient `R(x) = (xᵀ K x)/(xᵀ B_1 x)` with
/// `K = k_0² B − A`. `None` if the `B_1`-norm denominator vanishes.
fn rayleigh_beta_sq(k_op: &DMatrix<f64>, b1_re: &DMatrix<f64>, n: usize, x: &[f64]) -> Option<f64> {
    let xv = DMatrix::<f64>::from_column_slice(n, 1, x);
    let num = (xv.transpose() * k_op * &xv)[(0, 0)];
    let den = (xv.transpose() * b1_re * &xv)[(0, 0)];
    (den.abs() > 0.0).then_some(num / den)
}

/// A propagating cutoff-pencil candidate: its β-direct Rayleigh quotient
/// `β² = R(x_c)` (used both as a physics-informed shift and to order the
/// candidates) and its cutoff-pencil eigenvector `x_c` (interior-DoF
/// ordering, length `n`, used to seed the β-direct shift-and-invert).
type CutoffCandidate = (f64, Vec<f64>);

/// A raw cutoff eigenpair `(k_c², x_c)` (interior-DoF ordering, length `n`)
/// before the β-direct tagging / floor in [`cutoff_candidates`].
type CutoffEigenpair = (f64, Vec<f64>);

/// The output of a raw cutoff-eigenpair source ([`dense_cutoff_eigenpairs`] /
/// [`sparse_cutoff_eigenpairs`]): the eigenpairs plus the gradient-cluster
/// floor `k_c² ≤ floor` the caller drops.
type CutoffEigenpairs = (Vec<CutoffEigenpair>, f64);

/// Gather the propagating candidates of the **cutoff** pencil
/// `A x = k_c² B x`, each recovered by inverse iteration and tagged with
/// its β-direct Rayleigh quotient `β² = R(x_c)`, sorted by `R(x_c)`
/// **descending** (highest ε_eff / slowest wave first).
///
/// **Why no transverse-energy pre-filter here (Phase 1.3.1.1 step 5.6).**
/// The step-5.2/5.3 selection rejected candidates whose *cutoff-pencil*
/// eigenvector was not transverse-energy-dominated (`‖e_t‖²/‖x‖² < `
/// [`TRANSVERSE_ENERGY_FLOOR`]). At p=2 the curl-free gradient edge
/// functions `∇(λ_aλ_b)` enlarge the near-null cluster, so the **physical**
/// dominant mode's *cutoff-pencil* eigenvector mixes heavily with the
/// E_z / gradient block (its `‖e_t‖²/‖x‖²` drops to ≈0.03 on the ε_r=10.2
/// slab) even though its **true β-direct eigenvector is fully transverse**
/// (`‖e_t‖²/‖x‖² ≈ 1.0`). Pre-rejecting on the cutoff-pencil fraction
/// therefore discards the physical mode's shift and the selection locks
/// onto a non-dominant (lower-ε_eff) mode. The transverse screen is
/// instead applied to the **converged β-direct eigenvector** in
/// [`beta_direct_shift_invert`], where it is reliable; this gather keeps
/// every real, above-floor, positive-`R` cutoff candidate as a shift.
///
/// `floor` for the curl-free gradient null-space is the same
/// `k_c² ≤ max|k_c²| · 1e-6` cutoff the prior selection used; only the
/// transverse pre-filter is removed.
fn cutoff_candidates(
    a_re: &DMatrix<f64>,
    b_re: &DMatrix<f64>,
    k_op: &DMatrix<f64>,
    b1_re: &DMatrix<f64>,
    n: usize,
    k0_sq: f64,
) -> Result<Vec<CutoffCandidate>, yee_core::Error> {
    // Source the raw `(k_c², x_c)` cutoff eigenpairs either from the dense
    // `O(n³)` `complex_eigenvalues` (small `n`, the reference path) or from
    // the sparse shift-invert (large `n`, the Phase 1.3.1.1 step 5.7
    // production path; mesh scaling past the dense cap). Both yield the same
    // *physical* low-cutoff pairs; the tagging + sort below is identical.
    //
    // Each source also returns the **gradient null-space floor**: candidates
    // with `k_c² ≤ floor` are the curl-free gradient cluster and are dropped.
    // The dense path uses its original full-spectrum `max|k_c²|·1e-6`; the
    // sparse path uses a `k_0²`-relative floor (it does not form the full
    // spectrum) — see [`sparse_cutoff_eigenpairs`]. The floor is load-bearing
    // here because a homogeneous-guide gradient mode is *purely transverse*
    // (`E_z ≡ 0`, t-frac = 1), so the downstream transverse screen cannot
    // distinguish it — only its near-zero cutoff does.
    let (raw, spurious_floor): (Vec<(f64, Vec<f64>)>, f64) = if n <= DENSE_CUTOFF_DOF_THRESHOLD {
        dense_cutoff_eigenpairs(a_re, b_re)?
    } else {
        sparse_cutoff_eigenpairs(a_re, b_re, b1_re, n, k0_sq)?
    };

    let mut cands: Vec<CutoffCandidate> = Vec::new();
    for (k_c_sq, x) in raw {
        if !k_c_sq.is_finite() || k_c_sq <= spurious_floor {
            continue;
        }
        let total: f64 = x.iter().map(|&v| v * v).sum();
        if total <= 0.0 {
            continue;
        }
        // β-direct Rayleigh quotient on the cutoff-pencil eigenvector — the
        // physics-informed shift estimate. Only propagating candidates
        // (β² > 0) target a physical guided mode.
        let Some(rq) = rayleigh_beta_sq(k_op, b1_re, n, &x) else {
            continue;
        };
        if rq <= 0.0 || !rq.is_finite() {
            continue;
        }
        cands.push((rq, x));
    }
    // Sort by the β-direct Rayleigh quotient descending: the dominant
    // quasi-TEM mode is the slowest wave (largest β² / highest ε_eff), so
    // its shift is tried first.
    cands.sort_by(|p, q| q.0.partial_cmp(&p.0).unwrap_or(std::cmp::Ordering::Equal));
    Ok(cands)
}

/// Interior-DoF count at or below which [`cutoff_candidates`] uses the dense
/// `O(n³)` `complex_eigenvalues` cutoff eigendecomposition (the step-5.2
/// reference / fallback). Above it, the sparse shift-invert
/// ([`sparse_cutoff_eigenpairs`]) is used instead (Phase 1.3.1.1 step 5.7).
///
/// The threshold is set just above the validation-anchor sizes that DoD-1
/// pins bit-identical (the WR-90 6×6 / 8×8 homogeneous, vertical-slab, and
/// uniform-fill meshes land at `n ≤ 225` interior DoF) so those gates keep
/// running through the dense reference, while the finer ε_r=10.2 meshes
/// (DoD-2) and any realistic cross-section (DoD-3) take the sparse path.
const DENSE_CUTOFF_DOF_THRESHOLD: usize = 260;

/// Dense reference: every real cutoff eigenpair `(k_c², x_c)` of
/// `A x = k_c² B x`, recovered by `B⁻¹A` `complex_eigenvalues` + inverse
/// iteration, plus the **gradient-cluster floor** `max|k_c²|·1e-6` (the
/// original step-5.2 pencil filter, preserved verbatim so this path stays
/// bit-identical). `O(n³)`; the small-`n` reference / fallback path of
/// [`cutoff_candidates`]. Returns *all* real eigenpairs (gradient cluster
/// included); the caller floors out the non-physical ones.
fn dense_cutoff_eigenpairs(
    a_re: &DMatrix<f64>,
    b_re: &DMatrix<f64>,
) -> Result<CutoffEigenpairs, yee_core::Error> {
    let b_lu = b_re.clone().lu();
    let binv_a = b_lu.solve(a_re).ok_or_else(|| {
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
    let floor = max_abs * 1e-6;
    let mut out: Vec<(f64, Vec<f64>)> = Vec::new();
    for ev in kc_evals.iter() {
        if ev.im.abs() > 1e-6 * ev.re.abs().max(1.0) {
            continue;
        }
        let k_c_sq = ev.re;
        if !k_c_sq.is_finite() {
            continue;
        }
        let Some(x) = inverse_iterate(a_re, b_re, k_c_sq) else {
            continue;
        };
        out.push((k_c_sq, x));
    }
    Ok((out, floor))
}

/// Sparse production path (Phase 1.3.1.1 step 5.7): the few **lowest
/// strictly-positive** cutoff eigenpairs `(k_c², x_c)` of `A x = k_c² B x`,
/// found by a **block shift-and-invert** (faer sparse LU of `(A − σB)`, then
/// subspace inverse iteration with a `B`-inner-product Rayleigh-Ritz), with
/// **no** dense `O(n³)` eigendecomposition. This lifts the ~457-DoF dense cap
/// that limited the whole cross-section eigensolver.
///
/// **The gradient null cluster.** `A` (curl-curl + nodal stiffness) has a
/// large curl-free gradient null space; on this pencil that cluster lands at
/// k_c² ≤ 0 (slightly negative / ≈0), while the physical guided modes sit at
/// strictly positive k_c² (the dominant TE10 of the air guide at `(π/a)²`;
/// the dielectric-loaded modes lower). A shift placed *at a small positive
/// value* does not cleanly separate the cluster (it is the nearest spectrum
/// to any small σ), so instead we sweep a **σ ladder spanning the physical
/// k_c² window `(0, k_0²·ε_r,max)`** — the upper bound a propagating mode can
/// have (`β² > 0 ⇒ k_c² < k_0²·ε_eff ≤ k_0²·ε_r,max`). Each rung's block
/// shift-invert returns the eigenpairs nearest that σ; we keep only the
/// strictly-positive-k_c² ones (gradient cluster excluded by sign) and union
/// across rungs, deduplicating by k_c². At least one rung lands near the
/// dominant mode regardless of where it sits in the positive window, so the
/// physical dominant candidate is always captured — its β-direct screen +
/// highest-β² pick happen downstream, unchanged from the dense path.
///
/// `ε_r,max` is recovered from the diagonal ratio of the ε_r-weighted mass
/// `B` to the unweighted `B_1` over the transverse block (no geometry input).
fn sparse_cutoff_eigenpairs(
    a_re: &DMatrix<f64>,
    b_re: &DMatrix<f64>,
    b1_re: &DMatrix<f64>,
    n: usize,
    k0_sq: f64,
) -> Result<CutoffEigenpairs, yee_core::Error> {
    // ε_r,max from the diagonal mass ratio B_ii / B1_ii (B = ε_r·mass on the
    // diagonal blocks, B_1 = unweighted): the max over DoF is ε_r,max. Floor
    // at 1.0 (vacuum) so the window is never degenerate.
    let mut eps_max = 1.0_f64;
    for i in 0..n {
        let d1 = b1_re[(i, i)];
        if d1.abs() > 0.0 {
            let ratio = b_re[(i, i)] / d1;
            if ratio.is_finite() && ratio > eps_max {
                eps_max = ratio;
            }
        }
    }
    // Upper edge of the physical k_c² window. A small margin (1.05) keeps the
    // top rung from sitting exactly on the bound.
    let kc_window = (k0_sq * eps_max * 1.05).max(k0_sq);

    // σ ladder across the positive window. Geometric-ish spacing favours the
    // low end (where the dielectric-loaded dominant modes live) while still
    // reaching the air-cutoff fraction (~0.43·k_0²) and above. Block width
    // grabs a handful of eigenpairs nearest each σ; positive-k_c² survivors
    // are unioned. The block + ladder are intentionally modest — the
    // downstream β-direct shift-invert does the precise mode lock-on.
    const LADDER_FRACS: &[f64] = &[0.02, 0.06, 0.15, 0.3, 0.5, 0.75, 1.0];
    let block = (12usize).min(n);

    let mut pairs: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut last_err: Option<yee_core::Error> = None;
    for &frac in LADDER_FRACS {
        let sigma = (frac * kc_window).max(1.0);
        match block_shift_invert_cutoff(a_re, b_re, sigma, block, n) {
            Ok(found) => {
                for (k_c_sq, x) in found {
                    if k_c_sq.is_finite() && k_c_sq > 0.0 {
                        pairs.push((k_c_sq, x));
                    }
                }
            }
            // A rung whose `(A − σB)` is singular (σ on an eigenvalue) or that
            // fails to converge is skipped; other rungs cover the window.
            Err(e) => last_err = Some(e),
        }
    }

    if pairs.is_empty() {
        return Err(last_err.unwrap_or_else(|| {
            yee_core::Error::Numerical(
                "eigensolver(mixed): sparse cutoff shift-invert found no positive-k_c² mode \
                 across the σ ladder (gradient cluster only?)"
                    .into(),
            )
        }));
    }

    // Deduplicate by k_c² (the same physical mode is found by adjacent rungs):
    // sort ascending, drop near-duplicates within a tight relative tolerance.
    pairs.sort_by(|p, q| p.0.partial_cmp(&q.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut deduped: Vec<(f64, Vec<f64>)> = Vec::with_capacity(pairs.len());
    for (k_c_sq, x) in pairs {
        let dup = deduped
            .last()
            .is_some_and(|(prev, _)| (k_c_sq - prev).abs() <= 1e-6 * prev.abs().max(1.0));
        if !dup {
            deduped.push((k_c_sq, x));
        }
    }
    // Gradient null-space floor for the caller. The sparse path does not form
    // the full spectrum, so it cannot use the dense `max|k_c²|·1e-6`; instead
    // it floors relative to `k_0²`. The curl-free gradient cluster sits at
    // k_c² ≈ 0 (the homogeneous-guide leak measured at k_c² ≈ 2 ≈ 5e-5·k_0²),
    // while the physical dominant cutoff is `O((π/a)²) ≈ 0.43·k_0²` (lower for
    // dielectric loading, but still ≳ 0.07·k_0² on the validated slabs). A
    // floor of `k_0²·1e-3` therefore sits safely between the cluster and the
    // lowest physical mode.
    let floor = k0_sq * SPARSE_GRADIENT_FLOOR_FRAC;
    Ok((deduped, floor))
}

/// Sparse-path gradient-cluster floor as a fraction of `k_0²`: candidates
/// with `k_c² ≤ k_0²·this` are dropped as curl-free gradient noise. Set
/// well above the measured homogeneous-guide gradient leak (≈5e-5·k_0²) and
/// well below the lowest physical dominant cutoff (≳0.07·k_0² on the
/// validation slabs, `(π/a)² ≈ 0.43·k_0²` for the air guide). See
/// [`sparse_cutoff_eigenpairs`].
const SPARSE_GRADIENT_FLOOR_FRAC: f64 = 1e-3;

/// Block shift-and-invert for the **few cutoff eigenpairs nearest `sigma`** of
/// `A x = k_c² B x`. Factors `(A − σ B)` once via faer `sp_lu` (the yee-fem
/// `build_shifted` + step-5.3 sparse-LU pattern, re-implemented here — no
/// cross-crate coupling), carries a `block`-wide subspace through
/// inverse-iteration `S ← (A − σB)⁻¹ B S`, and runs a dense `B`-inner-product
/// Rayleigh-Ritz each sweep to extract the `block` Ritz pairs nearest `σ`
/// (largest shift-invert eigenvalue `θ = 1/(k_c² − σ)`). Returns
/// `(k_c², eigenvector)` for the converged Ritz pairs (interior-DoF ordering,
/// length `n`).
///
/// This is the block analogue of the single-vector [`inverse_iterate`]; the
/// block resolves the clustered low cutoff modes (and any near-degeneracy)
/// that a one-vector iteration would smear, mirroring why yee-fem's
/// [`LobpcgEigen`] carries a block. `B` is SPD-on-the-physical-subspace here
/// (it is indefinite globally because of the edge↔node coupling, but the
/// Rayleigh-Ritz uses the symmetric-definite reduction only on the small
/// projected pencil, which is well-conditioned for the recovered block).
fn block_shift_invert_cutoff(
    a_re: &DMatrix<f64>,
    b_re: &DMatrix<f64>,
    sigma: f64,
    block: usize,
    n: usize,
) -> Result<Vec<(f64, Vec<f64>)>, yee_core::Error> {
    // (A − σB) as a sparse column matrix from its nonzeros (entry-summing
    // triplets; exact-zero entries dropped — the step-5.3 / yee-fem pattern).
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::with_capacity(n * 8);
    for j in 0..n {
        for i in 0..n {
            let v = a_re[(i, j)] - sigma * b_re[(i, j)];
            if v != 0.0 {
                triplets.push(Triplet::new(i, j, v));
            }
        }
    }
    let shifted =
        SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).map_err(|e| {
            yee_core::Error::Numerical(format!(
                "eigensolver(mixed): failed to build sparse (A − σB) for cutoff shift-invert: {e:?}"
            ))
        })?;
    let lu: Lu<usize, f64> = shifted.sp_lu().map_err(|e| {
        yee_core::Error::Numerical(format!(
            "eigensolver(mixed): sparse LU of (A − σB) failed: {e:?} (cutoff shift σ={sigma})"
        ))
    })?;

    // Apply T = (A − σB)⁻¹ B to a single column.
    let apply_t = |col: &DMatrix<f64>| -> DMatrix<f64> {
        let rhs = b_re * col;
        let mut y = faer::Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = rhs[(i, 0)];
        }
        lu.solve_in_place_with_conj(faer::Conj::No, y.as_mut());
        DMatrix::<f64>::from_fn(n, 1, |i, _| y[(i, 0)])
    };

    // Deterministic block seed (distinct dominant spatial frequency per
    // column + a column-dependent spike, mirroring yee-fem `block_seed`), so
    // the subspace spans a genuinely `block`-dimensional space and the
    // eigensolve is bit-reproducible.
    let mut s: Vec<DMatrix<f64>> = (0..block)
        .map(|c| {
            let freq = (c as f64 + 1.0) * std::f64::consts::PI;
            let phase = 0.37 * c as f64;
            let mut v = DMatrix::<f64>::from_fn(n, 1, |i, _| {
                let t = (i as f64 + 0.5) / (n as f64);
                (freq * t + phase).cos() + 0.25 * (2.0 * freq * t).sin()
            });
            v[(c % n, 0)] += 1.0;
            v
        })
        .collect();

    let mut prev_kc: Vec<f64> = vec![f64::NAN; block];
    let mut ritz: Vec<(f64, DMatrix<f64>)> = Vec::new();
    let max_sweeps = 100usize;
    for _sweep in 0..max_sweeps {
        // Inverse-iterate every column: S ← T S.
        for col in s.iter_mut() {
            *col = apply_t(col);
        }
        // B-orthonormalize the block by modified Gram-Schmidt in the B-inner
        // product, dropping columns that collapse (soft-locking). `B` is
        // SPD on the recovered physical subspace; a non-positive B-norm
        // signals a column that fell into the indefinite coupling/gradient
        // direction — drop it.
        b_orthonormalize_block(&mut s, b_re, n);
        if s.is_empty() {
            break;
        }
        // Dense Rayleigh-Ritz on the B-orthonormal block: solve the small
        // projected standard pencil (SᵀAS) c = k_c² (SᵀBS) c. With S
        // B-orthonormal, SᵀBS ≈ I, but we form it for the Cholesky reduction
        // so rounding does not bias the Ritz values.
        let bb = s.len();
        let sa: Vec<DMatrix<f64>> = s.iter().map(|c| a_re * c).collect();
        let sb: Vec<DMatrix<f64>> = s.iter().map(|c| b_re * c).collect();
        let mut g_a = DMatrix::<f64>::zeros(bb, bb);
        let mut g_b = DMatrix::<f64>::zeros(bb, bb);
        for i in 0..bb {
            for j in i..bb {
                let va = (s[i].transpose() * &sa[j])[(0, 0)];
                let vb = (s[i].transpose() * &sb[j])[(0, 0)];
                g_a[(i, j)] = va;
                g_a[(j, i)] = va;
                g_b[(i, j)] = vb;
                g_b[(j, i)] = vb;
            }
        }
        let Some((vals, vecs)) = dense_small_gen_sym_eigen(&g_a, &g_b) else {
            break;
        };
        // Rotate the block into the Ritz basis: new S column k = Σ_r S[r]·C[r,k].
        let mut new_s: Vec<DMatrix<f64>> = Vec::with_capacity(bb);
        for k in 0..bb {
            let mut v = DMatrix::<f64>::zeros(n, 1);
            for (r, sr) in s.iter().enumerate() {
                let c = vecs[(r, k)];
                if c != 0.0 {
                    v += c * sr;
                }
            }
            new_s.push(v);
        }
        ritz = vals.iter().cloned().zip(new_s.iter().cloned()).collect();
        s = new_s;

        // Convergence: the (up to `block`) Ritz k_c² nearest σ have stopped
        // moving. Compare the sorted-by-|k_c²−σ| leading values.
        let mut idx: Vec<usize> = (0..bb).collect();
        idx.sort_by(|&p, &q| {
            (vals[p] - sigma)
                .abs()
                .partial_cmp(&(vals[q] - sigma).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut moved = false;
        for (rank, &i) in idx.iter().enumerate().take(prev_kc.len().min(bb)) {
            let pv = prev_kc[rank];
            if !pv.is_finite() || (vals[i] - pv).abs() > 1e-9 * vals[i].abs().max(1.0) {
                moved = true;
            }
            prev_kc[rank] = vals[i];
        }
        if !moved {
            break;
        }
    }

    if ritz.is_empty() {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver(mixed): cutoff block shift-invert collapsed to an empty subspace \
             at σ={sigma}"
        )));
    }
    Ok(ritz
        .into_iter()
        .map(|(kc, v)| (kc, v.column(0).iter().copied().collect::<Vec<f64>>()))
        .collect())
}

/// In-place block `B`-orthonormalization by modified Gram-Schmidt in the
/// `B`-inner product (`⟨u,v⟩_B = uᵀ B v`), dropping columns whose
/// post-orthogonalization `B`-norm is non-positive or below `√ε` (soft
/// locking — and, here, also the guard that rejects columns that drifted into
/// the indefinite coupling/gradient directions where `uᵀ B u ≤ 0`). Peer of
/// yee-fem's `block_m_orthonormalize`, on `DMatrix` columns.
fn b_orthonormalize_block(block: &mut Vec<DMatrix<f64>>, b_re: &DMatrix<f64>, n: usize) {
    let drop_tol = f64::EPSILON.sqrt();
    let mut accepted: Vec<DMatrix<f64>> = Vec::with_capacity(block.len());
    for col in block.drain(..) {
        let mut v = col;
        for ej in &accepted {
            let bv = b_re * &v;
            let coeff = (ej.transpose() * &bv)[(0, 0)];
            v -= coeff * ej;
        }
        let bv = b_re * &v;
        let norm_sq = (v.transpose() * &bv)[(0, 0)];
        if !norm_sq.is_finite() || norm_sq <= drop_tol * drop_tol {
            continue; // collapsed into the accepted span, or non-B-positive
        }
        let inv = 1.0 / norm_sq.sqrt();
        accepted.push(DMatrix::<f64>::from_fn(n, 1, |i, _| v[(i, 0)] * inv));
    }
    *block = accepted;
}

/// Solve the small dense **generalized symmetric** pencil `G_a c = λ G_b c`
/// (`G_a = SᵀAS`, `G_b = SᵀBS`) by the Cholesky-reduction path (peer of
/// yee-fem's `dense_gen_sym_eigen`): `G_b = L Lᵀ`, standard problem
/// `(L⁻¹ G_a L⁻ᵀ) y = λ y`, `nalgebra` symmetric eigensolve, back-transform
/// `c = L⁻ᵀ y`. Returns `(eigenvalues, eigenvectors)` (eigenvectors
/// column-stacked in the `S` basis), or `None` if `G_b` is not SPD (the block
/// is rank-deficient despite the soft-lock drop — the caller stops the sweep).
fn dense_small_gen_sym_eigen(
    g_a: &DMatrix<f64>,
    g_b: &DMatrix<f64>,
) -> Option<(Vec<f64>, DMatrix<f64>)> {
    let bb = g_a.nrows();
    let chol = g_b.clone().cholesky()?;
    let l = chol.l();
    let y = l.solve_lower_triangular(g_a)?;
    let z = l.solve_lower_triangular(&y.transpose())?;
    let a_tilde = z.transpose();
    let a_sym = 0.5 * (&a_tilde + a_tilde.transpose());
    let eig = SymmetricEigen::new(a_sym);
    let lt = l.transpose();
    let mut c = DMatrix::<f64>::zeros(bb, bb);
    let mut vals = vec![0.0_f64; bb];
    for col in 0..bb {
        vals[col] = eig.eigenvalues[col];
        let yvec = eig.eigenvectors.column(col).clone_owned();
        let cy = lt.solve_upper_triangular(&yvec)?;
        for row in 0..bb {
            c[(row, col)] = cy[row];
        }
    }
    Some((vals, c))
}

/// Recover `k_0²` from the β-direct operator `k_op = k_0² B − A` given `A` and
/// `B`, via the least-squares ratio `Σ k_op_ij·B_ij / Σ B_ij²` of
/// `(k_op + A) = k_0² B`. Used only to thread a consistent σ scale into
/// [`cutoff_candidates`] from the dense-reference caller; the result is
/// exact up to floating-point since `k_op + A = k_0² B` holds entrywise.
fn recover_k0_sq(k_op: &DMatrix<f64>, a_re: &DMatrix<f64>, b_re: &DMatrix<f64>) -> f64 {
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    let n = b_re.nrows();
    for j in 0..n {
        for i in 0..n {
            let kb = k_op[(i, j)] + a_re[(i, j)]; // = k_0² B_ij
            let b = b_re[(i, j)];
            num += kb * b;
            den += b * b;
        }
    }
    // `den = Σ B_ij²` is zero only if B is identically zero (a degenerate /
    // empty pencil that no physical cross-section produces at nonzero
    // frequency). The 0.0 fallback degrades gracefully — it collapses the
    // σ-ladder, the sparse solve finds no candidate, and the caller surfaces
    // the (counted) "no candidate converged" error rather than a silent
    // wrong β. The debug-assert documents the invariant (step-5.7 review P1-2).
    debug_assert!(
        den > 0.0,
        "recover_k0_sq: B is identically zero (degenerate pencil)"
    );
    if den > 0.0 { num / den } else { 0.0 }
}

/// Select the dominant mode of the **cutoff** pencil whose *cutoff-pencil*
/// eigenvector is transverse-energy-dominated and has the largest β-direct
/// Rayleigh quotient `β² = R(x_c)`. Returns the eigenvector and its `R(x_c)`.
///
/// This is the **dense reference** selection ([`solve_dense_mixed_rq`]),
/// which intentionally evaluates β² as the Rayleigh quotient *on the
/// cutoff-pencil eigenvector* (the step-5.2 hybrid). It is exercised only
/// on the small-`n` homogeneous / uniform-fill anchors, where the
/// cutoff-pencil and β-direct eigenvectors coincide and the dominant mode
/// is unambiguously transverse-dominated; the `TRANSVERSE_ENERGY_FLOOR`
/// pre-filter is therefore reliable here (unlike the p=2 inhomogeneous
/// production path, which uses [`cutoff_candidates`] + the converged-
/// eigenvector screen instead — see that function's docstring).
fn select_dominant_cutoff_rq(
    a_re: &DMatrix<f64>,
    b_re: &DMatrix<f64>,
    k_op: &DMatrix<f64>,
    b1_re: &DMatrix<f64>,
    n: usize,
    n_t: usize,
) -> Result<(Vec<f64>, f64), yee_core::Error> {
    let mut best: Option<(f64, Vec<f64>)> = None; // (R(x_c), eigenvector)
    // The dense-reference path is exercised on small-`n` anchors only; pass
    // `k0_sq` recovered from the supplied `k_op = k0² B − A` so the shared
    // [`cutoff_candidates`] sparse-or-dense dispatch sees a consistent scale.
    // (`select_dominant_cutoff_rq` is itself only reached on small `n`, so the
    // dense branch is taken regardless; `k0_sq` is needed only for the sparse
    // branch's σ ladder.)
    let k0_sq = recover_k0_sq(k_op, a_re, b_re);
    for (rq, x) in cutoff_candidates(a_re, b_re, k_op, b1_re, n, k0_sq)? {
        let total: f64 = x.iter().map(|&v| v * v).sum();
        let trans: f64 = (0..n_t).map(|i| x[i] * x[i]).sum();
        if total <= 0.0 || trans / total < TRANSVERSE_ENERGY_FLOOR {
            continue;
        }
        let take = match &best {
            None => true,
            Some((curr, _)) => rq > *curr,
        };
        if take {
            best = Some((rq, x));
        }
    }
    let (rq, x_sel) = best.ok_or_else(|| {
        yee_core::Error::Numerical(
            "eigensolver(mixed): no transverse-energy-dominated propagating cutoff mode above the \
             spurious floor"
                .into(),
        )
    })?;
    Ok((x_sel, rq))
}

/// Direct β-direct **sparse shift-and-invert** for the physical mode
/// nearest the shift `sigma0` of `K x = β² B_1 x` (`K = k_0² B − A`).
///
/// Builds `(K − σ₀ B_1)` as a faer [`SparseColMat`] (entry-summing
/// triplets, near-zeros dropped — the yee-fem `build_shifted` pattern),
/// factors once via `sp_lu`, and runs inverse iteration
/// `z ← (K − σ₀B_1)⁻¹ B_1 z` to the eigenpair nearest `σ₀`. Convergence is
/// declared when the β-direct Rayleigh quotient `R(z)` stops moving.
/// Returns `(β², x_true)` with `x_true` the converged true β-direct
/// eigenvector (interior-DoF ordering, length `n`).
///
/// `seed` optionally biases the inverse-iteration start (the cutoff-mode
/// eigenvector is a good seed — it already overlaps the physical mode);
/// `None` falls back to a deterministic non-symmetric ramp.
///
/// The converged mode is **screened**: if it is not transverse-energy-
/// dominated (`‖e_t‖²/‖x‖² <` [`TRANSVERSE_ENERGY_FLOOR`]) the solve
/// surfaces a [`yee_core::Error::Numerical`] (a spurious capture — should
/// not happen at a physics-informed shift, but guarded).
fn beta_direct_shift_invert(
    k_op: &DMatrix<f64>,
    b1_re: &DMatrix<f64>,
    n: usize,
    n_t: usize,
    sigma0: f64,
    seed: Option<&[f64]>,
) -> Result<(f64, Vec<f64>), yee_core::Error> {
    // Build (K − σ₀ B_1) as a sparse column matrix from its nonzeros.
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::with_capacity(n * 8);
    for j in 0..n {
        for i in 0..n {
            let v = k_op[(i, j)] - sigma0 * b1_re[(i, j)];
            if v != 0.0 {
                triplets.push(Triplet::new(i, j, v));
            }
        }
    }
    let shifted =
        SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets).map_err(|e| {
            yee_core::Error::Numerical(format!(
                "eigensolver(mixed): failed to build sparse (K − σ₀B_1): {e:?}"
            ))
        })?;
    let lu: Lu<usize, f64> = shifted.sp_lu().map_err(|e| {
        yee_core::Error::Numerical(format!(
            "eigensolver(mixed): sparse LU of (K − σ₀B_1) failed: {e:?} \
             (shift σ₀={sigma0} may sit on an eigenvalue — re-shift)"
        ))
    })?;

    // Inverse iteration z ← (K − σ₀B_1)⁻¹ (B_1 z), normalized each step.
    // Seed with the cutoff-mode eigenvector when available (it overlaps
    // the physical mode strongly), else a deterministic non-symmetric ramp
    // (avoids accidental orthogonality to the target eigenvector).
    let mut z: Vec<f64> = match seed {
        Some(s) if s.len() == n => s.to_vec(),
        _ => (0..n).map(|i| 1.0 + (i as f64) * 0.001).collect(),
    };
    {
        let nrm = z.iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm > 0.0 {
            for v in z.iter_mut() {
                *v /= nrm;
            }
        }
    }

    let mut beta_sq = sigma0;
    let mut beta_sq_prev = f64::NAN;
    let mut converged = false;
    for _ in 0..200 {
        // rhs = B_1 z
        let zm = DMatrix::<f64>::from_column_slice(n, 1, &z);
        let rhs_m = b1_re * &zm;
        // y = (K − σ₀B_1)⁻¹ rhs  via the sparse LU.
        let mut y = faer::Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = rhs_m[(i, 0)];
        }
        lu.solve_in_place_with_conj(faer::Conj::No, y.as_mut());
        let nrm = (0..n).map(|i| y[(i, 0)] * y[(i, 0)]).sum::<f64>().sqrt();
        if !nrm.is_finite() || nrm == 0.0 {
            return Err(yee_core::Error::Numerical(
                "eigensolver(mixed): β-direct inverse iteration collapsed to a zero iterate".into(),
            ));
        }
        for i in 0..n {
            z[i] = y[(i, 0)] / nrm;
        }
        beta_sq = rayleigh_beta_sq(k_op, b1_re, n, &z).ok_or_else(|| {
            yee_core::Error::Numerical(
                "eigensolver(mixed): β-direct Rayleigh quotient has a zero B_1-norm denominator"
                    .into(),
            )
        })?;
        // Converged when the Rayleigh quotient stops moving (the iterate
        // has locked onto the eigenvector nearest σ₀).
        if beta_sq_prev.is_finite()
            && (beta_sq - beta_sq_prev).abs() <= 1e-12 * beta_sq.abs().max(1.0)
        {
            converged = true;
            break;
        }
        beta_sq_prev = beta_sq;
    }
    if !converged {
        return Err(yee_core::Error::Numerical(
            "eigensolver(mixed): β-direct inverse iteration did not converge within 200 \
             iterations (shift σ₀ may be poorly aligned with the physical mode)"
                .into(),
        ));
    }

    // Screen the converged mode: it must be transverse-energy-dominated.
    let total: f64 = z.iter().map(|v| v * v).sum();
    let trans: f64 = (0..n_t).map(|i| z[i] * z[i]).sum();
    if total <= 0.0 || trans / total < TRANSVERSE_ENERGY_FLOOR {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver(mixed): β-direct shift-and-invert converged to a non-transverse \
             (spurious E_z / gradient) mode at σ₀={sigma0} (‖e_t‖²/‖x‖²={:.3}); re-shift",
            trans / total.max(f64::MIN_POSITIVE)
        )));
    }

    Ok((beta_sq, z))
}

/// Pack a converged interior-DoF eigenvector `x` (length `n = n_t + n_z`)
/// and its `β²` into a [`MixedEigenSolution`], splitting `E_t` (first
/// `n_t`) from `E_z` (the rest) and fixing the global sign deterministically
/// off the transverse block (largest-magnitude `E_t` component positive,
/// matching [`solve_dense`]). Downstream consumers renormalize against a
/// physical reference point.
fn pack_mixed_solution(beta_sq_re: f64, n_t: usize, n: usize, x: &[f64]) -> MixedEigenSolution {
    let argmax = (0..n_t)
        .max_by(|&a, &b| {
            x[a].abs()
                .partial_cmp(&x[b].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(0);
    let sign = if x[argmax] < 0.0 { -1.0 } else { 1.0 };
    let e_t: Vec<Complex64> = (0..n_t).map(|i| Complex64::new(sign * x[i], 0.0)).collect();
    let e_z: Vec<Complex64> = (n_t..n).map(|i| Complex64::new(sign * x[i], 0.0)).collect();
    MixedEigenSolution {
        beta_sq: Complex64::new(beta_sq_re, 0.0),
        e_t,
        e_z,
    }
}

/// Dense reference for the mixed β-direct solve: cutoff-pencil select +
/// β-direct **Rayleigh quotient on the cutoff-pencil eigenvector** (the
/// step-5.2 hybrid). Retained as the small-`n` reference / fallback and
/// exercised by the unit tests; the production path
/// ([`solve_dense_mixed`]) is the step-5.3 sparse direct solve, which
/// recovers the *true* β-direct eigenvector and so removes this path's
/// mesh-stable eigenvector-mismatch bias on inhomogeneous fills. On the
/// uniform-fill anchor and the homogeneous canary the two paths agree (the
/// cutoff and β-direct eigenvectors coincide there).
#[allow(dead_code)]
fn solve_dense_mixed_rq(
    asm: &AssembledMixed,
    freq_hz: f64,
) -> Result<MixedEigenSolution, yee_core::Error> {
    let n_t = asm.n_t;
    let n = n_t + asm.n_z;
    if n == 0 || n_t == 0 {
        return Err(yee_core::Error::Numerical(
            "eigensolver(mixed): empty DoF set (no interior edges?)".into(),
        ));
    }
    let (a_re, b_re, b1_re) = mixed_real_blocks(asm)?;
    let omega = std::f64::consts::TAU * freq_hz;
    let k0 = omega / yee_core::units::C0;
    let k_op = k0 * k0 * &b_re - &a_re;
    // Phase 1.3.1.1 step 5.6: select the transverse-dominated cutoff
    // candidate with the highest β-direct Rayleigh quotient (= the dominant
    // quasi-TEM mode); its RQ is β² directly — the dense reference's β² is
    // the RQ on the cutoff-pencil eigenvector (the documented
    // step-5.2-hybrid mismatch the production sparse path removes).
    let (x_sel, beta_sq_re) = select_dominant_cutoff_rq(&a_re, &b_re, &k_op, &b1_re, n, n_t)?;
    if beta_sq_re <= 0.0 {
        return Err(yee_core::Error::Numerical(format!(
            "eigensolver(mixed): dominant mode is evanescent at {freq_hz} Hz (β² = {beta_sq_re})"
        )));
    }
    Ok(pack_mixed_solution(beta_sq_re, n_t, n, &x_sel))
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

    /// Vertical-slab WR-90 mesh: lower-x half tagged 1, rest tagged 0
    /// (mirrors `eigensolver_inhomogeneous::vertical_slab_mesh`).
    fn vertical_slab_mesh(a: f64, b: f64, nx: usize, ny: usize) -> TriMesh2D {
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
                let xc = a * ((i as f64) + 0.5) / (nx as f64);
                let tag = if xc < a / 2.0 { 1u32 } else { 0u32 };
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
        // Step-5-review P1-1 load-bearing guard (the coupling-block coverage
        // for the production sparse β-direct path, step 5.3): on a
        // horizontal-slab guide, zeroing ONLY the off-diagonal coupling
        // block B_tz (= B_ztᵀ) must measurably change the dominant
        // eigenpair — proving the coupling participates and pinning that it
        // is non-trivially placed. (On a homogeneous or vertical-slab guide
        // this delta is zero because the dominant mode has E_z = 0; the
        // horizontal slab is exactly where the coupling bites.) In the
        // β-direct pencil the coupling enters both K (via B) and the −β² B_1
        // RHS metric, so the delta is large (≈49 %) — this is the
        // load-bearing signal the integration test
        // `coupling_block_loadbearing_horizontal_slab` defers to (the
        // recovered E_z fraction on the *true* β-direct eigenvector is small,
        // ≈2e-5, so it cannot itself anchor the guard at step 5.3).
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
            "coupling delta (horizontal slab, β-direct): β with coupling {beta_full:.4}, \
             β without {beta_nc:.4}, rel Δ {rel:.3e}"
        );
        // In the β-direct pencil the coupling is strongly load-bearing
        // (≈49 % shift). Require a large delta — well above floating-point
        // noise and above the ~1 % the step-5.2 hybrid showed — so the
        // guard fails loudly if the coupling block is ever inert/misplaced.
        assert!(
            rel > 0.1,
            "zeroing B_tz must substantially change the β-direct β (rel Δ {rel:.3e}); \
             a small delta means the coupling block is inert/misplaced"
        );
    }

    #[test]
    fn sparse_direct_lands_physical_mode_on_loaded_slab() {
        // Phase 1.3.1.1 step 5.3 (DoD-1, unit level): the production
        // `solve_dense_mixed` (now the sparse β-direct shift-and-invert)
        // must land the PHYSICAL (transverse-energy-dominated) mode on the
        // high-contrast horizontal slab — not the spurious E_z / curl-free
        // gradient cluster at β² ≈ k_0²⟨ε_r⟩. We assert (a) a strongly
        // transverse-dominated eigenvector, (b) a physically-sensible
        // ε_eff well above the area-average (field-concentrated in the
        // dielectric), and (c) a positive, finite β.
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

        let sol = solve_dense_mixed(&asm, freq_hz).unwrap();
        let beta_sq = sol.beta_sq.re;
        assert!(beta_sq > 0.0 && beta_sq.is_finite(), "β² = {beta_sq}");

        // Transverse-energy fraction of the recovered eigenvector.
        let et2: f64 = sol.e_t.iter().map(|z| z.norm_sqr()).sum();
        let ez2: f64 = sol.e_z.iter().map(|z| z.norm_sqr()).sum();
        let tfrac = et2 / (et2 + ez2);
        eprintln!(
            "sparse-direct loaded slab: β={:.3}, t-frac={tfrac:.4}",
            beta_sq.sqrt()
        );
        assert!(
            tfrac >= TRANSVERSE_ENERGY_FLOOR,
            "sparse direct must land a transverse-dominated mode (t-frac {tfrac:.4}); \
             a low fraction means it captured the spurious E_z / gradient branch"
        );

        // ε_eff = (β² + (π/a)²)/k_0² must be field-concentrated (≫ air, and
        // above the area-average ≈ 5.6) for the dominant LSM-to-y mode of a
        // half-ε_r=10.2-filled guide. The spurious E_z≈0 branch would sit
        // near ε_eff ≈ ⟨ε_r⟩ at the top of the β-direct spectrum.
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = omega / yee_core::units::C0;
        let kx = std::f64::consts::PI / a;
        let eps_eff = (beta_sq + kx * kx) / (k0 * k0);
        assert!(
            eps_eff > 4.0,
            "dominant mode ε_eff {eps_eff:.3} should be field-concentrated in the dielectric"
        );
    }

    #[test]
    fn sparse_direct_matches_dense_rq_on_uniform_fill() {
        // Phase 1.3.1.1 step 5.3: on a UNIFORMLY-filled guide the cutoff-
        // pencil eigenvector and the true β-direct eigenvector coincide
        // (B_1 = ε_r⁻¹ B up to the ε-independent coupling; the dominant
        // mode is purely transverse), so the sparse direct solve and the
        // dense Rayleigh-quotient fallback must agree to tight precision.
        // This pins that the new sparse path did not perturb the analytic-
        // anchor case (DoD-4 no-regression at the unit level).
        let a = 22.86e-3;
        let b = 10.16e-3;
        let freq_hz = 10e9;
        let mesh = rectangular_mesh(a, b, 6, 6);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(2.55, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        let table = EdgeTable::build(&mesh);
        let asm = assemble_mixed(&mesh, &eps, &mu, &table);

        let beta_sparse = solve_dense_mixed(&asm, freq_hz).unwrap().beta_sq.re.sqrt();
        let beta_rq = solve_dense_mixed_rq(&asm, freq_hz)
            .unwrap()
            .beta_sq
            .re
            .sqrt();
        let rel = (beta_sparse - beta_rq).abs() / beta_rq;
        eprintln!(
            "uniform fill ε_r=2.55: sparse-direct β={beta_sparse:.5}, dense-RQ β={beta_rq:.5}, \
             rel {rel:.3e}"
        );
        assert!(
            rel < 1e-3,
            "sparse direct β {beta_sparse} must match dense-RQ β {beta_rq} on a uniform fill \
             (rel {rel:.3e}); the eigenvectors coincide there"
        );
    }

    /// Run the production downstream selection (β-direct shift-invert from
    /// each candidate shift + transverse screen + highest-β²) given an
    /// already-gathered candidate list, returning the dominant β. Mirrors the
    /// step-2 loop of [`solve_dense_mixed`]; used by the M2 dense-vs-sparse
    /// agreement guard to isolate the **candidate source** from the rest of
    /// the pipeline.
    fn dominant_beta_from_candidates(
        cands: &[CutoffCandidate],
        k_op: &DMatrix<f64>,
        b1_re: &DMatrix<f64>,
        n: usize,
        n_t: usize,
    ) -> f64 {
        let mut best: Option<f64> = None;
        for (sigma0, x_c) in cands {
            if let Ok((beta_sq, _)) =
                beta_direct_shift_invert(k_op, b1_re, n, n_t, *sigma0, Some(x_c))
                && beta_sq > 0.0
                && beta_sq.is_finite()
            {
                best = Some(match best {
                    None => beta_sq,
                    Some(c) => c.max(beta_sq),
                });
            }
        }
        best.expect("a transverse-dominated propagating mode")
            .sqrt()
    }

    /// Tag a raw `(k_c², x_c)` list with the β-direct Rayleigh quotient and
    /// the positive-k_c² floor, then sort descending — the same post-
    /// processing [`cutoff_candidates`] applies, factored for the M2 test so
    /// it can compare the dense and sparse *sources* through identical
    /// tagging.
    fn tag_candidates(
        raw: Vec<(f64, Vec<f64>)>,
        spurious_floor: f64,
        k_op: &DMatrix<f64>,
        b1_re: &DMatrix<f64>,
        n: usize,
    ) -> Vec<CutoffCandidate> {
        let mut cands: Vec<CutoffCandidate> = Vec::new();
        for (k_c_sq, x) in raw {
            if !k_c_sq.is_finite() || k_c_sq <= spurious_floor {
                continue;
            }
            let total: f64 = x.iter().map(|&v| v * v).sum();
            if total <= 0.0 {
                continue;
            }
            let Some(rq) = rayleigh_beta_sq(k_op, b1_re, n, &x) else {
                continue;
            };
            if rq <= 0.0 || !rq.is_finite() {
                continue;
            }
            cands.push((rq, x));
        }
        cands.sort_by(|p, q| q.0.partial_cmp(&p.0).unwrap_or(std::cmp::Ordering::Equal));
        cands
    }

    /// M2 (Phase 1.3.1.1 step 5.7, DoD-1 — THE GUARD). The **sparse** cutoff
    /// shift-invert (`sparse_cutoff_eigenpairs`) must yield the SAME dominant
    /// mode as the **dense** `complex_eigenvalues` path
    /// (`dense_cutoff_eigenpairs`) at the existing validation meshes, when
    /// fed through the identical downstream selection. This pins that
    /// switching the candidate source to sparse does not change the selected
    /// physical mode — the non-negotiable agreement the escape-hatch forbids
    /// shipping without. Covered cases: homogeneous (ε_r=1), uniform fill
    /// (ε_r=2.55), vertical slab (ε_r=2.2), FR-4 horizontal slab (ε_r=4.4),
    /// and the high-contrast horizontal slab (ε_r=10.2).
    #[test]
    fn sparse_cutoff_agrees_with_dense_dominant_beta() {
        let a = 22.86e-3;
        let b = 10.16e-3;
        let freq_hz = 10e9;
        let omega = std::f64::consts::TAU * freq_hz;
        let k0 = omega / yee_core::units::C0;
        let k0_sq = k0 * k0;

        // (label, mesh, ε_r map, agreement tolerance). The homogeneous /
        // uniform cases are essentially exact (the dominant cutoff is well
        // isolated); the loaded slabs allow a small tol for the sparse
        // subspace's benign rounding vs the dense QR. All are far tighter
        // than the gate tolerances (≤5 % FR-4, etc.).
        let air = || {
            let mut e = HashMap::new();
            e.insert(0u32, Complex64::new(1.0, 0.0));
            e
        };
        let uniform = || {
            let mut e = HashMap::new();
            e.insert(0u32, Complex64::new(2.55, 0.0));
            e
        };
        let loaded = |eps_fill: f64| {
            let mut e = HashMap::new();
            e.insert(0u32, Complex64::new(1.0, 0.0));
            e.insert(1u32, Complex64::new(eps_fill, 0.0));
            e
        };

        let cases: Vec<(&str, TriMesh2D, HashMap<u32, Complex64>, f64)> = vec![
            (
                "homogeneous 6x6 ε_r=1",
                rectangular_mesh(a, b, 6, 6),
                air(),
                1e-6,
            ),
            (
                "uniform 6x6 ε_r=2.55",
                rectangular_mesh(a, b, 6, 6),
                uniform(),
                1e-6,
            ),
            (
                "vertical-slab 8x8 ε_r=2.2",
                vertical_slab_mesh(a, b, 8, 8),
                loaded(2.2),
                5e-4,
            ),
            (
                "FR-4 horiz 8x8 ε_r=4.4",
                horizontal_slab_mesh(a, b, 8, 8),
                loaded(4.4),
                5e-4,
            ),
            (
                "hi-contrast horiz 8x8 ε_r=10.2",
                horizontal_slab_mesh(a, b, 8, 8),
                loaded(10.2),
                5e-4,
            ),
            // The n > DENSE_CUTOFF_DOF_THRESHOLD (260) dense-vs-sparse
            // agreement (step-5.7 review P1-1) lives in the separate
            // `#[ignore]`'d `sparse_cutoff_agrees_with_dense_at_n_above_threshold`
            // — the dense O(n³) reference at n>260 is ~40 s, kept out of the
            // routine fast path (the sparse production path at large n is also
            // exercised end-to-end by the step-5.7 mesh-scaling / finer-mesh
            // tests).
        ];

        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        mu.insert(1u32, Complex64::new(1.0, 0.0));

        for (label, mesh, eps, tol) in cases {
            let table = EdgeTable::build(&mesh);
            let asm = assemble_mixed(&mesh, &eps, &mu, &table);
            let (a_re, b_re, b1_re) = mixed_real_blocks(&asm).unwrap();
            let n = asm.n_t + asm.n_z;
            let n_t = asm.n_t;
            let k_op = k0_sq * &b_re - &a_re;

            let (dense_raw, dense_floor) = dense_cutoff_eigenpairs(&a_re, &b_re).unwrap();
            let dense_cands = tag_candidates(dense_raw, dense_floor, &k_op, &b1_re, n);
            let beta_dense = dominant_beta_from_candidates(&dense_cands, &k_op, &b1_re, n, n_t);

            let (sparse_raw, sparse_floor) =
                sparse_cutoff_eigenpairs(&a_re, &b_re, &b1_re, n, k0_sq).unwrap();
            let sparse_cands = tag_candidates(sparse_raw, sparse_floor, &k_op, &b1_re, n);
            let beta_sparse = dominant_beta_from_candidates(&sparse_cands, &k_op, &b1_re, n, n_t);

            let rel = (beta_sparse - beta_dense).abs() / beta_dense.abs().max(1e-30);
            eprintln!(
                "M2 {label} (n={n}): dense β={beta_dense:.6}, sparse β={beta_sparse:.6}, \
                 rel {rel:.3e}  (#dense_cands={}, #sparse_cands={})",
                dense_cands.len(),
                sparse_cands.len()
            );
            assert!(
                rel <= tol,
                "{label}: sparse cutoff dominant β {beta_sparse} must agree with dense \
                 {beta_dense} (rel {rel:.3e} > tol {tol:.1e}) — DoD-1 is non-negotiable"
            );
        }
    }

    /// DoD-1 agreement at `n > DENSE_CUTOFF_DOF_THRESHOLD` (260) — the regime
    /// where production dispatch actually takes the sparse branch (step-5.7
    /// review P1-1). `#[ignore]`'d because the dense O(n³) reference at this
    /// size (~40 s at n≈529) is too slow for the routine fast path; run via
    /// `cargo test -p yee-mom --lib -- --include-ignored`. Verified manually:
    /// 12×12 FR-4 gives dense β = sparse β = 328.300878 (rel 5.2e-16).
    #[test]
    #[ignore = "dense O(n³) reference at n>260 (~40s); run with --include-ignored"]
    fn sparse_cutoff_agrees_with_dense_at_n_above_threshold() {
        let (a, b, freq_hz) = (22.86e-3, 10.16e-3, 10e9);
        let k0_sq = (std::f64::consts::TAU * freq_hz / yee_core::units::C0).powi(2);
        let mesh = horizontal_slab_mesh(a, b, 12, 12);
        let mut eps = HashMap::new();
        eps.insert(0u32, Complex64::new(1.0, 0.0));
        eps.insert(1u32, Complex64::new(4.4, 0.0));
        let mut mu = HashMap::new();
        mu.insert(0u32, Complex64::new(1.0, 0.0));
        mu.insert(1u32, Complex64::new(1.0, 0.0));

        let table = EdgeTable::build(&mesh);
        let asm = assemble_mixed(&mesh, &eps, &mu, &table);
        let (a_re, b_re, b1_re) = mixed_real_blocks(&asm).unwrap();
        let n = asm.n_t + asm.n_z;
        let n_t = asm.n_t;
        assert!(
            n > DENSE_CUTOFF_DOF_THRESHOLD,
            "mesh must exceed the threshold, got n={n}"
        );
        let k_op = k0_sq * &b_re - &a_re;

        let (dr, df) = dense_cutoff_eigenpairs(&a_re, &b_re).unwrap();
        let beta_dense = dominant_beta_from_candidates(
            &tag_candidates(dr, df, &k_op, &b1_re, n),
            &k_op,
            &b1_re,
            n,
            n_t,
        );
        let (sr, sf) = sparse_cutoff_eigenpairs(&a_re, &b_re, &b1_re, n, k0_sq).unwrap();
        let beta_sparse = dominant_beta_from_candidates(
            &tag_candidates(sr, sf, &k_op, &b1_re, n),
            &k_op,
            &b1_re,
            n,
            n_t,
        );
        let rel = (beta_sparse - beta_dense).abs() / beta_dense.abs().max(1e-30);
        assert!(
            rel <= 5e-3,
            "n={n}: sparse β {beta_sparse} must agree with dense {beta_dense} (rel {rel:.3e})"
        );
    }
}
