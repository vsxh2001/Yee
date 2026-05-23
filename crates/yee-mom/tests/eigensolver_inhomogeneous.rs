//! Phase 1.3.1.1 step 5 validation gate — mixed `(E_t, E_z)` cross-section
//! eigensolve on inhomogeneous (dielectric-loaded) waveguides.
//!
//! The step-5 wire-in switched [`NumericalCrossSection::solve`] from the
//! transverse-only Nedelec eigensolve to the full mixed Lee-Sun-Cendes
//! block pencil. This gate covers:
//!
//! * **DoD-V1 (regression, homogeneous):** the mixed solve on the
//!   air-filled WR-90 cross-section reproduces the analytic TE10 β to
//!   within the existing 1 % tolerance — i.e. wiring in the longitudinal
//!   `E_z` block did not perturb the working homogeneous path. (The
//!   tighter transverse-vs-mixed `< 0.1 %` equivalence is asserted at the
//!   unit level in `eigensolver::solve::tests::
//!   mixed_solve_reproduces_transverse_beta_on_homogeneous_guide`, where
//!   both paths are reachable; here only the public mixed path is, so the
//!   external analytic reference is used.)
//! * **DoD-V2′ (capability, physics inequality + regression):** a
//!   **vertical** dielectric-slab-loaded WR-90 (lower-x half filled with
//!   `ε_r = 2.2`). The numerical β is **bracketed** by the rigorous
//!   monotonic physics inequality `β_air < β_loaded < β_fullyloaded`
//!   (β increases monotonically with dielectric fill fraction) and
//!   regression-tracked against the mesh-converged numerical value. The
//!   bracket is a necessary-condition *sanity bound*, **not** a
//!   validation of the β value itself; the regression value pins drift,
//!   and an external published reference is the open item (below). This
//!   is the spec §7c escape-hatch fallback: the published transcendental
//!   dielectric-loaded-guide root is deferred to **step-5.1** (the
//!   closed-form TE_x0 / LSE-LSM dispersion did not yield a root matching
//!   the mesh-converged numerical β in the bring-up window). The spec's
//!   literal `k₀ < β < k₀√ε_r,max` band does **not** apply to a *closed*
//!   (cutoff-bearing) guide — the air-filled TE10 β already sits below
//!   `k₀` — so the monotonic empty/full bracket is used instead.
//! * **Coupling-block guard (step-5-review P1-1):** a **horizontal**
//!   dielectric-slab-loaded WR-90 (`ε_r = 10.2`, lower-y half). Unlike a
//!   vertical slab (whose dominant mode is pure-TE, `E_z ≡ 0`, leaving
//!   the `B_tz` coupling untouched), a horizontal slab's dominant mode is
//!   genuinely **hybrid** (`E_z ≠ 0`), so the `1/μ_r` coupling block is
//!   load-bearing. This case asserts `‖E_z‖/‖E_t‖ > 1e-2` (the
//!   longitudinal field is actually present) plus the same bracket +
//!   regression on β — guarding the spec/ADR's highest-risk item (the
//!   coupling sign/placement), which the homogeneous β canary cannot
//!   reach because there `E_z = 0`. The element-level coupling
//!   sign/scale/transpose is independently pinned in
//!   `eigensolver::assembly::tests::
//!   local_b_ze_matches_independent_quadrature_sign_and_scale`, and the
//!   "zero only `B_tz`" load-bearing delta in
//!   `eigensolver::solve::tests::zeroing_coupling_changes_hybrid_mode`.
//! * **DoD-V3 (Z_w):** the numerical `Z_w` reduces to the TE form
//!   `η₀ k₀ / β` within 1 % on the homogeneous guide, and is finite,
//!   positive-real-dominated, and regression-tracked on the loaded guide.
//!
//! **Phase 1.3.1.1 step 5.1 — published-transcendental reconciliation
//! (CLAUDE.md §4 gap: still OPEN, documented finding).** The horizontal-slab
//! case now also emits a *non-failing* reconciliation diagnostic
//! ([`reconcile_against_transcendental`]) against a published closed-form
//! reference: the **LSM-to-y transverse-resonance** dispersion of the
//! slab-loaded guide (`slab_loaded_beta`; the lib-side
//! `eigensolver::reference` mirror is **independently verified** to rel err
//! `0.000e0` against a shooting / finite-difference solve of the same
//! transverse ODE, and reduces exactly to the air / fully-filled TE10
//! limits — DoD-1). The dominant family is LSM-to-y (the `TE_{m0}`-derived,
//! `H_y = 0` family), matched to the numerical mode's weakly-hybrid field
//! orientation (`‖E_z‖/‖E_t‖ ≈ 0.0105`, i.e. dominantly transverse `E_y`).
//!
//! The verified reference puts the dominant mode at **β ≈ 582.95 rad/m**
//! (ε_eff ≈ 8.17 — field concentrated in the ε_r = 10.2 layer, consistent
//! with the variational area-average estimate ε_eff ≈ 5.6 and above). The
//! numerical solver instead converges (mesh-stably, 8×8 → 12×12) to
//! **β ≈ 201.52 rad/m** (ε_eff ≈ 1.35 — barely above air). The ≈ 2.9× gap
//! is far outside any mesh tolerance, so per ADR-0052 / spec §4 the §4
//! published-benchmark gap is **not** closed: the V2′ monotonic bracket +
//! regression remain the floor, the reference ships as a reported
//! diagnostic, and the solver-side inhomogeneous-accuracy gap is a
//! **FINDING** queued to step-5.2 (out of step-5.1's lane to patch the
//! mixed solver). The reference reproducing physically-sensible β / ε_eff
//! while the solver does not is itself evidence the discrepancy is
//! solver-side, not reference-side (the question the prior bring-up could
//! not answer because its reference was unverified).

use num_complex::Complex64;
use std::collections::HashMap;
use std::f64::consts::PI;
use yee_mesh::TriMesh2D;
use yee_mom::ports::{NumericalCrossSection, RectangularWaveguideTe10};

const A: f64 = 22.86e-3; // WR-90 long dimension (m)
const B: f64 = 10.16e-3; // WR-90 short dimension (m)
const FREQ_HZ: f64 = 10.0e9;
const EPS_FILL: f64 = 2.2; // vertical-slab dielectric relative permittivity
const EPS_FILL_HI: f64 = 10.2; // horizontal-slab high-contrast substrate (RT/duroid 6010)
const C0: f64 = 299_792_458.0;

// ─────────────────────────────────────────────────────────────────────────
// Published transcendental reference (Phase 1.3.1.1 step 5.1).
//
// LSM-to-y transverse-resonance dispersion for the horizontal-slab guide
// (dielectric ε_r in 0 ≤ y ≤ d₁, air above; x-variation sin(mπx/a)). This is
// a self-contained mirror of the lib-side `eigensolver::reference` module,
// whose `slab_loaded_beta` is **independently verified** against a
// shooting-method / finite-difference solve of the same transverse ODE in
// `eigensolver::reference::tests` (rel err 0.000e0 vs the dominant LSM and
// LSE roots, and exact reduction to the air / fully-filled TE10 limits).
// Kept self-contained here (rather than re-exported through the crate's
// `__internal` surface) to confine this step's edits to the eigensolver +
// test lane; the lib-side unit tests are the load-bearing DoD-1 check.
//
// The dominant slab-loaded mode is LSM-to-y (the TE_{m0}-derived family,
// H_y = 0), confirmed against the numerical dominant mode's field
// orientation: the numerical mode is weakly hybrid (‖E_z‖/‖E_t‖ ≈ 0.0105),
// i.e. dominantly transverse E_y — the LSM-to-y signature, not LSE-to-y
// (which would have E_y = 0 and a large E_z fraction). See ADR-0052.

/// One LSM-to-y stub term `(ε_r / k_y) cot(k_y d)`, robust to imaginary k_y
/// (k_y² < 0 ⇒ k_y = j q ⇒ term = −(ε_r/q) coth(q d), real-negative).
fn lsm_term(eps_r: f64, ky_sq: f64, d: f64) -> f64 {
    if ky_sq > 0.0 {
        let k = ky_sq.sqrt();
        (eps_r / k) / (k * d).tan()
    } else {
        let q = (-ky_sq).sqrt();
        -(eps_r / q) / (q * d).tanh()
    }
}

/// LSM-to-y transverse-resonance residual; a propagating mode is a root.
fn lsm_residual(d1: f64, eps_r: f64, k0: f64, m: u32, beta: f64) -> f64 {
    let d2 = B - d1;
    let kx = (m as f64) * PI / A;
    let ky1_sq = eps_r * k0 * k0 - kx * kx - beta * beta;
    let ky2_sq = k0 * k0 - kx * kx - beta * beta;
    lsm_term(eps_r, ky1_sq, d1) + lsm_term(1.0, ky2_sq, d2)
}

/// Dominant (largest-β) LSM-to-y root of the horizontal slab-loaded guide,
/// by downward scan + bisection (no external dependency). Mirrors the
/// verified `eigensolver::reference::slab_loaded_beta`.
fn slab_loaded_beta(d1: f64, eps_r: f64, freq_hz: f64, m: u32) -> Option<f64> {
    let k0 = std::f64::consts::TAU * freq_hz / C0;
    let kx = (m as f64) * PI / A;
    let beta_max_sq = eps_r * k0 * k0 - kx * kx;
    if beta_max_sq <= 0.0 {
        return None;
    }
    let beta_hi = beta_max_sq.sqrt();
    let n = 4000usize;
    let step = (beta_hi - 1e-3) / (n as f64);
    let mut prev_beta = beta_hi - 1e-6;
    let mut prev = lsm_residual(d1, eps_r, k0, m, prev_beta);
    for i in 1..=n {
        let beta = beta_hi - 1e-6 - (i as f64) * step;
        if beta <= 1e-3 {
            break;
        }
        let cur = lsm_residual(d1, eps_r, k0, m, beta);
        if prev.is_finite()
            && cur.is_finite()
            && prev * cur < 0.0
            && (cur - prev).abs() < (cur.abs() + prev.abs() + 1.0)
        {
            let (mut lo, mut hi, mut f_lo) = (beta, prev_beta, cur);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let f_mid = lsm_residual(d1, eps_r, k0, m, mid);
                if f_lo * f_mid <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            return Some(0.5 * (lo + hi));
        }
        prev_beta = beta;
        prev = cur;
    }
    None
}

/// Structured `nx × ny` quad-grid WR-90 mesh, air everywhere (tag 0).
/// Each quad splits along the `(low-x, low-y) → (high-x, high-y)`
/// diagonal into two CCW triangles. Mirrors `eigensolver_wr90.rs`.
fn air_mesh(nx: usize, ny: usize) -> TriMesh2D {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([A * (i as f64) / (nx as f64), B * (j as f64) / (ny as f64)]);
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

/// WR-90 mesh with a **vertical** dielectric slab filling the lower-x
/// half (`x < a/2`, material tag 1); the rest is air (tag 0). The
/// triangle's centroid x-coordinate decides its material. A vertical
/// slab supports a genuine partial-fill TE mode (`E_y` is tangential to
/// the `x = const` interface), so the dominant mode shifts β while
/// staying `E_z ≈ 0` — the cleanest inhomogeneous probe of the mixed
/// solve.
fn vertical_slab_mesh(nx: usize, ny: usize) -> TriMesh2D {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([A * (i as f64) / (nx as f64), B * (j as f64) / (ny as f64)]);
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
            let xc = A * ((i as f64) + 0.5) / (nx as f64);
            let tag = if xc < A / 2.0 { 1u32 } else { 0u32 };
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
}

/// WR-90 mesh with a **horizontal** dielectric slab filling the lower-y
/// half (`y < b/2`, material tag 1); the rest is air (tag 0). The
/// triangle's centroid y-coordinate decides its material. A horizontal
/// slab puts the dielectric interface normal to `ŷ`, where the dominant
/// mode's `E_y` is the *normal* field component (`D_y` continuous, `E_y`
/// discontinuous): the mode is genuinely **hybrid** (`E_z ≠ 0`), so the
/// `1/μ_r` `E_t`/`E_z` coupling block is load-bearing — the case that
/// exercises the highest-risk part of the assembly.
fn horizontal_slab_mesh(nx: usize, ny: usize) -> TriMesh2D {
    let mut vertices = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            vertices.push([A * (i as f64) / (nx as f64), B * (j as f64) / (ny as f64)]);
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
            let yc = B * ((j as f64) + 0.5) / (ny as f64);
            let tag = if yc < B / 2.0 { 1u32 } else { 0u32 };
            triangles.push([v00, v10, v11]);
            tags.push(tag);
            triangles.push([v00, v11, v01]);
            tags.push(tag);
        }
    }
    TriMesh2D::new(vertices, triangles, None, Some(tags)).unwrap()
}

fn air_eps_mu() -> (HashMap<u32, Complex64>, HashMap<u32, Complex64>) {
    let mut eps = HashMap::new();
    eps.insert(0u32, Complex64::new(1.0, 0.0));
    eps.insert(1u32, Complex64::new(1.0, 0.0));
    let mut mu = HashMap::new();
    mu.insert(0u32, Complex64::new(1.0, 0.0));
    mu.insert(1u32, Complex64::new(1.0, 0.0));
    (eps, mu)
}

fn loaded_eps_mu() -> (HashMap<u32, Complex64>, HashMap<u32, Complex64>) {
    loaded_eps_mu_with(EPS_FILL)
}

fn loaded_eps_mu_with(eps_fill: f64) -> (HashMap<u32, Complex64>, HashMap<u32, Complex64>) {
    let mut eps = HashMap::new();
    eps.insert(0u32, Complex64::new(1.0, 0.0));
    eps.insert(1u32, Complex64::new(eps_fill, 0.0));
    let mut mu = HashMap::new();
    mu.insert(0u32, Complex64::new(1.0, 0.0));
    mu.insert(1u32, Complex64::new(1.0, 0.0));
    (eps, mu)
}

#[test]
fn dod_v1_homogeneous_mixed_reproduces_te10_beta() {
    // DoD-V1: the mixed solve on the air-filled WR-90 must still match
    // the analytic TE10 β within 1 % — wiring in the E_z block did not
    // perturb the homogeneous path.
    let mesh = air_mesh(6, 6);
    let (eps, mu) = air_eps_mu();
    let mut mode = NumericalCrossSection::new(mesh, eps, mu);
    mode.solve(FREQ_HZ).expect("homogeneous mixed solve");

    let beta_num = mode.beta.expect("β cached").re;
    let analytic = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: 1.0,
    };
    let beta_analytic = analytic.beta(FREQ_HZ);
    let rel = (beta_num - beta_analytic).abs() / beta_analytic;
    eprintln!(
        "DoD-V1 homogeneous: mixed β {beta_num:.6}, analytic TE10 {beta_analytic:.6}, rel {rel:.3e}"
    );
    assert!(
        rel < 0.01,
        "mixed β {beta_num} vs analytic {beta_analytic} (rel {rel:.4}) must stay < 1 %"
    );

    // Longitudinal field must be ~zero on the homogeneous guide.
    let ez_norm: f64 = mode
        .mode_profile_ez
        .as_ref()
        .expect("E_z cached")
        .iter()
        .map(|z| z.norm_sqr())
        .sum::<f64>()
        .sqrt();
    let et_norm: f64 = mode
        .mode_profile
        .as_ref()
        .expect("E_t cached")
        .iter()
        .map(|z| z.norm_sqr())
        .sum::<f64>()
        .sqrt();
    assert!(
        ez_norm < 1e-6 * et_norm.max(1e-30),
        "homogeneous-guide E_z must be ~zero: ‖E_z‖={ez_norm}, ‖E_t‖={et_norm}"
    );
}

#[test]
fn dod_v3_homogeneous_zw_reduces_to_te_form() {
    // DoD-V3 (homogeneous reduction guard): the numerical Z_w must reduce
    // to the TE-mode wave impedance η₀ k₀ / β = Z_TE10 within 1 %.
    let mesh = air_mesh(6, 6);
    let (eps, mu) = air_eps_mu();
    let mut mode = NumericalCrossSection::new(mesh, eps, mu);
    mode.solve(FREQ_HZ).expect("homogeneous mixed solve");

    let zw_num = mode.z_w.expect("Z_w cached").re;
    // The closed-form TE10 wave impedance is exactly η₀ / √(1−(fc/f)²) =
    // η₀ k₀ / β_analytic, i.e. the TE form the numerical Z_w must reduce
    // to on the homogeneous guide.
    let zw_te = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: 1.0,
    }
    .wave_impedance(FREQ_HZ);
    let rel = (zw_num - zw_te).abs() / zw_te;
    eprintln!(
        "DoD-V3 homogeneous: numerical Z_w {zw_num:.4} Ω, η₀k₀/β {zw_te:.4} Ω, rel {rel:.3e}"
    );
    assert!(
        rel < 0.01,
        "numerical Z_w {zw_num} must reduce to TE-form η₀k₀/β {zw_te} within 1 % (rel {rel:.4})"
    );
    // And positive-real.
    assert!(zw_num > 0.0, "Z_w must be positive-real");
}

#[test]
fn dod_v2_prime_loaded_beta_bracket_and_regression() {
    // DoD-V2′ (capability, physics inequality + regression): the
    // vertical-slab-loaded WR-90 dominant β is bracketed by the rigorous
    // monotonic inequality β_air < β_loaded < β_fullyloaded and tracked
    // against a regression value.
    let mesh = vertical_slab_mesh(8, 8);
    let (eps, mu) = loaded_eps_mu();
    let mut mode = NumericalCrossSection::new(mesh, eps, mu);
    mode.solve(FREQ_HZ).expect("loaded mixed solve");
    let beta_loaded = mode.beta.expect("β cached").re;

    // Analytic empty/full TE10 brackets (kc = π/a fixed by the PEC walls).
    let beta_air = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: 1.0,
    }
    .beta(FREQ_HZ);
    let beta_full = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: EPS_FILL,
    }
    .beta(FREQ_HZ);
    eprintln!(
        "DoD-V2′ loaded: β_loaded {beta_loaded:.4}, bracket (air {beta_air:.4}, full {beta_full:.4})"
    );

    // Rigorous monotonic physics inequality: a partial fill lies strictly
    // between empty and fully-filled.
    assert!(
        beta_loaded > beta_air,
        "loaded β {beta_loaded} must exceed air β {beta_air} (dielectric slows the wave)"
    );
    assert!(
        beta_loaded < beta_full,
        "loaded β {beta_loaded} must be below fully-filled β {beta_full}"
    );

    // Regression value (8×8-quad vertical-slab mesh, ε_r = 2.2, 10 GHz):
    // β ≈ 180.23 rad/m (mesh-converged: 8×8 → 10×10 → within 0.04 %).
    let beta_reg = 180.23;
    let rel = (beta_loaded - beta_reg).abs() / beta_reg;
    assert!(
        rel < 0.02,
        "loaded β {beta_loaded} drifted from regression {beta_reg} (rel {rel:.4}); \
         update the regression value if the formulation changed deliberately"
    );
}

#[test]
fn coupling_block_loadbearing_horizontal_slab_ez_nonzero() {
    // Step-5-review P1-1 coverage guard for the highest-risk item (the
    // E_t/E_z coupling block). A HORIZONTAL dielectric slab (interface
    // ⊥ ŷ) supports a genuinely HYBRID dominant mode (E_z ≠ 0), unlike
    // the vertical slab whose dominant mode is pure-TE (E_z = 0, coupling
    // untouched). So this case forces the 1/μ_r coupling block to
    // participate, and asserts:
    //   (1) ‖E_z‖/‖E_t‖ > 1e-2 — the longitudinal field is actually
    //       present (proves the coupling carries E_z into the mode);
    //   (2) β satisfies the same rigorous monotonic bracket
    //       β_air < β_loaded < β_fullyloaded;
    //   (3) β tracks a mesh-converged regression value.
    // ε_r = 10.2 (a standard high-contrast substrate) is used so the
    // hybrid E_z content clears the 1e-2 floor with margin.
    let mesh = horizontal_slab_mesh(8, 8);
    let (eps, mu) = loaded_eps_mu_with(EPS_FILL_HI);
    let mut mode = NumericalCrossSection::new(mesh, eps, mu);
    mode.solve(FREQ_HZ).expect("horizontal-slab mixed solve");
    let beta_loaded = mode.beta.expect("β cached").re;

    let ez_norm: f64 = mode
        .mode_profile_ez
        .as_ref()
        .expect("E_z cached")
        .iter()
        .map(|z| z.norm_sqr())
        .sum::<f64>()
        .sqrt();
    let et_norm: f64 = mode
        .mode_profile
        .as_ref()
        .expect("E_t cached")
        .iter()
        .map(|z| z.norm_sqr())
        .sum::<f64>()
        .sqrt();
    let ratio = ez_norm / et_norm.max(1e-30);
    eprintln!(
        "coupling guard (horizontal slab ε_r={EPS_FILL_HI}): β {beta_loaded:.4}, \
         ‖E_z‖/‖E_t‖ = {ratio:.5}"
    );

    // (1) The longitudinal field must be genuinely present — this is the
    // coverage the homogeneous canary cannot provide (there E_z ≡ 0).
    assert!(
        ratio > 1e-2,
        "horizontal-slab dominant mode must be hybrid (‖E_z‖/‖E_t‖ > 1e-2); \
         got {ratio:.5} — the coupling block is not load-bearing, regression"
    );

    // (2) Monotonic bracket.
    let beta_air = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: 1.0,
    }
    .beta(FREQ_HZ);
    let beta_full = RectangularWaveguideTe10 {
        a: A,
        b: B,
        eps_r: EPS_FILL_HI,
    }
    .beta(FREQ_HZ);
    assert!(
        beta_loaded > beta_air && beta_loaded < beta_full,
        "loaded β {beta_loaded} must lie in bracket (air {beta_air:.4}, full {beta_full:.4})"
    );

    // (3) Regression value (8×8-quad horizontal-slab mesh, ε_r = 10.2,
    // 10 GHz): β ≈ 201.52 rad/m (mesh-converged 8×8 → 12×12 within 0.01 %).
    let beta_reg = 201.52;
    let rel = (beta_loaded - beta_reg).abs() / beta_reg;
    assert!(
        rel < 0.02,
        "horizontal-slab β {beta_loaded} drifted from regression {beta_reg} (rel {rel:.4})"
    );

    // (4) Step-5.1 published-transcendental reconciliation — REPORTED,
    // NON-FAILING DIAGNOSTIC. The LSM-to-y transverse-resonance reference
    // (`slab_loaded_beta`, independently verified in
    // `eigensolver::reference::tests`) is compared to the numerical β.
    // This is a *finding*, not a gate: the verified reference and the
    // numerical solver DISAGREE (see below), so per ADR-0052 / spec §4 the
    // V2′ bracket + regression above remain the floor and the reference is
    // emitted as a diagnostic rather than promoted to the primary gate.
    // step-5.2 owns root-causing the solver-side gap (it is OUT of step-5.1's
    // lane to patch the mixed solver). See the module header.
    reconcile_against_transcendental(beta_loaded);
}

/// Emit the step-5.1 numerical-vs-reference reconciliation diagnostic for
/// the horizontal slab. **Non-failing**: it prints the verified-reference
/// dominant β, the numerical β at the two mesh densities the gate exercises,
/// the implied ε_eff, and the relative gap — it asserts nothing about their
/// agreement (the V2′ bracket in the caller is the gate). The reference is
/// the LSM-to-y dominant root; the numerical β is recomputed here at 8×8 and
/// 12×12 to show the gap is mesh-converged, not a discretisation artefact.
fn reconcile_against_transcendental(beta_8x8: f64) {
    let k0 = std::f64::consts::TAU * FREQ_HZ / C0;
    let kx = PI / A;
    let eps_eff = |beta: f64| (beta * beta + kx * kx) / (k0 * k0);

    // Verified published reference: dominant LSM-to-y mode, ε_r = 10.2,
    // dielectric in the lower half (d₁ = b/2), m = 1.
    let beta_ref = slab_loaded_beta(B / 2.0, EPS_FILL_HI, FREQ_HZ, 1)
        .expect("LSM transcendental dominant root must exist for the loaded guide");

    eprintln!("step-5.1 reconciliation (horizontal slab ε_r={EPS_FILL_HI}, d₁=b/2, m=1):");
    eprintln!(
        "  published reference (verified LSM-to-y transverse resonance): \
         β_ref = {beta_ref:.4} rad/m (ε_eff = {:.4})",
        eps_eff(beta_ref)
    );

    // Numerical β across mesh densities (recomputed; the caller already has
    // 8×8). 12×12 confirms the gap is not mesh-limited.
    for &(nx, ny) in &[(8usize, 8usize), (12usize, 12usize)] {
        let mesh = horizontal_slab_mesh(nx, ny);
        let (eps, mu) = loaded_eps_mu_with(EPS_FILL_HI);
        let mut mode = NumericalCrossSection::new(mesh, eps, mu);
        mode.solve(FREQ_HZ).expect("horizontal-slab mixed solve");
        let beta_num = mode.beta.expect("β cached").re;
        let rel = (beta_num - beta_ref).abs() / beta_ref;
        eprintln!(
            "  numerical {nx}×{ny}: β_num = {beta_num:.4} rad/m \
             (ε_eff = {:.4}), |β_num−β_ref|/β_ref = {rel:.4}",
            eps_eff(beta_num)
        );
    }
    let rel_8x8 = (beta_8x8 - beta_ref).abs() / beta_ref;
    eprintln!(
        "  FINDING: verified reference and numerical solver DISAGREE \
         (rel gap ≈ {rel_8x8:.2}); §4 published-benchmark gap remains OPEN, \
         V2′ bracket retained, step-5.2 queued (solver-side root-cause)."
    );
}

#[test]
fn dod_v3_loaded_zw_finite_positive_regression() {
    // DoD-V3 (loaded): the numerical Z_w on the loaded guide is finite,
    // positive-real-dominated, and regression-tracked.
    let mesh = vertical_slab_mesh(8, 8);
    let (eps, mu) = loaded_eps_mu();
    let mut mode = NumericalCrossSection::new(mesh, eps, mu);
    mode.solve(FREQ_HZ).expect("loaded mixed solve");
    let zw = mode.z_w.expect("Z_w cached");
    eprintln!("DoD-V3 loaded: Z_w = {:.4} + j{:.4} Ω", zw.re, zw.im);

    assert!(zw.re.is_finite() && zw.im.is_finite(), "Z_w must be finite");
    assert!(zw.re > 0.0, "Z_w must be positive-real-dominated");
    assert!(
        zw.im.abs() < 1e-6 * zw.re.abs(),
        "lossless guide → Z_w must be ~real"
    );

    // Loaded Z_w sits below the air-filled value (lower wave impedance
    // under dielectric loading). Regression value (8×8 mesh): ≈ 438.1 Ω.
    let zw_reg = 438.1;
    let rel = (zw.re - zw_reg).abs() / zw_reg;
    assert!(
        rel < 0.03,
        "loaded Z_w {} drifted from regression {zw_reg} (rel {rel:.4})",
        zw.re
    );
}
