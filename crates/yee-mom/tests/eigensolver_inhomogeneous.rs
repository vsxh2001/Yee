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

use num_complex::Complex64;
use std::collections::HashMap;
use yee_mesh::TriMesh2D;
use yee_mom::ports::{NumericalCrossSection, RectangularWaveguideTe10};

const A: f64 = 22.86e-3; // WR-90 long dimension (m)
const B: f64 = 10.16e-3; // WR-90 short dimension (m)
const FREQ_HZ: f64 = 10.0e9;
const EPS_FILL: f64 = 2.2; // vertical-slab dielectric relative permittivity
const EPS_FILL_HI: f64 = 10.2; // horizontal-slab high-contrast substrate (RT/duroid 6010)

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
