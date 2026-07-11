//! FS.5a yield gates (ADR-0211):
//!
//! - `yield-mc-001`: the estimator reproduces the analytic normal CDF —
//!   pass iff x < z with x ~ N(0,1) gives yield Φ(z); the MC estimate
//!   must bracket Φ(z) within 1.5× its own Wilson 95 % CI at n = 100 000
//!   for z ∈ {−1, 0, 0.5, 1, 2}. (Deterministic seed ⇒ not flaky.)
//! - `yield-mc-002`: bit-identical determinism in the seed.
//! - `surrogate-yield-001` (the FULL-SUITE-ROADMAP FS.5 gate): yield
//!   through a trained GP surrogate vs brute-force MC on the closed form
//!   itself, patch-resonance testcase — same seed, so the sample streams
//!   are identical and the difference isolates surrogate error.
//!
//! All instant (< 1 s), non-ignored.

use nalgebra::{DMatrix, DVector};
use yee_surrogate::{GaussianProcess, ToleranceSpec, yield_estimate};

/// Standard normal CDF, Abramowitz–Stegun 7.1.26 (|ε| < 1.5e-7 — three
/// orders below the MC noise floor at n = 1e5).
fn phi(z: f64) -> f64 {
    let x = z / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    let erf_abs = 1.0 - poly * (-x * x).exp();
    let erf = if x < 0.0 { -erf_abs } else { erf_abs };
    0.5 * (1.0 + erf)
}

#[test]
fn yield_mc_001_matches_analytic_normal_cdf() {
    let spec = ToleranceSpec {
        nominal: vec![0.0],
        sigma: vec![1.0],
    };
    for &z in &[-1.0, 0.0, 0.5, 1.0, 2.0] {
        let est = yield_estimate(|p| p[0] < z, &spec, 100_000, 20260711);
        let err = (est.yield_frac - phi(z)).abs();
        assert!(
            err <= 1.5 * est.ci95_half_width,
            "z={z}: MC yield {} vs Phi(z) {} (err {err:.5}, ci {:.5})",
            est.yield_frac,
            phi(z),
            est.ci95_half_width
        );
    }
}

#[test]
fn yield_mc_002_deterministic_in_seed() {
    let spec = ToleranceSpec {
        nominal: vec![1.0, 2.0],
        sigma: vec![0.3, 0.1],
    };
    let pass = |p: &[f64]| p[0] + p[1] < 3.2;
    let a = yield_estimate(pass, &spec, 10_000, 99);
    let b = yield_estimate(pass, &spec, 10_000, 99);
    assert_eq!(a, b, "same seed must reproduce bit-identically");
    // Different seeds explore different sample streams.
    let others: Vec<usize> = [1u64, 2, 3]
        .iter()
        .map(|&s| yield_estimate(pass, &spec, 10_000, s).n_pass)
        .collect();
    assert!(
        others.iter().any(|&n| n != a.n_pass),
        "different seeds must not all reproduce n_pass={}",
        a.n_pass
    );
}

/// Patch resonance closed form, GHz: f = c / (2 L √ε_eff) with the
/// zeroth-order ε_eff = (ε_r + 1)/2. Nonlinear in both parameters.
fn patch_f_ghz(l_m: f64, eps_r: f64) -> f64 {
    let c = 299_792_458.0;
    c / (2.0 * l_m * ((eps_r + 1.0) / 2.0).sqrt()) / 1.0e9
}

#[test]
fn surrogate_yield_001_gp_matches_brute_force() {
    // FR-4 patch: L = 29 mm ± 0.1 mm (etch), ε_r = 4.4 ± 0.05 (batch
    // spread). Spec: resonance within ±40 MHz of nominal. Linearized
    // σ_f ≈ 18 MHz ⇒ analytic yield ≈ 2Φ(2.2)−1 ≈ 0.97 — deliberately
    // away from both 1.0 (trivial) and 0.5 (spec-meaningless).
    let (l0, er0) = (29.0e-3, 4.4);
    let (sl, ser) = (0.1e-3, 0.05);
    let f0 = patch_f_ghz(l0, er0);
    let band_ghz = 0.040;

    // Train the GP in σ-normalized coordinates (u = (p − nominal)/σ, so
    // both axes are O(1) and one RBF length scale fits both) on a 9×9
    // grid spanning ±3σ, target centred at f0.
    let grid: Vec<f64> = (0..9).map(|i| -3.0 + 0.75 * i as f64).collect();
    let n = grid.len() * grid.len();
    let mut x = DMatrix::<f64>::zeros(n, 2);
    let mut y = DVector::<f64>::zeros(n);
    for (row, (i, j)) in (0..9).flat_map(|i| (0..9).map(move |j| (i, j))).enumerate() {
        let (ul, uer) = (grid[i], grid[j]);
        x[(row, 0)] = ul;
        x[(row, 1)] = uer;
        y[row] = patch_f_ghz(l0 + sl * ul, er0 + ser * uer) - f0;
    }
    let gp = GaussianProcess::fit(x, y, 1.5, 0.1, 1e-6).expect("GP fit");

    let spec = ToleranceSpec {
        nominal: vec![l0, er0],
        sigma: vec![sl, ser],
    };
    let n_mc = 20_000;
    let seed = 20260711;

    // Brute force: the closed form itself.
    let brute = yield_estimate(
        |p| (patch_f_ghz(p[0], p[1]) - f0).abs() <= band_ghz,
        &spec,
        n_mc,
        seed,
    );
    // Surrogate: identical sample stream (same seed), GP posterior mean.
    let surr = yield_estimate(
        |p| {
            let u = DVector::from_row_slice(&[(p[0] - l0) / sl, (p[1] - er0) / ser]);
            gp.predict_mean(&u).abs() <= band_ghz
        },
        &spec,
        n_mc,
        seed,
    );

    // Sanity: the brute-force yield sits in the designed-for regime.
    assert!(
        brute.yield_frac > 0.90 && brute.yield_frac < 0.999,
        "brute-force yield {} outside the designed regime",
        brute.yield_frac
    );
    // The FS.5 roadmap gate: surrogate yield ≈ brute-force yield. Same
    // seed ⇒ the difference is pure surrogate error at the spec boundary.
    // Measured 2026-07-11: brute 0.9721, surrogate 0.9720 (Δ = 1e-4, and
    // both on the analytic linearization 2Φ(2.2)−1 ≈ 0.972); tolerance
    // 0.02 leaves room for future kernel/hyperparameter churn only.
    let delta = (surr.yield_frac - brute.yield_frac).abs();
    assert!(
        delta <= 0.02,
        "surrogate yield {} vs brute-force {} (delta {delta:.5})",
        surr.yield_frac,
        brute.yield_frac
    );
}
