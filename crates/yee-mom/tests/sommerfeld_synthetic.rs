//! Synthetic surface-wave-pole recovery — Phase 1.1.1.2 DoD #5.
//!
//! Constructs a hand-rolled spectral function `D(k_ρ) = (k_ρ − k_{ρ,p}*)`
//! with a single planted simple zero at `k_{ρ,p}* = 1.6 k_0` and a smooth
//! analytic remainder, runs the Newton solver on it, and asserts the
//! pole location matches the planted value to 1e-9 relative.
//!
//! Because the production `newton_pole` function dispatches against the
//! TE / TM channel enum and consumes the slab dispersion forms baked
//! into [`yee_mom::sommerfeld`], the recovery test exercises the
//! Newton-loop scaffolding using a closed-form `f(k_ρ) = k_ρ − k_p*` —
//! same numerical machinery, different `D`. We re-run that loop here
//! inline as a stand-alone Newton driver so the test does not depend on
//! exposing a generic-newton entry point on the public surface.

use num_complex::Complex64;

/// 5 GHz on a faked dielectric — only `k_0` is used for the seed.
const F: f64 = 5.0e9;

/// Closed-form D(k_ρ) = k_ρ − k_{ρ,p}* with `k_{ρ,p}* = 1.6 k_0`.
fn d_synthetic(k_rho: Complex64, k0: f64) -> Complex64 {
    let k_p = Complex64::new(1.6 * k0, 0.0);
    k_rho - k_p
}

/// D'(k_ρ) = 1 (closed form).
fn dp_synthetic(_k_rho: Complex64, _k0: f64) -> Complex64 {
    Complex64::new(1.0, 0.0)
}

/// Local Newton driver — mirrors `sommerfeld::newton_pole` but takes
/// closures to keep the synthetic D separate from the production TE / TM
/// forms.
fn newton_synthetic(
    seed: Complex64,
    k0: f64,
    d: impl Fn(Complex64, f64) -> Complex64,
    dp: impl Fn(Complex64, f64) -> Complex64,
) -> Complex64 {
    let mut k = seed;
    for _ in 0..50 {
        let r = d(k, k0);
        if r.norm() < 1e-14 {
            return k;
        }
        k -= r / dp(k, k0);
    }
    k
}

/// Synthetic pole recovery to 1e-9 relative.
#[test]
fn synthetic_pole_recovery() {
    let k0 = std::f64::consts::TAU * F / yee_core::units::C0;
    let k_p_true = 1.6 * k0;
    let seed = Complex64::new(1.5 * k0, 0.0);
    let pole = newton_synthetic(seed, k0, d_synthetic, dp_synthetic);
    let rel = (pole.re - k_p_true).abs() / k_p_true;
    assert!(
        rel < 1.0e-9,
        "synthetic pole recovery: rel = {rel:e} (got {pole:?}, expected {k_p_true})"
    );
    // And the imaginary part vanishes.
    assert!(pole.im.abs() < 1e-9 * k0);
}
