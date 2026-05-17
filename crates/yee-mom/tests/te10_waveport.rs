//! TE10 wave-port closed-form sanity checks.
//!
//! Phase 1.3.1.0 — `RectangularWaveguideTe10` analytic mode.
//! WR-90 (X-band) is the standard rectangular waveguide for 8.2 – 12.4 GHz
//! service; its TE10 cutoff is `c / (2 × 22.86 mm) ≈ 6.557 GHz`. This test
//! file checks the closed-form `f_c`, `β`, `Z_TE10`, and `E_y(x, y)`
//! profile against textbook values (Pozar, *Microwave Engineering* §3.3).

use yee_mom::ports::RectangularWaveguideTe10;

#[test]
fn wr90_cutoff_is_6_5_ghz() {
    // WR-90 (X-band) standard: a = 22.86 mm, b = 10.16 mm, air-filled.
    let wg = RectangularWaveguideTe10 {
        a: 22.86e-3,
        b: 10.16e-3,
        eps_r: 1.0,
    };
    let fc = wg.cutoff_hz();
    let expected = 2.99792458e8 / (2.0 * 22.86e-3); // ≈ 6.5575 GHz
    assert!(
        (fc - expected).abs() < 1.0,
        "WR-90 cutoff: got {} expected {}",
        fc,
        expected
    );
    assert!(
        (fc - 6.5575e9).abs() < 5e6,
        "WR-90 cutoff should be ~6.56 GHz, got {}",
        fc
    );
}

#[test]
fn wr90_at_10ghz_propagates() {
    let wg = RectangularWaveguideTe10 {
        a: 22.86e-3,
        b: 10.16e-3,
        eps_r: 1.0,
    };
    let f = 10.0e9;
    let beta = wg.beta(f);
    assert!(
        beta.is_finite() && beta > 0.0,
        "WR-90 at 10 GHz should propagate (β = {})",
        beta
    );
    // β should be less than free-space k (slower phase velocity perpendicular
    // to cutoff than free space, by the standard TE-mode dispersion relation).
    let k0 = 2.0 * std::f64::consts::PI * f / 2.99792458e8;
    assert!(beta < k0, "β_10 should be < k0 above cutoff");
}

#[test]
fn wr90_at_5ghz_is_below_cutoff() {
    let wg = RectangularWaveguideTe10 {
        a: 22.86e-3,
        b: 10.16e-3,
        eps_r: 1.0,
    };
    assert!(wg.beta(5.0e9).is_nan());
    assert!(wg.wave_impedance(5.0e9).is_nan());
}

#[test]
fn e_y_profile_is_zero_at_walls_and_peak_at_center() {
    let wg = RectangularWaveguideTe10 {
        a: 22.86e-3,
        b: 10.16e-3,
        eps_r: 1.0,
    };
    let center_x = 0.5 * wg.a;
    let center_y = 0.5 * wg.b;
    assert!((wg.e_y_profile(0.0, center_y) - 0.0).abs() < 1e-12);
    assert!((wg.e_y_profile(wg.a, center_y) - 0.0).abs() < 1e-12);
    assert!((wg.e_y_profile(center_x, center_y) - 1.0).abs() < 1e-12);
}

#[test]
fn wave_impedance_above_cutoff_exceeds_intrinsic() {
    let wg = RectangularWaveguideTe10 {
        a: 22.86e-3,
        b: 10.16e-3,
        eps_r: 1.0,
    };
    let eta0 = 376.730313668;
    let z = wg.wave_impedance(10.0e9);
    assert!(z > eta0, "Z_TE10 > η_0 above cutoff (got {})", z);
}
