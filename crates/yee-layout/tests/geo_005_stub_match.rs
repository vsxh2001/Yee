//! Gate `stub-match-001` (FS.6.2a, ADR-0219): the single-stub matching
//! synthesis is verified two ways — the published Pozar Example 5.2
//! stub position, and the machine contract that the synthesized
//! (d, l_open) pair nulls the combined reflection to numerical
//! precision for a whole grid of passive loads.

use std::f64::consts::{PI, TAU};

use yee_layout::single_stub_match;

/// Normalized line admittance a distance `d` toward the generator from
/// a plane with reflection `gamma`, plus an open stub of length `l`.
fn total_gamma(gamma: (f64, f64), beta: f64, d: f64, l: f64) -> f64 {
    let rot = -2.0 * beta * d;
    let gd = (
        gamma.0 * rot.cos() - gamma.1 * rot.sin(),
        gamma.0 * rot.sin() + gamma.1 * rot.cos(),
    );
    let denom = (1.0 + gd.0) * (1.0 + gd.0) + gd.1 * gd.1;
    let y_line = (
        (1.0 - gd.0 * gd.0 - gd.1 * gd.1) / denom,
        -2.0 * gd.1 / denom,
    );
    let y_tot = (y_line.0, y_line.1 + (beta * l).tan());
    let dn = (1.0 + y_tot.0) * (1.0 + y_tot.0) + y_tot.1 * y_tot.1;
    let g = (
        (1.0 - y_tot.0 * y_tot.0 - y_tot.1 * y_tot.1) / dn,
        2.0 * y_tot.1 / dn,
    );
    g.0.hypot(g.1)
}

#[test]
fn pozar_example_5_2_stub_position() {
    // Z_L = 60 − j80 Ω on 50 Ω: the first shunt-stub crossing is at
    // d = 0.110 λ (Pozar, Microwave Engineering, Example 5.2).
    let zl = (60.0 / 50.0, -80.0 / 50.0);
    let dn = (zl.0 + 1.0) * (zl.0 + 1.0) + zl.1 * zl.1;
    let gamma = ((zl.0 * zl.0 + zl.1 * zl.1 - 1.0) / dn, 2.0 * zl.1 / dn);
    let beta = TAU; // per λ, so d_m is in λ units
    let m = single_stub_match(gamma, beta);
    assert!(
        (m.d_m - 0.1104).abs() < 2e-3,
        "d = {:.4} λ vs Pozar 0.110 λ",
        m.d_m
    );
    // The pair nulls the reflection regardless of which textbook branch
    // the length lands on.
    assert!(total_gamma(gamma, beta, m.d_m, m.l_open_m) < 1e-9);
}

#[test]
fn synthesized_match_nulls_reflection_for_a_load_grid() {
    let beta = TAU / 0.0672; // a realistic FR-4 λg ≈ 67.2 mm at 2.45 GHz
    for i in 0..12 {
        for j in 1..9 {
            let mag = j as f64 * 0.1; // 0.1 ..= 0.8
            let phi = i as f64 * TAU / 12.0 - PI;
            let gamma = (mag * phi.cos(), mag * phi.sin());
            let m = single_stub_match(gamma, beta);
            let resid = total_gamma(gamma, beta, m.d_m, m.l_open_m);
            assert!(
                resid < 1e-9,
                "load |Γ| = {mag:.1}, φ = {phi:.2}: residual {resid:.2e} \
                 (d = {:.4} mm, l = {:.4} mm)",
                m.d_m * 1e3,
                m.l_open_m * 1e3
            );
            assert!(m.d_m >= 0.0 && m.d_m < PI / beta);
            assert!(m.l_open_m >= 0.0 && m.l_open_m < PI / beta);
        }
    }
}
