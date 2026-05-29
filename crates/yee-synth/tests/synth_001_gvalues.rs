//! `synth-001` — published lowpass-prototype g-values gate.
//!
//! Computed Butterworth and Chebyshev element values must match the published
//! tables (Pozar Tables 8.3/8.4; Matthaei-Young-Jones Table 4.05-2) to
//! `≤ 1e-3` absolute. These numbers are ground truth — if a check misses,
//! the formula (β/γ constants or the even-order load) is wrong; do NOT loosen
//! the tolerance.

use yee_synth::{Approximation, prototype};

const TOL: f64 = 1e-3;

/// Assert `g[1..=expected.len()]` matches `expected` and `g[N+1] == g_load`.
fn assert_gvalues(approx: Approximation, n: usize, expected: &[f64], g_load: f64) {
    let proto = prototype(approx, n);
    assert_eq!(
        proto.g.len(),
        n + 2,
        "prototype g-vector length should be N+2 = {} for {approx:?} N={n}",
        n + 2
    );
    for (k, &want) in expected.iter().enumerate() {
        let got = proto.g[k + 1]; // g[0]=g0, g1 at index 1
        assert!(
            (got - want).abs() < TOL,
            "{approx:?} N={n}: g{} = {got:.4}, expected {want:.4} (Δ={:.2e})",
            k + 1,
            (got - want).abs()
        );
    }
    let got_load = proto.g[n + 1];
    assert!(
        (got_load - g_load).abs() < TOL,
        "{approx:?} N={n}: g{} (load) = {got_load:.4}, expected {g_load:.4} (Δ={:.2e})",
        n + 1,
        (got_load - g_load).abs()
    );
}

#[test]
fn butterworth_n3() {
    assert_gvalues(Approximation::Butterworth, 3, &[1.0, 2.0, 1.0], 1.0);
}

#[test]
fn butterworth_n5() {
    assert_gvalues(
        Approximation::Butterworth,
        5,
        &[0.6180, 1.6180, 2.0000, 1.6180, 0.6180],
        1.0,
    );
}

#[test]
fn chebyshev_0p5db_n3() {
    assert_gvalues(
        Approximation::Chebyshev { ripple_db: 0.5 },
        3,
        &[1.5963, 1.0967, 1.5963],
        1.0,
    );
}

#[test]
fn chebyshev_0p5db_n5() {
    assert_gvalues(
        Approximation::Chebyshev { ripple_db: 0.5 },
        5,
        &[1.7058, 1.2296, 2.5408, 1.2296, 1.7058],
        1.0,
    );
}

#[test]
fn chebyshev_3db_n3() {
    assert_gvalues(
        Approximation::Chebyshev { ripple_db: 3.0 },
        3,
        &[3.3487, 0.7117, 3.3487],
        1.0,
    );
}

/// Even-order Chebyshev load check: `g_{N+1} = coth²(β/4) ≈ 1.9841` for the
/// 0.5 dB N=4 prototype. This is the formula most likely to be wrong if the
/// even-order branch is missing.
#[test]
fn chebyshev_0p5db_n4_even_order_load() {
    let proto = prototype(Approximation::Chebyshev { ripple_db: 0.5 }, 4);
    let g5 = proto.g[5]; // g_{N+1} for N=4
    assert!(
        (g5 - 1.9841).abs() < TOL,
        "Chebyshev 0.5 dB N=4: g5 (even-order load) = {g5:.4}, expected 1.9841 (Δ={:.2e})",
        (g5 - 1.9841).abs()
    );
}
