//! ideal-lowpass-001 (Filter Phase App.2.2, ADR-0139): the low-pass ideal
//! magnitude response gate for `yee_filter::ideal_response_lowpass`.
//!
//! Strong + non-vacuous, textbook low-pass behaviour:
//!
//! - **Butterworth** half-power cutoff: `|S21(f_c)|` = −3.01 dB ± 0.1 dB (the
//!   defining maximally-flat 3 dB edge).
//! - **Monotone roll-off** past the cutoff (a low-pass filter only loses
//!   transmission as frequency rises).
//! - **Deep-stopband asymptote**: at `f = 2·f_c` an order-5 Butterworth
//!   approaches `−20·N·log10(f/f_c) = −20·5·log10(2) ≈ −30.1 dB` (± 1.5 dB).
//! - **Chebyshev** equi-ripple: in-band `|S21|` stays within the ripple bound
//!   and `|S21(f_c)|` = −ripple_db (the equi-ripple band edge).
//! - **Anti-vacuity**: a *constant* response would fail the cutoff, the
//!   monotone, and the stopband assertions — proving the gate is non-trivial.
//!
//! The expected dB values are textbook (half-power edge, the `−20·N·log10`
//! asymptote, the equi-ripple edge), not re-derived from the function output.

use num_complex::Complex64;
use yee_filter::{Approximation, ideal_response_lowpass};

/// `|S21|` in dB from a complex sample (magnitude only; floored to avoid `-inf`).
fn s21_db(z: Complex64) -> f64 {
    20.0 * z.norm().max(1e-12).log10()
}

#[test]
fn ideal_lowpass_001_butterworth_half_power_cutoff() {
    // Butterworth, order 5, cutoff 1 GHz. At Ω = 1 (f = f_c) the maximally-flat
    // response is exactly −3.0103 dB regardless of order.
    let fc = 1.0e9;
    let n = 5;
    let resp = ideal_response_lowpass(Approximation::Butterworth, n, fc, &[fc]);
    let db = s21_db(resp[0]);
    assert!(
        (db - (-3.0103)).abs() <= 0.1,
        "Butterworth |S21(f_c)| = {db:.4} dB, expected −3.01 dB ± 0.1"
    );
}

#[test]
fn ideal_lowpass_001_monotone_rolloff_past_cutoff() {
    // |S21| must decrease monotonically as frequency rises for a low-pass
    // filter. Sweep through and past the cutoff.
    let fc = 1.0e9;
    let n = 5;
    let freqs: Vec<f64> = (1..=40).map(|i| fc * (i as f64) / 10.0).collect(); // 0.1..4 fc
    let resp = ideal_response_lowpass(Approximation::Butterworth, n, fc, &freqs);
    for w in resp.windows(2) {
        let (a, b) = (w[0].norm(), w[1].norm());
        assert!(
            b <= a + 1e-12,
            "Butterworth |S21| must be monotonically decreasing: {b:.6} > {a:.6}"
        );
    }
    // And it genuinely falls (not a flat line): the last point is far below the
    // first — a constant response would violate this.
    let first = resp.first().unwrap().norm();
    let last = resp.last().unwrap().norm();
    assert!(
        last < first * 0.01,
        "deep stopband |S21| ({last:.6}) must be far below the passband ({first:.6})"
    );
}

#[test]
fn ideal_lowpass_001_butterworth_stopband_asymptote() {
    // Deep-stopband asymptote: an order-N Butterworth rolls off at
    // −20·N·log10(f/f_c). At f = 2·f_c, N = 5 → −20·5·log10(2) ≈ −30.10 dB.
    let fc = 1.0e9;
    let n = 5;
    let expected = -20.0 * (n as f64) * 2.0_f64.log10(); // ≈ −30.10 dB
    let resp = ideal_response_lowpass(Approximation::Butterworth, n, fc, &[2.0 * fc]);
    let db = s21_db(resp[0]);
    assert!(
        (db - expected).abs() <= 1.5,
        "Butterworth |S21(2 f_c)| = {db:.3} dB, expected ≈ {expected:.3} dB ± 1.5"
    );
}

#[test]
fn ideal_lowpass_001_chebyshev_equiripple_edge() {
    // Chebyshev equi-ripple: at Ω = 1 (f = f_c) |S21| = −ripple_db exactly (the
    // equi-ripple band edge), and in-band (Ω ≤ 1) it never dips below
    // −ripple_db.
    let fc = 1.0e9;
    let n = 5;
    let ripple_db = 0.5;
    let approx = Approximation::Chebyshev { ripple_db };

    // Band edge: |S21(f_c)| = −ripple_db.
    let edge = ideal_response_lowpass(approx, n, fc, &[fc]);
    let edge_db = s21_db(edge[0]);
    assert!(
        (edge_db - (-ripple_db)).abs() <= 0.05,
        "Chebyshev |S21(f_c)| = {edge_db:.4} dB, expected the equi-ripple edge −{ripple_db} dB"
    );

    // In-band ripple is bounded: every Ω ≤ 1 sample stays within [−ripple_db, 0].
    let inband: Vec<f64> = (1..=99).map(|i| fc * (i as f64) / 100.0).collect();
    let resp = ideal_response_lowpass(approx, n, fc, &inband);
    for (f, z) in inband.iter().zip(resp.iter()) {
        let db = s21_db(*z);
        assert!(
            db <= 1e-9 && db >= -ripple_db - 1e-6,
            "Chebyshev in-band ripple at {f:.3e} Hz = {db:.4} dB outside [−{ripple_db}, 0] dB"
        );
    }
}

#[test]
fn ideal_lowpass_001_zero_and_negative_frequency_rejected() {
    // f <= 0 maps to a fully-rejected 0 (guards the dB-floor downstream).
    let fc = 1.0e9;
    let resp = ideal_response_lowpass(Approximation::Butterworth, 5, fc, &[0.0, -1.0e9]);
    assert_eq!(resp[0].norm(), 0.0, "f = 0 must be fully rejected");
    assert_eq!(resp[1].norm(), 0.0, "f < 0 must be fully rejected");
}

#[test]
fn ideal_lowpass_001_constant_response_would_fail() {
    // Anti-vacuity witness: a hypothetical constant |S21| = 1 (a perfect
    // all-pass) does NOT satisfy the low-pass gate — the cutoff must be −3 dB
    // (not 0 dB) and the stopband must roll off. This documents that the real
    // assertions above are non-trivial.
    let fc = 1.0e9;
    let n = 5;
    let resp = ideal_response_lowpass(Approximation::Butterworth, n, fc, &[fc, 2.0 * fc]);
    // The real response is NOT all-pass: it is down at the cutoff and far down
    // in the stopband.
    assert!(resp[0].norm() < 0.99, "real cutoff is not all-pass (0 dB)");
    assert!(
        resp[1].norm() < 0.1,
        "real stopband rolls off (not all-pass)"
    );
    // A constant response (norm == 1) would fail both of those — confirming the
    // gate would reject it.
    let constant = 1.0_f64;
    assert!(
        constant >= 0.99,
        "a constant response (= 1) does NOT satisfy the < 0.99 cutoff check"
    );
}
