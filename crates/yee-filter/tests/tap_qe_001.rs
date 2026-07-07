//! Gate `tap-qe-001` (R.4/F1.2.1, ADR-0197): the tapped half-wave-resonator
//! qe→tap inversion `t = (L/π)·acos(√((π/2)(Z0/Zr)/qe))` round-trips the
//! forward relation `Qe(t) = (π/2)(Z0/Zr)/cos²(πt/L)`, matches a
//! hand-computed value, is monotone (weaker coupling → tap toward the
//! fold), and rejects unrealizable external Qs.

use std::f64::consts::PI;

use yee_filter::{DimError, tap_offset_from_qe};

#[test]
fn tap_matches_hand_computed_value_and_round_trips() {
    // qe = 6.67 (e.g. g0·g1/FBW = 1/0.15), Z0 = Zr = 50 Ω, arm = 8 mm.
    let arm = 8.0e-3;
    let l = 2.0 * arm;
    let qe = 6.67;
    let t = tap_offset_from_qe(qe, 50.0, 50.0, l).unwrap();
    // Hand: t/L = acos(sqrt((π/2)/6.67))/π = 0.33868.
    assert!(
        (t / l - 0.33868).abs() < 1e-4,
        "t/L = {} vs hand-computed 0.33868",
        t / l
    );
    // Forward relation round-trips.
    let qe_back = (PI / 2.0) / (PI * t / l).cos().powi(2);
    assert!((qe_back - qe).abs() / qe < 1e-12, "round-trip qe {qe_back}");
}

#[test]
fn weaker_coupling_moves_the_tap_toward_the_fold() {
    let l = 16.0e-3;
    let t_strong = tap_offset_from_qe(2.0, 50.0, 50.0, l).unwrap();
    let t_weak = tap_offset_from_qe(20.0, 50.0, 50.0, l).unwrap();
    assert!(t_strong < t_weak, "{t_strong} vs {t_weak}");
    // Both on the physical half of the resonator (below the fold at L/2).
    assert!(t_strong > 0.0 && t_weak < l / 2.0);
}

#[test]
fn unrealizable_qe_is_rejected() {
    // qe below (π/2)(Z0/Zr): stronger than a tap at the antinode can couple.
    let err = tap_offset_from_qe(1.0, 50.0, 50.0, 16.0e-3).unwrap_err();
    assert!(matches!(err, DimError::TapNotRealizable { .. }), "{err}");
}
