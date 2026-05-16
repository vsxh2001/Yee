//! Integration tests for CODATA 2018 constants exposed by `yee_core::units`.

use yee_core::units::{C0, EPS0, ETA0, MU0};

/// Relative tolerance for CODATA reference comparisons.
const REL_TOL: f64 = 1e-12;

fn approx_rel(actual: f64, expected: f64, rel_tol: f64) -> bool {
    if expected == 0.0 {
        actual.abs() <= rel_tol
    } else {
        ((actual - expected) / expected).abs() <= rel_tol
    }
}

#[test]
fn c0_matches_codata_2018() {
    // c is defined exactly.
    assert_eq!(C0, 299_792_458.0);
}

#[test]
fn eps0_matches_codata_2018() {
    assert!(
        approx_rel(EPS0, 8.854_187_812_8e-12, REL_TOL),
        "EPS0 = {EPS0} not within {REL_TOL} of CODATA 2018"
    );
}

#[test]
fn mu0_matches_codata_2018() {
    assert!(
        approx_rel(MU0, 1.256_637_062_12e-6, REL_TOL),
        "MU0 = {MU0} not within {REL_TOL} of CODATA 2018"
    );
}

#[test]
fn eta0_matches_codata_2018() {
    assert!(
        approx_rel(ETA0, 376.730_313_668, REL_TOL),
        "ETA0 = {ETA0} not within {REL_TOL} of CODATA 2018"
    );
}

#[test]
fn eta0_is_consistent_with_sqrt_mu0_over_eps0() {
    // The three CODATA 2018 constants are each published to 12 significant
    // figures, but they are independently rounded — so the cross-consistency
    // η₀ ≟ √(μ₀/ε₀) only holds to roughly 5e-12 relative, slightly looser
    // than each constant's own 1e-12 reference tolerance.
    const CONSISTENCY_TOL: f64 = 5e-12;
    let derived = (MU0 / EPS0).sqrt();
    assert!(
        approx_rel(ETA0, derived, CONSISTENCY_TOL),
        "ETA0 = {ETA0} inconsistent with sqrt(MU0/EPS0) = {derived}"
    );
}
