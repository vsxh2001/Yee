//! Integration tests for `yee_core::FreqRange`.

use yee_core::{Error, FreqRange};

#[test]
fn new_accepts_valid_input() {
    let r = FreqRange::new(1.0e9, 2.0e9, 11).expect("valid range");
    assert_eq!(r.start_hz, 1.0e9);
    assert_eq!(r.stop_hz, 2.0e9);
    assert_eq!(r.n_points, 11);
}

#[test]
fn new_rejects_non_increasing_band() {
    let err = FreqRange::new(2.0e9, 1.0e9, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_equal_endpoints() {
    let err = FreqRange::new(1.0e9, 1.0e9, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_zero_points() {
    let err = FreqRange::new(1.0e9, 2.0e9, 0).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_positive_infinite_stop() {
    let err = FreqRange::new(1.0e9, f64::INFINITY, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_negative_infinite_start() {
    let err = FreqRange::new(f64::NEG_INFINITY, 2.0e9, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_nan_start() {
    let err = FreqRange::new(f64::NAN, 2.0e9, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn new_rejects_nan_stop() {
    let err = FreqRange::new(1.0e9, f64::NAN, 11).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)));
}

#[test]
fn iter_yields_single_point() {
    let r = FreqRange::new(1.5e9, 2.5e9, 1).expect("valid range");
    let pts: Vec<f64> = r.iter().collect();
    assert_eq!(pts, vec![1.5e9]);
}

#[test]
fn iter_yields_two_exact_endpoints() {
    let r = FreqRange::new(1.0e9, 2.0e9, 2).expect("valid range");
    let pts: Vec<f64> = r.iter().collect();
    assert_eq!(pts.len(), 2);
    assert_eq!(pts[0], 1.0e9);
    assert_eq!(pts[1], 2.0e9);
}

#[test]
fn iter_yields_evenly_spaced_points_with_exact_endpoints() {
    let r = FreqRange::new(1.0e9, 2.0e9, 5).expect("valid range");
    let pts: Vec<f64> = r.iter().collect();
    assert_eq!(pts.len(), 5);
    // Endpoints exact.
    assert_eq!(pts[0], 1.0e9);
    assert_eq!(pts[4], 2.0e9);
    // Step size consistency: differences should match (1e9 - 0) / 4 within fp.
    // f64 lerp on these magnitudes is well under 1 µHz of absolute error.
    let step = 0.25e9;
    for (k, pt) in pts.iter().enumerate().skip(1) {
        let expected = 1.0e9 + (k as f64) * step;
        let diff = (pt - expected).abs();
        assert!(
            diff <= 1.0e-6,
            "pts[{k}] = {pt} not within 1e-6 Hz of {expected}"
        );
    }
}

#[test]
fn iter_count_matches_n_points() {
    let r = FreqRange::new(2.4e9, 5.8e9, 101).expect("valid range");
    assert_eq!(r.iter().count(), 101);
}

#[test]
fn iter_size_hint_reports_exact_remaining() {
    let r = FreqRange::new(1.0e9, 2.0e9, 7).expect("valid range");
    let mut it = r.iter();
    assert_eq!(it.size_hint(), (7, Some(7)));
    let _ = it.next();
    assert_eq!(it.size_hint(), (6, Some(6)));
    // ExactSizeIterator agrees with size_hint.
    assert_eq!(it.len(), 6);
}
