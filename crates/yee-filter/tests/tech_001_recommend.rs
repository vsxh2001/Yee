//! `tech-001` — guided technique-recommender canonical gate (App.2.0, ADR-0136).
//!
//! Pins the deterministic decision tree of [`yee_filter::recommend_technique`]
//! to the canonical spec→technique table (spec DoD §1) so the thresholds (the
//! ≈500 MHz distributed floor and the 5 % / 20 % fractional-bandwidth bands)
//! cannot silently drift. The gate is **non-vacuous**: the table spans every
//! leaf of the tree and four distinct techniques, so a recommender that returns
//! a constant cannot pass it.

use yee_filter::{
    Approximation, FilterSpec, RealizationTechnique, Response, SpecMask, recommend_technique,
};

/// Build a spec for the given response / centre-or-cutoff frequency /
/// fractional bandwidth (the only fields the recommender reads). The mask /
/// order / approximation are filled with representative values exactly as the
/// in-crate tests and the studio `demo_spec` construct a [`FilterSpec`].
fn spec(response: Response, f0_hz: f64, fbw: f64) -> FilterSpec {
    FilterSpec {
        response,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz,
        fbw,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 10.0,
            stopband: vec![],
        },
    }
}

/// The canonical table: each `(response, f0, fbw)` maps to exactly this primary.
fn canonical_table() -> Vec<(FilterSpec, RealizationTechnique)> {
    use RealizationTechnique::*;
    vec![
        (spec(Response::Bandpass, 100e6, 0.05), LumpedLc),
        (spec(Response::Bandpass, 2.4e9, 0.05), EdgeCoupled),
        (spec(Response::Bandpass, 2.4e9, 0.25), EdgeCoupled),
        (spec(Response::Bandpass, 5e9, 0.02), Interdigital),
        (spec(Response::Lowpass, 1e9, 0.0), SteppedImpedance),
        (spec(Response::Lowpass, 50e6, 0.0), LumpedLc),
        (spec(Response::Highpass, 1e9, 0.0), LumpedLc),
    ]
}

#[test]
fn canonical_spec_to_technique_table() {
    for (s, expected) in canonical_table() {
        let rec = recommend_technique(&s);
        assert_eq!(
            rec.primary, expected,
            "({:?}, {:.3e} Hz, fbw={}) should recommend {:?}, got {:?} — rationale: {}",
            s.response, s.f0_hz, s.fbw, expected, rec.primary, rec.rationale
        );
    }
}

#[test]
fn every_recommendation_has_a_nonempty_rationale() {
    for (s, _) in canonical_table() {
        let rec = recommend_technique(&s);
        assert!(
            !rec.rationale.trim().is_empty(),
            "({:?}, {:.3e} Hz) produced an empty rationale",
            s.response,
            s.f0_hz
        );
    }
}

#[test]
fn primary_is_not_in_alternatives() {
    for (s, _) in canonical_table() {
        let rec = recommend_technique(&s);
        assert!(
            !rec.alternatives.iter().any(|(t, _)| *t == rec.primary),
            "({:?}, {:.3e} Hz): primary {:?} appears in its own alternatives",
            s.response,
            s.f0_hz,
            rec.primary
        );
    }
}

/// Non-vacuity guard: the table must exercise at least four distinct techniques,
/// so a constant recommender provably cannot pass `canonical_spec_to_technique_table`.
#[test]
fn table_is_non_vacuous() {
    let mut seen: Vec<RealizationTechnique> = Vec::new();
    for (_, expected) in canonical_table() {
        if !seen.contains(&expected) {
            seen.push(expected);
        }
    }
    assert!(
        seen.len() >= 4,
        "the canonical table must span >= 4 distinct techniques to be non-vacuous, saw {}: {:?}",
        seen.len(),
        seen
    );
}
