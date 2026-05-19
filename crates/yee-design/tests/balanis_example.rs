//! Phase 3.nl.0 R2 — integration tests for the Balanis Ch. 14 initial-estimate
//! calculator.
//!
//! Tests are kept as straight-line `#[test]` functions (no fixture
//! abstraction) so a future Balanis-Example-14.2 row can be appended without
//! re-reading boilerplate. The crate is a leaf in the workspace; no other
//! crate's tests will pull these in.
//!
//! References:
//! - Balanis, *Antenna Theory: Analysis and Design*, 4th ed.,
//!   §14.2 Example 14.1 (10 GHz, ε_r = 2.2, h = 1.588 mm; W ≈ 11.86 mm,
//!   L ≈ 9.06 mm).
//! - CLAUDE.md §4 — published-benchmark validation case requirement.

use yee_design::{
    DesignIntent, GeometryFamily, InitialEstimate, Provenance, Substrate, substrate_library,
};

/// Build a minimal [`DesignIntent`] for a rectangular patch on an explicit
/// substrate. Centralised here so each test row stays one assertion block.
fn patch_intent(f_hz: f64, eps_r: f64, h_mm: f64) -> DesignIntent {
    DesignIntent {
        family: GeometryFamily::RectangularPatch,
        target_frequency_hz: f_hz,
        substrate: Substrate::Explicit {
            eps_r,
            h_mm,
            loss_tangent: 0.0009,
        },
        gain_target_dbi: None,
        bandwidth_target_mhz: None,
        source_prompt: format!("test patch f={f_hz} eps_r={eps_r} h={h_mm}mm"),
        provenance: Provenance {
            source: "offline".to_string(),
            model: None,
            temperature: None,
            schema_version: "1".to_string(),
            substrate_library_version: substrate_library().version.clone(),
        },
    }
}

/// **Balanis Example 14.1** (`Antenna Theory: Analysis and Design`, 4th ed.,
/// §14.2): design a rectangular microstrip antenna for `f = 10 GHz` on a
/// substrate with `ε_r = 2.2` and height `h = 1.588 mm` (= 0.1588 cm).
///
/// Published outputs (Balanis Table 14.2 / worked-example body):
///   - `W ≈ 11.86 mm` (1.186 cm)
///   - `L ≈ 9.06 mm`  (0.906 cm)
///
/// Tolerance: ±0.5% per the R2 brief (CLAUDE.md §4 published-benchmark
/// requirement).
#[test]
fn balanis_example_14_1_matches_published() {
    let intent = patch_intent(10.0e9, 2.2, 1.588);
    let est = InitialEstimate::from_intent(&intent).expect("estimate must succeed");

    let w_mm = est.width_m * 1.0e3;
    let l_mm = est.length_m * 1.0e3;

    let w_err = (w_mm - 11.86).abs() / 11.86;
    let l_err = (l_mm - 9.06).abs() / 9.06;

    assert!(
        w_err < 0.005,
        "W = {w_mm:.4} mm vs published 11.86 mm (rel err {:.4}% > 0.5%)",
        w_err * 100.0
    );
    assert!(
        l_err < 0.005,
        "L = {l_mm:.4} mm vs published 9.06 mm (rel err {:.4}% > 0.5%)",
        l_err * 100.0
    );

    // Bonus assertions: y₀ and w_feed should be finite + plausible. We do
    // *not* gate on a precise value because Balanis Example 14.1 itself does
    // not specify a numeric inset offset — the inset feed is computed in a
    // follow-on example.
    assert!(est.inset_offset_m.is_finite() && est.inset_offset_m >= 0.0);
    assert!(est.feed_width_m.is_finite() && est.feed_width_m > 0.0);
    // ΔL and ε_reff should sit in textbook ranges.
    let eps_reff_expected = 1.972; // Balanis Eq. 14-1 at this substrate.
    assert!(
        (est.eps_reff - eps_reff_expected).abs() / eps_reff_expected < 0.005,
        "ε_reff = {} vs expected {eps_reff_expected}",
        est.eps_reff
    );
}

/// 2.4 GHz on FR-4 (ε_r ≈ 4.4, h = 1.6 mm). No published reference value to
/// gate on; the brief asks only that the result be finite + positive +
/// plausible-magnitude. This is the canonical Phase-3 `mom-003` target so
/// having this row in the test suite doubles as a smoke test for the
/// downstream solver gate (Phase 3.nl.0 R6, `nl-001`).
#[test]
fn fr4_2g4_patch_is_plausible_and_finite() {
    let intent = patch_intent(2.4e9, 4.4, 1.6);
    let est = InitialEstimate::from_intent(&intent).expect("estimate must succeed");

    // Every dimension is finite and strictly positive.
    for (name, v) in [
        ("W", est.width_m),
        ("L", est.length_m),
        ("y_0", est.inset_offset_m),
        ("w_feed", est.feed_width_m),
        ("eps_reff", est.eps_reff),
        ("delta_L", est.delta_l_m),
    ] {
        assert!(v.is_finite(), "{name} not finite: {v}");
        assert!(v > 0.0 || name == "y_0", "{name} not positive: {v}");
    }
    // y₀ specifically may be 0 in the degenerate R_in > R_edge case; for
    // FR-4 at 2.4 GHz with a 50 Ω feed we expect a strictly positive offset.
    assert!(est.inset_offset_m > 0.0);

    // Plausibility: at 2.4 GHz a half-wavelength in air is ~62.5 mm; on FR-4
    // (ε_reff ≈ 3.3) the guided half-wavelength is ~34 mm. The patch length
    // is the half-guided-wavelength minus 2·ΔL, so a ballpark of 25-35 mm is
    // expected. W on this substrate is typically ~25-45 mm.
    let w_mm = est.width_m * 1.0e3;
    let l_mm = est.length_m * 1.0e3;
    assert!(
        (15.0..60.0).contains(&w_mm),
        "W = {w_mm} mm not in plausible 15-60 mm range"
    );
    assert!(
        (20.0..45.0).contains(&l_mm),
        "L = {l_mm} mm not in plausible 20-45 mm range"
    );

    // Inset offset must be inside the patch.
    assert!(
        est.inset_offset_m < est.length_m,
        "y_0 = {} m ≥ L = {} m",
        est.inset_offset_m,
        est.length_m
    );

    // Feed width: 50 Ω on FR-4 / 1.6 mm is canonically ≈ 3 mm
    // (Hammerstad–Jensen). A wider 1-6 mm window keeps the assertion away
    // from formula-form quibbles.
    let wf_mm = est.feed_width_m * 1.0e3;
    assert!(
        (1.0..6.0).contains(&wf_mm),
        "w_feed = {wf_mm} mm not in plausible 1-6 mm range for 50 Ω on FR-4"
    );
}

/// Pure-function property: [`InitialEstimate::from_intent`] called twice with
/// the same input must return bit-identical output. The R2 brief makes this
/// explicit — the spec §8 byte-identity guarantee for the downstream emit
/// stage relies on the estimator being deterministic.
#[test]
fn from_intent_is_pure_idempotent() {
    let intent = patch_intent(5.8e9, 3.55, 0.508); // 5.8 GHz on RO4003C-like.
    let a = InitialEstimate::from_intent(&intent).expect("a");
    let b = InitialEstimate::from_intent(&intent).expect("b");

    // PartialEq compares every f64 field with `==`. The function does no
    // I/O, no randomness, and no global mutable state, so the bits must
    // match exactly — not approximately.
    assert_eq!(a, b, "from_intent is not deterministic");
}
