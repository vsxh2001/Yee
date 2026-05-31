//! dim-stepped-001 (Filter Phase F1.2.3): Pozar Example 8.6 published-benchmark.
//!
//! Synthesize the stepped-impedance low-pass filter of Pozar *Microwave
//! Engineering* §8.6, Example 8.6 — maximally-flat (Butterworth) order N = 6,
//! f_c = 2.5 GHz, Z₀ = 50 Ω, Z_high = 120 Ω, Z_low = 20 Ω — and assert the six
//! section **electrical lengths** (degrees, source → load) reproduce Pozar's
//! published table within ±1.0°:
//! `[11.85°, 33.76°, 44.28°, 46.12°, 32.41°, 12.34°]`.
//!
//! The test derives `βl` from `dimension_stepped_impedance` (it does NOT
//! hardcode the computed value as the expected), asserts the low-Z-first
//! alternation (`sections[0].high_z == false` — the standard low-pass prototype
//! begins with a shunt capacitor), and that every physical length is positive
//! and finite. Non-vacuous: six distinct published values, so a constant
//! synthesizer fails. Patterned on `hairpin_dim_001` / `dim_002_sanity`.

use yee_filter::dimension_stepped_impedance;
use yee_layout::Substrate;
use yee_synth::{Approximation, prototype};

/// Pozar Example 8.6 published section electrical lengths, degrees, source→load.
const POZAR_8_6_BETAL_DEG: [f64; 6] = [11.85, 33.76, 44.28, 46.12, 32.41, 12.34];

#[test]
fn dim_stepped_001_pozar_example_8_6() {
    // Maximally-flat (Butterworth) order N = 6 low-pass prototype.
    let proto = prototype(Approximation::Butterworth, 6);
    assert_eq!(proto.order(), 6, "Pozar Ex 8.6 is an order-6 prototype");

    // A representative microstrip substrate (same shape the edge-coupled / hairpin
    // gates build). The Pozar electrical lengths are substrate-independent — they
    // depend only on the g-values and the impedance ratios — so any physical
    // substrate exercises the same βl formula.
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };

    let dims = dimension_stepped_impedance(&proto, 2.5e9, 50.0, 120.0, 20.0, &substrate)
        .expect("Pozar Ex 8.6 fixture should dimension without error");

    assert_eq!(
        dims.sections.len(),
        6,
        "order-6 prototype → six line sections"
    );

    // Low-Z-first alternation: the standard low-pass prototype begins with a
    // shunt capacitor, so section 0 is the low-Z line.
    assert!(
        !dims.sections[0].high_z,
        "section 0 must be the low-Z (shunt-capacitor) line"
    );
    for (i, sec) in dims.sections.iter().enumerate() {
        let expect_high_z = i % 2 == 1;
        assert_eq!(
            sec.high_z, expect_high_z,
            "section {i} alternation wrong: high_z = {} but expected {expect_high_z}",
            sec.high_z
        );
    }

    // Each section's electrical length (derived from the function output, NOT
    // hardcoded) must match Pozar's published degrees within ±1.0°.
    for (i, sec) in dims.sections.iter().enumerate() {
        let betal_deg = sec.electrical_length_rad.to_degrees();
        let expected = POZAR_8_6_BETAL_DEG[i];
        let err = (betal_deg - expected).abs();
        assert!(
            err <= 1.0,
            "section {i} βl = {betal_deg:.3}° but Pozar §8.6 published {expected:.2}° \
             (error {err:.3}° > 1.0°)"
        );
    }

    // Every physical length positive and finite.
    for (i, sec) in dims.sections.iter().enumerate() {
        assert!(
            sec.length_m.is_finite() && sec.length_m > 0.0,
            "section {i} length_m = {:.6e} m must be finite and > 0",
            sec.length_m
        );
        assert!(
            sec.width_m.is_finite() && sec.width_m > 0.0,
            "section {i} width_m = {:.6e} m must be finite and > 0",
            sec.width_m
        );
    }

    // Non-vacuous: the published targets are six distinct values, so a constant
    // synthesizer would fail. Confirm the realized lengths are not all equal.
    let first = dims.sections[0].electrical_length_rad;
    assert!(
        dims.sections
            .iter()
            .any(|s| (s.electrical_length_rad - first).abs() > 1e-6),
        "section electrical lengths must not all be equal (constant fails)"
    );
}
