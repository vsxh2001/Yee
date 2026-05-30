//! hairpin-dim-001 (Filter Phase F1.2.2): bisection-inversion round-trip.
//!
//! Mirror of `dim_001_inversion_roundtrip` for the hairpin topology. Synthesize
//! the committed Chebyshev 0.5 dB N=5 BPF (f0 = 2 GHz, FBW = 0.10, Z0 = 50 Ω) on
//! FR-4 (εr = 4.4, h = 1.6 mm), dimension it as a hairpin, then re-evaluate the
//! coupled-microstrip model on each solved gap and assert the recovered coupling
//! coefficient reproduces the target `FBW · m_{i,i+1}` to < 1 % relative — i.e.
//! the gap bisection (shared with edge-coupled, since adjacent hairpins couple
//! through the edge gap between their arms) inverted `coupling_coefficient`
//! correctly. Also asserts `gaps_m.len() == N − 1` and `arm_length_m ≈ λ_g/4`.

use yee_filter::{Approximation, FilterSpec, Response, SpecMask, dimension_hairpin, synthesize};
use yee_layout::{Substrate, coupled_microstrip, coupling_coefficient, eps_eff, microstrip_width};

/// Speed of light in vacuum, m/s (exact, SI definition).
const C: f64 = 299_792_458.0;

/// FR-4 Chebyshev 0.5 dB N=5 bandpass spec, matching the existing yee-layout
/// test substrate (εr = 4.4, h = 1.6 mm).
fn fixture() -> (FilterSpec, Substrate) {
    let spec = FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 10.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    };
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    (spec, substrate)
}

#[test]
fn hairpin_dim_001_inversion_roundtrip() {
    let (spec, substrate) = fixture();
    let proj = synthesize(&spec);
    let dims = dimension_hairpin(&proj, &substrate)
        .expect("N=5 coupled-resonator fixture should dimension as a hairpin without error");

    assert_eq!(dims.gaps_m.len(), 4, "N=5 → 4 inter-resonator gaps");
    assert_eq!(dims.target_k.len(), 4);

    // Round-trip: each solved gap must re-evaluate to its target coupling < 1 %.
    for (i, (&gap, &target)) in dims.gaps_m.iter().zip(dims.target_k.iter()).enumerate() {
        let realized = coupling_coefficient(&coupled_microstrip(
            dims.line_width_m,
            gap,
            substrate.height_m,
            substrate.eps_r,
        ));
        let rel = (realized - target).abs() / target.abs();
        assert!(
            rel < 0.01,
            "gap[{i}] = {gap:.6e} m realizes k = {realized:.6} but target_k = {target:.6} \
             (rel error {rel:.4} >= 1%)"
        );
    }

    // Arm length must be a quarter guided wavelength at f0 (the U-folded
    // half-wave is two ≈λ/4 arms): arm_length_m == c / (4·f0·√ε_eff).
    let e_eff = eps_eff(dims.line_width_m, substrate.height_m, substrate.eps_r);
    let expected_arm = C / (4.0 * spec.f0_hz * e_eff.sqrt());
    let arm_rel = (dims.arm_length_m - expected_arm).abs() / expected_arm;
    assert!(
        arm_rel < 0.02,
        "arm_length_m = {:.6e} m but λ_g/4 = {expected_arm:.6e} m (rel error {arm_rel:.4})",
        dims.arm_length_m
    );

    // Line width matches the Hammerstad-Jensen Z0 synthesis (same as the feed).
    let expected_w = microstrip_width(spec.z0_ohm, substrate.eps_r, substrate.height_m);
    assert!(
        (dims.line_width_m - expected_w).abs() / expected_w < 1e-12,
        "line_width_m should equal microstrip_width(z0)"
    );
}
