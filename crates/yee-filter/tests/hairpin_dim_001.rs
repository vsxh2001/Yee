//! hairpin-dim-001 (Filter Phase F1.2.2): bisection-inversion round-trip.
//!
//! Mirror of `dim_001_inversion_roundtrip` for the hairpin topology. Synthesize
//! the committed Chebyshev 0.5 dB N=5 BPF (f0 = 2 GHz, FBW = 0.10, Z0 = 50 Ω) on
//! FR-4 (εr = 4.4, h = 1.6 mm), dimension it as a hairpin, then re-evaluate the
//! coupled-microstrip model on each solved gap and assert the recovered coupling
//! coefficient reproduces the target `FBW · m_{i,i+1}` to < 1 % relative — i.e.
//! the gap bisection (shared with edge-coupled, since adjacent hairpins couple
//! through the edge gap between their arms) inverted `coupling_coefficient`
//! correctly. Also asserts `gaps_m.len() == N − 1` and the **fold-corrected**
//! arm length `arm = (λ_g/2 − fold_spacing)/2` (R.4a: the U's midline — arm +
//! bend + arm — is the half-wave; the original `λ_g/4` form left every
//! resonator electrically long by the bend path, measured as a wrecked
//! passband by the first `engine-bpf-verify-001` run).

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

    // Fold- and corner-corrected arm length (R.4a + R.6): the resonator
    // midline — arm + fold + arm — is the half-wave, and each 90° corner
    // shortens the electrical path by κ·w (κ = 0.85, calibrated against
    // the R.4 measured seed detune), so
    // arm == (λ_g/2 − fold_spacing)/2 + 0.85·w.
    let e_eff = eps_eff(dims.line_width_m, substrate.height_m, substrate.eps_r);
    let halfwave = C / (2.0 * spec.f0_hz * e_eff.sqrt());
    let expected_arm = (halfwave - dims.fold_spacing_m) / 2.0 + 0.85 * dims.line_width_m;
    let arm_rel = (dims.arm_length_m - expected_arm).abs() / expected_arm;
    assert!(
        arm_rel < 0.02,
        "arm_length_m = {:.6e} m but (λ_g/2 − fold)/2 = {expected_arm:.6e} m \
         (rel error {arm_rel:.4})",
        dims.arm_length_m
    );

    // Line width matches the Hammerstad-Jensen Z0 synthesis (same as the feed).
    let expected_w = microstrip_width(spec.z0_ohm, substrate.eps_r, substrate.height_m);
    assert!(
        (dims.line_width_m - expected_w).abs() / expected_w < 1e-12,
        "line_width_m should equal microstrip_width(z0)"
    );
}

/// R.6: thinner (higher-impedance) resonator lines make the previously
/// TapNotRealizable thick stack dimension cleanly — the tap minimum scales
/// with (Z0/Zr) and the fold consumes fewer millimetres.
#[test]
fn hairpin_dim_zr_unlocks_the_thick_stack() {
    let (spec, substrate) = fixture();
    let project = synthesize(&spec);

    let opts = yee_filter::HairpinOptions {
        resonator_z_ohm: Some(70.0),
        ..yee_filter::HairpinOptions::default()
    };
    let dims = yee_filter::dimension_hairpin_opts(&project, &substrate, &opts)
        .expect("Zr = 70 should dimension the thick stack");

    // Thinner resonator line than the Z0 feed.
    let w_z0 = microstrip_width(spec.z0_ohm, substrate.eps_r, substrate.height_m);
    assert!(dims.line_width_m < w_z0, "resonator line must be thinner");
    assert!(
        (dims.feed_width_m - w_z0).abs() / w_z0 < 1e-12,
        "feed stays a Z0 line"
    );
    // The tap sits on the arm with the feed half-width clear of the bend.
    assert!(dims.tap_offset_m + dims.feed_width_m / 2.0 <= dims.arm_length_m);
    // Defaults (Zr = None) reproduce dimension_hairpin_with_fold exactly.
    let d_default = yee_filter::dimension_hairpin_opts(
        &project,
        &substrate,
        &yee_filter::HairpinOptions::default(),
    );
    let d_fold = yee_filter::dimension_hairpin_with_fold(&project, &substrate, 2.0);
    assert_eq!(format!("{d_default:?}"), format!("{d_fold:?}"));
}
