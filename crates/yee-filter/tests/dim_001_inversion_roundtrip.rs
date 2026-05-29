//! dim-001 (Filter Phase F1.2.0): bisection-inversion round-trip.
//!
//! Synthesize the committed Chebyshev 0.5 dB N=5 BPF (f0 = 2 GHz, FBW = 0.10,
//! Z0 = 50 Ω) on FR-4 (εr = 4.4, h = 1.6 mm), dimension it, then re-evaluate the
//! coupled-microstrip model on each solved gap and assert the recovered coupling
//! coefficient reproduces the target `FBW · m_{i,i+1}` to < 1 % relative — i.e.
//! the gap bisection inverted the (validated) `coupling_coefficient` correctly.

use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, dimension_edge_coupled, synthesize,
};
use yee_layout::{Substrate, coupled_microstrip, coupling_coefficient};

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
fn dim_001_inversion_roundtrip() {
    let (spec, substrate) = fixture();
    let proj = synthesize(&spec);
    let dims = dimension_edge_coupled(&proj, &substrate)
        .expect("N=5 coupled-resonator fixture should dimension without error");

    assert_eq!(dims.gaps_m.len(), 4, "N=5 → 4 inter-resonator gaps");
    assert_eq!(dims.target_k.len(), 4);

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
}
