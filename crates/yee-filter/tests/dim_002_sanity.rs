//! dim-002 (Filter Phase F1.2.0): physical-sanity checks on the synthesized
//! edge-coupled dimensions.
//!
//! For the same Chebyshev 0.5 dB N=5 / FR-4 fixture as dim-001: every gap > 0;
//! gaps strictly decrease as their target coupling increases (tighter coupling
//! → smaller gap, from the monotonic coupled-line model); the line width is
//! exactly `microstrip_width(z0, εr, h)`; the resonator length is within ±2 % of
//! `c/(2·f0·√ε_eff)`; and all dimensions sit in a physically sane µm–mm range.

use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, dimension_edge_coupled, synthesize,
};
use yee_layout::{Substrate, eps_eff, microstrip_width};

const C: f64 = 299_792_458.0;

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
fn dim_002_sanity() {
    let (spec, substrate) = fixture();
    let proj = synthesize(&spec);
    let dims = dimension_edge_coupled(&proj, &substrate).expect("fixture dimensions");

    // Every gap strictly positive.
    for (i, &g) in dims.gaps_m.iter().enumerate() {
        assert!(g > 0.0, "gap[{i}] = {g:.6e} m must be > 0");
    }

    // Gaps decrease as target_k increases (tighter coupling → smaller gap).
    // Sort the (target_k, gap) pairs by ascending target_k; the gaps must come
    // out monotonically non-increasing, and *strictly* decreasing across
    // distinct target_k. (A symmetric prototype — like this Chebyshev N=5 — has
    // mirror-equal couplings `[k0, k1, k1, k0]`, so equal targets must give
    // equal gaps; only the bisection inverse being well-defined makes that hold.)
    let mut pairs: Vec<(f64, f64)> = dims
        .target_k
        .iter()
        .copied()
        .zip(dims.gaps_m.iter().copied())
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut saw_distinct = false;
    for w in pairs.windows(2) {
        let (k_lo, g_lo) = w[0];
        let (k_hi, g_hi) = w[1];
        if k_hi > k_lo {
            // Distinct targets: strictly larger k ⇒ strictly smaller gap.
            saw_distinct = true;
            assert!(
                g_hi < g_lo,
                "tighter coupling must give a smaller gap: k {k_lo:.6}->{k_hi:.6} but \
                 gap {g_lo:.6e}->{g_hi:.6e} did not strictly decrease"
            );
        } else {
            // Equal targets: equal gaps (the inverse is single-valued). Allow a
            // bisection-tolerance slack (gaps agree to the convergence band).
            let rel = (g_hi - g_lo).abs() / g_lo;
            assert!(
                rel < 1e-3,
                "equal target_k ({k_lo:.6} == {k_hi:.6}) must give equal gaps but \
                 {g_lo:.6e} vs {g_hi:.6e} (rel {rel:.4})"
            );
        }
    }
    assert!(
        saw_distinct,
        "fixture must contain at least two distinct couplings to test monotonicity"
    );

    // Line width is exactly the HJ synthesis width.
    let w_expected = microstrip_width(spec.z0_ohm, substrate.eps_r, substrate.height_m);
    assert_eq!(
        dims.line_width_m, w_expected,
        "line_width_m must equal microstrip_width(z0, εr, h) exactly"
    );

    // Resonator length within ±2 % of c/(2·f0·√ε_eff).
    let e_eff = eps_eff(dims.line_width_m, substrate.height_m, substrate.eps_r);
    let len_expected = C / (2.0 * spec.f0_hz * e_eff.sqrt());
    let rel = (dims.resonator_length_m - len_expected).abs() / len_expected;
    assert!(
        rel < 0.02,
        "resonator_length_m = {:.6e} m vs λ_g/2 = {len_expected:.6e} m (rel {rel:.4} >= 2%)",
        dims.resonator_length_m
    );

    // Physically-sane ranges (ADR-0097: "a physically sane µm–mm range").
    //
    // The in-plane *coupling features* — line width and inter-resonator gaps —
    // sit in the µm–mm window the DoD names: [1 µm, 20 mm]. The resonator
    // length is a half guided wavelength, however, and a λ_g/2 resonator at
    // f0 = 2 GHz on FR-4 (ε_eff ≈ 3.3) is intrinsically ≈ 41 mm — longer than a
    // 20 mm coupling-feature cap can be. (The DoD's literal "all dimensions in
    // [1 µm, 20 mm]" is dimensionally inconsistent with a 2 GHz half-wave
    // resonator; see the F1.2.0 report finding.) The length is therefore graded
    // against its own physically-sane half-wave bound: a 1–10 GHz half-wave
    // microstrip resonator is tens of mm, so [1 mm, 200 mm] is the sane window.
    let feature_in_range = |x: f64| (1e-6..=20e-3).contains(&x);
    assert!(
        feature_in_range(dims.line_width_m),
        "line_width_m {:.6e} out of [1µm, 20mm]",
        dims.line_width_m
    );
    for (i, &g) in dims.gaps_m.iter().enumerate() {
        assert!(feature_in_range(g), "gap[{i}] {g:.6e} out of [1µm, 20mm]");
    }
    assert!(
        (1e-3..=200e-3).contains(&dims.resonator_length_m),
        "resonator_length_m {:.6e} out of the sane half-wave window [1mm, 200mm]",
        dims.resonator_length_m
    );
}
