//! dim-combline-001 (Filter Phase F1.2.5): combline dimensional synthesis gate.
//!
//! The PROPER, non-tautological gate for [`yee_filter::dimension_combline`]
//! (Hong & Lancaster *Microstrip Filters for RF/Microwave Applications* §5.2.5,
//! eqs 5.42–5.46), in three parts:
//!
//! 1. **Published benchmark (H&L eq 5.46) — the non-tautological core.**
//!    Synthesize the 5-pole 0.1 dB Chebyshev band-pass at FBW = 0.10 (the book's
//!    worked combline example, g = [1.1468, 1.3712, 1.9750, 1.3712, 1.1468]) and
//!    assert the synthesized external Q and adjacent couplings reproduce H&L's
//!    *published* design numbers — Qe ≈ 11.468, M₁₂ ≈ 0.07975, M₂₃ ≈ 0.06077 —
//!    plus the FBW = 0.15 pseudocombline point (Qe ≈ 7.645, M₁₂ ≈ 0.11962,
//!    M₂₃ ≈ 0.09115). These are specific published external numbers exercising the
//!    full g → Qe/M chain; a constant or wrong synthesis fails by ≫ tol. We do
//!    NOT hardcode the synthesizer's own output as expected — we compare it to the
//!    book.
//! 2. **Combline-distinct resonance consistency (first-principles, NOT the cap
//!    tautology).** From a `dimension_combline` result (θ0 = 45° = π/4) build the
//!    loaded short-circuited-stub input susceptance *independently*,
//!    `B(f) = −(1/Z0)·cot(θ0·f/f0) + 2π·f·C_L`, root-find `B(f) = 0` over a band
//!    around f0, and assert the root equals f0 within ±1%. This re-derives
//!    resonance from the admittance (catching a wrong cap / length / sign) — it
//!    does **not** assert `C_L == cot(θ0)/(2π·f0·Z0)` back (which would be
//!    tautological — the engine's own emit formula).
//! 3. The coupling gaps are solved/bracketed (no clamping), `θ0 ≥ π/2` → a
//!    `DimError`, and all physical dims + the loading cap are positive and finite.

use std::f64::consts::{FRAC_PI_4, PI};

use yee_filter::{
    Approximation, DimError, FilterProject, FilterSpec, Response, SpecMask, dimension_combline,
    synthesize,
};
use yee_layout::Substrate;

/// FR-4 representative substrate (matches the hairpin / stepped-Z gates).
fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// Build the 5-pole 0.1 dB Chebyshev band-pass at the requested fractional
/// bandwidth (H&L §5.2.5 worked combline example).
fn spec_5pole_cheb_01db(fbw: f64) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.1 },
        f0_hz: 2.0e9,
        fbw,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.1,
            return_loss_db: 16.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    }
}

/// External Q + the two distinct adjacent couplings (M₁₂, M₂₃) of a synthesized
/// 5-pole project: `M_{i,i+1} = FBW · m[i][i+1]` (= H&L's `M_{i,i+1}`; see
/// `yee_filter::dimension` module docs).
fn qe_m12_m23(project: &FilterProject) -> (f64, f64, f64) {
    let fbw = project.spec.fbw;
    let m12 = fbw * project.coupling.m[0][1];
    let m23 = fbw * project.coupling.m[1][2];
    (project.coupling.qe_in, m12, m23)
}

#[test]
fn dim_combline_001_hl_eq_5_46_published_benchmark() {
    // --- FBW = 0.10 (H&L eq 5.46 worked combline design) -------------------
    let proj_010 = synthesize(&spec_5pole_cheb_01db(0.10));
    let (qe, m12, m23) = qe_m12_m23(&proj_010);

    // Published H&L numbers (g = [1.1468, 1.3712, 1.9750, 1.3712, 1.1468]).
    const HL_QE_010: f64 = 11.468;
    const HL_M12_010: f64 = 0.07975;
    const HL_M23_010: f64 = 0.06077;

    let rel = |got: f64, want: f64| (got - want).abs() / want.abs();
    assert!(
        rel(qe, HL_QE_010) < 0.01,
        "FBW=0.10 Qe = {qe:.6} vs H&L 11.468 (rel {:.5} >= 1%)",
        rel(qe, HL_QE_010)
    );
    assert!(
        rel(m12, HL_M12_010) < 0.01,
        "FBW=0.10 M12 = {m12:.6} vs H&L 0.07975 (rel {:.5} >= 1%)",
        rel(m12, HL_M12_010)
    );
    assert!(
        rel(m23, HL_M23_010) < 0.01,
        "FBW=0.10 M23 = {m23:.6} vs H&L 0.06077 (rel {:.5} >= 1%)",
        rel(m23, HL_M23_010)
    );
    // yee-synth's Chebyshev g-values reproduce H&L's to ~5 decimals, so the
    // numbers also clear the much tighter ±1e-3 absolute band — pin that too so a
    // future g-value drift is caught immediately, not silently masked by the ±1%.
    assert!(
        (qe - HL_QE_010).abs() < 1e-3,
        "FBW=0.10 Qe = {qe:.6} not within 1e-3 of H&L 11.468"
    );
    assert!(
        (m12 - HL_M12_010).abs() < 1e-3 && (m23 - HL_M23_010).abs() < 1e-3,
        "FBW=0.10 M12 = {m12:.6} / M23 = {m23:.6} not within 1e-3 of H&L 0.07975 / 0.06077"
    );

    // --- FBW = 0.15 (pseudocombline point) ---------------------------------
    let proj_015 = synthesize(&spec_5pole_cheb_01db(0.15));
    let (qe2, m12_2, m23_2) = qe_m12_m23(&proj_015);

    const HL_QE_015: f64 = 7.645;
    const HL_M12_015: f64 = 0.11962;
    const HL_M23_015: f64 = 0.09115;

    assert!(
        rel(qe2, HL_QE_015) < 0.01,
        "FBW=0.15 Qe = {qe2:.6} vs H&L 7.645 (rel {:.5} >= 1%)",
        rel(qe2, HL_QE_015)
    );
    assert!(
        rel(m12_2, HL_M12_015) < 0.01,
        "FBW=0.15 M12 = {m12_2:.6} vs H&L 0.11962 (rel {:.5} >= 1%)",
        rel(m12_2, HL_M12_015)
    );
    assert!(
        rel(m23_2, HL_M23_015) < 0.01,
        "FBW=0.15 M23 = {m23_2:.6} vs H&L 0.09115 (rel {:.5} >= 1%)",
        rel(m23_2, HL_M23_015)
    );

    // Surface the synthesized vs published numbers for the verification log.
    println!(
        "H&L eq 5.46  FBW=0.10: Qe={qe:.6} (H&L 11.468)  M12={m12:.6} (0.07975)  \
         M23={m23:.6} (0.06077)"
    );
    println!(
        "H&L eq 5.46  FBW=0.15: Qe={qe2:.6} (H&L 7.645)  M12={m12_2:.6} (0.11962)  \
         M23={m23_2:.6} (0.09115)"
    );
}

#[test]
fn dim_combline_001_resonance_consistency() {
    // Demo band-pass spec at FBW = 0.10, θ0 = 45° = π/4.
    let spec = spec_5pole_cheb_01db(0.10);
    let z0 = spec.z0_ohm;
    let f0 = spec.f0_hz;
    let proj = synthesize(&spec);
    let dims = dimension_combline(&proj, FRAC_PI_4, &substrate())
        .expect("N=5 coupled-resonator combline fixture should dimension without error");

    let theta0 = dims.theta0_rad;
    let c_l = dims.loading_cap_f;

    // INDEPENDENT loaded short-circuited-stub susceptance:
    //   B(f) = −(1/Z0)·cot(θ0·f/f0) + 2π·f·C_L    (β(f)·L = θ0·(f/f0)).
    // This does NOT invert the engine's C_L formula; it re-derives resonance from
    // the admittance, so a wrong cap / length / dispersion / sign is caught.
    let b_of = |f: f64| {
        let arg = theta0 * f / f0;
        let cot = arg.cos() / arg.sin();
        -(1.0 / z0) * cot + 2.0 * PI * f * c_l
    };

    // Root-find B(f) = 0 by bisection over [0.5·f0, 1.5·f0]. The short-circuited
    // stub susceptance −(1/Z0)·cot rises monotonically with f over (0, π/2) and
    // the cap term is monotone increasing, so B is monotone increasing across the
    // bracket → a single sign change at the resonance.
    let mut lo = 0.5 * f0;
    let mut hi = 1.5 * f0;
    let (b_lo, b_hi) = (b_of(lo), b_of(hi));
    assert!(
        b_lo < 0.0 && b_hi > 0.0,
        "B(f) must change sign across [0.5 f0, 1.5 f0]: B(lo) = {b_lo:.6}, B(hi) = {b_hi:.6}"
    );
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if b_of(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let f_res = 0.5 * (lo + hi);
    let rel = (f_res - f0).abs() / f0;
    println!(
        "combline resonance: independent B(f)=0 root f_res = {f_res:.6e} Hz vs f0 = {f0:.6e} Hz \
         (rel {rel:.6}); theta0 = {theta0:.6} rad, C_L = {c_l:.6e} F"
    );
    assert!(
        rel < 0.01,
        "loaded-stub resonance root f_res = {f_res:.6e} Hz is not within 1% of f0 = {f0:.6e} Hz \
         (rel {rel:.6}) — cap/length/sign bug"
    );
}

#[test]
fn dim_combline_001_dims_solved_and_positive() {
    let spec = spec_5pole_cheb_01db(0.10);
    let proj = synthesize(&spec);
    let sub = substrate();
    let dims = dimension_combline(&proj, FRAC_PI_4, &sub)
        .expect("combline fixture should dimension without error");

    // N = 5 → 4 inter-resonator gaps, index-aligned with target_k.
    assert_eq!(dims.gaps_m.len(), 4, "N=5 → 4 inter-resonator gaps");
    assert_eq!(dims.target_k.len(), 4);

    // Every solved gap re-evaluates to its target coupling (< 1% relative) — the
    // shared solve_gap bisection bracketed the target (no clamping).
    for (i, (&gap, &target)) in dims.gaps_m.iter().zip(dims.target_k.iter()).enumerate() {
        let realized = yee_layout::coupling_coefficient(&yee_layout::coupled_microstrip(
            dims.line_width_m,
            gap,
            sub.height_m,
            sub.eps_r,
        ));
        let rel = (realized - target).abs() / target.abs();
        assert!(
            rel < 0.01,
            "gap[{i}] = {gap:.6e} m realizes k = {realized:.6} but target_k = {target:.6} \
             (rel {rel:.4} >= 1%)"
        );
        assert!(
            gap.is_finite() && gap > 0.0,
            "gap[{i}] = {gap:.6e} m must be finite and > 0"
        );
    }

    // All physical dims + the loading cap positive and finite.
    assert!(
        dims.line_width_m.is_finite() && dims.line_width_m > 0.0,
        "line_width_m = {:.6e} must be finite and > 0",
        dims.line_width_m
    );
    assert!(
        dims.resonator_length_m.is_finite() && dims.resonator_length_m > 0.0,
        "resonator_length_m = {:.6e} must be finite and > 0",
        dims.resonator_length_m
    );
    assert!(
        dims.loading_cap_f.is_finite() && dims.loading_cap_f > 0.0,
        "loading_cap_f = {:.6e} must be finite and > 0",
        dims.loading_cap_f
    );
    assert!(
        (dims.theta0_rad - FRAC_PI_4).abs() < 1e-12,
        "theta0_rad should round-trip the requested π/4"
    );

    // Sanity: at θ0 = π/4, L = (π/4)/β(f0) = λ_g/8 — a quarter of the edge-coupled
    // λ_g/2 straight resonator (compact), confirming the short-line geometry.
    let e_eff = yee_layout::eps_eff(dims.line_width_m, sub.height_m, sub.eps_r);
    let beta0 = 2.0 * PI * spec.f0_hz * e_eff.sqrt() / 299_792_458.0;
    let expected_len = FRAC_PI_4 / beta0;
    assert!(
        (dims.resonator_length_m - expected_len).abs() / expected_len < 1e-9,
        "resonator_length_m = {:.6e} m should equal θ0/β(f0) = {expected_len:.6e} m",
        dims.resonator_length_m
    );
}

#[test]
fn dim_combline_001_theta0_out_of_range_errors() {
    let proj = synthesize(&spec_5pole_cheb_01db(0.10));
    let sub = substrate();

    // θ0 = π/2 (exactly the upper bound) → InvalidTheta0 (cot = 0 → C_L = 0).
    match dimension_combline(&proj, std::f64::consts::FRAC_PI_2, &sub) {
        Err(DimError::InvalidTheta0(t)) => {
            assert!((t - std::f64::consts::FRAC_PI_2).abs() < 1e-12)
        }
        other => panic!("θ0 = π/2 must give InvalidTheta0, got {other:?}"),
    }
    // θ0 > π/2 → InvalidTheta0 (cot < 0 → non-physical negative cap).
    assert!(matches!(
        dimension_combline(&proj, 2.0, &sub),
        Err(DimError::InvalidTheta0(_))
    ));
    // θ0 = 0 → InvalidTheta0.
    assert!(matches!(
        dimension_combline(&proj, 0.0, &sub),
        Err(DimError::InvalidTheta0(_))
    ));
    // A valid θ0 in (0, π/2) succeeds.
    assert!(dimension_combline(&proj, FRAC_PI_4, &sub).is_ok());
}
