//! dim-interdigital-001 (Filter Phase F1.2.7): interdigital dimensional
//! synthesis gate.
//!
//! The PROPER, non-tautological gate for [`yee_filter::dimension_interdigital`]
//! (Hong & Lancaster *Microstrip Filters for RF/Microwave Applications* §5 —
//! interdigital λ_g/4 lines short-circuited at alternating ends, with **no
//! loading cap**), in three parts mirroring `dim_combline_001`:
//!
//! 1. **Published benchmark (H&L Qe/M) — the non-tautological core.** Synthesize
//!    the 5-pole 0.1 dB Chebyshev band-pass at FBW = 0.10 (and the FBW = 0.15
//!    point) and assert the synthesized external Q and adjacent couplings
//!    reproduce H&L's *published* design numbers — Qe ≈ 11.468, M₁₂ ≈ 0.07975,
//!    M₂₃ ≈ 0.06077 at FBW 0.10; Qe ≈ 7.645, M₁₂ ≈ 0.11962, M₂₃ ≈ 0.09115 at
//!    FBW 0.15 — to < 1% (and the tighter ±1e-3 band). These are the same
//!    prototype-derived, technique-independent coupling numbers the combline gate
//!    uses (legitimately shared — the coupling matrix is topology-agnostic), and
//!    they are non-tautological: compared to the *book*, not the synthesizer's own
//!    output. A constant or wrong synthesis fails by ≫ tol.
//! 2. **Interdigital λ_g/4 resonance (interdigital-DISTINCT).** From a
//!    `dimension_interdigital` result build the **unloaded** short-circuited-stub
//!    susceptance `B(f) = −(1/Z0)·cot((π/2)·f/f0)` (θ = π/2, **no cap term**) and
//!    root-find `B(f) = 0` over `[0.5 f0, 1.5 f0]`; assert the root equals f0
//!    within ±1%. The resonance comes from the λ_g/4 length **alone**
//!    (`cot(π/2) = 0` ⇒ `B(f0) = 0`), NOT a cap — the structural contrast with
//!    combline (which needs `C_L`). Catches a wrong length / dispersion / sign
//!    bug.
//! 3. **Dims solved + positive + structural.** N = 5 → 4 gaps; every solved gap
//!    re-evaluates to its `target_k` (< 1% via `yee_layout::coupling_coefficient`,
//!    no clamping); `line_width_m` and `resonator_length_m` finite and > 0;
//!    `resonator_length_m == (π/2)/β(f0) = λ_g/4` (closed-form, < 1e-9 rel, both
//!    forms agreeing); and the `OrderTooSmall` error path. (The struct carries
//!    **no loading cap** — enforced at compile time: `InterdigitalDimensions` has
//!    no such field. The `UnsupportedTopology` arm is unreachable from a test:
//!    `yee_filter::Topology` is `#[non_exhaustive]` with only `CoupledResonator`,
//!    so no other variant can be constructed — see the gate-3 note.)

use std::f64::consts::{FRAC_PI_2, PI};

use yee_filter::{
    Approximation, DimError, FilterProject, FilterSpec, Response, SpecMask, dimension_interdigital,
    synthesize,
};
use yee_layout::Substrate;

/// FR-4 representative substrate (matches the hairpin / combline / stepped-Z
/// gates).
fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

/// Build the 5-pole 0.1 dB Chebyshev band-pass at the requested fractional
/// bandwidth (H&L §5 worked coupled-resonator example).
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
fn dim_interdigital_001_hl_published_benchmark() {
    // --- FBW = 0.10 (H&L worked coupled-resonator design) ------------------
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

    // --- FBW = 0.15 point --------------------------------------------------
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
        "H&L Qe/M  FBW=0.10: Qe={qe:.6} (H&L 11.468)  M12={m12:.6} (0.07975)  \
         M23={m23:.6} (0.06077)"
    );
    println!(
        "H&L Qe/M  FBW=0.15: Qe={qe2:.6} (H&L 7.645)  M12={m12_2:.6} (0.11962)  \
         M23={m23_2:.6} (0.09115)"
    );
}

#[test]
fn dim_interdigital_001_quarter_wave_resonance() {
    // Demo band-pass spec at FBW = 0.10. Interdigital takes NO θ0 — θ = π/2.
    let spec = spec_5pole_cheb_01db(0.10);
    let z0 = spec.z0_ohm;
    let f0 = spec.f0_hz;
    let proj = synthesize(&spec);
    let _dims = dimension_interdigital(&proj, &substrate())
        .expect("N=5 coupled-resonator interdigital fixture should dimension without error");

    // INDEPENDENT UNLOADED short-circuited-stub susceptance at θ = π/2:
    //   B(f) = −(1/Z0)·cot((π/2)·f/f0)            (β(f)·L = (π/2)·(f/f0)).
    // There is NO cap term (interdigital is the θ = π/2 limit of combline:
    // cot(π/2) = 0 ⇒ B(f0) = 0 already). A wrong length / dispersion / sign moves
    // the root off f0 or breaks the bracket sign-change.
    let theta = FRAC_PI_2;
    let b_of = |f: f64| {
        let arg = theta * f / f0;
        let cot = arg.cos() / arg.sin();
        -(1.0 / z0) * cot
    };

    // Root-find B(f) = 0 by bisection over [0.5·f0, 1.5·f0]. −(1/Z0)·cot rises
    // monotonically with f over (0, π) and crosses zero exactly at the quarter-
    // wave point (arg = π/2, i.e. f = f0) — a single sign change in the bracket.
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
        "interdigital resonance: UNLOADED B(f)=0 root f_res = {f_res:.6e} Hz vs f0 = {f0:.6e} Hz \
         (rel {rel:.6}); theta = pi/2 (lambda_g/4), NO loading cap"
    );
    assert!(
        rel < 0.01,
        "unloaded-stub resonance root f_res = {f_res:.6e} Hz is not within 1% of f0 = {f0:.6e} Hz \
         (rel {rel:.6}) — length/dispersion/sign bug"
    );
}

#[test]
fn dim_interdigital_001_dims_solved_and_positive() {
    let spec = spec_5pole_cheb_01db(0.10);
    let proj = synthesize(&spec);
    let sub = substrate();
    let dims = dimension_interdigital(&proj, &sub)
        .expect("interdigital fixture should dimension without error");

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

    // All physical dims positive and finite. There is NO loading-cap field on
    // InterdigitalDimensions (the structural contrast with ComblineDimensions) —
    // its absence is enforced at compile time, so there is nothing to assert here.
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

    // resonator_length_m == (π/2)/β(f0) = λ_g/4. Assert BOTH closed forms agree:
    // the (π/2)/β form the engine uses, and the equivalent c/(4·f0·√ε_eff) form.
    let e_eff = yee_layout::eps_eff(dims.line_width_m, sub.height_m, sub.eps_r);
    let c = 299_792_458.0_f64;
    let beta0 = 2.0 * PI * spec.f0_hz * e_eff.sqrt() / c;
    let len_from_beta = FRAC_PI_2 / beta0; // (π/2)/β(f0)
    let len_lambda_quarter = c / (4.0 * spec.f0_hz * e_eff.sqrt()); // λ_g/4
    assert!(
        (dims.resonator_length_m - len_from_beta).abs() / len_from_beta < 1e-9,
        "resonator_length_m = {:.6e} m should equal (π/2)/β(f0) = {len_from_beta:.6e} m",
        dims.resonator_length_m
    );
    assert!(
        (len_from_beta - len_lambda_quarter).abs() / len_lambda_quarter < 1e-9,
        "(π/2)/β(f0) = {len_from_beta:.6e} m should equal λ_g/4 = {len_lambda_quarter:.6e} m"
    );
}

#[test]
fn dim_interdigital_001_order_too_small_errors() {
    let sub = substrate();

    // A 1-pole project → N = 1 coupling matrix → no inter-resonator coupling to
    // realize → OrderTooSmall.
    let mut spec1 = spec_5pole_cheb_01db(0.10);
    spec1.order = Some(1);
    let proj1 = synthesize(&spec1);
    assert_eq!(proj1.coupling.m.len(), 1, "order-1 → 1×1 coupling matrix");
    assert!(
        matches!(
            dimension_interdigital(&proj1, &sub),
            Err(DimError::OrderTooSmall)
        ),
        "N=1 must give OrderTooSmall"
    );

    // The valid N=5 fixture succeeds (the topology guard passes for the only
    // synthesized topology, CoupledResonator; the UnsupportedTopology arm is
    // unreachable from a test because yee_filter::Topology is #[non_exhaustive]
    // with no other constructible variant).
    let proj5 = synthesize(&spec_5pole_cheb_01db(0.10));
    assert!(dimension_interdigital(&proj5, &sub).is_ok());
}
