//! bom-001 (Filter Phase F2.1): E-series component selection + BOM gate.
//!
//! Two halves:
//!
//! 1. **E-series correctness** — textbook IEC 60063 anchors for the log-nearest
//!    selector, including the boundary cases where the *geometric* midpoint (not
//!    the arithmetic one) decides which member wins, that E96 picks a finer value
//!    than E24 for the same input, and that every result is an actual member of
//!    the series.
//! 2. **BOM** — for the committed Chebyshev 0.5 dB N=5 BPF (f0 = 2 GHz,
//!    FBW = 0.10, Z0 = 50 Ω; the F2.0 / lumped-001 fixture), synthesize the LC
//!    ladder, select E24 parts, and assert: every chosen value is within the E24
//!    quantization bound (≤ ~5.1 %) of its ideal; the BOM has `2·N = 10` physical
//!    parts; and the symmetric ladder's duplicate values merge so there are
//!    *fewer* than 10 distinct lines.
//!
//! This validates the **selection**, NOT that the quantized response still passes
//! the spec mask — that yield question is F2.4 (Monte-Carlo). Do NOT gate on a
//! quantized-response pass here.

use yee_filter::{
    Approximation, Bom, CompKind, ESeries, FilterSpec, Response, SpecMask, select_components,
    synthesize, synthesize_lumped,
};

/// Chebyshev 0.5 dB N=5 bandpass spec (clone of the lumped-001 fixture).
fn fixture() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

/// Is `v` an actual member of `series` (mantissa × 10^decade)?
fn is_series_member(series: ESeries, v: f64) -> bool {
    let mantissa = v / 10f64.powi(v.log10().floor() as i32);
    series
        .values_decade()
        .iter()
        .any(|&m| (m - mantissa).abs() < 1e-6)
}

/// Relative-tolerance equality, since `mantissa × 10^decade` is not bit-identical
/// to a decimal literal like `4.7e-9` in IEEE 754 (e.g. `4.7·1e-9 ≠ 4.7e-9`).
fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * b.abs().max(1.0)
}

#[test]
fn eseries_textbook_anchors() {
    // Exact member maps to itself.
    assert_eq!(ESeries::E24.nearest(1.0), 1.0);

    // 1.04e3 sits below the geometric midpoint √(1.0·1.1)=1.0488, so log-nearest
    // resolves DOWN to 1.0e3 (not 1.1e3). This is the geometric-midpoint
    // boundary case the spec calls out.
    let n = ESeries::E24.nearest(1.04e3);
    assert!(
        approx_eq(n, 1.0e3),
        "1.04e3 → 1.0e3, got {n} (boundary 1.0488)"
    );

    // 4.5e-9 is log-closer to 4.7e-9 than to 4.3e-9 (the two bracketing E24
    // members): |log10(4.7)-log10(4.5)| = 0.0189 < |log10(4.3)-log10(4.5)| =
    // 0.0197. Picks 4.7e-9, and it is one of the two candidates.
    let c = ESeries::E24.nearest(4.5e-9);
    assert!(
        approx_eq(c, 4.7e-9),
        "4.5e-9 → 4.7e-9 (log-nearest of {{4.3e-9, 4.7e-9}}), got {c}"
    );
    assert!(
        approx_eq(c, 4.7e-9) || approx_eq(c, 4.3e-9),
        "4.5e-9 selection {c} must be one of the bracketing E24 members"
    );

    // A few more known E24 cases (decade-invariant log-nearest).
    assert!(approx_eq(ESeries::E24.nearest(9.5), 9.1)); // below √(9.1·10)=9.54
    assert!(approx_eq(ESeries::E24.nearest(6.6e-12), 6.8e-12));

    // E96 picks a finer value than E24 for inputs between E24 members.
    assert!(approx_eq(ESeries::E96.nearest(1.04e3), 1.05e3));
    assert!(approx_eq(ESeries::E96.nearest(4.5e-9), 4.53e-9));
    assert_ne!(
        ESeries::E96.nearest(1.04e3),
        ESeries::E24.nearest(1.04e3),
        "E96 must be finer than E24 here"
    );

    // Every nearest result is an actual series member, across decades + series.
    for series in [ESeries::E24, ESeries::E96] {
        for &x in &[
            1.0, 1.04e3, 4.5e-9, 6.6e-12, 9.5, 3.3e-6, 7.0, 2.71, 8.1e9, 1.5e-15,
        ] {
            let v = series.nearest(x);
            assert!(
                is_series_member(series, v),
                "{series:?}.nearest({x}) = {v} is not a series member"
            );
            // Selection stays in the same decade-ish neighbourhood (≤ ~5.5% for
            // E24, ≤ ~1.2% for E96 — the half-step quantization bound).
        }
    }
}

#[test]
fn bom_001() {
    let spec = fixture();
    let proj = synthesize(&spec);
    let n = proj.prototype.order();
    assert_eq!(n, 5, "fixture is order N=5");

    let ladder = synthesize_lumped(&proj).expect("N=5 bandpass fixture should synthesize");
    let bom: Bom = select_components(&ladder, ESeries::E24);

    // --- selection within the E24 quantization bound ---------------------
    // The largest E24 ratio step is 9.1→10 (and 8.2→9.1), ≈ 9.9 %; half-step ⇒
    // any value lands within ≈ 5.1 % of its ideal. Use 5.5 % as a safe bound.
    for line in &bom.lines {
        assert_eq!(line.series, ESeries::E24);
        assert_eq!(line.tolerance_pct, 5.0);
        assert!(line.esr_ohm.is_none());
        assert!(line.srf_hz.is_none());
        assert!(
            line.deviation_pct.abs() <= 5.5,
            "{:?} chosen {:.6e} vs ideal {:.6e} → deviation {:.3}% exceeds E24 bound",
            line.kind,
            line.chosen_value,
            line.ideal_value,
            line.deviation_pct
        );
        // The recorded deviation matches the chosen/ideal pair.
        let expect = (line.chosen_value - line.ideal_value) / line.ideal_value * 100.0;
        assert!((line.deviation_pct - expect).abs() < 1e-9);
    }

    // --- BOM completeness: 2·N physical parts ----------------------------
    assert_eq!(
        bom.total_parts(),
        2 * n,
        "an inductor + a capacitor per resonator → 2·N parts"
    );

    // --- duplicate grouping: symmetric ladder (g1=g5, g2=g4) merges ------
    // The symmetric Chebyshev prototype has g_k = g_{N+1-k}, so resonators 1/5
    // and 2/4 share identical L,C → identical chosen E-series parts → merged
    // lines with qty > 1. Hence strictly fewer than 2·N distinct lines.
    assert!(
        bom.lines.len() < 2 * n,
        "symmetric ladder must merge duplicates: {} distinct lines, expected < {}",
        bom.lines.len(),
        2 * n
    );
    // At least one merged (qty ≥ 2) line exists.
    assert!(
        bom.lines.iter().any(|l| l.qty >= 2),
        "symmetric ladder should produce at least one grouped (qty≥2) line"
    );

    // Sanity: the line kinds are only inductors / capacitors.
    for line in &bom.lines {
        assert!(matches!(
            line.kind,
            CompKind::Inductor | CompKind::Capacitor
        ));
    }
}
