//! Oracle ground-truth reference for the FEM-EM filter-S21 gate (FEM-EM brick 6,
//! ADR-0153).
//!
//! 3-pole Chebyshev 0.5 dB band-pass, f0 = 2 GHz, FBW = 10%, Z0 = 50 Ω.
//! Sweeps the realized lumped LC ladder S21 (the curve every EM method must
//! reproduce) AND the closed-form ideal response as a cross-check, and emits the
//! reference grid.
//!
//! The `main` is **assertion-backed**: it panics (non-zero exit) unless the
//! emitted reference grid has passband insertion-loss ripple ≤ 0.5 dB and its
//! −3 dB edges sit at ≈1.9 / 2.1 GHz (the 10% FBW band edges). The same checks
//! run as `#[test]`s so `cargo test -p yee-filter --examples` exercises them.

use yee_filter::{
    Approximation, FilterSpec, LumpedLadder, Response, SpecMask, ideal_response, ladder_s21,
    synthesize, synthesize_lumped,
};
use yee_synth::lowpass_to_bandpass;

fn db(mag: f64) -> f64 {
    20.0 * mag.log10()
}

/// Realization margin (dB) on the in-band ripple bound, mirroring the
/// [`lumped_001`](../../tests/lumped_001.rs) gate. The narrow-band LC transform
/// (Pozar §8.3) is *geometrically* symmetric (`f1·f2 = f0²`) so its ripple peaks
/// sit at the geometric edges, fractionally offset from the arithmetic `|Ω| = 1`
/// sample points — the realized N=3 ripple is therefore ≈ 0.5001 dB, a ~1e-4 dB
/// arithmetic-vs-geometric edge mismatch, **not** a synthesis error. `1e-3 dB`
/// is three orders of magnitude tighter than any component tolerance, so a
/// genuinely broken realization still fails. This is documented, not a
/// weakening.
const RIPPLE_REALIZATION_TOL_DB: f64 = 1e-3;

fn spec_for(fbw: f64) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw,
        order: Some(3),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(1.6e9, 0.0), (2.4e9, 0.0)],
        },
    }
}

/// The canonical reference ladder: 3-pole Cheb 0.5 dB BPF, f0=2 GHz, 10% FBW.
/// (Used by the `#[test]`s; `main` builds the same ladder alongside its project.)
#[cfg(test)]
fn reference_ladder() -> LumpedLadder {
    synthesize_lumped(&synthesize(&spec_for(0.10))).expect("bandpass N=3 synthesizes")
}

/// The realized-ladder −3 dB edges over `[1.5, 2.5] GHz`, found by fine scan.
///
/// Returns `(f_lo_hz, f_hi_hz)` — the first up-crossing and the following
/// down-crossing of the −3 dB level. `None` if the band is not resolved.
fn minus3db_edges(ladder: &LumpedLadder) -> Option<(f64, f64)> {
    let fine = 200_001;
    let mut lo_edge = None;
    let mut hi_edge = None;
    let mut prev_below = true;
    for i in 0..fine {
        let f = 1.5e9 + (2.5e9 - 1.5e9) * (i as f64) / ((fine - 1) as f64);
        let d = db(ladder_s21(ladder, f).norm());
        let below = d < -3.0;
        if prev_below && !below && lo_edge.is_none() {
            lo_edge = Some(f);
        }
        if !prev_below && below && lo_edge.is_some() {
            hi_edge = Some(f);
        }
        prev_below = below;
    }
    match (lo_edge, hi_edge) {
        (Some(lo), Some(hi)) => Some((lo, hi)),
        _ => None,
    }
}

/// Worst-case **equiripple** passband insertion-loss ripple (dB) of the realized
/// ladder. "In-band" is the Chebyshev equiripple region `|Ω| ≤ 1` under the
/// low-pass→band-pass map (mirroring [`yee_filter::mask_verdict`]) — **not** the
/// wider −3 dB band, whose edges sit at the −3 dB skirt where IL has already
/// climbed well past the 0.5 dB ripple floor. Ripple is `max_IL − min_IL` over
/// that band; for a 0.5 dB Chebyshev it should be ≈ 0.5 dB.
fn passband_ripple_db(ladder: &LumpedLadder, f0_hz: f64, fbw: f64) -> f64 {
    let n = 20_001;
    let mut min_il = f64::INFINITY;
    let mut max_il = f64::NEG_INFINITY;
    // Scan generously around the band and keep only the |Ω| <= 1 samples.
    for i in 0..n {
        let f = f0_hz * 0.75 + (f0_hz * 1.25 - f0_hz * 0.75) * (i as f64) / ((n - 1) as f64);
        if f <= 0.0 || lowpass_to_bandpass(f, f0_hz, fbw).abs() > 1.0 {
            continue;
        }
        let il = -db(ladder_s21(ladder, f).norm()); // insertion loss, dB (>= 0)
        min_il = min_il.min(il);
        max_il = max_il.max(il);
    }
    (max_il - min_il).max(0.0)
}

/// Print the -3 dB edges of a ladder over [1.5,2.5] GHz by fine scan.
fn report_edges(ladder: &LumpedLadder, label: &str) {
    if let Some((lo, hi)) = minus3db_edges(ladder) {
        println!(
            "{label} -3 dB edges: {:.4} / {:.4} GHz  BW={:.1} MHz ({:.2}% FBW)  Q~{:.1}",
            lo / 1e9,
            hi / 1e9,
            (hi - lo) / 1e6,
            (hi - lo) / 2.0e9 * 100.0,
            2.0e9 / (hi - lo)
        );
    }
}

fn main() {
    let spec = spec_for(0.10);

    let project = synthesize(&spec);
    let ladder = synthesize_lumped(&project).expect("bandpass N=3 synthesizes");

    println!("=== 3-pole Chebyshev 0.5 dB BPF: f0=2GHz FBW=10% Z0=50ohm ===");
    println!("prototype g-values: {:?}", project.prototype.g);
    println!("ladder resonators ({}):", ladder.resonators.len());
    for (i, r) in ladder.resonators.iter().enumerate() {
        println!(
            "  res[{i}] {:?}  L={:.6e} H  C={:.6e} F  (LC*w0^2 = {:.9})",
            r.branch,
            r.l_henry,
            r.c_farad,
            r.l_henry * r.c_farad * (std::f64::consts::TAU * 2.0e9).powi(2)
        );
    }

    // f1/f2 band edges (arithmetic): f0*(1 +/- FBW/2)
    let f1 = 2.0e9 * (1.0 - 0.10 / 2.0);
    let f2 = 2.0e9 * (1.0 + 0.10 / 2.0);
    println!(
        "\narithmetic band edges: f1={:.4} GHz  f2={:.4} GHz",
        f1 / 1e9,
        f2 / 1e9
    );

    // Dense table 1.5 - 2.5 GHz, plus the key callout points.
    println!(
        "\n{:>10}  {:>14}  {:>14}",
        "f (GHz)", "|S21| ladder dB", "|S21| ideal dB"
    );
    let n = 101;
    for i in 0..n {
        let f = 1.5e9 + (2.5e9 - 1.5e9) * (i as f64) / ((n - 1) as f64);
        let s21_lad = ladder_s21(&ladder, f).norm();
        let s21_ideal = ideal_response(&project, &[f])[0].norm();
        println!(
            "{:>10.4}  {:>14.4}  {:>14.4}",
            f / 1e9,
            db(s21_lad),
            db(s21_ideal)
        );
    }

    // Key callout points.
    println!("\n=== KEY REFERENCE VALUES ===");
    let pts = [
        ("stopband lo  1.60 GHz", 1.6e9),
        ("band edge f1 1.90 GHz", f1),
        ("center f0    2.00 GHz", 2.0e9),
        ("band edge f2 2.10 GHz", f2),
        ("stopband hi  2.40 GHz", 2.4e9),
    ];
    for (label, f) in pts {
        let s21_lad = db(ladder_s21(&ladder, f).norm());
        let s21_ideal = db(ideal_response(&project, &[f])[0].norm());
        println!(
            "  {label}:  ladder {:>9.4} dB   ideal {:>9.4} dB",
            s21_lad, s21_ideal
        );
    }

    // Find the -3 dB edges of the LADDER response by fine scan, and report.
    let edges = minus3db_edges(&ladder);
    if let Some((lo, hi)) = edges {
        println!(
            "\n-3 dB edges (ladder): f_lo={:.4} GHz  f_hi={:.4} GHz  BW={:.1} MHz  ({:.2}% FBW)",
            lo / 1e9,
            hi / 1e9,
            (hi - lo) / 1e6,
            (hi - lo) / 2.0e9 * 100.0
        );
    }

    // Robustness reference: a higher-Q (5% FBW) variant to challenge winners.
    println!("\n=== ROBUSTNESS REFERENCE: 5% FBW (higher-Q) ===");
    let hq_proj = synthesize(&spec_for(0.05));
    let hq_ladder = synthesize_lumped(&hq_proj).expect("hq synthesizes");
    report_edges(&hq_ladder, "5% FBW");
    report_edges(&ladder, "10% FBW");
    for (label, f) in [
        ("1.80 GHz", 1.80e9),
        ("1.95 GHz", 1.95e9),
        ("2.00 GHz", 2.00e9),
        ("2.05 GHz", 2.05e9),
        ("2.20 GHz", 2.20e9),
    ] {
        println!(
            "  5% FBW {label}: ladder {:>9.4} dB",
            db(ladder_s21(&hq_ladder, f).norm())
        );
    }

    // ---------------------------------------------------------------------
    // GATE (also covered by the #[test]s below): the emitted reference grid
    // must be a real 0.5 dB / 10% FBW band-pass. Panic (non-zero exit) if not,
    // so `cargo run --example oracle_reference` is itself a check.
    // ---------------------------------------------------------------------
    let (f_lo, f_hi) = edges.expect("ladder -3 dB band must resolve over 1.5-2.5 GHz");
    let ripple = passband_ripple_db(&ladder, spec.f0_hz, spec.fbw);
    println!(
        "\n[gate] equiripple passband ripple = {:.4} dB (<= 0.5)   -3dB edges = {:.4}/{:.4} GHz (~1.9/2.1)",
        ripple,
        f_lo / 1e9,
        f_hi / 1e9
    );
    assert!(
        ripple <= 0.5 + RIPPLE_REALIZATION_TOL_DB,
        "passband ripple {ripple:.4} dB exceeds the 0.5 dB Chebyshev spec \
         (+{RIPPLE_REALIZATION_TOL_DB:.0e} dB realization margin)"
    );
    // 10% FBW band-pass: arithmetic edges f1=1.9, f2=2.1 GHz. The realized
    // narrow-band LC ladder lands within ~1% of these (the slight asymmetry is
    // the arithmetic-vs-geometric band-edge mismatch of the LC transform).
    assert!(
        (f_lo - 1.90e9).abs() <= 0.03e9,
        "low -3dB edge {:.4} GHz is not ~1.90 GHz",
        f_lo / 1e9
    );
    assert!(
        (f_hi - 2.10e9).abs() <= 0.03e9,
        "high -3dB edge {:.4} GHz is not ~2.10 GHz",
        f_hi / 1e9
    );
    println!("[gate] oracle_reference OK");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_passband_ripple_within_half_db() {
        let ladder = reference_ladder();
        let ripple = passband_ripple_db(&ladder, 2.0e9, 0.10);
        assert!(
            ripple <= 0.5 + RIPPLE_REALIZATION_TOL_DB,
            "equiripple passband ripple {ripple:.4} dB exceeds 0.5 dB \
             (+{RIPPLE_REALIZATION_TOL_DB:.0e} dB realization margin)"
        );
        // Sanity: a real 0.5 dB Chebyshev *uses* its ripple budget — an
        // odd-order Cheb has equiripple maxima down to the ripple floor, so the
        // observed ripple should be close to (not far below) 0.5 dB. If it were
        // ~0 the band detection or the ladder would be wrong.
        assert!(
            ripple > 0.3,
            "ripple {ripple:.4} dB is suspiciously flat for a 0.5 dB Chebyshev"
        );
    }

    #[test]
    fn reference_minus3db_edges_at_band_edges() {
        let ladder = reference_ladder();
        let (f_lo, f_hi) = minus3db_edges(&ladder).expect("band resolves");
        // 10% FBW arithmetic band edges are 1.9 / 2.1 GHz.
        assert!(
            (f_lo - 1.90e9).abs() <= 0.03e9,
            "low edge {:.4} GHz not ~1.90 GHz",
            f_lo / 1e9
        );
        assert!(
            (f_hi - 2.10e9).abs() <= 0.03e9,
            "high edge {:.4} GHz not ~2.10 GHz",
            f_hi / 1e9
        );
        // The fractional bandwidth of the realized passband is ~10%.
        let fbw = (f_hi - f_lo) / 2.0e9;
        assert!(
            (fbw - 0.10).abs() <= 0.02,
            "realized FBW {:.3} not ~0.10",
            fbw
        );
    }

    #[test]
    fn reference_center_is_near_zero_loss() {
        // At f0 a lossless 0.5 dB Chebyshev BPF returns near 0 dB (one of the
        // equiripple maxima of an odd-order Chebyshev sits at band center).
        let ladder = reference_ladder();
        let il_center = -db(ladder_s21(&ladder, 2.0e9).norm());
        assert!(
            il_center >= -1e-6 && il_center <= 0.5 + 1e-6,
            "center insertion loss {il_center:.4} dB out of [0, 0.5]"
        );
    }
}
