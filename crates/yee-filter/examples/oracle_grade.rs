//! Oracle grader for the FEM-EM filter-S21 gate (FEM-EM brick 6, ADR-0153).
//!
//! Feed extracted `(freq_GHz, |S21|_dB)` points on the command line as `f:db`
//! pairs; the grader prints the reference value and the error at each, a
//! PASS/FAIL summary against the skeptic tolerances, AND the
//! **geometric-asymmetry discriminator**: this band-pass has its lower stopband
//! notch (1.6 GHz) genuinely DEEPER than its upper notch (2.4 GHz). The
//! asymmetry is real, not a fitted artifact: under the arithmetic-bandwidth
//! low-pass→band-pass map `Ω = (1/Δ)(f/f0 − f0/f)`, the geometric images of f0
//! are unequal — `|Ω(1.6 GHz)| = 4.5 > |Ω(2.4 GHz)| ≈ 3.67` — so the lower notch
//! sits farther out on the rejection skirt. The realized lumped ladder
//! reproduces this (≈ −41.8 dB at 1.6 GHz vs ≈ −36.3 dB at 2.4 GHz, margin
//! ≈ +5.5 dB), so the *correct* reference PASSES the discriminator. A response
//! that came back **symmetric** (margin ≈ 0) or **inverted** (upper notch
//! deeper, margin < 0) is FLAGGED: it has lost the band-pass-mapping asymmetry
//! and is not a faithful realization (the honest anti-fitted-artifact check — a
//! flat/symmetric curve is not evidence the EM solve captured the real
//! response).
//!
//! Reference: 3-pole Cheb 0.5 dB BPF, f0=2 GHz, FBW=10%, Z0=50 Ω.
//! Tolerances: passband |err| <= 2 dB, rejection |err| <= 5 dB.

use yee_filter::{
    Approximation, FilterSpec, LumpedLadder, Response, SpecMask, ladder_s21, synthesize,
    synthesize_lumped,
};

fn db(mag: f64) -> f64 {
    20.0 * mag.log10()
}

/// The two stopband notch frequencies the asymmetry discriminator compares.
const F_LO_NOTCH_GHZ: f64 = 1.6;
const F_HI_NOTCH_GHZ: f64 = 2.4;

fn reference_spec() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(3),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![],
        },
    }
}

fn reference_ladder() -> LumpedLadder {
    synthesize_lumped(&synthesize(&reference_spec())).expect("bandpass N=3 synthesizes")
}

/// Result of the geometric-asymmetry discriminator over an extracted |S21|
/// curve.
#[derive(Debug, Clone, Copy, PartialEq)]
struct AsymmetryVerdict {
    /// Rejection depth (positive dB down) at the lower 1.6 GHz notch.
    depth_lo_db: f64,
    /// Rejection depth (positive dB down) at the upper 2.4 GHz notch.
    depth_hi_db: f64,
    /// `depth_lo − depth_hi`: positive means the lower notch is deeper (the
    /// physical microstrip signature); ≤ 0 means symmetric / inverted.
    margin_db: f64,
    /// `true` iff the lower notch is deeper than the upper by at least the
    /// required margin.
    pass: bool,
}

/// Linearly interpolate the extracted curve `pts` (sorted `(f_ghz, db)`) at
/// `f_ghz`. Clamps to the endpoints outside the sampled range. `None` if the
/// curve is empty.
fn interp_db(pts: &[(f64, f64)], f_ghz: f64) -> Option<f64> {
    if pts.is_empty() {
        return None;
    }
    if f_ghz <= pts[0].0 {
        return Some(pts[0].1);
    }
    if f_ghz >= pts[pts.len() - 1].0 {
        return Some(pts[pts.len() - 1].1);
    }
    for w in pts.windows(2) {
        let (f0, d0) = w[0];
        let (f1, d1) = w[1];
        if (f0..=f1).contains(&f_ghz) {
            let t = if (f1 - f0).abs() < 1e-15 {
                0.0
            } else {
                (f_ghz - f0) / (f1 - f0)
            };
            return Some(d0 + t * (d1 - d0));
        }
    }
    Some(pts[pts.len() - 1].1)
}

/// Run the geometric-asymmetry discriminator over an extracted curve.
///
/// Computes the rejection depth (positive dB = how far |S21| is below 0 dB) at
/// the 1.6 GHz and 2.4 GHz notches by interpolating the supplied curve, then
/// checks the lower notch is DEEPER by at least `margin_db`. This is a real
/// computation over the curve — it never hardcodes the verdict.
fn asymmetry_discriminator(pts: &[(f64, f64)], margin_db: f64) -> Option<AsymmetryVerdict> {
    let d_lo = interp_db(pts, F_LO_NOTCH_GHZ)?;
    let d_hi = interp_db(pts, F_HI_NOTCH_GHZ)?;
    // Depth = how far DOWN from 0 dB (a -40 dB point has depth +40).
    let depth_lo_db = -d_lo;
    let depth_hi_db = -d_hi;
    let margin = depth_lo_db - depth_hi_db;
    Some(AsymmetryVerdict {
        depth_lo_db,
        depth_hi_db,
        margin_db: margin,
        pass: margin >= margin_db,
    })
}

/// Parse `f:db` CLI args into a sorted `(f_ghz, db)` curve.
fn parse_pairs<I: IntoIterator<Item = String>>(args: I) -> Vec<(f64, f64)> {
    let mut pts: Vec<(f64, f64)> = args
        .into_iter()
        .filter_map(|arg| {
            let (fs, ds) = arg.split_once(':')?;
            let f = fs.parse::<f64>().ok()?;
            let d = ds.parse::<f64>().ok()?;
            Some((f, d))
        })
        .collect();
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    pts
}

/// Passband (near-band) tolerance in dB: `|extracted − reference|` must stay
/// within this over ~[1.85, 2.15] GHz.
const PASSBAND_TOL_DB: f64 = 2.0;
/// Stopband / rejection-skirt tolerance in dB (looser — the skirt is steep, so
/// a small frequency offset reads as a large dB error there).
const REJECTION_TOL_DB: f64 = 5.0;

/// Worst-case band errors + asymmetry result of grading an extracted curve
/// against the reference ladder. `asymmetry_pass` is `None` when the curve
/// lacks the points the discriminator needs.
struct GradeOutcome {
    worst_pass_db: f64,
    worst_rej_db: f64,
    asymmetry_pass: Option<bool>,
}

/// Pure grading decision: compare each extracted point to the reference ladder
/// and run the geometric-asymmetry discriminator. This is the function the
/// exit-code gate AND its tests use, so the pass/fail decision can never drift
/// from what is tested.
fn evaluate(pairs: &[(f64, f64)]) -> GradeOutcome {
    let ladder = reference_ladder();
    let mut worst_pass_db = 0.0_f64;
    let mut worst_rej_db = 0.0_f64;
    for &(f_ghz, d_meas) in pairs {
        let d_ref = db(ladder_s21(&ladder, f_ghz * 1e9).norm());
        let err = (d_meas - d_ref).abs();
        if (1.85..=2.15).contains(&f_ghz) {
            worst_pass_db = worst_pass_db.max(err);
        } else {
            worst_rej_db = worst_rej_db.max(err);
        }
    }
    GradeOutcome {
        worst_pass_db,
        worst_rej_db,
        asymmetry_pass: asymmetry_discriminator(pairs, ASYMMETRY_MARGIN_DB).map(|v| v.pass),
    }
}

/// The gate: a submission FAILS if either band exceeds tolerance OR the
/// geometric-asymmetry discriminator flags it (a symmetric/inverted curve has
/// lost the band-pass-mapping asymmetry and is not a credible EM result).
fn outcome_failed(o: &GradeOutcome) -> bool {
    o.worst_pass_db > PASSBAND_TOL_DB
        || o.worst_rej_db > REJECTION_TOL_DB
        || o.asymmetry_pass == Some(false)
}

fn main() {
    let pairs = parse_pairs(std::env::args().skip(1));
    if pairs.is_empty() {
        eprintln!("usage: oracle_grade 1.6:-40.1 2.0:-0.2 2.4:-33.5 ...");
        std::process::exit(2);
    }
    let ladder = reference_ladder();

    // Passband is the -0.5 dB ripple band ~[1.9, 2.1] GHz; outside is stopband
    // skirt. Use 1.85-2.15 GHz as the "in/near band" tol zone.
    println!(
        "{:>10}  {:>10}  {:>10}  {:>8}  verdict",
        "f(GHz)", "extracted", "reference", "err"
    );
    for &(f_ghz, d_meas) in &pairs {
        let d_ref = db(ladder_s21(&ladder, f_ghz * 1e9).norm());
        let err = d_meas - d_ref;
        let tol = if (1.85..=2.15).contains(&f_ghz) {
            PASSBAND_TOL_DB
        } else {
            REJECTION_TOL_DB
        };
        println!(
            "{:>10.4}  {:>10.3}  {:>10.3}  {:>+8.3}  {}",
            f_ghz,
            d_meas,
            d_ref,
            err,
            if err.abs() <= tol { "ok" } else { "FAIL" }
        );
    }

    let outcome = evaluate(&pairs);
    println!(
        "\nworst passband err = {:.3} dB (tol {})   worst rejection err = {:.3} dB (tol {})",
        outcome.worst_pass_db, PASSBAND_TOL_DB, outcome.worst_rej_db, REJECTION_TOL_DB
    );
    if let Some(v) = asymmetry_discriminator(&pairs, ASYMMETRY_MARGIN_DB) {
        println!(
            "\n[asymmetry] depth(1.6 GHz)={:.2} dB  depth(2.4 GHz)={:.2} dB  margin={:+.2} dB  -> {}",
            v.depth_lo_db,
            v.depth_hi_db,
            v.margin_db,
            if v.pass {
                "PASS (lower notch deeper: physical microstrip signature)"
            } else {
                "FLAG (symmetric/inverted: not a geometry-aware EM result)"
            }
        );
    }

    // The exit code IS the gate: 0 only if every point is within tolerance AND
    // the geometric-asymmetry discriminator passed. (The prior version always
    // exited 0 on the CLI path, so a CI/B7 caller checking `$?` would rubber-
    // stamp a broken submission — review P0.)
    if outcome_failed(&outcome) {
        eprintln!(
            "\n[oracle_grade] FAIL: submission exceeds tolerance or fails the \
             geometric-asymmetry discriminator"
        );
        std::process::exit(1);
    }
    println!("\n[oracle_grade] PASS");
}

/// Minimum required depth excess (lower − upper notch) for the asymmetry
/// discriminator to PASS. The faithful reference clears this with ≈ +5.5 dB of
/// margin; the 1 dB bar cleanly separates that from a symmetric response
/// (margin ≈ 0, FLAGGED) and an inverted one (margin < 0, FLAGGED).
const ASYMMETRY_MARGIN_DB: f64 = 1.0;

#[cfg(test)]
mod tests {
    use super::*;

    /// Sample the realized ladder |S21| (dB) onto a curve over [1.5, 2.5] GHz.
    fn sample_ladder_curve(ladder: &LumpedLadder, n: usize) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let f_ghz = 1.5 + (2.5 - 1.5) * (i as f64) / ((n - 1) as f64);
                (f_ghz, db(ladder_s21(ladder, f_ghz * 1e9).norm()))
            })
            .collect()
    }

    #[test]
    fn gate_predicate_flags_each_failure_mode_independently() {
        // All-pass: within both tolerances, asymmetry passed -> NOT failed.
        assert!(!outcome_failed(&GradeOutcome {
            worst_pass_db: 0.1,
            worst_rej_db: 0.1,
            asymmetry_pass: Some(true),
        }));
        // Passband over tolerance alone -> failed.
        assert!(outcome_failed(&GradeOutcome {
            worst_pass_db: PASSBAND_TOL_DB + 0.5,
            worst_rej_db: 0.1,
            asymmetry_pass: Some(true),
        }));
        // Rejection over tolerance alone -> failed.
        assert!(outcome_failed(&GradeOutcome {
            worst_pass_db: 0.1,
            worst_rej_db: REJECTION_TOL_DB + 0.5,
            asymmetry_pass: Some(true),
        }));
        // Within tolerance but asymmetry FLAGGED -> failed (the anti-fitted-
        // artifact path: a symmetric/inverted curve must not pass the gate).
        assert!(outcome_failed(&GradeOutcome {
            worst_pass_db: 0.1,
            worst_rej_db: 0.1,
            asymmetry_pass: Some(false),
        }));
    }

    /// End-to-end: the faithful reference curve passes the FULL gate (tolerance
    /// against itself + asymmetry). If this failed, the exit-code gate would
    /// reject a perfect submission.
    #[test]
    fn gate_passes_faithful_reference_curve() {
        let ladder = reference_ladder();
        let curve = sample_ladder_curve(&ladder, 2001);
        assert!(
            !outcome_failed(&evaluate(&curve)),
            "the faithful reference must PASS the full grading gate"
        );
    }

    /// End-to-end: a submission matching the reference everywhere EXCEPT a +3 dB
    /// error at 2.0 GHz (passband tol 2.0) must fail the gate -> the CLI path now
    /// exits non-zero (review P0: it previously always exited 0).
    #[test]
    fn gate_fails_over_tolerance_submission() {
        let ladder = reference_ladder();
        let mut curve = sample_ladder_curve(&ladder, 2001);
        // Bump the point nearest 2.0 GHz by +3 dB.
        let idx = curve
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (a.0 - 2.0).abs().partial_cmp(&(b.0 - 2.0).abs()).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        curve[idx].1 += 3.0;
        let outcome = evaluate(&curve);
        assert!(
            outcome.worst_pass_db > PASSBAND_TOL_DB,
            "worst passband err {:.3} should exceed tol after +3 dB bump",
            outcome.worst_pass_db
        );
        assert!(
            outcome_failed(&outcome),
            "an in-band +3 dB error must fail the gate"
        );
    }

    /// (a) Grading the reference against a *perturbed copy of itself* must pass
    /// the tolerance mask: a small (sub-tolerance) perturbation stays "ok"
    /// everywhere, and a perturbation that blows past tolerance FAILs. This
    /// proves the grader actually compares against the reference rather than
    /// rubber-stamping.
    #[test]
    fn grading_perturbed_self_passes_within_tol_and_fails_beyond() {
        let ladder = reference_ladder();
        let key_pts_ghz = [1.6, 1.9, 2.0, 2.1, 2.4];

        // Sub-tolerance perturbation: +1.0 dB everywhere (passband tol is 2.0,
        // rejection tol is 5.0) -> every point stays within tolerance.
        for &f_ghz in &key_pts_ghz {
            let f = f_ghz * 1e9;
            let d_ref = db(ladder_s21(&ladder, f).norm());
            let d_meas = d_ref + 1.0;
            let err = (d_meas - d_ref).abs();
            let in_band = (1.85..=2.15).contains(&f_ghz);
            let tol = if in_band { 2.0 } else { 5.0 };
            assert!(
                err <= tol,
                "sub-tol perturbation flagged at {f_ghz} GHz: err {err:.3} > tol {tol}"
            );
        }

        // Over-tolerance perturbation in the passband: +3.0 dB at 2.0 GHz
        // (passband tol 2.0) MUST fail. If the grader passed this, it would be
        // tautological.
        let f = 2.0e9;
        let d_ref = db(ladder_s21(&ladder, f).norm());
        let d_meas = d_ref + 3.0;
        let err = (d_meas - d_ref).abs();
        assert!(
            err > 2.0,
            "an in-band +3 dB error ({err:.3}) should exceed the 2.0 dB passband tol"
        );
    }

    /// The faithful reference (realized lumped ladder) reproduces the real
    /// band-pass-mapping asymmetry — lower notch ≈ −41.8 dB, upper ≈ −36.3 dB,
    /// margin ≈ +5.5 dB — so it must PASS the discriminator. This is the
    /// ground-truth-passes-its-own-check guard: if the *reference* failed, the
    /// discriminator's sign/orientation would be wrong.
    #[test]
    fn discriminator_passes_faithful_reference() {
        let ladder = reference_ladder();
        let curve = sample_ladder_curve(&ladder, 2001);
        let v = asymmetry_discriminator(&curve, ASYMMETRY_MARGIN_DB).expect("curve non-empty");
        // Sanity-pin the measured depths to the known reference values.
        assert!(
            (v.depth_lo_db - 41.77).abs() < 0.5,
            "lower-notch depth {:.3} dB not ~41.8",
            v.depth_lo_db
        );
        assert!(
            (v.depth_hi_db - 36.26).abs() < 0.5,
            "upper-notch depth {:.3} dB not ~36.3",
            v.depth_hi_db
        );
        assert!(
            v.margin_db > 4.0,
            "reference asymmetry margin {:+.3} dB should be ~+5.5",
            v.margin_db
        );
        assert!(
            v.pass,
            "the faithful reference must PASS the asymmetry discriminator"
        );
    }

    /// A SYMMETRIC response (equal notch depths) is the anti-fitted-artifact
    /// target: margin ≈ 0, so it must be FLAGGED. A method whose curve came back
    /// symmetric has lost the band-pass-mapping asymmetry and is not credited.
    #[test]
    fn discriminator_flags_symmetric_input() {
        // Equal -38 dB notches either side of a ~0 dB passband.
        let curve = vec![
            (1.60, -38.0),
            (1.90, -3.0),
            (2.00, -0.2),
            (2.10, -3.0),
            (2.40, -38.0),
        ];
        let v = asymmetry_discriminator(&curve, ASYMMETRY_MARGIN_DB).expect("curve non-empty");
        assert!(
            v.margin_db.abs() < 1e-9,
            "symmetric curve should have ~0 margin, got {:+.3} dB",
            v.margin_db
        );
        assert!(
            !v.pass,
            "symmetric (equal-depth) response must be FLAGGED, not passed"
        );
    }

    /// On a synthetic ASYMMETRIC curve (lower notch deeper than upper) the
    /// discriminator must FIRE (pass): depth(1.6) > depth(2.4). This is a hand-
    /// built curve, so the depths are known and the verdict is checkable.
    #[test]
    fn discriminator_fires_on_asymmetric_input() {
        // Lower notch 1.6 GHz at -45 dB (depth 45); upper notch 2.4 GHz at
        // -30 dB (depth 30); passband ~0 dB. margin = 45 - 30 = +15 dB.
        let curve = vec![
            (1.50, -30.0),
            (1.60, -45.0),
            (1.70, -20.0),
            (1.90, -3.0),
            (2.00, -0.2),
            (2.10, -3.0),
            (2.30, -18.0),
            (2.40, -30.0),
            (2.50, -22.0),
        ];
        let v = asymmetry_discriminator(&curve, ASYMMETRY_MARGIN_DB).expect("curve non-empty");
        assert!(
            (v.depth_lo_db - 45.0).abs() < 1e-9,
            "depth_lo {:.3} != 45",
            v.depth_lo_db
        );
        assert!(
            (v.depth_hi_db - 30.0).abs() < 1e-9,
            "depth_hi {:.3} != 30",
            v.depth_hi_db
        );
        assert!(
            (v.margin_db - 15.0).abs() < 1e-9,
            "margin {:.3} != 15",
            v.margin_db
        );
        assert!(
            v.pass,
            "deeper-lower-notch curve must PASS the discriminator"
        );
    }

    /// An INVERTED curve (upper notch deeper than lower) must be FLAGGED:
    /// margin < 0, pass == false. Guards against the discriminator accidentally
    /// using `abs()` (which would credit either asymmetry direction).
    #[test]
    fn discriminator_flags_inverted_input() {
        let curve = vec![
            (1.60, -28.0), // depth 28
            (2.00, -0.2),
            (2.40, -50.0), // depth 50 (upper deeper -> wrong direction)
        ];
        let v = asymmetry_discriminator(&curve, ASYMMETRY_MARGIN_DB).expect("curve non-empty");
        assert!(
            v.margin_db < 0.0,
            "inverted curve should have negative margin, got {:+.3}",
            v.margin_db
        );
        assert!(
            !v.pass,
            "inverted (upper-deeper) curve must be FLAGGED, not passed"
        );
    }

    /// Interpolation sanity: midpoint of two points is their average, and the
    /// exact-frequency lookups land on the supplied values.
    #[test]
    fn interp_is_linear_and_exact_on_nodes() {
        let curve = vec![(1.6, -40.0), (2.0, 0.0), (2.4, -30.0)];
        assert!((interp_db(&curve, 1.6).unwrap() + 40.0).abs() < 1e-9);
        assert!((interp_db(&curve, 2.4).unwrap() + 30.0).abs() < 1e-9);
        // Midpoint 1.8 GHz between (1.6,-40) and (2.0,0) -> -20.
        assert!((interp_db(&curve, 1.8).unwrap() + 20.0).abs() < 1e-9);
        // Clamps outside range.
        assert!((interp_db(&curve, 1.0).unwrap() + 40.0).abs() < 1e-9);
        assert!((interp_db(&curve, 3.0).unwrap() + 30.0).abs() < 1e-9);
    }

    /// `parse_pairs` parses and sorts `f:db` tokens and ignores garbage.
    #[test]
    fn parse_pairs_sorts_and_filters() {
        let args = ["2.4:-30.0", "1.6:-40.0", "junk", "2.0:-0.2", "x:y"]
            .iter()
            .map(|s| s.to_string());
        let pts = parse_pairs(args);
        assert_eq!(pts.len(), 3, "two garbage tokens dropped");
        assert_eq!(pts[0].0, 1.6, "sorted ascending by frequency");
        assert_eq!(pts[2].0, 2.4);
    }
}
