//! coupling-matrix-s-001 (Filter Phase T9, ADR-0172): the complex
//! coupling-matrix → S-parameter gate for [`yee_filter::coupling_matrix_s_params`].
//!
//! This validates the complex `(S11, S21)` model (Hong & Lancaster *Microstrip
//! Filters* 2nd ed §8.1 general formulation, eq (8.30)–(8.31):
//! `[A] = [q] + jΩ·[U] − j·[m]`, `S21 = (2/√(qe_in·qe_out))·[A]⁻¹_{N1}`,
//! `S11 = 1 − (2/qe_in)·[A]⁻¹_{11}`) against an **independent** route, for
//! synthesized Chebyshev filters at **N = 3** and **N = 5**:
//!
//! 1. **Magnitude agreement (the load-bearing, NON-circular check):** `|S21|`
//!    from `coupling_matrix_s_params` matches `|S21|` from
//!    [`yee_filter::ideal_response`] — the closed-form characteristic-function
//!    magnitude (Chebyshev `1/(1+ε²T_N²(Ω))`), a *completely different*
//!    derivation — to within [`MAG_TOL`] across the band. Two independent
//!    synthesis routes agreeing is what validates the complex model: a wrong
//!    sign, scale, `qe`-normalization, or `Ω`-map breaks this by ~0.8–1.0 in
//!    `|S21|` (measured), so [`MAG_TOL`] is **far** below that gap and far above
//!    the observed ~2e-5 residual. **Do NOT weaken [`MAG_TOL`] to force green —
//!    the magnitude agreement IS the validation (ADR-0172 escape hatch).**
//! 2. **Losslessness:** the synthesized `[m]` is lossless, so
//!    `|S11|² + |S21|² ≈ 1` at every sweep point (tol [`LOSSLESS_TOL`]).
//! 3. **Phase is non-trivial + continuous:** `S21` (and `S11`) carry *varying*
//!    phase across the band — distinguishing the complex model from
//!    `ideal_response`'s flat (zero) phase — with no spurious ±2π jump between
//!    adjacent close samples (the phase is physically continuous).
//!
//! Pure-compute, non-`#[ignore]`'d: no FDTD, no `rand`, sub-millisecond.

use num_complex::Complex64;
use std::f64::consts::PI;
use yee_filter::{
    Approximation, CouplingMatrix, FilterSpec, Response, SpecMask, coupling_matrix_s_params,
    ideal_response, synthesize,
};

/// Magnitude-agreement tolerance on `|S21|` (`coupling_matrix_s_params` vs the
/// independent `ideal_response`). The observed worst-case deviation is ~2e-5; a
/// wrong sign / scale / `qe`-normalization / `Ω`-map deviates by ~0.8–1.0. This
/// `2e-3` bound is ~100× the real residual yet ~400× below the failure gap —
/// tight enough that any convention error fails. **Do NOT loosen this.**
const MAG_TOL: f64 = 2.0e-3;

/// Losslessness tolerance: `|S11|² + |S21|² ≈ 1`. The model is exactly lossless
/// (the residual is floating-point round-off, ~1e-15), so `1e-9` is generous
/// while still rejecting a non-passive S-matrix.
const LOSSLESS_TOL: f64 = 1.0e-9;

/// Build a Chebyshev 0.5 dB bandpass spec of the given order (f0 = 2 GHz,
/// FBW = 10 %, Z0 = 50 Ω) — the project's standard fixture shape.
fn cheb_spec(order: usize) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(order),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

/// In/near-band sweep grid: 0.85·f0 … 1.15·f0, `n` points. Dense enough that
/// adjacent points differ by a small phase step (so the continuity check is
/// meaningful), and wide enough to span the passband + skirts.
fn band_sweep(f0: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| f0 * (0.85 + 0.30 * i as f64 / (n - 1) as f64))
        .collect()
}

#[test]
fn coupling_matrix_s_001_magnitude_agreement_n3_n5() {
    // The load-bearing non-circular check: complex |S21| == characteristic-
    // function |S21| for BOTH N = 3 and N = 5.
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        let freqs = band_sweep(spec.f0_hz, 201);

        let ideal = ideal_response(&proj, &freqs);
        let cms = coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw);
        assert_eq!(cms.len(), freqs.len());

        let mut max_dev = 0.0_f64;
        for (k, &f) in freqs.iter().enumerate() {
            let mag_ideal = ideal[k].norm();
            let (_s11, s21) = cms[k];
            let dev = (mag_ideal - s21.norm()).abs();
            max_dev = max_dev.max(dev);
            assert!(
                dev <= MAG_TOL,
                "N={order} f={:.4} GHz: |S21| disagreement {dev:.3e} > {MAG_TOL:.1e} \
                 (coupling-matrix {:.6} vs ideal_response {mag_ideal:.6}) — a wrong \
                 sign/scale/Ω-map. Do NOT weaken MAG_TOL.",
                f / 1e9,
                s21.norm(),
            );
        }
        // Witness the achieved agreement is genuinely tight (not just under tol):
        // the real residual is ~2e-5, two orders below MAG_TOL.
        assert!(
            max_dev < MAG_TOL * 0.5,
            "N={order}: worst |S21| deviation {max_dev:.3e} should sit well inside \
             MAG_TOL ({MAG_TOL:.1e}); a near-tol value hints at a partial convention error"
        );
    }
}

#[test]
fn coupling_matrix_s_001_lossless() {
    // The synthesized coupling matrix is lossless: |S11|² + |S21|² == 1.
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        let freqs = band_sweep(spec.f0_hz, 201);
        let cms = coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw);

        for (k, &f) in freqs.iter().enumerate() {
            let (s11, s21) = cms[k];
            let power = s11.norm().powi(2) + s21.norm().powi(2);
            assert!(
                (power - 1.0).abs() <= LOSSLESS_TOL,
                "N={order} f={:.4} GHz: |S11|²+|S21|² = {power:.12} ≠ 1 (residual {:.3e} > {LOSSLESS_TOL:.0e}) \
                 — a lossless matrix must conserve power",
                f / 1e9,
                (power - 1.0).abs(),
            );
        }
    }
}

#[test]
fn coupling_matrix_s_001_phase_nontrivial_and_continuous() {
    // The complex model carries REAL phase (unlike ideal_response's flat zero
    // phase), and that phase is physically continuous (no spurious ±2π jumps).
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        // Dense grid so adjacent points differ by a SMALL phase step — this makes
        // the "no spurious ±2π jump" check non-vacuous.
        let freqs = band_sweep(spec.f0_hz, 401);
        let cms = coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw);

        // 1) ideal_response is flat-phase: every sample has ~zero argument.
        let ideal = ideal_response(&proj, &freqs);
        let ideal_phase_span = phase_span(ideal.iter().map(|z| z.arg()));
        assert!(
            ideal_phase_span < 1e-9,
            "ideal_response is supposed to be flat-phase, but its S21 phase spans \
             {ideal_phase_span:.3e} rad — the gate's premise (complex model adds phase) is wrong"
        );

        // 2) The complex S21 phase VARIES appreciably across the band.
        let s21_phases: Vec<f64> = cms.iter().map(|(_s11, s21)| s21.arg()).collect();
        let s21_span = phase_span(s21_phases.iter().copied());
        assert!(
            s21_span > 1.0,
            "N={order}: complex S21 phase span {s21_span:.3} rad is too flat — the model \
             should carry non-trivial phase (vs ideal_response's ~0)"
        );

        // 3) No spurious ±2π discontinuity: the smallest-arc step between adjacent
        // samples stays well below π (a real ±2π wrap would show ~2π raw / ~π arc).
        let max_jump = max_adjacent_phase_step(&s21_phases);
        assert!(
            max_jump < PI / 2.0,
            "N={order}: adjacent S21 phase step {max_jump:.3} rad ≥ π/2 — a spurious ±2π jump \
             (the phase must be continuous on this dense grid)"
        );

        // 4) S11 likewise carries real (varying) phase — not the lossless-quadrature
        // constant the old magnitude path used.
        let s11_phases: Vec<f64> = cms.iter().map(|(s11, _s21)| s11.arg()).collect();
        let s11_span = phase_span(s11_phases.iter().copied());
        assert!(
            s11_span > 1.0,
            "N={order}: complex S11 phase span {s11_span:.3} rad is too flat — S11 should carry \
             real phase too"
        );
    }
}

#[test]
fn coupling_matrix_s_001_degenerate_inputs() {
    // N = 0 (empty matrix) → one (0, 0) pair per frequency; f ≤ 0 → (1, 0)
    // (full reflection), matching ideal_response's f ≤ 0 → |S21| = 0 floor.
    let empty = CouplingMatrix {
        m: vec![],
        qe_in: 10.0,
        qe_out: 10.0,
    };
    let r = coupling_matrix_s_params(&empty, &[1.0e9, 2.0e9], 2.0e9, 0.10);
    assert_eq!(r.len(), 2);
    for (s11, s21) in r {
        assert_eq!(s11, Complex64::new(0.0, 0.0));
        assert_eq!(s21, Complex64::new(0.0, 0.0));
    }

    let proj = synthesize(&cheb_spec(3));
    let r = coupling_matrix_s_params(&proj.coupling, &[0.0, -1.0e9], 2.0e9, 0.10);
    assert_eq!(r.len(), 2);
    for (s11, s21) in r {
        assert_eq!(s21.norm(), 0.0, "f ≤ 0 must have no transmission");
        assert_eq!(s11, Complex64::new(1.0, 0.0), "f ≤ 0 must fully reflect");
    }
}

/// Span (max − min) of a sequence of phase angles, in radians.
fn phase_span(phases: impl Iterator<Item = f64>) -> f64 {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for p in phases {
        lo = lo.min(p);
        hi = hi.max(p);
    }
    hi - lo
}

/// Largest smallest-arc step between adjacent phase samples, in `[0, π]`. A
/// genuine ±2π `atan2` wrap appears as a raw `~2π` jump whose smallest-arc value
/// is `~0`; a real discontinuity would be `~π`. This collapses the harmless wrap
/// so only a *physical* jump trips the continuity check.
fn max_adjacent_phase_step(phases: &[f64]) -> f64 {
    let mut worst = 0.0_f64;
    for w in phases.windows(2) {
        let mut d = (w[1] - w[0]).abs();
        if d > PI {
            d = (2.0 * PI - d).abs(); // collapse the trivial atan2 wrap.
        }
        worst = worst.max(d);
    }
    worst
}
