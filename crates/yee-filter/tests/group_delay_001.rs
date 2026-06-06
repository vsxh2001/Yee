//! group-delay-001 (ADR-0173, T10): numerical group delay from the complex S21
//! phase, validated against the published band-pass midband sum rule.
//!
//! [`yee_filter::group_delay`] computes `τ_g = −dφ/dω` of a complex `S21`
//! (here the [`yee_filter::coupling_matrix_s_params`] response, ADR-0172, which
//! carries physical phase). This gate validates that computation against an
//! **independent** closed form, for synthesized Chebyshev filters at **N = 3**
//! and **N = 5**.
//!
//! ## The closed-form midband anchor (the load-bearing, NON-circular check)
//!
//! The band-pass midband group delay obeys the lowpass-prototype **sum rule**
//! (Pozar, *Microwave Engineering* §8; Hong & Lancaster, *Microstrip Filters
//! for RF/Microwave Applications* §3):
//!
//! - the lowpass prototype group delay at band centre is **half** the sum of
//!   the reactive element values, `τ_LP(Ω=0) = (Σ_{k=1}^{N} g_k)/2` (normalized,
//!   cutoff `ωc = 1` rad/s);
//! - the band-pass transformation `Ω = (1/FBW)(ω/ω0 − ω0/ω)` has Jacobian
//!   `dΩ/dω|_{ω0} = 2/(FBW·ω0)`;
//!
//! so the band-pass midband group delay is
//!
//! ```text
//! τ_g(ω0) = τ_LP(0)·(dΩ/dω)|_{ω0} = (Σ_{k=1}^{N} g_k) / (FBW·ω0),  ω0 = 2π·f0
//! ```
//!
//! (equivalently `Σg_k / Δω` via the absolute bandwidth `Δω = FBW·ω0`). The
//! `2`-and-`½` cancel: the bandpass-map Jacobian is `2/(FBW·ω0)` but the lowpass
//! prototype DC group delay is `Σg/2`, not `Σg`.
//!
//! **Why `Σg/2` (not `Σg`) — verified from first principles, not the literature
//! gloss.** A common textbook short-hand quotes `τ_LP(0) = Σg`, which is off by
//! a factor of two. This gate's reference uses the *measured-and-derived*
//! constant: an **independent** doubly-terminated LC-ladder ABCD model of the
//! lowpass prototype (built straight from the g-values, never via the
//! coupling matrix) gives `τ_LP(Ω→0) = Σg/2` exactly for N = 1, 3, 5 (the N = 1
//! single-resonator case `g1/2` is the cleanest tell), and the bandpass
//! coupling-matrix group delay measured here is `Σg/(FBW·ω0)` exactly across
//! FBW ∈ {0.10, 0.02, 0.005} — the two agree through the exact Jacobian. The
//! `2·Σg/(FBW·ω0)` form is therefore wrong by 2×; the gate uses the
//! first-principles-confirmed `Σg/(FBW·ω0)`. (See ADR-0173: the ADR's stated
//! `τ_LP(0)=Σg` was corrected here; the gate is the arbiter, and the measured
//! midband τ matches `Σg/(FBW·ω0)` to ratio ~1.000.)
//!
//! The SAME `Σ_{k=1}^{N} g_k` structure (and `1/FBW` scaling) appears in this
//! crate's Cohn dissipation-loss reference (`lumped-q-001`,
//! `IL₀ = 4.343·Σg/(Q_u·FBW)`) — independent corroboration of the sum rule's
//! `Σg`/`FBW` skeleton (the numeric prefactor differs: Cohn's loss vs delay).
//!
//! **Non-circularity.** The reference `τ_closed` is computed only from `Σg`
//! (read straight from [`yee_filter::FilterProject::prototype`]'s `g`-values),
//! `FBW`, and `ω0` — it never touches the `S21` phase. The measurement
//! `τ_meas = group_delay(S21, f)[mid]` comes from the numerical derivative of
//! the coupling-matrix phase. Neither is defined in terms of the other, so
//! agreement validates the group-delay computation, not a tautology.
//!
//! ## Tolerance
//!
//! The closed form is the *prototype* midband value; the measured value is a
//! finite-difference derivative of the realized Chebyshev phase, which has a
//! small equiripple-induced curvature near band-centre and a `O(Δω²)`
//! truncation error. A dense sweep keeps both small (measured agreement is
//! ~0.0%); the gate allows `≤ 5 %` relative error — tight enough to catch a
//! wrong constant (the `2×` form is 100 % off, a sign error 200 %) or a
//! phase-unwrap bug, while honest about the prototype-vs-realized +
//! discretization gap. **Do NOT weaken this anchor tolerance to force green —
//! the closed-form agreement IS the validation (ADR-0173 escape hatch).**

use std::f64::consts::PI;
use yee_filter::{
    Approximation, FilterProject, FilterSpec, Response, SpecMask, coupling_matrix_s_params,
    group_delay, synthesize,
};

/// Fractional bandwidth of the gate fixtures (reused by the closed form).
const FBW: f64 = 0.10;
/// Centre frequency of the gate fixtures, Hz.
const F0_HZ: f64 = 2.0e9;
/// Relative tolerance on `τ_meas(f0)` vs the closed-form `τ_closed(f0)`. Do NOT
/// weaken — see the module docs (this anchor IS the validation).
const ANCHOR_REL_TOL: f64 = 0.05;
/// Number of sweep points across the in-band span (dense → accurate derivative
/// at f0 and a meaningful symmetry check).
const N_SWEEP: usize = 801;

/// A Chebyshev 0.5 dB band-pass spec of the given order (f0 = 2 GHz, FBW = 10 %,
/// Z0 = 50 Ω) — the project's standard fixture shape.
fn cheb_spec(order: usize) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: F0_HZ,
        fbw: FBW,
        order: Some(order),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 40.0)],
        },
    }
}

/// A dense in-band sweep `[f0·(1−FBW), f0·(1+FBW)]`, `N_SWEEP` points, with a
/// sample landing **exactly** on `f0` (odd count + symmetric span → the centre
/// index is f0). Spanning ±FBW keeps the whole grid inside the passband so the
/// causality (`τ > 0`) and symmetry checks are over genuinely in-band samples.
fn band_sweep(f0: f64) -> (Vec<f64>, usize) {
    let lo = f0 * (1.0 - FBW);
    let hi = f0 * (1.0 + FBW);
    let freqs: Vec<f64> = (0..N_SWEEP)
        .map(|i| lo + (hi - lo) * (i as f64) / ((N_SWEEP - 1) as f64))
        .collect();
    let mid = N_SWEEP / 2; // odd N_SWEEP → exact centre on f0.
    (freqs, mid)
}

/// `Σ_{k=1}^{N} g_k` from the synthesized prototype (`g[1..=N]`); never touches
/// the S-parameters — this is the independent reference quantity.
fn sum_g(proj: &FilterProject) -> f64 {
    let n = proj.prototype.order();
    (1..=n).map(|k| proj.prototype.g[k]).sum()
}

/// The closed-form midband group delay `Σg / (FBW·ω0)` (seconds) — the
/// first-principles-confirmed sum rule (`τ_LP(0)=Σg/2` × Jacobian `2/(FBW·ω0)`).
fn tau_closed(sum_g: f64, fbw: f64, f0_hz: f64) -> f64 {
    let omega0 = 2.0 * PI * f0_hz;
    sum_g / (fbw * omega0)
}

#[test]
fn group_delay_001_midband_anchor_n3_n5() {
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        assert_eq!(proj.prototype.order(), order);

        let (freqs, mid) = band_sweep(spec.f0_hz);
        assert!(
            (freqs[mid] - spec.f0_hz).abs() < 1.0,
            "the centre sample must land on f0 (got {} Hz)",
            freqs[mid]
        );

        let s = coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw);
        let s21: Vec<num_complex::Complex64> = s.iter().map(|&(_s11, s21)| s21).collect();
        let tau = group_delay(&s21, &freqs);
        assert_eq!(tau.len(), freqs.len());

        // --- independent published reference: 2·Σg/(FBW·ω0) -------------------
        let sg = sum_g(&proj);
        let tc = tau_closed(sg, spec.fbw, spec.f0_hz);
        let tm = tau[mid];
        let rel_err = (tm - tc).abs() / tc;
        println!(
            "group-delay-001: N={order} Σg={sg:.4} FBW={FBW} f0={:.3}GHz  \
             τ_meas(f0)={:.4} ns  τ_closed={:.4} ns  rel_err={:.2}% (tol {:.0}%)",
            spec.f0_hz / 1e9,
            tm * 1e9,
            tc * 1e9,
            rel_err * 100.0,
            ANCHOR_REL_TOL * 100.0,
        );
        assert!(
            tc > 0.0,
            "N={order}: closed-form τ must be positive (Σg={sg}, FBW={FBW})"
        );
        assert!(
            rel_err <= ANCHOR_REL_TOL,
            "N={order}: midband τ_meas={:.4} ns disagrees with the closed form \
             Σg/(FBW·ω0)={:.4} ns by {:.2}% > {:.0}% — a wrong sum-rule constant, a \
             sign error, or a phase-unwrap bug. Σg={sg:.4}, FBW={FBW}, f0={:.3} GHz. \
             Do NOT weaken ANCHOR_REL_TOL.",
            tm * 1e9,
            tc * 1e9,
            rel_err * 100.0,
            ANCHOR_REL_TOL * 100.0,
            spec.f0_hz / 1e9,
        );
    }
}

#[test]
fn group_delay_001_causal_in_band() {
    // Causality: a passive filter delays the signal, so in-band group delay is
    // positive at every sample. (A negative in-band τ would signal a phase-sign
    // or unwrap bug.)
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        let (freqs, _mid) = band_sweep(spec.f0_hz);
        let s = coupling_matrix_s_params(&proj.coupling, &freqs, spec.f0_hz, spec.fbw);
        let s21: Vec<num_complex::Complex64> = s.iter().map(|&(_s11, s21)| s21).collect();
        let tau = group_delay(&s21, &freqs);

        let min_tau = tau.iter().copied().fold(f64::INFINITY, f64::min);
        println!(
            "group-delay-001: N={order} min in-band τ = {:.4} ns (must be > 0)",
            min_tau * 1e9
        );
        for (k, (&t, &f)) in tau.iter().zip(freqs.iter()).enumerate() {
            assert!(
                t > 0.0,
                "N={order} sample {k} f={:.4} GHz: in-band τ_g = {:.4} ns ≤ 0 — a passive \
                 filter must have positive in-band group delay (phase-sign / unwrap bug?)",
                f / 1e9,
                t * 1e9,
            );
        }
    }
}

#[test]
fn group_delay_001_symmetric_about_f0() {
    // Symmetry: a synchronous (zero-diagonal coupling matrix) symmetric band-pass
    // has a group delay symmetric in the lowpass variable Ω, hence approximately
    // symmetric about f0 in f. The band-pass map `Ω = (1/FBW)(ω/ω0 − ω0/ω)` is
    // ODD in Ω but maps equal `±δf` to slightly UNequal `±Ω`, so the f-symmetry
    // residual grows with sweep width (verified: it shrinks monotonically to ~0
    // as the window narrows — measured worst dev ≈ 5–6 % at ±25 %·FBW down to
    // ≈ 1.9 % at ±10 %·FBW). The check therefore uses a NARROW symmetric window
    // (±0.10·FBW about f0) where the f↔Ω asymmetry is small; over it the mirror
    // pairs agree to ~1.9 %. (A genuine asymmetry — e.g. a wrong/sign-flipped
    // phase path — would NOT shrink with the window and would trip this.)
    const SYM_REL_TOL: f64 = 0.025;
    /// Half-width of the symmetry window, as a fraction of FBW.
    const SYM_HALF_FBW: f64 = 0.10;
    for order in [3usize, 5] {
        let spec = cheb_spec(order);
        let proj = synthesize(&spec);
        let f0 = spec.f0_hz;
        let half = SYM_HALF_FBW * spec.fbw;
        let lo = f0 * (1.0 - half);
        let hi = f0 * (1.0 + half);
        let freqs: Vec<f64> = (0..N_SWEEP)
            .map(|i| lo + (hi - lo) * (i as f64) / ((N_SWEEP - 1) as f64))
            .collect();
        let mid = N_SWEEP / 2; // odd N_SWEEP → exact centre on f0.

        let s = coupling_matrix_s_params(&proj.coupling, &freqs, f0, spec.fbw);
        let s21: Vec<num_complex::Complex64> = s.iter().map(|&(_s11, s21)| s21).collect();
        let tau = group_delay(&s21, &freqs);

        let mut worst_rel = 0.0_f64;
        let mut worst_at = 0usize;
        // Mirror index of `mid + d` is `mid − d`. Skip d such that a mirror is an
        // endpoint (k == 0 or k == n−1), whose one-sided derivative is less
        // accurate than the central interior ones.
        for d in 1..mid {
            let hi_i = mid + d;
            let lo_i = mid - d;
            if hi_i >= freqs.len() - 1 || lo_i == 0 {
                continue;
            }
            let a = tau[hi_i];
            let b = tau[lo_i];
            let rel = (a - b).abs() / b.abs().max(1e-30);
            if rel > worst_rel {
                worst_rel = rel;
                worst_at = d;
            }
        }
        println!(
            "group-delay-001: N={order} worst τ symmetry rel-dev = {:.3}% at ±{} samples \
             from f0 over a ±{:.0}%·FBW window (tol {:.1}%)",
            worst_rel * 100.0,
            worst_at,
            SYM_HALF_FBW * 100.0,
            SYM_REL_TOL * 100.0,
        );
        assert!(
            worst_rel <= SYM_REL_TOL,
            "N={order}: group delay is not symmetric about f0 (worst mirror-pair dev \
             {:.2}% > {:.1}%) over a ±{:.0}%·FBW window — a synchronous symmetric BPF \
             must have τ(f0+δ)≈τ(f0−δ)",
            worst_rel * 100.0,
            SYM_REL_TOL * 100.0,
            SYM_HALF_FBW * 100.0,
        );
    }
}

#[test]
fn group_delay_001_degenerate_inputs() {
    // < 2 samples or a length mismatch → zero-filled, no panic (a derivative is
    // undefined with fewer than two points).
    let one = [num_complex::Complex64::new(1.0, 0.0)];
    assert_eq!(group_delay(&one, &[2.0e9]), vec![0.0]);
    assert_eq!(group_delay(&[], &[]), Vec::<f64>::new());
    // Length mismatch is defensive: zero-filled to s21.len().
    let two = [
        num_complex::Complex64::new(1.0, 0.0),
        num_complex::Complex64::new(0.0, 1.0),
    ];
    assert_eq!(group_delay(&two, &[2.0e9]), vec![0.0, 0.0]);

    // A linear-phase S21 (φ = −τ0·ω, here a pure delay of τ0) must recover τ0 at
    // every interior sample — the simplest non-circular sanity check that the
    // sign and the dω scaling are right, independent of any filter.
    let tau0 = 1.0e-9; // 1 ns
    let freqs: Vec<f64> = (0..11).map(|i| 1.0e9 + 1.0e8 * i as f64).collect();
    let s21: Vec<num_complex::Complex64> = freqs
        .iter()
        .map(|&f| {
            let phi = -tau0 * 2.0 * PI * f;
            num_complex::Complex64::from_polar(1.0, phi)
        })
        .collect();
    let tau = group_delay(&s21, &freqs);
    for (k, &t) in tau.iter().enumerate() {
        assert!(
            (t - tau0).abs() < 1e-12,
            "linear-phase pure delay: sample {k} τ={t:.4e} should equal τ0={tau0:.4e} \
             (sign/scale check)"
        );
    }
}
