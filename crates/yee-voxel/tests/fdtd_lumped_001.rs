//! Validation gate **fdtd-lumped-001** — full-wave FDTD EM simulation of a
//! synthesized lumped-LC band-pass filter board, cross-checked against the
//! analytic circuit response [`yee_filter::ladder_s21`] (Filter Phase F2.3,
//! ADR-0115).
//!
//! # What this gate proves
//!
//! It is the goal's named "EM simulation" of the lumped track — the full-wave
//! analogue of the distributed `fdtd-line-eeff-001` propagation gate. A
//! synthesized [`yee_filter::LumpedLadder`] (F2.0) is placed as SMD footprints
//! on a microstrip board (F2.2), voxelized (F1.1a), each L/C dropped on the
//! grid as a [`yee_fdtd::LumpedRlcPort`] (series branch → one series-RLC; shunt
//! branch → pure-L ‖ pure-C), driven / sensed across two ports, time-stepped,
//! and DFT'd into `|S21|(f)`. The gate confirms the FDTD response has the
//! correct *band-pass shape*: ≈ 0 dB in-band and a deep stopband — and that it
//! agrees with the analytic circuit `ladder_s21` of the *same* ladder.
//!
//! # Tolerance
//!
//! Deliberately **loose** (coarse-grid staircased FDTD with a one-way
//! circuit→field lumped coupling is approximate): the gate asserts
//!
//! - in-band (`f0 ≈ 2 GHz`) transmission within a few dB of 0 dB, and
//! - ≥ ~20 dB rejection at the 2.4 GHz stopband point,
//!
//! for *both* the FDTD `|S21|` and the analytic `ladder_s21` — i.e. the FDTD
//! reproduces the same pass/stop verdict as the circuit. This is a real
//! assertion (passband-vs-stopband behaviour), not a no-op; tightening to an
//! exact-match tolerance is a follow-on once the skeleton is green.
//!
//! # Why `#[ignore]`'d + CI-routed
//!
//! Two multi-minute FDTD solves (DUT + thru) → never runs in the default
//! `cargo test`. It runs in a dedicated CI `--release` job (`fdtd-lumped-gate`
//! in `.github/workflows/ci.yml`), mirroring the `fdtd-line-eeff` /
//! `mom-001` / GPU-nightly `--release -- --ignored` idiom. Per CLAUDE.md §4 the
//! feature must not merge until this gate is GREEN.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-voxel --release -- --ignored fdtd_lumped_001 --nocapture
//! ```

use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, ladder_s21, synthesize, synthesize_lumped,
};
use yee_layout::Substrate;
use yee_voxel::{LumpedSimConfig, simulate_lumped_board};

/// FR-4 substrate.
fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// Filter centre (Hz).
const F0_HZ: f64 = 2.0e9;
/// Stopband cross-check point (Hz).
const F_STOP_HZ: f64 = 2.4e9;

/// Insertion loss (dB, positive = loss) of a transmission magnitude.
fn il_db(mag: f64) -> f64 {
    -20.0 * mag.max(1e-300).log10()
}

/// Nearest swept sample to `f` in a `(freq, |S21|)` sweep.
fn s21_at(sweep: &[(f64, f64)], f: f64) -> f64 {
    sweep
        .iter()
        .min_by(|a, b| (a.0 - f).abs().total_cmp(&(b.0 - f).abs()))
        .map(|&(_, s)| s)
        .expect("non-empty sweep")
}

/// `fdtd-lumped-001`: FDTD `|S21|` of a Chebyshev N=5 lumped band-pass vs the
/// analytic `ladder_s21` — band-pass shape (in-band ≈ 0 dB, ≥ ~20 dB stopband).
#[test]
#[ignore = "slow: two multi-minute FDTD solves; fdtd-lumped-001 lumped-LC EM gate (F2.3, ADR-0115); run with --release --ignored"]
fn fdtd_lumped_001_matches_analytic_bandpass_shape() {
    // ------------------------------------------------------------------
    // Synthesize a Chebyshev N=5 lumped band-pass ladder at 2 GHz, 10 % FBW.
    // ------------------------------------------------------------------
    let spec = FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: F0_HZ,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 10.0,
            stopband: vec![(F_STOP_HZ, 30.0)],
        },
    };
    let project = synthesize(&spec);
    let ladder = synthesize_lumped(&project).expect("bandpass ladder synthesizes");
    assert_eq!(ladder.resonators.len(), 5, "Chebyshev N=5 → 5 resonators");

    let substrate = fr4();

    // ------------------------------------------------------------------
    // Analytic circuit cross-reference: ladder_s21 of the SAME ladder.
    // ------------------------------------------------------------------
    let ana_pass = ladder_s21(&ladder, F0_HZ).norm();
    let ana_stop = ladder_s21(&ladder, F_STOP_HZ).norm();
    let ana_pass_il = il_db(ana_pass);
    let ana_stop_rej = il_db(ana_stop);

    // ------------------------------------------------------------------
    // FDTD EM simulation of the placed board (DUT + thru normalization).
    // ------------------------------------------------------------------
    // CW per-frequency steady-state drive (F2.3-d, ADR-0128): a small frequency
    // set spanning the passband and stopband — the gate-check points (2.0 GHz
    // passband, 2.4 GHz stopband) plus a handful for the sweep shape. Each
    // frequency costs two full FDTD solves (DUT + thru), so this is NOT a fine
    // sweep.
    let cfg = LumpedSimConfig {
        cw_freqs_hz: vec![1.6e9, 1.8e9, F0_HZ, 2.2e9, F_STOP_HZ, 2.6e9],
        ..LumpedSimConfig::default()
    };
    let sweep = simulate_lumped_board(&ladder, &substrate, &cfg);
    assert_eq!(sweep.len(), cfg.cw_freqs_hz.len());

    let fdtd_pass = s21_at(&sweep, F0_HZ);
    let fdtd_stop = s21_at(&sweep, F_STOP_HZ);
    let fdtd_pass_il = il_db(fdtd_pass);
    let fdtd_stop_rej = il_db(fdtd_stop);

    eprintln!(
        "\nfdtd-lumped-001 lumped-LC EM gate (F2.3, ADR-0115)
  ladder:        Chebyshev N=5, f0 = {:.2} GHz, FBW = 10 %, Z0 = 50 Ω
  cross-check f: passband {:.2} GHz, stopband {:.2} GHz
  analytic:      |S21|(f0)  = {:.4}  ({:.2} dB IL),  |S21|(2.4G) = {:.4}  ({:.1} dB rej)
  FDTD:          |S21|(f0)  = {:.4}  ({:.2} dB IL),  |S21|(2.4G) = {:.4}  ({:.1} dB rej)",
        F0_HZ * 1e-9,
        F0_HZ * 1e-9,
        F_STOP_HZ * 1e-9,
        ana_pass,
        ana_pass_il,
        ana_stop,
        ana_stop_rej,
        fdtd_pass,
        fdtd_pass_il,
        fdtd_stop,
        fdtd_stop_rej,
    );
    eprintln!("  full sweep (GHz, |S21|, dB):");
    for &(f, s) in &sweep {
        eprintln!("    {:6.3}  {:.5}  {:7.2} dB", f * 1e-9, s, il_db(s));
    }

    // ------------------------------------------------------------------
    // Sanity-check the analytic reference IS a band-pass (it always is for
    // this synthesis — a guard against a broken cross-reference).
    // ------------------------------------------------------------------
    assert!(
        ana_pass_il < 3.0,
        "analytic ladder is not low-loss in-band: {ana_pass_il:.2} dB IL at f0"
    );
    assert!(
        ana_stop_rej >= 20.0,
        "analytic ladder does not reject 2.4 GHz: only {ana_stop_rej:.1} dB"
    );

    // ------------------------------------------------------------------
    // FDTD must reproduce the band-pass shape (loose, coarse-grid tolerances):
    //   - in-band within a few dB of 0 dB,
    //   - ≥ ~20 dB rejection at the stopband point.
    // ------------------------------------------------------------------
    assert!(
        fdtd_pass_il < 6.0,
        "fdtd-lumped-001 FAILED (passband): FDTD in-band IL = {fdtd_pass_il:.2} dB \
         (|S21| = {fdtd_pass:.4}); expected within a few dB of 0 dB"
    );
    assert!(
        fdtd_stop_rej >= 20.0,
        "fdtd-lumped-001 FAILED (stopband): FDTD rejection at {:.2} GHz = {fdtd_stop_rej:.1} dB \
         (|S21| = {fdtd_stop:.4}); expected ≥ 20 dB",
        F_STOP_HZ * 1e-9,
    );

    // ------------------------------------------------------------------
    // Cross-validation: FDTD and analytic agree on the pass/stop verdict
    // (the in-band transmission is far above the stopband transmission in
    // both, by a comparable margin).
    // ------------------------------------------------------------------
    assert!(
        fdtd_pass > fdtd_stop * 5.0,
        "fdtd-lumped-001 FAILED: FDTD passband/stopband contrast too small \
         (|S21|(f0) = {fdtd_pass:.4} vs |S21|(2.4G) = {fdtd_stop:.4})"
    );
}
