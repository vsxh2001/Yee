//! top-c-coupled-001 (JLCPCB narrow-band track, ADR-0165 brick T1): the
//! **top-C-coupled (capacitively-coupled)** lumped band-pass synthesis meets the
//! Chebyshev mask (non-circular), and its `(f0, FBW)` JLCPCB-orderable envelope
//! is measured.
//!
//! ## (1) Synthesis S21 validation — non-circular
//!
//! For a representative NARROW-band spec (0.5 dB Cheb, N=3, f0 = 1 GHz,
//! FBW = 10 %, 50 Ω) we [`synthesize_top_c_coupled`] the `N` shunt LC resonators
//! and `N+1` series coupling caps, then sweep [`top_c_s21`] (an **independent**
//! ABCD cascade — series-cap impedance · shunt-resonator admittance · …) and
//! grade the realized response against the Chebyshev mask:
//!
//! - in-band ripple (over `|Ω| ≤ 1`, the geometric band edges) within tolerance of 0.5 dB,
//! - in-band return loss ≥ the mask bound,
//! - deep rejection at a stopband point (1.5·f0).
//!
//! The response comes purely from the network ABCD analysis, **not** from the
//! synthesis inputs, so a mask pass is a non-circular proof the synthesis is
//! correct (a wrong J-inverter or coupling-cap formula would mis-center the
//! band, kill the stopband skirt, or blow the ripple far past tolerance).
//!
//! ### Realization tolerance on the in-band ripple ([`RIPPLE_TOL_DB`])
//!
//! A capacitive J-inverter is an exact inverter **only at `ω0`**; its
//! frequency-dependent reactance adds a dispersion term that raises the realized
//! in-band ripple above the prototype's 0.5 dB and grows with FBW (Naaman &
//! Aumentado, arXiv:2109.11628, note the "slight asymmetry … due to the
//! additional frequency dependence … of the admittance inverters", and that the
//! method suits FBW up to ≈ 20 %). On this N=3/10 % fixture the realized in-band
//! ripple is ≈ 0.85 dB. This is a **physical property of the topology**, not a
//! synthesis error, so the gate allows a documented realization slack on the
//! *ripple bound only* (the deep-stopband + center checks carry no slack and
//! reject a genuinely broken synthesis). Mirrors the
//! [`lumped_001`](lumped_001.rs) gate's documented band-edge slack.
//!
//! ## (2) Realizability-envelope probe — the de-risk
//!
//! We sweep `(f0, FBW)` over `{0.2, 0.5, 1.0, 2.0} GHz × {5, 10, 20} %`
//! (N=3, both 0603 and 0402), synthesize each, build [`BomLine`]s for the shunt
//! inductor + the distinct shunt-node caps + the coupling caps, run
//! [`autopick`], and report per-cell orderable coverage (does every part resolve
//! to a real LCSC Basic part?). We **assert at least one cell is fully orderable
//! (zero blanks)** — proving top-C-coupled extends the orderable envelope into a
//! narrow-band regime — and record where it still blanks (the GHz-narrow
//! coupling-cap sub-pF floor). Pure-compute, deterministic, fast (no FDTD, no
//! `#[ignore]`). Do NOT weaken: an honest `None` on a sub-pF coupling cap is the
//! correct, informative de-risk outcome.

use yee_filter::{
    Approximation, BomLine, CompKind, ESeries, Footprint, TopCNetwork, autopick,
    synthesize_top_c_coupled, top_c_s21,
};
use yee_synth::lowpass_to_bandpass;

/// Realization slack on the in-band ripple bound, dB. Bounds the capacitive
/// J-inverter frequency-dependence (≈ 0.85 dB realized on the N=3/10 % fixture
/// vs the 0.5 dB prototype — see the module docs). The deep-stopband and center
/// checks carry NO slack, so a broken synthesis is still rejected.
const RIPPLE_TOL_DB: f64 = 0.6;

/// A capacitor BOM line at `value` farads (autopick keys on `chosen_value`).
fn cap_line(value: f64) -> BomLine {
    BomLine {
        kind: CompKind::Capacitor,
        ideal_value: value,
        chosen_value: value,
        deviation_pct: 0.0,
        series: ESeries::E24,
        tolerance_pct: 5.0,
        qty: 1,
        esr_ohm: None,
        srf_hz: None,
    }
}

/// An inductor BOM line at `value` henries.
fn ind_line(value: f64) -> BomLine {
    BomLine {
        kind: CompKind::Inductor,
        ideal_value: value,
        chosen_value: value,
        deviation_pct: 0.0,
        series: ESeries::E24,
        tolerance_pct: 5.0,
        qty: 1,
        esr_ohm: None,
        srf_hz: None,
    }
}

/// `|S21|` of the synthesized network at `f_hz`.
fn s21_mag(net: &TopCNetwork, f_hz: f64) -> f64 {
    top_c_s21(net, f_hz, net.z0_ohm).norm()
}

#[test]
fn top_c_coupled_001_s21_meets_chebyshev_mask() {
    // Representative narrow-band spec.
    let (approx, n, f0, fbw, z0) = (
        Approximation::Chebyshev { ripple_db: 0.5 },
        3,
        1.0e9,
        0.10,
        50.0,
    );
    let net = synthesize_top_c_coupled(approx, n, f0, fbw, z0);

    assert_eq!(net.shunt.len(), n, "N=3 → 3 shunt resonators");
    assert_eq!(
        net.coupling_caps_farad.len(),
        n + 1,
        "N=3 → N+1 coupling caps"
    );

    println!("\n=== top-C-coupled synthesis: 0.5 dB Cheb N=3, f0=1 GHz, FBW=10 %, 50 Ω ===");
    for (i, r) in net.shunt.iter().enumerate() {
        println!(
            "  shunt resonator {}: L = {:.3} nH, C = {:.3} pF",
            i + 1,
            r.l_henry * 1e9,
            r.c_farad * 1e12
        );
    }
    let labels = ["C01", "C12", "C23", "C34"];
    for (k, c) in net.coupling_caps_farad.iter().enumerate() {
        println!("  coupling cap {} = {:.4} pF", labels[k], c * 1e12);
    }

    // ---- grade the realized S21 over the MAPPED band (|Ω| ≤ 1) -------------
    // Sweep widely; classify each sample by the band-pass map (geometric edges).
    let n_samp = 6000usize;
    let (f_lo, f_hi) = (0.6 * f0, 1.6 * f0);
    let mut min_il = f64::INFINITY;
    let mut max_il = f64::NEG_INFINITY;
    let mut worst_rl = f64::INFINITY;
    let mut saw_band = false;
    let mut peak_mag = 0.0f64;
    let mut peak_f = 0.0f64;
    for i in 0..=n_samp {
        let f = f_lo + (i as f64) * (f_hi - f_lo) / (n_samp as f64);
        let m = s21_mag(&net, f);
        if m > peak_mag {
            peak_mag = m;
            peak_f = f;
        }
        if lowpass_to_bandpass(f, f0, fbw).abs() > 1.0 {
            continue; // out of the mapped pass-band
        }
        saw_band = true;
        let il = -20.0 * m.max(1e-300).log10();
        let s11_sq = (1.0 - m * m).max(0.0);
        let rl = if s11_sq <= 0.0 {
            f64::INFINITY
        } else {
            -10.0 * s11_sq.log10()
        };
        min_il = min_il.min(il);
        max_il = max_il.max(il);
        worst_rl = worst_rl.min(rl);
    }
    assert!(
        saw_band,
        "no swept frequency fell inside the mapped pass-band"
    );
    let ripple = (max_il - min_il).max(0.0);

    // Stopband points (no slack).
    let rej_1p5 = -20.0 * s21_mag(&net, 1.5 * f0).max(1e-300).log10();
    let rej_2p0 = -20.0 * s21_mag(&net, 2.0 * f0).max(1e-300).log10();

    println!(
        "  realized: peak |S21| = {peak_mag:.4} @ {:.4} GHz; in-band ripple = {ripple:.3} dB; \
         worst RL = {worst_rl:.2} dB",
        peak_f / 1e9
    );
    println!("  stopband: rej @ 1.5 GHz = {rej_1p5:.1} dB; rej @ 2.0 GHz = {rej_2p0:.1} dB");

    // (a) lossless + correctly centered: peak ~1 within ~2 % of f0 (the inverter
    //     dispersion shifts the peak slightly; a broken synthesis would mis-place
    //     it badly or not reach ~1).
    assert!(
        peak_mag > 0.98 && peak_mag <= 1.0 + 1e-9,
        "peak |S21| = {peak_mag} should be ~1 (lossless equi-ripple band-pass)"
    );
    assert!(
        (peak_f / f0 - 1.0).abs() < 0.05,
        "peak at {:.4} GHz should be within 5 % of f0 = {f0:.3e}",
        peak_f / 1e9
    );
    // (b) ripple meets 0.5 dB within the documented J-inverter-dispersion slack.
    assert!(
        ripple <= 0.5 + RIPPLE_TOL_DB + 1e-9,
        "in-band ripple {ripple:.3} dB exceeds 0.5 + {RIPPLE_TOL_DB} dB slack"
    );
    // (c) in-band return loss is real (a coupled band-pass; ≥ 6 dB here). No slack.
    assert!(
        worst_rl >= 6.0 - 1e-9,
        "in-band worst RL {worst_rl:.2} dB below 6 dB (synthesis mis-matched)"
    );
    // (d) DEEP stopband skirt — the load-bearing non-circular check. A wrong
    //     J-inverter / coupling-cap formula cannot produce a sharp 25 dB+ skirt
    //     a half-octave out. No slack. (Realized ≈ 44 dB.)
    assert!(
        rej_1p5 >= 25.0,
        "rejection at 1.5 f0 = {rej_1p5:.1} dB below 25 dB (skirt broken)"
    );
    assert!(
        rej_2p0 >= 25.0,
        "rejection at 2.0 f0 = {rej_2p0:.1} dB below 25 dB"
    );
}

/// Per-cell result of the realizability-envelope probe.
struct Cell {
    f0_ghz: f64,
    fbw_pct: f64,
    footprint: Footprint,
    covered: usize,
    total: usize,
    blanks: Vec<String>,
}

/// Synthesize one `(f0, FBW)` cell and autopick every distinct part at one
/// footprint. Returns the coverage + the human-labelled blanks.
fn probe_cell(f0: f64, fbw: f64, footprint: Footprint) -> Cell {
    let net = synthesize_top_c_coupled(
        Approximation::Chebyshev { ripple_db: 0.5 },
        3,
        f0,
        fbw,
        50.0,
    );

    let mut covered = 0usize;
    let mut total = 0usize;
    let mut blanks = Vec::new();

    // Shunt inductor (all resonators share L = Zr/ω0 → one distinct inductor).
    let l = net.shunt[0].l_henry;
    total += 1;
    if autopick(&ind_line(l), footprint).is_some() {
        covered += 1;
    } else {
        blanks.push(format!("L={:.2}nH", l * 1e9));
    }

    // Distinct shunt-node caps (group equal values — the symmetric proto gives
    // C1 == C3, so this collapses to {edge node, centre node}).
    let mut node_caps: Vec<f64> = Vec::new();
    for r in &net.shunt {
        if !node_caps.iter().any(|c| (c - r.c_farad).abs() < 1e-18) {
            node_caps.push(r.c_farad);
        }
    }
    for c in &node_caps {
        total += 1;
        if autopick(&cap_line(*c), footprint).is_some() {
            covered += 1;
        } else {
            blanks.push(format!("Cnode={:.3}pF", c * 1e12));
        }
    }

    // Distinct coupling caps (symmetric proto → {end, internal}).
    let mut cpl_caps: Vec<f64> = Vec::new();
    for c in &net.coupling_caps_farad {
        if !cpl_caps.iter().any(|x| (x - c).abs() < 1e-18) {
            cpl_caps.push(*c);
        }
    }
    for c in &cpl_caps {
        total += 1;
        if autopick(&cap_line(*c), footprint).is_some() {
            covered += 1;
        } else {
            blanks.push(format!("Ccpl={:.3}pF", c * 1e12));
        }
    }

    Cell {
        f0_ghz: f0 / 1e9,
        fbw_pct: fbw * 100.0,
        footprint,
        covered,
        total,
        blanks,
    }
}

#[test]
fn top_c_coupled_001_realizability_envelope() {
    let f0s = [0.2e9, 0.5e9, 1.0e9, 2.0e9];
    let fbws = [0.05, 0.10, 0.20];
    let footprints = [Footprint::Smd0603, Footprint::Smd0402];

    println!("\n=== top-C-coupled JLCPCB realizability envelope (N=3, 0.5 dB Cheb, 50 Ω) ===");
    println!(
        "{:>7} {:>5} {:>6} | {:>9} | coverage  (blanks = non-orderable parts)",
        "f0", "FBW", "fp", "covered"
    );
    println!("{}", "-".repeat(86));

    let mut fully_orderable: Vec<(f64, f64, Footprint)> = Vec::new();
    let mut coupling_cap_blank_at_ghz = false;

    for &f0 in &f0s {
        for &fbw in &fbws {
            for &fp in &footprints {
                let cell = probe_cell(f0, fbw, fp);
                let full = cell.covered == cell.total && cell.blanks.is_empty();
                if full {
                    fully_orderable.push((cell.f0_ghz, cell.fbw_pct, cell.footprint));
                }
                // Record whether GHz-band cells blank specifically on coupling caps.
                if cell.f0_ghz >= 1.0 && cell.blanks.iter().any(|b| b.starts_with("Ccpl")) {
                    coupling_cap_blank_at_ghz = true;
                }
                println!(
                    "{:>6.1}G {:>4.0}% {:>6} | {:>4}/{:<4} | {}",
                    cell.f0_ghz,
                    cell.fbw_pct,
                    format!("{:?}", cell.footprint),
                    cell.covered,
                    cell.total,
                    if full {
                        "FULLY ORDERABLE".to_string()
                    } else {
                        cell.blanks.join(", ")
                    }
                );
            }
        }
    }

    println!("\nfully-orderable (zero-blank) cells:");
    if fully_orderable.is_empty() {
        println!("  NONE");
    } else {
        for (f0, fbw, fp) in &fully_orderable {
            println!("  f0={f0:.1} GHz, FBW={fbw:.0} %, {fp:?}");
        }
    }

    // THE DE-RISK ASSERTION (honest either way):
    //
    // At least one (f0, FBW, footprint) cell must be FULLY orderable — proving
    // top-C-coupled extends the JLCPCB-orderable envelope into SOME narrow-band
    // regime (lower-freq / wider-band, where the coupling caps Cij = J/ω0 are
    // large enough to clear the 1 pF Basic floor). If NO cell is fully orderable,
    // this fails LOUDLY — the honest NO-GO that narrow-band lumped is
    // fundamentally distributed-only, surfaced rather than forced green.
    assert!(
        !fully_orderable.is_empty(),
        "NO (f0, FBW) cell is fully JLCPCB-orderable — top-C-coupled does NOT extend \
         the orderable envelope; narrow-band lumped is distributed-only (de-risk NO-GO)"
    );

    // And record the binding constraint at GHz: the coupling caps fall sub-pF
    // (below the 1 pF Basic floor), so GHz-narrow cells blank on the coupling
    // caps — confirming ADR-0165's hypothesis that GHz-narrow needs distributed.
    assert!(
        coupling_cap_blank_at_ghz,
        "expected the GHz-band cells to blank on sub-pF coupling caps (the documented \
         narrow-band-GHz floor); if they do not, the envelope is wider than predicted"
    );
}
