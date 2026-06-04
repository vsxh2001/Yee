//! lumped-q-001 (ADR-0160): finite-Q lumped response vs Cohn's dissipation loss.
//!
//! The lossless [`ladder_s21`](yee_filter::ladder_s21) shows the *ideal*
//! Chebyshev response (zero midband insertion loss). Real filters degrade
//! because the L/C components have finite unloaded quality factor `Q_u`. The
//! finite-Q response [`ladder_s21_lossy`](yee_filter::ladder_s21_lossy) adds a
//! per-resonator dissipation term; this gate validates the **midband insertion
//! loss** it produces against a published closed form:
//!
//! > Cohn's narrowband dissipation-loss formula (S. B. Cohn, "Dissipation Loss
//! > in Multiple-Coupled-Resonator Filters," Proc. IRE, 1959; Hong & Lancaster,
//! > §3.2):
//! >
//! >     IL₀ (dB) ≈ 4.343 · (Σ_{k=1}^{n} g_k) / (Q_u · FBW)
//!
//! where `g_k` are the low-pass prototype values and `FBW` the fractional
//! bandwidth.
//!
//! ## Non-circularity
//!
//! The two sides are independently derived. The **reference** `IL_cohn` is the
//! published closed form, computed only from `Σg` (read straight from
//! [`yee_synth::prototype`]), `Q_u`, and `FBW` — it never touches the ladder's
//! S-parameters. The **measurement** `IL_meas = −20·log10|S21(f₀)|` comes from
//! the lossy ABCD cascade over the realized L/C element values. Neither is
//! defined in terms of the other; agreement is evidence the lossy element model
//! is physically right, not a tautology.
//!
//! ## Tolerance
//!
//! Cohn's formula is a **narrowband** approximation (first order in `1/Q_u`,
//! valid for small `FBW`). At `FBW = 0.10` we are at the edge of "narrowband",
//! so the gate allows `≤ 15 %` relative error — tight enough to reject a wrong
//! loss model (a factor-of-two element-model bug is ~100 % off) but honest about
//! the approximation. Do NOT weaken this to force green; if the model and Cohn
//! disagree by more, that is a real finding.

use yee_filter::{
    Approximation, FilterSpec, LcBranch, Response, SpecMask, ladder_s21, ladder_s21_lossy,
    synthesize, synthesize_lumped,
};
use yee_synth::prototype;

/// Fractional bandwidth of the gate fixture (reused by the Cohn reference).
const FBW: f64 = 0.10;
/// Centre frequency of the gate fixture, Hz.
const F0_HZ: f64 = 2.0e9;
/// 4.343 = 10/ln(10) = 20·log10(e)/2 = the dB-per-neper constant in Cohn's form.
const DB_PER_NEPER: f64 = 4.343;
/// Narrowband-approximation relative tolerance on IL_meas vs IL_cohn (FBW=10 %
/// is at the edge of the narrowband regime). Do NOT weaken — see the module docs.
const COHN_REL_TOL: f64 = 0.15;

/// The 3-pole 0.5 dB Chebyshev band-pass spec the lumped track synthesizes
/// (same family / FBW / f0 as the F2.x lumped path; order 3 so `Σg ≈ 4.289`).
fn fixture() -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: F0_HZ,
        fbw: FBW,
        order: Some(3),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 9.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    }
}

/// `IL_meas` in dB at `f_hz` for the lossy ladder at unloaded `q_unloaded`.
fn insertion_loss_db(ladder: &yee_filter::LumpedLadder, f_hz: f64, q_unloaded: f64) -> f64 {
    let mag = ladder_s21_lossy(ladder, f_hz, q_unloaded).norm();
    -20.0 * mag.max(1e-300).log10()
}

#[test]
fn lumped_q_001() {
    let spec = fixture();
    let proj = synthesize(&spec);
    let n = proj.prototype.order();
    assert_eq!(n, 3, "fixture is order N=3");

    // --- independent published reference: Cohn's Σg --------------------------
    // Σg = Σ_{k=1}^{n} g_k, read straight from yee_synth (g[1..=n]); never from
    // the ladder's S-parameters.
    let proto = prototype(spec.approximation, n);
    let sum_g: f64 = (1..=n).map(|k| proto.g[k]).sum();
    assert!(
        (sum_g - 4.289).abs() < 5e-3,
        "Σg = {sum_g:.4}, expected ≈ 4.289 for the 3-pole 0.5 dB Cheb"
    );

    let ladder = synthesize_lumped(&proj).expect("N=3 bandpass fixture should synthesize");
    assert_eq!(ladder.resonators.len(), 3, "N=3 → 3 LC resonators");

    let q_u = 100.0;

    // --- independent measurement: lossy ABCD S21 at band-centre -------------
    let il_meas = insertion_loss_db(&ladder, spec.f0_hz, q_u);
    // --- Cohn's closed-form midband loss ------------------------------------
    let il_cohn = DB_PER_NEPER * sum_g / (q_u * FBW);

    let rel_err = (il_meas - il_cohn).abs() / il_cohn;
    println!(
        "lumped-q-001: Q_u={q_u} FBW={FBW} Σg={sum_g:.4}  IL_meas={il_meas:.4} dB  \
         IL_cohn={il_cohn:.4} dB  rel_err={:.2}% (tol {:.0}%)",
        rel_err * 100.0,
        COHN_REL_TOL * 100.0
    );

    assert!(
        rel_err <= COHN_REL_TOL,
        "midband IL_meas={il_meas:.4} dB disagrees with Cohn IL_cohn={il_cohn:.4} dB by \
         {:.2}% > {:.0}% (narrowband tol). Σg={sum_g:.4}, Q_u={q_u}, FBW={FBW}.",
        rel_err * 100.0,
        COHN_REL_TOL * 100.0
    );

    // --- tripwire (a): lossless (Q=∞) is ~0 dB AND bit-equals ladder_s21 -----
    let s21_inf = ladder_s21_lossy(&ladder, spec.f0_hz, f64::INFINITY);
    let il_lossless = -20.0 * s21_inf.norm().max(1e-300).log10();
    println!("lumped-q-001: lossless IL₀ = {il_lossless:.6} dB (expect ≈ 0)");
    assert!(
        il_lossless.abs() <= 0.01,
        "lossless midband IL₀ = {il_lossless:.6} dB should be ≈ 0 (≤ 0.01)"
    );
    assert_eq!(
        s21_inf,
        ladder_s21(&ladder, spec.f0_hz),
        "Q=∞ lossy S21 must be bit-identical to lossless ladder_s21"
    );

    // --- tripwire (b): finite-Q loss is present and material ----------------
    assert!(
        il_meas > 0.5,
        "finite-Q IL_meas={il_meas:.4} dB should be > 0.5 dB (loss present + material)"
    );

    // --- tripwire (c): 1/Q_u scaling — IL(2·Q) ≈ 0.5·IL(Q) within ~10% ------
    let il_q200 = insertion_loss_db(&ladder, spec.f0_hz, 2.0 * q_u);
    let scaling_ratio = il_q200 / il_meas; // Cohn predicts ≈ 0.5
    println!(
        "lumped-q-001: IL(Q={q_u})={il_meas:.4} dB  IL(Q={})={il_q200:.4} dB  \
         ratio={scaling_ratio:.4} (expect ≈ 0.5)",
        2.0 * q_u
    );
    assert!(
        (scaling_ratio - 0.5).abs() <= 0.10 * 0.5,
        "1/Q scaling broken: IL(2Q)/IL(Q) = {scaling_ratio:.4}, expected 0.5 ± 10%"
    );

    // A loss term lives on every resonator branch (both Series and Shunt arms),
    // so the gate exercises both R = ω₀L/Q (series) and G = ω₀C/Q (shunt) paths.
    let has_shunt = ladder
        .resonators
        .iter()
        .any(|r| r.branch == LcBranch::Shunt);
    let has_series = ladder
        .resonators
        .iter()
        .any(|r| r.branch == LcBranch::Series);
    assert!(
        has_shunt && has_series,
        "N=3 shunt-first ladder should have both shunt and series resonators"
    );
}
