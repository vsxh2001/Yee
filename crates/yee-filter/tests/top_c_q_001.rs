//! top-c-q-001 (ADR-0170 brick T6): finite-Q top-C response vs Cohn's loss.
//!
//! The lossless [`top_c_s21`](yee_filter::top_c_s21) shows the *ideal* top-C-
//! coupled band-pass response (zero midband insertion loss). Real filters
//! degrade because the shunt-resonator L/C components have finite unloaded
//! quality factor `Q_u`. The finite-Q response
//! [`top_c_s21_lossy`](yee_filter::top_c_s21_lossy) adds a per-resonator
//! dissipation conductance (`G = ω₀·C/Q_u`); this gate validates the **midband
//! insertion loss** it produces against the same published closed form the
//! lumped-ladder gate (`lumped-q-001`, ADR-0160) uses:
//!
//! > Cohn's narrowband dissipation-loss formula (S. B. Cohn, "Dissipation Loss
//! > in Multiple-Coupled-Resonator Filters," Proc. IRE, 1959; Hong & Lancaster,
//! > §3.2):
//! >
//! >     IL₀ (dB) ≈ 4.343 · (Σ_{k=1}^{n} g_k) / (Q_u · FBW)
//!
//! where `g_k` are the low-pass prototype values and `FBW` the fractional
//! bandwidth. Cohn's formula is **prototype-based / topology-independent** (it
//! depends only on the g-values, `Q_u`, and FBW), so it is the correct
//! independent reference for the top-C topology too — exactly as ADR-0170 §2
//! notes.
//!
//! ## Non-circularity
//!
//! The two sides are independently derived. The **reference** `IL_cohn` is the
//! published closed form, computed only from `Σg` (read straight from
//! [`yee_synth::prototype`]), `Q_u`, and `FBW` — it never touches the network's
//! S-parameters. The **measurement** `IL_meas = −20·log10|S21(f₀)|` comes from
//! the lossy ABCD cascade over the realized shunt-resonator + coupling-cap
//! values. Neither is defined in terms of the other; agreement is evidence the
//! lossy shunt-resonator model is physically right, not a tautology.
//!
//! ## Loss model (bare resonating cap)
//!
//! [`top_c_s21_lossy`] keys each shunt node's loss conductance to the **bare
//! resonating** capacitance `C_bare = 1/(Zr·ω₀)`, not the synthesis-reduced
//! physical node cap — because a top-C node *resonates* at `ω₀` with `C_bare`
//! (the negative-leg-absorbed coupling caps add back to it), so the resonator's
//! stored energy lives in `C_bare` and its unloaded Q is `ω₀·C_bare/G`. This is
//! the faithful physical mirror of the lumped ladder's resonator loss (where the
//! branch cap *is* the full resonating cap). Keying the loss to the reduced node
//! cap instead would model an effective Q higher than `Q_u` and undershoot Cohn
//! by ~28 % at FBW = 10 %; the bare-cap model matches Cohn to ≈ 2.2 % — see the
//! [`top_c_s21_lossy`] doc comment.
//!
//! ## Tolerance
//!
//! Cohn's formula is a **narrowband** approximation (first order in `1/Q_u`,
//! valid for small `FBW`). At `FBW = 0.10` we are at the edge of "narrowband",
//! so the gate allows `≤ 15 %` relative error — the SAME tolerance/spirit as
//! `lumped-q-001`: tight enough to reject a wrong loss model (a factor-of-two
//! element-model bug is ~100 % off) but honest about the approximation. Do NOT
//! weaken this to force green; if the model and Cohn disagree by more, that is a
//! real finding (try a narrower-FBW fixture first — Cohn is narrowband — before
//! suspecting the loss model).

use yee_filter::{Approximation, synthesize_top_c_coupled, top_c_s21, top_c_s21_lossy};
use yee_synth::prototype;

/// Fractional bandwidth of the gate fixture (reused by the Cohn reference).
const FBW: f64 = 0.10;
/// Centre frequency of the gate fixture, Hz.
const F0_HZ: f64 = 1.0e9;
/// System reference impedance, Ω (also the chosen resonator `Zr`).
const Z0_OHM: f64 = 50.0;
/// Filter order (N=3 shunt resonators; same family as lumped-q-001, Σg ≈ 4.289).
const ORDER: usize = 3;
/// 4.343 = 10/ln(10) = 20·log10(e)/2 = the dB-per-neper constant in Cohn's form.
const DB_PER_NEPER: f64 = 4.343;
/// Narrowband-approximation relative tolerance on IL_meas vs IL_cohn (FBW=10 %
/// is at the edge of the narrowband regime). Do NOT weaken — see the module docs.
const COHN_REL_TOL: f64 = 0.15;

/// The 3-pole 0.5 dB Chebyshev approximation the top-C track synthesizes (same
/// family / FBW / f0 as the `top_c_coupled_001` gate; order 3 so `Σg ≈ 4.289`).
const APPROX: Approximation = Approximation::Chebyshev { ripple_db: 0.5 };

/// `IL_meas` in dB at `f_hz` for the lossy top-C network at unloaded `q_unloaded`.
fn insertion_loss_db(net: &yee_filter::TopCNetwork, f_hz: f64, q_unloaded: f64) -> f64 {
    let mag = top_c_s21_lossy(net, f_hz, Z0_OHM, q_unloaded).norm();
    -20.0 * mag.max(1e-300).log10()
}

#[test]
fn top_c_q_001() {
    let net = synthesize_top_c_coupled(APPROX, ORDER, F0_HZ, FBW, Z0_OHM);
    assert_eq!(net.shunt.len(), ORDER, "N=3 → 3 shunt resonators");
    assert_eq!(
        net.coupling_caps_farad.len(),
        ORDER + 1,
        "N=3 → N+1=4 coupling caps"
    );

    // --- independent published reference: Cohn's Σg --------------------------
    // Σg = Σ_{k=1}^{n} g_k, read straight from yee_synth (g[1..=n]); never from
    // the network's S-parameters.
    let proto = prototype(APPROX, ORDER);
    let sum_g: f64 = (1..=ORDER).map(|k| proto.g[k]).sum();
    assert!(
        (sum_g - 4.289).abs() < 5e-3,
        "Σg = {sum_g:.4}, expected ≈ 4.289 for the 3-pole 0.5 dB Cheb"
    );

    let q_u = 100.0;

    // --- independent measurement: lossy ABCD S21 at band-centre -------------
    let il_meas = insertion_loss_db(&net, F0_HZ, q_u);
    // --- Cohn's closed-form midband loss ------------------------------------
    let il_cohn = DB_PER_NEPER * sum_g / (q_u * FBW);

    let rel_err = (il_meas - il_cohn).abs() / il_cohn;
    println!(
        "top-c-q-001: Q_u={q_u} FBW={FBW} Σg={sum_g:.4}  IL_meas={il_meas:.4} dB  \
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

    // --- tripwire (a): lossless (Q=∞) bit-equals top_c_s21 over a sweep ------
    // top_c_s21 delegates to top_c_s21_lossy(.., ∞), so this is by construction;
    // assert it bit-identically over 0.5–1.5·f0 (50 pts) to pin the delegation.
    for i in 0..=50 {
        let f = 0.5 * F0_HZ + (i as f64) * (1.0 * F0_HZ / 50.0); // 0.5·f0 .. 1.5·f0
        let lossless = top_c_s21(&net, f, Z0_OHM);
        let lossy_inf = top_c_s21_lossy(&net, f, Z0_OHM, f64::INFINITY);
        assert_eq!(
            lossy_inf, lossless,
            "Q=∞ lossy S21 must be bit-identical to lossless top_c_s21 at f={f:e}"
        );
        // A non-positive Q means "no loss" too — same bit-identical result.
        let lossy_zero = top_c_s21_lossy(&net, f, Z0_OHM, 0.0);
        let lossy_neg = top_c_s21_lossy(&net, f, Z0_OHM, -10.0);
        assert_eq!(
            lossy_zero, lossless,
            "Q_u=0 must be the bit-identical lossless response at f={f:e}"
        );
        assert_eq!(
            lossy_neg, lossless,
            "Q_u<0 must be the bit-identical lossless response at f={f:e}"
        );
    }
    // Lossless midband IL₀ is ≈ 0 dB (sanity on the bit-identical lossless path).
    let il_lossless = -20.0 * top_c_s21(&net, F0_HZ, Z0_OHM).norm().max(1e-300).log10();
    println!("top-c-q-001: lossless IL₀ = {il_lossless:.6} dB (expect ≈ 0)");
    assert!(
        il_lossless.abs() <= 0.01,
        "lossless midband IL₀ = {il_lossless:.6} dB should be ≈ 0 (≤ 0.01)"
    );

    // --- tripwire (b): finite-Q loss is present and material ----------------
    assert!(
        il_meas > 0.5,
        "finite-Q IL_meas={il_meas:.4} dB should be > 0.5 dB (loss present + material)"
    );

    // --- tripwire (c): 1/Q_u scaling — IL(2·Q) ≈ 0.5·IL(Q) within ~10% ------
    let il_q200 = insertion_loss_db(&net, F0_HZ, 2.0 * q_u);
    let scaling_ratio = il_q200 / il_meas; // Cohn predicts ≈ 0.5
    println!(
        "top-c-q-001: IL(Q={q_u})={il_meas:.4} dB  IL(Q={})={il_q200:.4} dB  \
         ratio={scaling_ratio:.4} (expect ≈ 0.5)",
        2.0 * q_u
    );
    assert!(
        (scaling_ratio - 0.5).abs() <= 0.10 * 0.5,
        "1/Q scaling broken: IL(2Q)/IL(Q) = {scaling_ratio:.4}, expected 0.5 ± 10%"
    );
}
