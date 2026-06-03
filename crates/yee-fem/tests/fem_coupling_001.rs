//! `fem-coupling-001` — FEM coupled-microstrip resonator-pair coupling
//! coefficient `k` gate (FEM-EM brick K1, ADR-0155).
//!
//! This is the **production PASS/FAIL gate** promoted from the GO de-risk probe
//! (`spike/fem-coupled-k-probe`, `933940f`). The extraction logic itself lives
//! in `src` ([`yee_fem::coupled_resonator_k`]); this test fixes the probe
//! geometry, runs the heavy FEM driven sweep, and asserts the two
//! measurement-driven tripwires that make it a real gate:
//!
//! 1. **Two cleanly resolvable peaks** — `peaks_resolvable` true AND the valley
//!    between them sits a real margin (≥ 6 dB) below the SHALLOWER peak. The
//!    probe's valley was 19.4 dB below, so 6 dB has honest headroom; this is a
//!    *re-smearing tripwire* — a future regression that over-couples the feeds
//!    or coarsens the mesh until the two peaks merge into one bump fails here.
//! 2. **`|k_fem − k_eps| / k_eps ≤ 0.30`** where `k_eps` is the even/odd
//!    ε_eff-split `(f_e²−f_o²)/(f_e²+f_o²)` from the Kirschning-Jansen
//!    `coupled_microstrip` even/odd ε_eff — the **like-for-like** analytic
//!    predictor of a resonant split (the FEM measures a resonant split, so this
//!    is the matching reference). The probe measured 17.2 % vs `k_eps`, so 30 %
//!    has headroom AND catches a real regression (a floored / smeared k). The
//!    impedance-k `k_imp = coupling_coefficient(...)` is **reported for
//!    traceability** but is NOT the gate ref: K2 (ADR-0155 Update) found `k_imp`
//!    and `k_eps` diverge at strong coupling (`k_imp/k_eps = 1.375` at
//!    S = 1.5 mm, outside the `[1.0,1.3]` comparability band the src encodes),
//!    so grading a resonant-split measurement vs `k_imp` is apples-to-oranges at
//!    tight gaps.
//!
//! ## Non-circular
//!
//! The reference is the Kirschning-Jansen quasi-static closed-form — the gate
//! grades `k_fem` against `k_eps`, derived from the even/odd ε_eff of
//! [`yee_layout::coupled_microstrip`] (the same closed-form family as
//! [`yee_layout::coupling_coefficient`]); the FEM is a full-wave Maxwell solve
//! on the meshed geometry — the analytic model does not set up the FEM mesh. The
//! two-peaks tripwire + the k-tolerance are real measured quantities (not a
//! tautology), so a regression that smears the peaks or floors `k` cannot pass.
//!
//! ## Tolerance honesty
//!
//! The 30 % tolerance and 6 dB valley margin are NOT to be weakened to force
//! green. The probe measured 17.2 % vs the gate ref `k_eps` (and 25.6 % vs the
//! impedance-k `k_imp`) and a 19.4 dB valley; there is a KNOWN systematic (both
//! FEM peaks pulled low by finite mesh dispersion + the feed-gap capacitive
//! load) that a finer mesh / a Qe-de-embedded peak would tighten (documented,
//! not hidden). The gate reports BOTH analytic references (`k_eps` graded,
//! `k_imp` for traceability).
//!
//! ## GATING — CRITICAL
//!
//! Multi-minute driven SWEEP (one per-ω sparse LU per frequency point; the
//! probe was ~280 s at ~63 k tets). `#[ignore]`'d so the debug
//! `cargo test --workspace` never runs it; run only in `--release`, boxed:
//!
//! ```text
//! YEE_BOX_DIR=$(pwd) YEE_BOX_MEM=14g YEE_BOX_CPUS=3 scripts/yee-box.sh bash -c '\
//!   cargo test -p yee-fem --release --test fem_coupling_001 \
//!   -- --ignored --nocapture'
//! ```

use yee_fem::{CoupledResonatorGeom, coupled_resonator_k};
use yee_layout::{coupled_microstrip, coupling_coefficient};

/// Sweep resolution across the split (matches the probe's 61 points → 10 MHz
/// step over 2.10–2.70 GHz, ~14 points across a ~140 MHz split).
const N_PTS: usize = 61;

/// k-tolerance vs the like-for-like ε_eff-split (`k_eps`). The probe measured
/// 17.2 % vs `k_eps`; 30 % has headroom and catches a real regression. Do NOT
/// weaken to force green (ADR-0155).
const K_TOL_FRAC: f64 = 0.30;

/// Required valley depth below the shallower peak (dB) for "two resolvable
/// peaks". The probe valley was 19.4 dB below; 6 dB is a real re-smearing
/// tripwire with headroom. Do NOT weaken (ADR-0155).
const VALLEY_MARGIN_DB: f64 = 6.0;

/// FEM coupled-microstrip resonator-pair coupling `k` gate.
///
/// Calls [`yee_fem::coupled_resonator_k`] at the ADR-0155 probe geometry
/// (W = 1 mm, S = 2 mm, h = 1 mm, ε_r = 4.4, f0 = 2.4 GHz), prints the full
/// auditable |S21|(f) spectrum + the k references + peak/valley structure, then
/// asserts the two real tripwires (resolvable peaks with a ≥ 6 dB valley margin;
/// `|k_fem − k_eps| / k_eps ≤ 0.30` vs the like-for-like ε_eff-split, with
/// `k_imp` reported for traceability).
#[test]
#[ignore = "K1 gate: multi-minute driven SWEEP (one per-ω sparse LU per point); run only in --release, boxed"]
fn fem_coupling_001() {
    let geom = CoupledResonatorGeom::probe();

    // Independent analytic reference (the gate's named `k_imp`), recomputed here
    // so the assertion does not rely on the `src` result's copy of it.
    let cm = coupled_microstrip(geom.trace_w, geom.gap_s, geom.sub_h, geom.eps_r);
    let k_imp = coupling_coefficient(&cm);

    eprintln!(
        "[fem-coupling-001] geometry: W={:.3}mm S={:.3}mm h={:.3}mm eps_r={} f0={:.2}GHz \
         L(λg/2)={:.3}mm  box=({:.2},?,{:.2})mm",
        geom.trace_w * 1e3,
        geom.gap_s * 1e3,
        geom.sub_h * 1e3,
        geom.eps_r,
        geom.f0_hz / 1e9,
        geom.resonator_length_m() * 1e3,
        geom.box_w * 1e3,
        geom.box_h * 1e3,
    );

    let t0 = std::time::Instant::now();
    let res = coupled_resonator_k(&geom, N_PTS).expect("coupled_resonator_k driven sweep must run");
    let wall = t0.elapsed().as_secs_f64();

    // ---- Auditable spectrum --------------------------------------------------
    eprintln!("\n{:>8}  {:>10}", "f(GHz)", "S21 dB");
    for &(f_ghz, d) in &res.s21_db {
        eprintln!("{f_ghz:>8.3}  {d:>10.2}");
    }

    let err_imp = (res.k_fem - k_imp).abs() / k_imp * 100.0;
    let err_eps = (res.k_fem - res.k_eps_ref).abs() / res.k_eps_ref * 100.0;
    let valley_margin = res.peak_lo_db.min(res.peak_hi_db) - res.valley_db;

    eprintln!(
        "\n==== FEM-COUPLING-001 GATE (K1, ADR-0155) ====\n\
         solve wall         : {:.1} s ({} pts)\n\
         peaks_resolvable   : {}\n\
         f_lo / f_hi (FEM)  : {:.4} / {:.4} GHz\n\
         peak dB (lo / hi)  : {:.2} / {:.2}\n\
         valley             : {:.2} dB  (margin below shallower peak = {:.2} dB; need ≥ {:.1})\n\
         k_fem  (split)     : {:.4}\n\
         k_imp  (analytic)  : {:.4}  -> err {:.1}%   <- coupled-line `coupling_coefficient` (traceability; diverges at strong coupling)\n\
         k_eps  (analytic)  : {:.4}  -> err {:.1}%   <- even/odd ε_eff resonant-split (GATE ref, like-for-like)\n\
         k tolerance        : {:.1}%\n\
         ==============================================",
        wall,
        N_PTS,
        res.peaks_resolvable,
        res.f_lo_hz / 1e9,
        res.f_hi_hz / 1e9,
        res.peak_lo_db,
        res.peak_hi_db,
        res.valley_db,
        valley_margin,
        VALLEY_MARGIN_DB,
        res.k_fem,
        k_imp,
        err_imp,
        res.k_eps_ref,
        err_eps,
        K_TOL_FRAC * 100.0,
    );

    // ---- Tripwire (1): two cleanly resolvable peaks --------------------------
    // A re-smearing regression (over-coupled feeds / coarsened mesh merging the
    // two peaks into one bump) fails here. The probe valley was 19.4 dB below
    // both peaks; the 6 dB bar has headroom and is a real measured threshold.
    assert!(
        res.peaks_resolvable,
        "fem-coupling-001: the two transmission peaks did NOT resolve — \
         `peaks_resolvable` is false (only one maximum, or no finite valley \
         between two peaks). The even/odd modes did not split into two clean \
         peaks (coupling too weak/strong, or both feeds not tapping). Full \
         spectrum printed above. Do NOT lower the resolution to force this."
    );
    assert!(
        valley_margin >= VALLEY_MARGIN_DB,
        "fem-coupling-001: re-smearing tripwire — the valley ({:.2} dB) sits only \
         {:.2} dB below the shallower peak ({:.2} dB), need ≥ {:.1} dB. The two \
         peaks have smeared together (over-coupling / coarse mesh across the \
         split). The probe valley was 19.4 dB below; do NOT weaken the margin. \
         Full spectrum printed above.",
        res.valley_db,
        valley_margin,
        res.peak_lo_db.min(res.peak_hi_db),
        VALLEY_MARGIN_DB,
    );

    // ---- Tripwire (2): k within 30 % of the like-for-like ε_eff-split --------
    // GRADED vs the resonant-split's natural analytic predictor k_eps (the FEM
    // measures a resonant split; k_eps is the same split formula on the analytic
    // even/odd ε_eff). The probe measured 17.2 % vs k_eps; 30 % has headroom AND
    // catches a real regression (a floored / smeared k). k_imp
    // (`coupling_coefficient`) is reported for traceability but is NOT the gate
    // ref: K2 (ADR-0155 Update) found it diverges from the resonant split at
    // strong coupling (k_imp/k_eps = 1.375 at S = 1.5 mm, outside the [1.0,1.3]
    // comparability band the src encodes), so grading a resonant split vs k_imp
    // is apples-to-oranges at tight gaps. Non-circular: KJ closed-form ε_eff vs
    // full-wave FEM. Do NOT weaken to force green (ADR-0155).
    assert!(
        err_eps <= K_TOL_FRAC * 100.0,
        "fem-coupling-001: k off — k_fem = {:.4} vs the like-for-like ε_eff-split \
         k_eps = {:.4} is {:.1}% off, exceeding the {:.1}% gate. (vs the \
         impedance-k k_imp = {:.4}, reported for traceability: {:.1}% — note k_imp \
         diverges from the resonant split at strong coupling, K2 finding.) The \
         probe measured 17.2% vs k_eps; a materially larger error means a \
         promotion bug (wrong geometry/sweep/peak-finder) — report the number, \
         do NOT lower the tolerance. Full spectrum printed above.",
        res.k_fem,
        res.k_eps_ref,
        err_eps,
        K_TOL_FRAC * 100.0,
        k_imp,
        err_imp,
    );
}
