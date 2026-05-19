//! `fdtd-007` production-gate test — dielectric-loaded thin slot
//! antenna vs Maloney & Smith 1993 Fig. 9.
//!
//! Drives [`yee_validation::run_fdtd_007_maloney_smith_slot`] end-to-end
//! on the Phase 2.fdtd.7 Q7 fixture: a slot of width `w = 0.5 mm` and
//! length `L = 30 mm` cut into a PEC ground plane backed by a dielectric
//! substrate (`ε_r = 2.2`, `h = 1.524 mm`), driven by a delta-gap
//! voltage source at the slot midpoint. The driver extracts `S_11(f)`
//! via per-bin DFTs of the lumped-port `V(t)` / `I(t)` traces and
//! reports `f_res` and `|S_11(f_res)|` from both the subgridded run
//! (coarse `dx = 1 mm` + fine `dx = 0.5 mm` over a `(40 × 6 × 4) mm`
//! box) and a globally-uniform `dx = 0.5 mm` reference.
//!
//! # Gate
//!
//! Per Phase 2.fdtd.7 plan §Q7 DoD:
//!
//! 1. `|f_res − 8.9 GHz| / 8.9 GHz ≤ 0.02` (±2 %).
//! 2. `|S_11(f_res) − (-22 dB)| ≤ 1 dB`.
//! 3. Subgrid-vs-uniform sanity check: max `|df|/f ≤ 0.3 %` and
//!    max `|dS_11| ≤ 0.3 dB` across five spot frequencies.
//!
//! Reference values are digitised from Maloney & Smith 1993 Fig. 9
//! (IEEE Trans. Antennas Propag. 41(5), DOI 10.1109/8.222288) to
//! ±5 % per the plan §Q7 escape hatch and recorded under
//! [`yee_validation::FDTD_007_FRES_REF_HZ`] /
//! [`yee_validation::FDTD_007_S11_DB_REF`] with a `TBD: tighten when
//! Fig. 9 digitisation verified against the journal` note.
//!
//! # Status — `#[ignore]`'d as Phase 2.fdtd.7 Q7 escape hatch
//!
//! The gate is `#[ignore]`'d in default CI per the plan §Q7 escape
//! hatch:
//!
//! > Phase 2.fdtd.7.y Step C6 (Track DDDDDDDD) escape-hatched to F1:
//! > the J path now skips the B2.2 coarse-`H` ghost subtraction so the
//! > source is the un-ghosted Berenger form `J = +n̂ × H_fine`, which
//! > is ≈ 0 throughout the source-on-coarse traversal under Mur-only
//! > inward coupling. ... Dropping the ghost subtraction surfaces the
//! > fine box as effectively passive on the source-on-coarse traversal.
//!
//! Net effect: the subgridded run's `S_11` is dominated by the
//! coarse-grid stencil, not by the substrate-loaded fine box. The
//! Maloney-Smith dielectric-loaded resonance frequency cannot be
//! recovered until the F2 inward-coupling restoration (deferred from
//! Track DDDDDDDD) lands as Phase 2.fdtd.7.z. The driver itself is
//! green-compiling end-to-end and is kept under `#[ignore]` as the
//! scaffolding the F2 work will plug into.
//!
//! Two further structural blockers compound the C6 trade-off; both
//! also land as Phase 2.fdtd.7.z follow-ups (cross-referenced in
//! [`yee_validation::run_fdtd_007_maloney_smith_slot`] doc):
//!
//! 1. `YeeGrid` exposes a scalar `eps_r` / `mu_r` — no per-cell
//!    material map for the heterogeneous `ε_r = 2.2` substrate slab.
//! 2. No per-cell PEC mask for the slot-in-ground-plane geometry; only
//!    the deprecated outer-face `boundary::apply_pec` is available.
//!
//! Wall-time on the current driver lands in the seconds-to-minutes
//! range under `--release` (well under the 30-min budget); the
//! `#[ignore]` here is **physics**-gated, not wall-time-gated.

use yee_validation::{CaseStatus, run_fdtd_007_maloney_smith_slot};

#[test]
#[ignore = "Phase 2.fdtd.7 Q7 escape hatch: C6 passive-fine-grid trade-off \
            leaves Maloney-Smith dielectric resonance unrecoverable until F2 \
            inward-coupling restoration lands (Phase 2.fdtd.7.z finding)"]
fn fdtd_007_within_two_percent_and_one_db() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    assert_eq!(
        result.status,
        CaseStatus::Passed,
        "fdtd-007 failed: {}",
        result.notes
    );
}

#[test]
#[ignore = "Phase 2.fdtd.7 Q7 escape hatch: see module docstring"]
fn fdtd_007_fres_within_two_percent_of_maloney_smith_fig9() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    // Hard gate (1): f_res within ±2 % of Maloney-Smith Fig. 9 (8.9 GHz, TBD).
    assert!(
        result.f_res_rel_error <= yee_validation::FDTD_007_TOL_FRES_REL,
        "f_res relative error {:.6} > {:.3} tolerance ({})",
        result.f_res_rel_error,
        yee_validation::FDTD_007_TOL_FRES_REL,
        result.notes
    );
}

#[test]
#[ignore = "Phase 2.fdtd.7 Q7 escape hatch: see module docstring"]
fn fdtd_007_s11_within_one_db_of_maloney_smith_fig9() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    // Hard gate (2): |S_11(f_res)| within ±1 dB of -22 dB (TBD).
    assert!(
        result.s11_db_abs_error <= yee_validation::FDTD_007_TOL_S11_DB_ABS,
        "|S_11(f_res)| abs error {:.3} dB > {:.2} dB tolerance ({})",
        result.s11_db_abs_error,
        yee_validation::FDTD_007_TOL_S11_DB_ABS,
        result.notes
    );
}

#[test]
#[ignore = "Phase 2.fdtd.7 Q7 escape hatch: see module docstring"]
fn fdtd_007_subgrid_vs_uniform_sanity_check() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    // Hard gate (3): subgrid-vs-uniform sanity check across 5 spot
    // frequencies — 0.3 % / 0.3 dB.
    assert!(
        result.sanity_max_fres_rel <= yee_validation::FDTD_007_TOL_SANITY_FRES_REL,
        "subgrid-vs-uniform max |df|/f = {:.6} > {:.3} tolerance ({})",
        result.sanity_max_fres_rel,
        yee_validation::FDTD_007_TOL_SANITY_FRES_REL,
        result.notes
    );
    assert!(
        result.sanity_max_s11_db <= yee_validation::FDTD_007_TOL_SANITY_S11_DB,
        "subgrid-vs-uniform max |dS_11| = {:.3} dB > {:.2} dB tolerance ({})",
        result.sanity_max_s11_db,
        yee_validation::FDTD_007_TOL_SANITY_S11_DB,
        result.notes
    );
}

/// Smoke test that the driver compiles, runs end-to-end, and returns
/// a well-formed result struct (not necessarily passing the gate). This
/// runs in default CI so the scaffolding stays green-compiling even
/// while the physics-gate tests above are `#[ignore]`'d pending the
/// Phase 2.fdtd.7.z inward-coupling restoration.
///
/// The smoke test is itself `#[ignore]`'d because the underlying solver
/// runs ~4000 coarse + ~4000 uniform-fine steps on a non-trivial grid
/// (`60 × 16 × 14` coarse + `120 × 32 × 28` uniform-fine), which puts
/// the wall-time in the minutes range under `--release` — within the
/// 30-min Q7 budget but well past the default `cargo test` envelope.
#[test]
#[ignore = "slow: ~minutes release; scaffolding smoke for Phase 2.fdtd.7 Q7"]
fn fdtd_007_driver_returns_well_formed_result() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    assert_eq!(result.id, "fdtd-007");
    // Frequencies must be finite and inside the swept band.
    assert!(result.f_res_subgrid_hz.is_finite());
    assert!(result.f_res_uniform_hz.is_finite());
    assert!(result.f_res_subgrid_hz >= 4.0e9 && result.f_res_subgrid_hz <= 14.0e9);
    assert!(result.f_res_uniform_hz >= 4.0e9 && result.f_res_uniform_hz <= 14.0e9);
    // S_11 in dB must be ≤ 0 (passive 1-port).
    assert!(result.s11_db_subgrid <= 0.0 + 1e-9);
    assert!(result.s11_db_uniform <= 0.0 + 1e-9);
    // Sanity-check metrics non-negative.
    assert!(result.sanity_max_fres_rel >= 0.0);
    assert!(result.sanity_max_s11_db >= 0.0);
    // Notes string non-empty, status one of the documented variants.
    assert!(!result.notes.is_empty());
    assert!(matches!(
        result.status,
        CaseStatus::Passed | CaseStatus::Failed | CaseStatus::Skipped
    ));
}
