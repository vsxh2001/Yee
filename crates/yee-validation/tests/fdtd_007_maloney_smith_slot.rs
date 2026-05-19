//! `fdtd-007` validation gate — dielectric-loaded thin slot antenna
//! vs Maloney & Smith 1993 Fig. 9.
//!
//! Drives [`yee_validation::run_fdtd_007_maloney_smith_slot`] end-to-end
//! on the rewired Phase 2.fdtd.7.z fixture: a slot of width `w = 0.5 mm`
//! and length `L = 30 mm` cut into a PEC ground plane backed by a
//! dielectric substrate (`ε_r = 2.2`, `h = 1.524 mm`), driven by a
//! delta-gap voltage source at the slot midpoint. The driver extracts
//! `S_11(f)` via per-bin DFTs of the lumped-port `V(t)` / `I(t)` traces
//! and reports `f_res` / `|S_11(f_res)|`.
//!
//! # Track UUUUUUUU rewire
//!
//! Track UUUUUUUU (Phase 2.fdtd.7.z) rewired the driver to consume the
//! per-cell ε map ([`yee_fdtd::YeeGrid::with_eps_r_cells`]; MMMMMMMM
//! `cb6f8ed`), per-component PEC mask
//! ([`yee_fdtd::YeeGrid::with_pec_mask_ey`] / `..._ez`; same commit),
//! and CPML-per-cell-ε coupling (PPPPPPPP `c57592f`) infrastructure
//! that landed against `main`. The substrate slab is now a true
//! per-cell `ε_r = 2.2` region; the ground plane is a per-component
//! PEC mask on the `i = i_gp` plane with the slot rectangle cut out;
//! both half-spaces are CPML-terminated for radiation.
//!
//! # Gate
//!
//! The original Phase 2.fdtd.7 Q7 DoD called for `±2 %` on `f_res`
//! against Maloney-Smith Fig. 9 (`8.9 GHz`, digitised TBD). Per the
//! brief's escape hatch:
//!
//! > f_res > 5 % off Fig. 9 → digitisation accuracy may need
//! > re-checking; report measured value and leave the gate
//! > `#[ignore]`'d with the actual measured value documented. Do NOT
//! > relax to > 5 %.
//!
//! The rewired uniform-fine driver lands at **`f_res ≈ 5.30 GHz`** with
//! `|S_11(f_res)| ≈ −6.26 dB` (see commit body for the full diagnostic
//! dump). That is `|df|/f ≈ 0.40` against the `8.9 GHz` reference —
//! well past the `±5 %` digitisation envelope — so the physics gates
//! against Fig. 9 stay `#[ignore]`'d pending one of:
//!
//! 1. **Reference-figure verification.** The `8.9 GHz` value carries
//!    a TBD flag (see [`yee_validation::FDTD_007_FRES_REF_HZ`]
//!    doc-comment); the only paper cited in the Phase 2.fdtd.7 spec
//!    (Maloney & Smith 1993 IEEE T-AP 41(5), "Wu-King resistive
//!    monopole") is on a *cylindrical monopole*, not a slot, so the
//!    Fig. 9 attribution itself may be wrong. The escape hatch reads
//!    this as "report measured value, leave gate `#[ignore]`'d" until
//!    the journal-figure scan is verified.
//! 2. **Geometry / mode interpretation.** The measured `5.3 GHz`
//!    resonance is consistent with the slot's *half-wave* mode in the
//!    `ε_eff ≈ 1.6` slab approximation (`c / (2 · 30 mm · √1.6)
//!    ≈ 3.95 GHz`, shifted upward by the finite slot width and the
//!    `(40 × 80 × 20)` mm grid's lateral mode loading). The published
//!    `8.9 GHz` may be a higher-order mode, a different slot length,
//!    or a different substrate thickness that the spec digitised
//!    incorrectly.
//!
//! Both follow-ups land as `fdtd-007.1` work. The driver itself is
//! green-compiling end-to-end and produces a well-formed result with
//! the expected `S_11 ≤ 0` passivity sanity check.
//!
//! # Subgridded variant — still `#[ignore]`'d
//!
//! The Phase 2.fdtd.7.y Step C6 un-ghosted-J Berenger closure leaves
//! the fine sub-grid effectively passive on source-on-coarse drive
//! (see LLLLLLLL commit body item 3). The subgridded variant of
//! `fdtd-007` remains `#[ignore]`'d pending the F2 inward-coupling
//! restoration (deferred from Track DDDDDDDD).
//!
//! Wall-time on the uniform-fine driver is < 10 s release on a
//! 44 × 90 × 31 grid; the `#[ignore]` on the smoke test below covers
//! debug-build / non-release-flag invocations where the run lands in
//! the minutes range.

use yee_validation::{CaseStatus, run_fdtd_007_maloney_smith_slot};

/// Smoke test — runs the rewired uniform-fine `fdtd-007` driver and
/// verifies it returns a well-formed [`yee_validation::Fdtd007ValidationResult`]
/// with finite measurements, the standard `id`, and a `S_11 ≤ 0`
/// passivity invariant.
///
/// **Track UUUUUUUU un-ignore:** this is the "uniform-fine sanity"
/// gate the brief asks to retire from the `#[ignore]` list. The driver
/// no longer runs a subgridded coarse + fine pair (the subgridded
/// variant is blocked on C6 inward coupling); the uniform-fine path
/// is fully functional after the MMMMMMMM + PPPPPPPP infrastructure
/// landings.
///
/// Wall-time is < 10 s release; debug builds land in the
/// minutes range, so the `--release` invocation is recommended.
#[test]
fn fdtd_007_uniform_fine_smoke() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");

    assert_eq!(result.id, "fdtd-007");
    // Frequencies must be finite and inside the swept band.
    assert!(result.f_res_uniform_hz.is_finite());
    assert!(
        result.f_res_uniform_hz >= 4.0e9 && result.f_res_uniform_hz <= 14.0e9,
        "f_res_uniform = {:.3} GHz outside the [4, 14] GHz sweep band",
        result.f_res_uniform_hz * 1e-9,
    );
    // For API compatibility, the subgridded fields are aliased to the
    // uniform-fine measurement when the subgridded path is skipped.
    assert!(result.f_res_subgrid_hz.is_finite());
    assert!(result.f_res_subgrid_hz >= 4.0e9 && result.f_res_subgrid_hz <= 14.0e9);

    // Passivity invariant: `|S_11|` of a passive 1-port must be ≤ 1
    // (i.e. `≤ 0 dB`).
    assert!(
        result.s11_db_uniform <= 0.0 + 1.0e-9,
        "|S_11(f_res)| = {:.3} dB exceeds passive-1-port bound (> 0 dB) — \
         the driver produced an active or numerically unstable response",
        result.s11_db_uniform,
    );

    // Sanity metrics non-negative.
    assert!(result.sanity_max_fres_rel >= 0.0);
    assert!(result.sanity_max_s11_db >= 0.0);

    // Diagnostic notes must be populated and the status must be one of
    // the documented variants.
    assert!(!result.notes.is_empty());
    assert!(matches!(
        result.status,
        CaseStatus::Passed | CaseStatus::Failed | CaseStatus::Skipped
    ));

    // Emit the notes string under `cargo test -- --nocapture` so the
    // operator can confirm the measured numbers without parsing the
    // result struct.
    eprintln!("fdtd-007 uniform-fine notes: {}", result.notes);
}

#[test]
#[ignore = "Phase 2.fdtd.7.z escape hatch (Track UUUUUUUU): measured \
            f_res ≈ 5.30 GHz on the uniform-fine driver is |df|/f ≈ 0.40 \
            off the Maloney-Smith 1993 Fig. 9 reference (8.9 GHz, TBD), \
            outside the ±5 % digitisation envelope. Gate stays \
            #[ignore]'d pending journal-figure verification or geometry \
            re-interpretation — see module docstring."]
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
#[ignore = "Phase 2.fdtd.7.z escape hatch (Track UUUUUUUU): measured \
            f_res ≈ 5.30 GHz vs Fig. 9 ref 8.9 GHz (TBD) — see module \
            docstring."]
fn fdtd_007_fres_within_two_percent_of_maloney_smith_fig9() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    assert!(
        result.f_res_rel_error <= yee_validation::FDTD_007_TOL_FRES_REL,
        "f_res relative error {:.6} > {:.3} tolerance ({})",
        result.f_res_rel_error,
        yee_validation::FDTD_007_TOL_FRES_REL,
        result.notes
    );
}

#[test]
#[ignore = "Phase 2.fdtd.7.z escape hatch (Track UUUUUUUU): measured \
            |S_11(f_res)| ≈ -6.26 dB vs Fig. 9 ref -22 dB (TBD); cavity-Q \
            sensitivity deferred to fdtd-007.1 radiation-CPML \
            calibration — see module docstring."]
fn fdtd_007_s11_within_one_db_of_maloney_smith_fig9() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
    assert!(
        result.s11_db_abs_error <= yee_validation::FDTD_007_TOL_S11_DB_ABS,
        "|S_11(f_res)| abs error {:.3} dB > {:.2} dB tolerance ({})",
        result.s11_db_abs_error,
        yee_validation::FDTD_007_TOL_S11_DB_ABS,
        result.notes
    );
}

#[test]
#[ignore = "Phase 2.fdtd.7.z (Track UUUUUUUU): the subgridded variant \
            remains blocked by the Phase 2.fdtd.7.y Step C6 un-ghosted-J \
            Berenger trade-off (F2 inward-coupling restoration deferred \
            from Track DDDDDDDD); the driver no longer runs a \
            subgrid-vs-uniform comparator, so this sanity gate cannot \
            execute its intended check."]
fn fdtd_007_subgrid_vs_uniform_sanity_check() {
    let result = run_fdtd_007_maloney_smith_slot().expect("fdtd-007 driver");
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
