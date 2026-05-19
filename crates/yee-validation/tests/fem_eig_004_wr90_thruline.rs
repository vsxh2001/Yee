//! `fem-eig-004` production-gate test — WR-90 two-port thru-line at
//! 10 GHz (Phase 4.fem.eig.3 step F6). Exercises the multi-port
//! `sweep_matrix` entry point introduced in F5 with F1+F2 coupled
//! exact-Whitney-1 modal RHS + projection.
//!
//! Drives [`yee_validation::run_fem_eig_004_wr90_thruline`] end-to-end
//! on a lossless air-filled WR-90 section (`a × b × d = 22.86 × 10.16
//! × 30 mm`) meshed with `(12, 6, 18)` Kuhn 6-tet bricks (~7.8 k tets).
//! Both end faces (`z = 0` and `z = d`) carry the TE_{10} mode; the
//! four sidewalls are PEC.
//!
//! ## Gate decomposition
//!
//! Three hard gates at the 10 GHz center of a five-point sweep
//! `{9.8, 9.9, 10.0, 10.1, 10.2} GHz`:
//!
//! * **(A)** `|S_{21}(10 GHz)|` within `±0.1 dB` of 0 dB — perfect
//!   transmission through a lossless thru-line.
//! * **(B)** `|S_{11}(10 GHz)| < -20 dB` — low matched-port
//!   reflection.
//! * **(C)** `|S_{12}(10 GHz) − S_{21}(10 GHz)| < 1e-3` — reciprocity
//!   invariant (Pozar §4.3 passive lossless network).
//!
//! ## References
//!
//! * Phase 4.fem.eig.3 design spec
//!   `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md`
//!   §8 (fem-eig-004 gate criteria).
//! * Phase 4.fem.eig.3 plan
//!   `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-3.md` step F6.
//! * Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §4.3
//!   (reciprocity `S_{p,q} = S_{q,p}` for lossless passive multi-ports).
//! * Sheen, D. M., Ali, S. M., Abouzahra, M. D., Katehi, P. B. L.,
//!   *IEEE Trans. MTT* 38(7) (1990), pp. 849-857 — eq. 7 multi-port
//!   column-extraction convention used by
//!   [`yee_fem::OpenBoundarySolver::sweep_matrix`].

use yee_validation::run_fem_eig_004_wr90_thruline;

/// Smoke gate — the driver runs to completion and emits a finite
/// five-point sweep of the 2 × 2 S-matrix. Default-CI.
#[test]
fn fem_eig_004_driver_runs_and_emits_finite_sweep() {
    let result = run_fem_eig_004_wr90_thruline().expect("fem-eig-004 driver");

    assert_eq!(
        result.frequencies_hz.len(),
        5,
        "fem-eig-004 should sweep 5 frequencies; got {}: {}",
        result.frequencies_hz.len(),
        result.notes
    );
    assert_eq!(
        result.s.len(),
        5,
        "fem-eig-004 should produce 5 S-matrices; got {}: {}",
        result.s.len(),
        result.notes
    );
    for (k, s_k) in result.s.iter().enumerate() {
        assert_eq!(
            s_k.shape(),
            (2, 2),
            "fem-eig-004 s[{k}] should be 2×2; got {:?}",
            s_k.shape()
        );
        for q in 0..2 {
            for p in 0..2 {
                let v = s_k[(q, p)];
                assert!(
                    v.norm().is_finite(),
                    "fem-eig-004 s[{k}][({q}, {p})] = {v} is non-finite: {}",
                    result.notes
                );
            }
        }
    }

    eprintln!(
        "fem-eig-004 smoke summary: at 10 GHz |S_21| = {:.3} dB, \
         |S_11| = {:.3} dB, |S_12 − S_21| = {:.3e}",
        result.s21_db_at_10ghz, result.s11_db_at_10ghz, result.reciprocity_residual_at_10ghz,
    );
}

/// Gate (A) — through-line transmission. `|S_{21}(10 GHz)|` within
/// `±0.1 dB` of 0 dB (i.e. linear magnitude in `[0.9886, 1.0116]`).
/// A lossless thru-line of length `d ≈ 30 mm` at 10 GHz on
/// WR-90 transmits the incident TE_{10} mode with unity magnitude
/// modulo the modal-projection discretisation budget — the gate is
/// the headline cross-port-projection sanity check.
///
/// **Default-CI** at Phase 4.fem.eig.3 F6 with F1+F2 coupled exact-
/// Whitney-1 enabled: the driver measures `|S_{21}(10 GHz)| ≈ -0.045
/// dB` on the spec-scale `(12, 6, 18) = 7.8 k tets` mesh, comfortably
/// inside the ±0.1 dB window.
#[test]
fn fem_eig_004_through_transmission_gate() {
    let result = run_fem_eig_004_wr90_thruline().expect("fem-eig-004 driver");
    assert!(
        result.gate_a_through_transmission_ok,
        "fem-eig-004 gate (A) FAILED: |S_21(10 GHz)| = {:.3} dB outside ±0.1 dB \
         of 0 dB: {}",
        result.s21_db_at_10ghz, result.notes
    );
}

/// Gate (B) — matched-port reflection. `|S_{11}(10 GHz)| < -20 dB`.
/// Both end faces are TE_{10} wave-ports of identical cross-section,
/// so the incident TE_{10} mode at port 0 sees a "matched" boundary
/// modulo modal-projection discretisation.
///
/// **Default-CI** at Phase 4.fem.eig.3 F6 with F1+F2 enabled: the
/// driver measures `|S_{11}(10 GHz)| ≈ -53 dB`, far below the
/// -20 dB headroom — the coupled exact-Whitney-1 wave-port pair
/// realises matching cleanly at this resolution.
#[test]
fn fem_eig_004_matched_port_reflection_gate() {
    let result = run_fem_eig_004_wr90_thruline().expect("fem-eig-004 driver");
    assert!(
        result.gate_b_matched_reflection_ok,
        "fem-eig-004 gate (B) FAILED: |S_11(10 GHz)| = {:.3} dB does not satisfy \
         < -20 dB: {}",
        result.s11_db_at_10ghz, result.notes
    );
}

/// Gate (C) — reciprocity. `|S_{12}(10 GHz) − S_{21}(10 GHz)| < 1e-3`.
/// The thru-line is a passive lossless reciprocal structure (Pozar
/// §4.3 continuum identity); reciprocity is the cleanest invariant
/// the multi-port sweep should satisfy regardless of the modal-source
/// discretisation, because both off-diagonal entries are projected
/// using the **same** Whitney-1 basis and the **same** LU factor
/// (only the excited-port RHS differs between the (q=1, p=0) and
/// (q=0, p=1) columns).
///
/// **Default-CI** because reciprocity does not require strict
/// matched-port physics — it is a structural symmetry of the driven
/// matrix induced by the wave-port bilinear form, and the F5 multi-
/// port plumbing must preserve it bit-for-bit modulo round-off. A
/// failure here points at a systemic asymmetry bug in `sweep_matrix`
/// (e.g. swapped row/column indexing).
#[test]
fn fem_eig_004_reciprocity_gate() {
    let result = run_fem_eig_004_wr90_thruline().expect("fem-eig-004 driver");
    assert!(
        result.gate_c_reciprocity_ok,
        "fem-eig-004 gate (C) FAILED: |S_12 − S_21| at 10 GHz = {:.3e} exceeds \
         1e-3 (passive lossless WR-90 thru-line should be reciprocal): {}",
        result.reciprocity_residual_at_10ghz, result.notes
    );
}
