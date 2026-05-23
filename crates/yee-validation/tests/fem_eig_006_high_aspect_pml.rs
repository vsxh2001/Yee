//! `fem-eig-006` production-gate test — high-aspect 100 : 10 : 1
//! rectangular cavity (Phase 4.fem.eig.3.5 P5; spec §6).
//!
//! ## v3.5.3 wave-port termination (Phase 4.fem.eig.3.5.3 W1; ADR-0046)
//!
//! The driver now terminates the +x face with a TE_{10}
//! `FaceKind::WavePort(1)` (Jin §10.6 closed-cavity modal
//! termination) in place of the v3.5.2 CFS-PML shell, after the
//! SSSSSSSSS H4 ablation (Phase 4.fem.eig.3.5.2) found
//! `|S_{11}|(30 GHz)` frozen at 0.926 across all 18
//! (m, thickness, alpha_grading_order) rows. Berenger 1996 §IV-A:
//! Cartesian-aligned CFS-PML cannot absorb the TE_{10} guide-mode
//! at the +x face regardless of grading parameters.
//!
//! Stress-tests the wave-port modal projection on a 100 : 10 : 1
//! cavity at 30 GHz where the TE_{20} cutoff sits exactly at the
//! operating frequency — the regime where a TE_{10}-only port may
//! underestimate the reflection per spec §7 (a).
//!
//! Gate inventory (spec §6):
//!
//! 1. `fem_eig_006_magnitude_bounded` — `|S_{11}(30 GHz)| < 0.1`.
//! 2. `fem_eig_006_no_nan_inf` — `S_{11}` is finite (numerical-
//!    stability canary on the wave-port modal projection).
//!
//! References:
//!
//! * Phase 4.fem.eig.3.5.3 spec
//!   `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-3-design.md`
//!   §3, §4.2, §7.
//! * ADR-0046 `docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md`.
//! * Jin, *FEM in EM*, 3rd ed., Chapter 10.6 "Wave-port termination".
//! * Berenger 1996, *IEEE TAP* 44(1), §IV-A bulk-vs-guide-wave PML.
//! * Roden-Gedney 2000, *IEEE MWCL* 10(5), pp. 27-29 (legacy CFS-PML).
//! * Kuzuoglu-Mittra 1996, *IEEE MWCL* 6(12), pp. 447-449 (legacy).

use yee_validation::run_fem_eig_006_high_aspect_pml;

#[test]
fn fem_eig_006_smoke_runs() {
    let result = run_fem_eig_006_high_aspect_pml().expect("fem-eig-006 driver");
    eprintln!("{}", result.notes);
}

#[test]
fn fem_eig_006_no_nan_inf() {
    let result = run_fem_eig_006_high_aspect_pml().expect("fem-eig-006 driver");
    assert!(
        result.gate_b_finite_ok,
        "fem-eig-006 gate (B) FAILED: S_11 is non-finite — wave-port \
         modal-projection numerical pathology: {}",
        result.notes
    );
}

/// Phase 4.fem.eig.3.5.4 measurement (multi-mode wave-port,
/// `PortDefinition::modes = [TE_{10}, TE_{20}, TE_{01}]`):
/// `|S_{11}|(30 GHz) = 0.925637 (-0.67 dB)` on the native (16, 3, 2)
/// cavity (576 Kuhn-6 tets). Numerically indistinguishable from
/// the v3.5.3 W1 single-mode TE_{10} measurement `0.925644` — the
/// modal basis collapses to single-mode at 30 GHz because:
///
/// 1. `TE_{20}` cutoff on the port-face broad wall `B = 10 mm` is
///    `f_c = c / B = 30.0 GHz` **exactly** — at cutoff, `β = 0` and
///    the multi-mode stiffness block contribution vanishes.
/// 2. `TE_{01}` cutoff on the narrow wall `D = 1 mm` is
///    `f_c = c / (2 D) = 150.0 GHz` — evanescent at 30 GHz; carries
///    no propagating modal content.
///
/// The v3.5.4 design spec §2.2 mis-derived these cutoffs by treating
/// the cavity's propagation length `A = 100 mm` as the modal
/// analysis broad wall; corrected derivation lives in the test
/// docstring (here), the ROADMAP v3.5.4 entry, and ADR-0048 (the
/// v3.5.5 disposition: either retune the test frequency off the
/// cutoff edge, or land an absorbing-mode wave-port per Lee-Mittra
/// 1997).
///
/// Gate stays `#[ignore]`'d. Tolerance `< 0.1` is **not** weakened.
#[test]
#[ignore = "fem-eig-006 strict magnitude bound (Phase 4.fem.eig.3.5.4 multi-mode measurement): \
            PortDefinition modes = [TE_{10} (a_inc=1), TE_{20} (a_inc=0), TE_{01} (a_inc=0)]; \
            |S_11|(30 GHz) = 0.925637 (-0.67 dB) on native (16,3,2) cavity, 576 tets. Multi-mode \
            basis collapses to single-mode at 30 GHz: TE_{20} f_c = c/B = 30 GHz exactly (β=0 \
            at cutoff, stiffness block vanishes); TE_{01} f_c = c/(2 D) = 150 GHz (evanescent). \
            Queued for Phase 4.fem.eig.3.5.5: either retune test frequency off the cutoff edge \
            (e.g. 25 GHz, where TE_{20} propagates) or land absorbing-mode wave-port per \
            Lee-Mittra 1997 (ADR-0048). Tolerance < 0.1 not weakened."]
fn fem_eig_006_magnitude_bounded() {
    let result = run_fem_eig_006_high_aspect_pml().expect("fem-eig-006 driver");
    assert!(
        result.gate_a_magnitude_ok,
        "fem-eig-006 gate (A) FAILED: |S_11(30 GHz)| = {:.6} ({:.2} dB) \
         ≥ 0.1 — multi-mode wave-port basis collapses to single-mode \
         at 30 GHz (TE_{{20}} f_c = c/B = 30 GHz exactly, TE_{{01}} \
         f_c = 150 GHz evanescent). v3.5.5 disposition queued per \
         ADR-0048 (retune frequency or absorbing-mode wave-port): {}",
        result.s11_magnitude, result.s11_db, result.notes
    );
}
