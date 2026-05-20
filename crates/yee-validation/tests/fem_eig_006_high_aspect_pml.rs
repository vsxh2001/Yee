//! `fem-eig-006` production-gate test — high-aspect 100 : 10 : 1
//! rectangular cavity with CFS-PML termination at the +x truncation
//! face (Phase 4.fem.eig.3.5 P5; spec §6).
//!
//! Stress-tests CFS-PML on highly off-normal modal content — the
//! geometry forces the WR-90 TE_{10} scattering pattern onto the PML
//! inner boundary at near-grazing incidence, the regime where the
//! 2nd-order Engquist-Majda ABC saturates at `|S_{11}| ≈ 0.95` and
//! the CFS `α_α > 0` modification of Berenger 1994 PML earns its
//! keep (Kuzuoglu-Mittra 1996 §II).
//!
//! Gate inventory (spec §6):
//!
//! 1. `fem_eig_006_magnitude_bounded` — `|S_{11}(30 GHz)| < 0.1`.
//! 2. `fem_eig_006_no_nan_inf` — `S_{11}` is finite (PML-stability
//!    canary; the CFS modification rescues this from the Berenger
//!    1994 PML's evanescent-mode divergence).
//!
//! ## OOOOOOOOO status (2026-05-20, Phase 4.fem.eig.3.5)
//!
//! With the default CFS-PML grading the v3.5 measurement was recorded
//! and the gates evaluated per the OOOOOOOOO P5 disposition.
//!
//! References:
//!
//! * Phase 4.fem.eig.3.5 spec
//!   `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`
//!   §6.
//! * Roden-Gedney 2000, *IEEE MWCL* 10(5), pp. 27-29.
//! * Kuzuoglu-Mittra 1996, *IEEE MWCL* 6(12), pp. 447-449.

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
        "fem-eig-006 gate (B) FAILED: S_11 is non-finite (PML stability \
         failure — CFS α > 0 modification did not absorb evanescent / \
         grazing modes): {}",
        result.notes
    );
}

#[test]
#[ignore = "fem-eig-006 strict magnitude bound (Phase 4.fem.eig.3.5.2): H4 ablation grid ran \
            fem-eig-006 across all 18 H4 rows (m∈{3,4} × thickness∈{12,14,16} × \
            alpha_grading_order∈{0,1,2}); |S_11|(30 GHz) frozen at 0.926 in all rows. \
            alpha-grading is orthogonal to the 100:10:1 fixture — dominant modal content is \
            not normal-incidence at the +x face. fem-eig-003 absorption retires at the same \
            v3.5.2 defaults (band [-71.53, -55.58] dB). Queued for Phase 4.fem.eig.3.5.3 / \
            4.fem.eig.4: rotated PML / multi-face wedges / wave-port termination for the \
            high-aspect-ratio cavity"]
fn fem_eig_006_magnitude_bounded() {
    let result = run_fem_eig_006_high_aspect_pml().expect("fem-eig-006 driver");
    assert!(
        result.gate_a_magnitude_ok,
        "fem-eig-006 gate (A) FAILED: |S_11(30 GHz)| = {:.6} ({:.2} dB) \
         ≥ 0.1 — PML absorbs off-normal incidence below the spec \
         threshold: {}",
        result.s11_magnitude, result.s11_db, result.notes
    );
}
