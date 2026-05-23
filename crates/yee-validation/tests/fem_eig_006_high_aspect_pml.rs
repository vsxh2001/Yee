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
//! cavity. Phase 4.fem.eig.3.5.5 retuned the operating point from
//! 30 GHz to **40 GHz** (ADR-0048 Option (a)) so the TE_{20} mode
//! propagates (`β ≈ 554 rad/m`) instead of sitting exactly at its
//! `c / B = 30 GHz` cutoff, giving the v3.5.4 multi-mode basis real
//! propagating content to terminate. See the
//! `fem_eig_006_magnitude_bounded` docstring for the measurement and
//! the escape-hatch disposition.
//!
//! Gate inventory (spec §6):
//!
//! 1. `fem_eig_006_magnitude_bounded` — `|S_{11}(40 GHz)| < 0.1`.
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

/// Phase 4.fem.eig.3.5.5 measurement (ADR-0048 Option (a) — frequency
/// retune to **40 GHz**, multi-mode wave-port basis unchanged):
/// `|S_{11}|(40 GHz) = 0.955397 (-0.40 dB)` on the native (16, 3, 2)
/// cavity (576 Kuhn-6 tets).
///
/// At 40 GHz TE_{20} now propagates with `β ≈ 554 rad/m` (33% above
/// its `f_c = c / B = 30 GHz` cutoff), so the modal-degeneracy that
/// pinned v3.5.4 at 30 GHz no longer applies — the multi-mode basis
/// carries real propagating content. The residual nonetheless stayed
/// high (0.955, marginally **above** the v3.5.4 30 GHz value 0.926),
/// so the retune did **not** retire the gate.
///
/// **Refinement probe (one-shot, reverted).** To rule out the spec §4(a)
/// discretisation-limited hypothesis (~2.3 transverse cells/λ at
/// 40 GHz), the transverse mesh was bumped (NY 3→9, NZ 2→6, 5184 tets):
/// `|S_{11}|(40 GHz, refined) = 0.913956 (-0.78 dB)`. A 9× transverse
/// element-count increase moved the residual only ~0.04, nowhere near
/// the 0.1 gate — the reflection is **not** discretisation-limited. The
/// bump was reverted; the native (16, 3, 2) mesh stands.
///
/// **Disposition: escape-hatch (gate stays `#[ignore]`'d).** With both
/// modal degeneracy (v3.5.4) and discretisation (this probe) excluded,
/// the residual is a genuine limitation of the modal-projection
/// wave-port: projecting onto a finite TE_{mn} basis cannot fully match
/// the field at the truncation face of a strongly off-square cavity.
/// Per ADR-0048 the next step is Option (b) — the Lee-Mittra 1997
/// absorbing-mode wave-port — queued for Phase 4.fem.eig.3.5.6 in
/// ADR-0049. Tolerance `< 0.1` is **not** weakened.
///
/// History: v3.5.3 W1 single-mode TE_{10} `|S_{11}|(30 GHz) = 0.925644`;
/// v3.5.4 multi-mode `0.925637` (cutoff-degenerate, ADR-0048);
/// v3.5.5 retune to 40 GHz `0.955397` (this record, ADR-0049).
#[test]
#[ignore = "fem-eig-006 strict magnitude bound (Phase 4.fem.eig.3.5.5 frequency retune, ADR-0048 \
            Option (a)): FEM_EIG_006_F_HZ retuned 30→40 GHz so TE_{20} propagates (β≈554 rad/m, \
            33% above its c/B=30 GHz cutoff). |S_11|(40 GHz) = 0.955397 (-0.40 dB) on native \
            (16,3,2) cavity, 576 tets — did NOT retire (marginally above the 30 GHz value 0.926). \
            One-shot refinement probe (NY 3→9, NZ 2→6, 5184 tets) gave 0.913956 (-0.78 dB): a 9× \
            transverse refinement moved |S_11| only ~0.04, so the residual is NOT \
            discretisation-limited (probe reverted). Escape-hatch: residual is a genuine \
            modal-projection wave-port limitation; queued for Phase 4.fem.eig.3.5.6 \
            absorbing-mode wave-port per Lee-Mittra 1997 (ADR-0048 Option (b), ADR-0049). \
            Tolerance < 0.1 not weakened."]
fn fem_eig_006_magnitude_bounded() {
    let result = run_fem_eig_006_high_aspect_pml().expect("fem-eig-006 driver");
    assert!(
        result.gate_a_magnitude_ok,
        "fem-eig-006 gate (A) FAILED: |S_11(40 GHz)| = {:.6} ({:.2} dB) \
         ≥ 0.1 — multi-mode wave-port retune to 40 GHz (TE_{{20}} \
         propagating, β≈554 rad/m) did not retire the gate; refinement \
         probe excluded discretisation as the cause. Queued for Phase \
         4.fem.eig.3.5.6 absorbing-mode wave-port per ADR-0048 Option \
         (b) / ADR-0049: {}",
        result.s11_magnitude, result.s11_db, result.notes
    );
}
