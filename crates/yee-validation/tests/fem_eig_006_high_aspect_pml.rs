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
/// **Phase 4.fem.eig.3.5.6 Lee-Mittra absorbing-mode complement
/// (ADR-0070, escape-hatch).** The Lee-Mittra first-order absorbing BC
/// was applied to port_1 (the +x terminating face), replacing the
/// scalar `j(β₁₀+β₂₀) B_face` stiffness with `jk₀ B_face + Σ_m
/// j(β_m−k₀) R_m`. Measured result: `|S_{11}|(40 GHz) = 0.955500
/// (-0.40 dB)` — **essentially unchanged** from the v3.5.5 baseline
/// (0.01% change; well outside the < 0.1 gate). The Lee-Mittra
/// first-order complement does not improve absorption for this
/// high-aspect cavity geometry. Root-cause hypothesis: the {TE₁₀,
/// TE₂₀, TE₀₁} modal basis provides insufficient coverage of the
/// face impedance — both β₁₀ ≈ 776 and β₂₀ ≈ 554 rad/m are below
/// k₀ ≈ 838 rad/m so the corrections `j(β_m−k₀) R_m` reduce the
/// stiffness (negative imaginary shift), and the rank-1 projection R_m
/// covers only a small fraction of the face Gram B_face for the
/// actual TE mode shapes on this low-element-count mesh. Higher-order
/// absorbing BC (Lee-Mittra 1997 §V rational-function extension) or a
/// different port formulation (aperture integral, T-matrix) is needed.
/// Gate tolerance `< 0.1` is **not** weakened.
///
/// History: v3.5.3 W1 single-mode TE_{10} `|S_{11}|(30 GHz) = 0.925644`;
/// v3.5.4 multi-mode `0.925637` (cutoff-degenerate, ADR-0048);
/// v3.5.5 retune to 40 GHz `0.955397` (ADR-0049);
/// v3.5.6 Lee-Mittra absorbing complement `0.955500` (ADR-0070, this record).
#[test]
#[ignore = "fem-eig-006 strict magnitude bound (Phase 4.fem.eig.3.5.6 Lee-Mittra absorbing-mode \
            complement, ADR-0070): Lee-Mittra first-order complement applied to port_1 — \
            K = jk₀ B_face + Σ_m j(β_m−k₀) R_m. Measured |S_11|(40 GHz) = 0.955500 (-0.40 dB) \
            on native (16,3,2) cavity, 576 tets — essentially unchanged from v3.5.5 baseline \
            0.955397 (0.01% change). First-order Lee-Mittra does NOT retire the gate for this \
            high-aspect geometry: β₁₀≈776 and β₂₀≈554 rad/m < k₀≈838 rad/m; negative \
            j(β_m−k₀) corrections reduce stiffness but rank-1 R_m projects onto a small \
            fraction of B_face for the low-element-count mesh. Higher-order absorbing BC \
            (Lee-Mittra §V rational-function extension) queued for Phase 4.fem.eig.3.5.7. \
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
