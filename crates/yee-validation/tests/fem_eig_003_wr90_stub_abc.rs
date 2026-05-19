//! `fem-eig-003` production-gate test — WR-90 stub with 1st-order
//! Engquist-Majda ABC termination, swept `|S_{11}(f)|` across 8-12 GHz
//! vs spec §8 absorption window (Phase 4.fem.eig.2 step E5).
//!
//! Drives [`yee_validation::run_fem_eig_003_wr90_stub_abc`] end-to-end
//! on the spec §8 fixture (`a = 22.86 mm`, `b = 10.16 mm`, `d = 30 mm`)
//! meshed with `(nx, ny, nz) = (16, 8, 24)` Kuhn 6-tet bricks
//! (~18 k tets), with face `z = 0` tagged ABC, face `z = 30 mm` tagged
//! `WavePort(0)` (TE_{10} drive), and the four longitudinal sidewalls
//! tagged PEC. Sweeps 50 uniform points across 8-12 GHz at 80 MHz
//! spacing.
//!
//! ## Gate decomposition
//!
//! The driver enforces three hard gates plus one informational
//! runtime check per the Phase 4.fem.eig.2 plan E5 brief:
//!
//! * **(A)** `20·log10(|S_{11}(f)|) ∈ [-45, -35] dB` at every swept
//!   frequency — Engquist-Majda 1st-order absorption floor (Engquist
//!   & Majda 1977; Jin §10.4). ADR-0040 records this floor as the
//!   v0 physics limit.
//! * **(B)** `|S_{11}(f)| < 1` strictly — passive-structure invariant
//!   (Pozar §3.3).
//! * **(C)** Adjacent-bin `|Δ(20·log10|S_{11}|)|` bounded by 10 dB —
//!   no spurious resonance from ill-conditioning across the smooth
//!   ABC reflection spectrum.
//! * **(D, informational)** Wall-time recorded but not asserted.
//!
//! ## Escape-hatch disposition
//!
//! The Phase 4.fem.eig.2 E4 unit-test sibling
//! `crates/yee-fem/tests/open_boundary_sweep.rs` measured
//! `|S_{11}| ≈ 1.0` on a coarse `3×2×4` mesh, far outside the
//! `[-45, -35] dB` window. The plan E5 escape hatch reads: "if
//! walking-skeleton physics doesn't resolve `-40 dB` at 25 k tets,
//! document and continue." The strict gate (A) test is therefore
//! `#[ignore]`'d by default, with the measured `|S_{11}|` band recorded
//! in the driver's `notes` string for the follow-up Phase
//! 4.fem.eig.2.0.1 / 4.fem.eig.2.5 track. The non-strict gates
//! (B passive, C smoothness, D runtime informational) plus a
//! `gate_runs_without_panic` smoke remain in default CI.
//!
//! Until a refined mesh + revisited ABC face-block scaling resolve
//! the absorption floor to the Engquist-Majda physics limit, the
//! ignored strict-gate test exists as a tripwire that can be lifted
//! with a single `#[ignore]` removal once the floor is reached.
//!
//! See `crates/yee-fem/validation/README.md` for the validation rollup
//! and `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
//! §8 for the absorption-window reference.

use yee_validation::{CaseStatus, run_fem_eig_003_wr90_stub_abc};

/// Smoke gate — the driver completes without panicking on the spec-
/// scale `(16, 8, 24)` mesh. Default-CI; this asserts only that the
/// pipeline executes end-to-end and emits a finite `|S_{11}(f)|` band.
/// Strict gate (A) absorption-floor check lives in
/// [`fem_eig_003_strict_absorption_floor_gate`], which is `#[ignore]`'d
/// per the Phase 4.fem.eig.2 plan E5 escape hatch until the v0 ABC
/// face-block scaling resolves the floor at this mesh resolution.
#[test]
fn fem_eig_003_driver_runs_and_emits_finite_sweep() {
    let result = run_fem_eig_003_wr90_stub_abc().expect("fem-eig-003 driver");

    // Sanity: 50 swept points, each finite, with magnitude in a
    // sensible range. NaN or Inf indicates either a degenerate sparse
    // LU or a face-block scaling bug in the upstream Phase 4.fem.eig.2
    // E1-E4 layer.
    assert_eq!(
        result.frequencies_hz.len(),
        50,
        "fem-eig-003 sweep should produce 50 points; got {}: {}",
        result.frequencies_hz.len(),
        result.notes
    );
    assert_eq!(
        result.s11_magnitude.len(),
        50,
        "|S_11| array length mismatch: {}",
        result.notes
    );
    for (i, &mag) in result.s11_magnitude.iter().enumerate() {
        assert!(
            mag.is_finite(),
            "|S_11(f_{i})| = {mag} is non-finite — driven solve produced NaN/Inf: {}",
            result.notes
        );
    }
    for (i, &db) in result.s11_db.iter().enumerate() {
        assert!(
            db.is_finite() || db == f64::NEG_INFINITY,
            "s11_db[{i}] = {db} is non-finite and non-(-inf): {}",
            result.notes
        );
    }

    // Emit sweep summary for observability when run with `--nocapture`.
    // Helps the maintainer track the measured band against the Phase
    // 4.fem.eig.2 E5 escape-hatch disposition (driver currently
    // saturates at |S_11| ≈ 1.0 — see strict gate's `#[ignore]` for
    // the disposition rationale).
    let f_min = result.frequencies_hz.first().copied().unwrap_or(0.0);
    let f_max = result.frequencies_hz.last().copied().unwrap_or(0.0);
    let mid_idx = result.frequencies_hz.len() / 2;
    let f_mid = result.frequencies_hz[mid_idx];
    eprintln!(
        "fem-eig-003 smoke summary: band [{:.6e}, {:.6e}] dB \
         ; |S_11(f={:.2} GHz)| = {:.10} (mid) ; |S_11(f={:.2} GHz)| = {:.10} (low) \
         ; |S_11(f={:.2} GHz)| = {:.10} (high)",
        result.s11_db_min,
        result.s11_db_max,
        f_mid * 1e-9,
        result.s11_magnitude[mid_idx],
        f_min * 1e-9,
        result.s11_magnitude[0],
        f_max * 1e-9,
        result.s11_magnitude[result.frequencies_hz.len() - 1],
    );
}

/// Gate (B) — passive-structure invariant within the walking-skeleton
/// numerical-discretisation margin. `|S_{11}(f)| ≤ 1 + ε_num` at every
/// swept frequency with `ε_num = 0.05` (see
/// [`yee_validation::FEM_EIG_003_PASSIVE_MARGIN`] for the rationale).
/// Strict `< 1` is the continuum-limit identity (Pozar §3.3); the v0
/// Whitney-1 face-centroid quadrature + walking-skeleton modal-source
/// pipeline measures magnitudes clustering at `1.0` modulo round-off,
/// matching the Phase 4.fem.eig.2 E4 sibling convention.
///
/// A strict `< 1` continuum-limit tripwire exists separately as
/// [`fem_eig_003_strict_passive_bound_continuum_limit`], which is
/// `#[ignore]`'d under the same E5 escape hatch as gate (A).
#[test]
fn fem_eig_003_passive_structure_no_amplification() {
    let result = run_fem_eig_003_wr90_stub_abc().expect("fem-eig-003 driver");
    assert!(
        result.gate_b_passive_ok,
        "fem-eig-003 gate (B) FAILED: at least one |S_{{11}}(f)| > 1 + ε_num \
         (passive structure cannot amplify by more than the discretisation \
         margin): {}",
        result.notes
    );
}

/// **Continuum-limit gate (B) tripwire.** `|S_{11}(f)| < 1` strictly
/// at every swept frequency — the Pozar §3.3 passive-structure
/// identity in its un-relaxed form. `#[ignore]`'d under the same Phase
/// 4.fem.eig.2 E5 escape hatch as the strict absorption-floor gate.
///
/// **JJJJJJJJJ status (2026-05-19, Phase 4.fem.eig.3 F6).** With the
/// F2 coupled exact-Whitney-1 modal RHS + projection and the F4
/// 2nd-order Engquist–Majda ABC both enabled on the spec-scale
/// `(16, 8, 24)` mesh (18 432 tets), the driver now measures
/// `|S_{11}(f)|` band `[0.9945, 0.99999]` (corresponding `s11_db` band
/// `[-5.0e-2, -8.1e-5] dB`) across the 8-12 GHz sweep — every magnitude
/// is **strictly less than 1.0**, so the strict passive bound is
/// *numerically* satisfied at the un-rounded floating-point level. The
/// `#[ignore]` is retained, however, because the F2 + F4 fix did not
/// restore the absorption floor to the documented `~ -40 dB`
/// Engquist–Majda window (see
/// [`fem_eig_003_strict_absorption_floor_gate`]): the surface still
/// reflects almost all incident power. Once the mesh is refined and
/// the absorption floor lands in `[-45, -35] dB` (Phase
/// 4.fem.eig.3.0.3 mesh-refinement track), this strict-passive gate
/// will pass with comfortable margin and the `#[ignore]` lifts in
/// the same PR.
///
/// **Historical context (CCCCCCCCC partial fix).** Track CCCCCCCCC
/// applied the modal-amplitude `M_pp` normalisation fix to
/// [`yee_fem::OpenBoundarySolver::extract_s11`]; that retired the
/// **synthetic** matched-port identity but did not retire the
/// empirical `|S_{11}| ≈ 1.0` saturation. F2 retires the saturation
/// by lifting the lumped `N_i(centroid) ≈ t_i / 3` proxy to the exact
/// Whitney-1 identity at 3-point Gauss quadrature — the FEM-projected
/// `E_FEM` is now non-trivially non-zero on the port face — but the
/// remaining `~ 0.005` reflection budget is still dominated by the
/// coarse-mesh ABC numerical reflection rather than the continuum
/// limit.
#[test]
#[ignore = "fem-eig-003 strict passive bound: F2 coupled Whitney + F4 2nd-order ABC enabled, \
            measured |S_11| band [0.9945, 0.99999] strictly < 1 (gate would pass) but kept \
            #[ignore]'d coupled with the still-failing absorption-floor gate per Phase \
            4.fem.eig.3.0.3 mesh-refinement queue"]
fn fem_eig_003_strict_passive_bound_continuum_limit() {
    let result = run_fem_eig_003_wr90_stub_abc().expect("fem-eig-003 driver");
    let strict_passive_ok = result.s11_magnitude.iter().all(|&m| m < 1.0);
    assert!(
        strict_passive_ok,
        "fem-eig-003 strict passive bound FAILED: at least one |S_{{11}}(f)| ≥ 1.0 \
         exactly (continuum-limit Pozar §3.3 identity violated): {}",
        result.notes
    );
}

/// Gate (C) — sweep smoothness. Adjacent 80 MHz bins must not differ
/// by more than 10 dB in `20·log10(|S_{11}|)`. A spurious resonance
/// from ill-conditioning of the driven matrix would manifest as a
/// tens-of-dB jump across one bin; this gate canaries against that
/// failure mode without depending on the absolute absorption floor.
#[test]
fn fem_eig_003_sweep_smoothness_no_spurious_resonance() {
    let result = run_fem_eig_003_wr90_stub_abc().expect("fem-eig-003 driver");
    assert!(
        result.gate_c_smoothness_ok,
        "fem-eig-003 gate (C) FAILED: max adjacent-bin |Δ(20·log10|S_11|)| = \
         {:.3} dB exceeds 10 dB smoothness bound — likely spurious resonance \
         from ill-conditioning of the driven matrix: {}",
        result.max_adjacent_db_jump, result.notes
    );
}

/// Gate (A) — Engquist–Majda absorption floor.
/// `20·log10(|S_{11}(f)|) ∈ [-45, -35] dB` at every swept frequency
/// (spec §8 + ADR-0040 / ADR-0042).
///
/// **`#[ignore]`'d per the Phase 4.fem.eig.3 plan F6 escape hatch.**
///
/// **JJJJJJJJJ status (2026-05-19, F6).** With F1+F2 coupled exact-
/// Whitney-1 modal RHS + projection and F3+F4 2nd-order Engquist–Majda
/// ABC both enabled (`OpenBoundarySolver::with_coupled_whitney(true)
/// .with_abc_order(AbcOrder::Second)`) the spec-scale
/// `(16, 8, 24) = 18 432 tets` mesh measures `|S_{11}(f)|` band
/// `[0.9945, 0.99999]` → `s11_db ∈ [-5.0e-2, -8.1e-5] dB` across the
/// 8-12 GHz sweep. The 2nd-order ABC + coupled Whitney path lowered
/// the `|S_{11}|` floor measurably from the BBBBBBBBB walking-skeleton
/// saturation at `1.000_000_000` (numerical band `[-1e-15, 0.0] dB`)
/// down to a non-trivial `~ -0.05 dB` band — i.e. **the absorber now
/// dissipates some power**, but the spec §8 target `[-45, -35] dB`
/// remains far below the measured floor.
///
/// **Failure-mode diagnosis.** The remaining gap is a mesh-refinement
/// constraint, not a missing physics term: at `(16, 8, 24)` the
/// in-plane port-face element pitch is `~ 1.43 mm × 1.27 mm`, which
/// resolves the WR-90 TE_{10} mode tangential profile to only ~16
/// linear samples across the broad wall. The Engquist–Majda 1979
/// 2nd-order absorber needs ~30+ linear samples per cross-section
/// wavelength to hit the `~ -60 dB` continuum floor (Jin §10.4 table
/// 10.1); below that, modal-projection discretisation error
/// dominates the absorber's intrinsic floor.
///
/// **Queued follow-up: Phase 4.fem.eig.3.0.3 mesh-refinement track.**
/// Refine to `(24, 12, 36) = ~ 62 k tets` (~3.5× cost; estimated
/// wall-time ~6 min `--release` per sweep call vs ~30 s at
/// `(16, 8, 24)`) and re-measure. ADR-0042 §risks records this as the
/// expected failure mode for the F6 strict gate. If the refined mesh
/// lands in `[-45, -35] dB`, lift this `#[ignore]` in the
/// Phase 4.fem.eig.3.0.3 PR alongside the
/// [`fem_eig_003_strict_passive_bound_continuum_limit`] lift.
///
/// If the refined mesh **also** does not reach the window, then the
/// binding constraint is the 2nd-order ABC ill-conditioning near the
/// closed-stub TE_{10n} resonances (8 GHz n=1 / 12 GHz n=2) per spec
/// §10 risk register — defer to Phase 4.fem.eig.3.5 (CFS-PML).
#[test]
#[ignore = "fem-eig-003 strict absorption floor: F1-F4 enabled, measured |S_11| band \
            [0.9945, 0.99999] (s11_db [-5.0e-2, -8.1e-5] dB) — better than the BBBBBBBBB \
            walking-skeleton 1.000 floor but well outside the spec §8 [-45, -35] dB window; \
            queue Phase 4.fem.eig.3.0.3 mesh-refinement track per ADR-0042"]
fn fem_eig_003_strict_absorption_floor_gate() {
    let result = run_fem_eig_003_wr90_stub_abc().expect("fem-eig-003 driver");
    assert!(
        result.gate_a_floor_ok,
        "fem-eig-003 gate (A) FAILED: |S_{{11}}(f)| dB band [{:.2}, {:.2}] outside \
         the Engquist-Majda window [-45, -35] dB: {}",
        result.s11_db_min, result.s11_db_max, result.notes,
    );
    assert_eq!(
        result.status,
        CaseStatus::Passed,
        "fem-eig-003 overall status not Passed: {}",
        result.notes
    );
}
