//! `fem-eig-002` production-gate test — lossy single-pole Drude SiO₂
//! rectangular metallic cavity TE_{101} resonance vs the analytic
//! continuum dispersion relation (spec §9 / ADR-0039 §9 / Pozar §3.1
//! lossy-dielectric extension).
//!
//! Drives [`yee_validation::run_fem_eig_002_lossy_sio2_cavity`] end-to-end
//! on the spec §9 cavity (`a = 10 mm`, `b = 5 mm`, `d = 20 mm`) meshed
//! with `(nx, ny, nz) = (8, 4, 16)` Kuhn 6-tet bricks (3072 tets,
//! ~12 k edges, ~2 k interior DoFs), populates the bulk filler with a
//! single-pole Drude material (`ε_∞ = 3.78`, `ω_p = 2π · 0.4 GHz`,
//! `γ = 2π · 2.0 GHz` — fused-silica `ε_∞` + exaggerated loss per
//! ADR-0039 §9), and asserts the four hard checks from the Phase
//! 4.fem.eig.1 spec §9:
//!
//! 1. **(A) Re(f) bound** — `|Re(f_FEM) − Re(f_analytic)| /
//!    Re(f_analytic) ≤ 0.5 %`.
//! 2. **(B) Im(f) / Q bound** — `|Im(f_FEM) − Im(f_analytic)| /
//!    |Im(f_analytic)| ≤ 5 %` (looser per spec §9 because Im(f)
//!    extraction is more sensitive to the Newton residual floor).
//! 3. **(C) Newton iteration cap** — the outer Newton tracker
//!    converges in ≤ 8 iterations from the lossless free-space
//!    warm-start ω₀ = c · √((π/a)² + (π/d)²) (≈ 2π · 16.77 GHz).
//! 4. **(D) No bisection / divergence fallback** — the driver does not
//!    return [`yee_fem::DispersiveError::NewtonDidNotConverge`].
//!
//! The analytic complex reference is computed inline by the driver via
//! inner Newton on the continuum dispersion relation
//! `ω² ε_Drude(ω) / c² = (π/a)² + (π/d)²` (spec §9.1) — the same
//! closed-form ε(ω) as the FEM material model, so any discrepancy is a
//! discretisation / Newton-tracker error, not a material-model mismatch.
//!
//! ## Test consolidation
//!
//! The four hard assertions land in a **single `#[test]` body** rather
//! than the four-test layout used by `fem_eig_001_*`. The reason is
//! the FEM solve itself (one full Newton-tracked dispersive eigenmode
//! on a 3072-tet mesh) dominates the per-test cost; splitting the
//! assertions across four tests would re-run the same expensive solve
//! four times. Each gate failure still produces a distinct, readable
//! assertion message because the `result.notes` string is folded into
//! every panic message, and the `result.status` field reports the
//! conjunction so an "overall pass" check is still cheap.
//!
//! See `crates/yee-fem/validation/README.md` for the validation rollup
//! and `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//! §9 for the Drude reference parameters. Wall-time budget is
//! documented as `< 60 s` in `--release` (spec §9 informational gate);
//! the test itself does **not** enforce a wall-time bound — that would
//! be environment-dependent — but is recorded in the validation
//! README's "wall-time" column.

use yee_validation::{CaseStatus, run_fem_eig_002_lossy_sio2_cavity};

#[test]
fn fem_eig_002_lossy_sio2_cavity_all_four_hard_gates() {
    let result = run_fem_eig_002_lossy_sio2_cavity()
        .expect("fem-eig-002 driver: gate (D) — no divergence — must hold");

    // Sanity: the measured complex frequency must be finite. Guards
    // against NaN propagation from a degenerate Newton step that does
    // not technically fail the convergence test but yields garbage.
    assert!(
        result.f_measured_hz.re.is_finite() && result.f_measured_hz.im.is_finite(),
        "f_measured contains non-finite component: {}",
        result.notes
    );
    assert!(
        result.f_analytic_hz.re.is_finite() && result.f_analytic_hz.im.is_finite(),
        "f_analytic contains non-finite component: {}",
        result.notes
    );

    // Spec §9 hard gate (A): Re(f_FEM) within ±0.5 % of Re(f_analytic).
    assert!(
        result.re_f_rel_error <= 0.005,
        "gate (A) FAILED: Re(f) relative error {:.6} > 0.5 % tolerance: {}",
        result.re_f_rel_error,
        result.notes
    );

    // Spec §9 hard gate (B): Im(f_FEM) within ±5 % of Im(f_analytic).
    assert!(
        result.im_f_rel_error <= 0.05,
        "gate (B) FAILED: Im(f) relative error {:.6} > 5 % tolerance: {}",
        result.im_f_rel_error,
        result.notes
    );

    // Spec §9 hard gate (C): outer Newton converges in ≤ 8 iterations.
    assert!(
        result.newton_iterations <= 8,
        "gate (C) FAILED: Newton iterations {} > 8 cap: {}",
        result.newton_iterations,
        result.notes
    );

    // Spec §9 hard gate (D): the driver returned `Ok(_)` rather than
    // surfacing a `NewtonDidNotConverge` from the underlying
    // `DispersiveSolver::solve_with_newton`. The `.expect` at the top
    // of the body is the actual gate (D) check; this assertion makes
    // it explicit for the reader and pins the overall status field as
    // a redundant cross-check.
    assert_eq!(
        result.status,
        CaseStatus::Passed,
        "fem-eig-002 overall status not Passed despite individual gates ok: {}",
        result.notes
    );
}
