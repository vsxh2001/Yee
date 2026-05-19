//! Phase 2.fdtd.7.x B4 — Q6 round-trip energy-balance gate.
//!
//! Discrete energy-balance gate for the Berenger 2006 Huygens-surface
//! fine ↔ coarse closure wired into [`SubgriddedSolver::step`] by Phase
//! 2.fdtd.7.x B2 ([`crates/yee-fdtd/tests/berenger_traversal.rs`] is the
//! companion 100-step stability canary; this file's gate is the
//! long-time energy-conservation canary).
//!
//! ## What this verifies
//!
//! Initialises a Gaussian-modulated sinusoidal `E_z` pulse inside the
//! fine subregion of a closed PEC cavity (`(64, 64, 64)` coarse cells +
//! a centred `(16, 16, 16)`-coarse-cell fine nest at 2× refinement →
//! 32×32×32 fine cells). No CPML — the outer boundary is hard PEC on
//! every face, so the total electromagnetic energy `W(t)` is exactly
//! conserved in the continuum, and the discrete Yee leapfrog conserves
//! `W` up to truncation drift that grows polynomially with `N`. The
//! pulse propagates outward through the Huygens interface into the
//! coarse region, reflects off the PEC walls, and re-traverses the
//! interface back into the fine region. Round-trip energy drift over
//! `N` coarse steps is the gate.
//!
//! ## Energy integral
//!
//! ```text
//! W(t) = (eps_0 / 2) · sum |E|^2 · dV  +  (mu_0 / 2) · sum |H|^2 · dV
//! ```
//!
//! Computed as the **sum of two disjoint contributions** so the fine
//! interior is not double-counted with the coarse storage at the same
//! footprint:
//!
//! - **Coarse contribution:** every coarse cell strictly outside the
//!   fine box footprint (`(i, j, k) ∉ [lo, hi)`) contributes
//!   `(eps_0/2) (E_x^2 + E_y^2 + E_z^2) · dV_c + (mu_0/2)(H_x^2 +
//!   H_y^2 + H_z^2) · dV_c`.
//! - **Fine contribution:** every fine cell contributes the same with
//!   `dV_f = dV_c / 8` (2× refinement on each axis).
//!
//! The half-integer Yee staggering means each component has a slightly
//! different storage footprint; the disjoint-cell rule below uses cell
//! indices `(i, j, k)` matching the **primary-cell index** to decide
//! membership. This drops a thin layer of half-cell-staggered
//! contributions at the box boundary; the Berenger closure handles the
//! interface plane exactly via equivalent currents (so the small
//! boundary-cell discrepancy is bounded by the Berenger closure error,
//! which is what the gate measures).
//!
//! ## Status (Phase 2.fdtd.7.x B4)
//!
//! Both the 1000-step smoke variant and the 10 000-step strict gate
//! are **`#[ignore]`'d** per the AAAAAAA plan B4 escape hatch — the
//! Berenger 2006 closure as landed by B2 surfaces an energy drift far
//! above the spec §5 0.5% bound on this fine-seeded cavity geometry
//! (≈ 1.10× per-coarse-step growth, ≈ 1.3 × 10³ at `N = 30`,
//! catastrophic divergence well before `N = 1000`). The dominant
//! failure mode matches spec §6 risk 3 (Berenger `J` / `M`
//! time-centering mismatch — silent on the Q5 source-driven traversal
//! gate, surfaces on a fine-initial-state cavity). The bound is
//! preserved per the escape hatch's "do NOT widen the bound"
//! directive; resolution is deferred to a future spec amendment
//! (`Phase 2.fdtd.7.y`).
//!
//! The `q6_energy_accounting_initial_state` test (always-on) verifies
//! that the energy-accounting machinery is wired correctly and
//! reports a finite, strictly positive `W(0)`. The drift-asserting
//! tests are kept compiled but `#[ignore]`'d so they can be exercised
//! on demand with `--include-ignored`.
//!
//! Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
//! §5 (Q6 row), §6 risk 3 (time-centering canary).
//! Plan: `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-x-berenger-huygens.md`
//! step B4.
//! ADR:  `docs/src/decisions/0035-berenger-huygens-subgridding.md`.

use ndarray::Array3;
use yee_core::units::{C0, EPS0, MU0};
use yee_fdtd::{SubgridRegion, SubgriddedSolver, WalkingSkeletonSolver, YeeGrid};

/// Coarse grid extent (cells) along each axis. Cubic so the fine region
/// can be centred.
const N_C: usize = 64;
/// Coarse cell size (m). Vacuum.
const DX_C: f64 = 1.0e-3;
/// Fine subregion extent (coarse cells) along each axis. Centred in
/// the parent so the round-trip distance to every PEC wall is equal.
const SG_EXTENT: usize = 16;
/// Fine subregion lower corner (coarse-cell indices, inclusive).
const SG_LO: (usize, usize, usize) = (
    (N_C - SG_EXTENT) / 2,
    (N_C - SG_EXTENT) / 2,
    (N_C - SG_EXTENT) / 2,
);
/// Fine subregion upper corner (coarse-cell indices, exclusive).
const SG_HI: (usize, usize, usize) = (
    SG_LO.0 + SG_EXTENT,
    SG_LO.1 + SG_EXTENT,
    SG_LO.2 + SG_EXTENT,
);
/// Drift tolerance per spec §5. Round-trip `|W(N) − W(0)| / W(0)`.
/// Preserved at 0.5% per the AAAAAAA plan B4 escape-hatch instruction
/// to *not* widen the bound when the closure surfaces a drift > 0.5%.
const DRIFT_BOUND: f64 = 5.0e-3;
/// Short-N smoke-variant step count (kept for completeness; the test
/// is `#[ignore]`'d alongside the 10 000-step gate — see the note on
/// `q6_round_trip_smoke_1000_steps`).
const N_SHORT: usize = 1000;

/// Centre of the Gaussian-modulated sinusoid initial pulse, in **fine**
/// cells. Fine grid is `2 · SG_EXTENT = 32` cells per axis; the centre
/// is `(16, 16, 16)`.
const PULSE_CENTRE_F: (f64, f64, f64) = (SG_EXTENT as f64, SG_EXTENT as f64, SG_EXTENT as f64);
/// Gaussian envelope width (fine cells). Chosen so the pulse fits well
/// inside the 32-fine-cell box with negligible amplitude at the box
/// boundary at `t = 0` (`exp(-((SG_EXTENT/2)/sigma)^2) ≈ 3e-3` for
/// `sigma = 4`).
const PULSE_SIGMA_F: f64 = 4.0;
/// Carrier wavelength of the modulation (fine cells). 8 fine cells per
/// wavelength = 16 coarse-cells-equivalent, well-resolved on the Yee
/// stencil.
const PULSE_LAMBDA_F: f64 = 8.0;

/// Seed a Gaussian-modulated sinusoidal `E_z` initial condition on the
/// fine grid only. Returns nothing — mutates `region.fine_grid_mut()`.
///
/// `E_z(x, y, z, t = 0) = exp(-r²/(2σ²)) · sin(k · (x - x_c))`
///
/// with `r² = (x - x_c)² + (y - y_c)² + (z - z_c)²`, all coordinates in
/// **fine-cell** units measured from the fine-grid origin. The fine
/// grid's `E_z` storage has shape `[fine_nx + 1, fine_ny + 1, fine_nz]`;
/// we seed every interior cell, leaving the outermost layer at zero
/// (those are Dirichlet-set by the Q3 coarse→fine interpolation each
/// step and would be overwritten anyway).
fn seed_gaussian_pulse(region: &mut SubgridRegion) {
    let fine = region.fine_grid_mut();
    let (fnx_p1, fny_p1, fnz) = fine.ez.dim();
    let k_wave = 2.0 * std::f64::consts::PI / PULSE_LAMBDA_F;
    for i in 1..fnx_p1 - 1 {
        for j in 1..fny_p1 - 1 {
            for k in 0..fnz {
                let dx = (i as f64) - PULSE_CENTRE_F.0;
                let dy = (j as f64) - PULSE_CENTRE_F.1;
                let dz = (k as f64 + 0.5) - PULSE_CENTRE_F.2;
                let r2 = dx * dx + dy * dy + dz * dz;
                let envelope = (-0.5 * r2 / (PULSE_SIGMA_F * PULSE_SIGMA_F)).exp();
                fine.ez[(i, j, k)] = envelope * (k_wave * dx).sin();
            }
        }
    }
}

/// `0.5 · eps_0 · eps_r · sum |E|^2 · dV  +  0.5 · mu_0 · mu_r · sum |H|^2 · dV`
/// over the coarse grid, **excluding** cells whose primary-cell index
/// `(i, j, k)` lies in the fine-box footprint `[lo, hi)`.
///
/// The Yee components have slightly different shapes (`E_x` on edges
/// parallel to x has shape `[nx, ny+1, nz+1]`, etc.); we iterate each
/// component independently and use the primary-cell membership of the
/// underlying index. The half-integer staggered boundary cells of each
/// component are treated as "inside the box" iff their primary-cell
/// index lies in `[lo, hi)` along every axis. This is conservative
/// (excludes one row of half-cell-staggered contributions on each face)
/// and is symmetric for `t = 0` and `t = N`, so it cancels in the drift
/// ratio `|W(N) − W(0)| / W(0)`.
fn coarse_energy_outside_fine_box(
    grid: &YeeGrid,
    lo: (usize, usize, usize),
    hi: (usize, usize, usize),
) -> f64 {
    let dv = grid.dx * grid.dy * grid.dz;
    let mut e_energy = 0.0;
    let mut h_energy = 0.0;

    let inside = |i: usize, j: usize, k: usize| -> bool {
        i >= lo.0 && i < hi.0 && j >= lo.1 && j < hi.1 && k >= lo.2 && k < hi.2
    };

    accumulate_component(&grid.ex, &mut e_energy, &inside);
    accumulate_component(&grid.ey, &mut e_energy, &inside);
    accumulate_component(&grid.ez, &mut e_energy, &inside);
    accumulate_component(&grid.hx, &mut h_energy, &inside);
    accumulate_component(&grid.hy, &mut h_energy, &inside);
    accumulate_component(&grid.hz, &mut h_energy, &inside);

    0.5 * grid.eps_r * EPS0 * e_energy * dv + 0.5 * grid.mu_r * MU0 * h_energy * dv
}

/// `0.5 · eps_0 · eps_r · sum |E|^2 · dV_f  +  0.5 · mu_0 · mu_r · sum |H|^2 · dV_f`
/// over the entire fine grid.
fn fine_energy(fine: &YeeGrid) -> f64 {
    let dv = fine.dx * fine.dy * fine.dz;
    let mut e2 = 0.0;
    let mut h2 = 0.0;
    for &v in fine.ex.iter() {
        e2 += v * v;
    }
    for &v in fine.ey.iter() {
        e2 += v * v;
    }
    for &v in fine.ez.iter() {
        e2 += v * v;
    }
    for &v in fine.hx.iter() {
        h2 += v * v;
    }
    for &v in fine.hy.iter() {
        h2 += v * v;
    }
    for &v in fine.hz.iter() {
        h2 += v * v;
    }
    0.5 * fine.eps_r * EPS0 * e2 * dv + 0.5 * fine.mu_r * MU0 * h2 * dv
}

/// Sum `v²` over every cell of `arr` whose index `(i, j, k)` is
/// **not** in the fine-box footprint (per `inside`). `arr` is one of
/// the six Yee component arrays — its index naming follows the
/// primary-cell convention (the array index is the primary-cell index
/// modulo a half-cell offset baked into the staggering, which we treat
/// as integer).
fn accumulate_component<F>(arr: &Array3<f64>, acc: &mut f64, inside: &F)
where
    F: Fn(usize, usize, usize) -> bool,
{
    let (nx, ny, nz) = arr.dim();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if inside(i, j, k) {
                    continue;
                }
                let v = arr[(i, j, k)];
                *acc += v * v;
            }
        }
    }
}

/// Build a fresh subgridded closed-PEC cavity and seed the
/// Gaussian-modulated sinusoid initial pulse. The coarse solver uses
/// hard PEC (no CPML) so the round-trip is well-defined; the fine grid
/// inherits the Q3 Dirichlet boundary which is fed by the coarse
/// snapshot pair every coarse step.
fn build_cavity() -> SubgriddedSolver {
    let coarse_grid = YeeGrid::vacuum(N_C, N_C, N_C, DX_C);
    let inner = WalkingSkeletonSolver::new(coarse_grid);
    let region = SubgridRegion::new(inner.grid(), SG_LO, SG_HI)
        .expect("SubgridRegion::new accepts a centred 16-coarse-cell nest in a 64³ parent");
    let mut sub = SubgriddedSolver::new(inner).with_region(region);
    seed_gaussian_pulse(sub.region_mut().expect("region attached"));
    sub
}

/// Drive the closed-PEC subgridded cavity for `n_steps` coarse steps
/// and return `(W_0, W_N)`. The fine pulse propagates outward through
/// the Huygens interface, reflects off the PEC walls of the coarse
/// cavity, and re-traverses the interface back into the fine region.
fn run_and_measure(n_steps: usize) -> (f64, f64) {
    let mut sub = build_cavity();
    let w0 = total_energy(&sub);
    for _ in 0..n_steps {
        sub.step();
    }
    let wn = total_energy(&sub);
    (w0, wn)
}

/// Total stored EM energy = coarse-outside-box + fine.
fn total_energy(sub: &SubgriddedSolver) -> f64 {
    let coarse = coarse_energy_outside_fine_box(sub.inner().grid(), SG_LO, SG_HI);
    let fine = fine_energy(sub.region().expect("region attached").fine_grid());
    coarse + fine
}

/// Verify that the energy-accounting machinery is wired correctly: at
/// `t = 0` after seeding the fine-grid Gaussian-modulated sinusoid,
/// `W(0)` is finite and strictly positive. This is the *only*
/// always-on test in this file — see the `#[ignore]` notes on
/// `q6_round_trip_smoke_1000_steps` and `q6_round_trip_10000_steps`
/// below for why the assertion variants are gated.
#[test]
fn q6_energy_accounting_initial_state() {
    let sub = build_cavity();
    let w0 = total_energy(&sub);
    assert!(w0.is_finite(), "W(0) must be finite, got {w0}");
    assert!(
        w0 > 0.0,
        "W(0) must be strictly positive after seeding, got {w0}"
    );
    let transit_steps = (2.0 * (N_C as f64) * DX_C / C0 / sub_dt_coarse()) as usize;
    eprintln!(
        "Q6 initial state: W(0) = {w0:.6e} J on a {N_C}^3 coarse + 2x32^3 fine \
         closed-PEC cavity; one round-trip ≈ {transit_steps} coarse steps"
    );
}

/// Smoke variant of the Q6 gate at `N = 1000` coarse steps with the
/// strict 0.5% drift bound preserved. **`#[ignore]`'d** alongside the
/// 10 000-step gate: per the AAAAAAA plan B4 escape hatch, when the
/// Berenger 2006 closure surfaces a drift `> 0.5%`, the test is
/// marked `#[ignore]` and the bound is preserved (not widened) — the
/// fix is a future spec amendment (`Phase 2.fdtd.7.y`), not a
/// regression-tracker tolerance bump.
///
/// Empirical drift on the spec geometry under the B2-landed Berenger
/// closure: at `N = 30` coarse steps the drift is already ≈ 1.3 ×
/// 10³, and grows roughly exponentially (≈ 1.10× per coarse step) —
/// far past the 0.5% bound. The dominant divergence mode under this
/// initial condition (Gaussian-modulated sinusoid seeded **inside the
/// fine subregion** with zero coarse field at `t = 0`) is the
/// Berenger time-centering mismatch flagged in spec §6 risk 3 — the
/// `J = +n̂ × H` source is sampled at `t = n + 1/2` and the
/// `M = -n̂ × E` source at `t = n + 1`, but the fine grid's outer
/// `E_t` Dirichlet boundary is held at the **coarse**-side
/// interpolation rather than re-evolved through the Berenger surface
/// (so the round-trip energy ledger does not close at the discrete
/// level for fine-seeded pulses). This was silent on the Q5 source-
/// driven traversal gate (no fine-initial-state) and surfaces here as
/// the spec §6 risk 3 canary predicted.
///
/// Run with `cargo test -p yee-fdtd --release --test
/// subgrid_energy_balance -- --include-ignored` to reproduce the
/// regression-tracked drift values.
#[test]
#[ignore = "Phase 2.fdtd.7.y C6 (Track DDDDDDDD escape-hatched to F1 — drop J-side coarse-H \
            ghost subtraction): finite (≈ 75% drift at N = 1000, ≈ 79% at N = 10000) but \
            still above the 0.5% bound. Q6 seeds an initial pulse on the fine grid; under \
            Mur-only inward coupling + un-ghosted J source, the fine-side pulse couples into \
            the coarse via J = +n̂ × H_fine ≠ 0 but the magnetic equivalent current \
            M = -n̂ × (E_post - E_pre) is now the compensating-source form (not the canonical \
            -n̂ × E_fine), so the energy ledger does not close at the discrete level on this \
            fine-seeded geometry. Resolution requires restoring a stable inward coupling \
            channel together with the matching M-source re-balance — the F2 candidate that \
            re-introduced positive-feedback instability during C6 bring-up. Deferred to a \
            Phase 2.fdtd.7.y.α spec amendment."]
fn q6_round_trip_smoke_1000_steps() {
    let (w0, wn) = run_and_measure(N_SHORT);
    assert!(w0.is_finite() && w0 > 0.0, "W(0) must be finite, got {w0}");
    assert!(
        wn.is_finite(),
        "W({N_SHORT}) non-finite: energy balance diverged"
    );
    let drift = (wn - w0).abs() / w0;
    eprintln!(
        "Q6 smoke @ N = {N_SHORT}: W(0) = {w0:.6e} J, W(N) = {wn:.6e} J, \
         |dW|/W(0) = {:.4}% (bound {:.1}%)",
        100.0 * drift,
        100.0 * DRIFT_BOUND,
    );
    assert!(
        drift <= DRIFT_BOUND,
        "Q6 smoke energy drift {:.4}% exceeds {:.1}% bound at N = {N_SHORT}",
        100.0 * drift,
        100.0 * DRIFT_BOUND,
    );
}

/// Q6 spec gate: 10 000 coarse steps in a closed PEC cavity, drift
/// bound `|W(N) − W(0)| / W(0) ≤ 0.5%`. **`#[ignore]`'d** per the
/// AAAAAAA plan B4 escape hatch — the Berenger closure currently
/// surfaces drift far above the 0.5% bound (see the
/// `q6_round_trip_smoke_1000_steps` doc-comment for the underlying
/// failure mode). The bound is preserved (NOT widened); resolution is
/// a future spec amendment, not a tolerance change.
///
/// Run with `cargo test -p yee-fdtd --release --test
/// subgrid_energy_balance -- --include-ignored` to exercise.
#[test]
#[ignore = "Phase 2.fdtd.7.y C6 (Track DDDDDDDD escape-hatched to F1): drift ≈ 79% at \
            N = 10000 with W(N) finite throughout. See the matching note on \
            q6_round_trip_smoke_1000_steps for the un-ghosted-J + compensating-M asymmetry \
            rationale on this fine-seeded cavity. Resolution deferred to a Phase 2.fdtd.7.y.α \
            spec amendment that restores stable inward coupling + the matching M-source \
            re-balance."]
fn q6_round_trip_10000_steps() {
    const N: usize = 10_000;
    let (w0, wn) = run_and_measure(N);
    assert!(w0.is_finite() && w0 > 0.0, "W(0) must be finite, got {w0}");
    assert!(
        wn.is_finite(),
        "W({N}) non-finite: energy balance diverged at the 10000-step gate"
    );
    let drift = (wn - w0).abs() / w0;
    eprintln!(
        "Q6 gate @ N = {N}: W(0) = {w0:.6e} J, W(N) = {wn:.6e} J, \
         |dW|/W(0) = {:.4}% (bound {:.1}%)",
        100.0 * drift,
        100.0 * DRIFT_BOUND,
    );
    assert!(
        drift <= DRIFT_BOUND,
        "Q6 energy drift {:.4}% exceeds {:.1}% bound at N = {N}",
        100.0 * drift,
        100.0 * DRIFT_BOUND,
    );
}

/// Helper: the coarse-grid `dt`. Built once from `YeeGrid::vacuum`'s
/// 0.9× Courant factor, which is deterministic and avoids leaking the
/// `YeeGrid` constructor into the smoke-test sanity check.
fn sub_dt_coarse() -> f64 {
    let g = YeeGrid::vacuum(N_C, N_C, N_C, DX_C);
    g.dt
}
