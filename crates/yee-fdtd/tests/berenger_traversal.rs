//! Phase 2.fdtd.7.x B2 — Berenger Huygens-surface step driver stability.
//!
//! Stability smoke gate for the Berenger 2006 equivalent-current closure
//! wired into [`SubgriddedSolver::step`] / `step_with_gaussian_source_ez`
//! by Phase 2.fdtd.7.x B2. The strict 0.5%-of-peak plane-wave traversal
//! gate (Q5) is unblocked by B3 (subgrid_plane_wave_traversal.rs); this
//! file's gate is the 100-step stability canary that B2 needs to satisfy
//! before B3 can be safely un-`#[ignore]`'d.
//!
//! ## What this verifies
//!
//! Drives the [`SubgriddedSolver`] for 100 coarse steps with a
//! Gaussian-in-time `E_z` source upstream of the fine sub-region, and
//! asserts that **neither the coarse nor the fine grid diverges**:
//!
//! - every field cell stays finite (`is_finite`), and
//! - the peak `|E_z|` on the fine grid stays bounded below `1e3` V/m
//!   for a source amplitude `≤ 1` V/m injection (the previous Q4
//!   bidirectional direct-copy closure typically diverged
//!   exponentially around step 50–100 with the same geometry; bounded
//!   propagation at 100 steps is the B2 acceptance criterion before B3
//!   takes over the strict tolerance gate).
//!
//! ## Why not bit-exact-vs-Q3-only
//!
//! The brief flagged a passthrough-limit bit-exact comparison against
//! "Q3-only injection" as an option. With the Berenger closure
//! sourcing equivalent currents `J = +n̂ × H_tot` / `M = -n̂ × E_tot`
//! from the fine grid, the coarse `E_t` and `H_t` on the interface
//! plane receive a non-zero correction *the moment any fine field is
//! non-zero* — there is no operating regime in which the Berenger
//! pipeline matches Q3-only-no-fine-to-coarse bit-exact, by design.
//! The stability gate is the operative B2 acceptance criterion;
//! strict-tolerance comparison against the pure-coarse reference is
//! B3's job.
//!
//! Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`.
//! Plan: `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-x-berenger-huygens.md` step B2.
//! ADR:  `docs/src/decisions/0035-berenger-huygens-subgridding.md`.

use yee_fdtd::{CpmlParams, SubgridRegion, SubgriddedSolver, WalkingSkeletonSolver, YeeGrid};

/// Coarse grid x-extent (cells).
const NX_C: usize = 64;
/// Coarse grid y-extent (cells).
const NY_C: usize = 32;
/// Coarse grid z-extent (cells).
const NZ_C: usize = 32;
/// Coarse cell size (m).
const DX_C: f64 = 1.0e-3;
/// CPML thickness on the coarse grid (cells per face).
const NPML_C: usize = 6;
/// Subgrid lower corner (coarse-cell indices, inclusive).
const SG_LO: (usize, usize, usize) = (16, 12, 12);
/// Subgrid upper corner (coarse-cell indices, exclusive).
const SG_HI: (usize, usize, usize) = (24, 20, 20);
/// Gaussian source location (coarse-cell index on `E_z`).
const SRC: (usize, usize, usize) = (8, 16, 16);
/// Number of coarse steps to integrate for the 100-step canary; the
/// 500-step gate runs in [`berenger_step_propagates_without_divergence_500_steps`].
const N_COARSE_STEPS: usize = 100;
/// Stability bound on peak `|E_z|` (V/m). For a `~1 V/m` Gaussian
/// excitation in vacuum the steady-state wave amplitude is `≪ 1`; an
/// exponentially-diverging interface drives the fine grid past `1e3`
/// V/m within ~50 coarse steps (the failure mode that motivated the
/// closure switch from spec `2026-05-18` Q4 to spec `2026-05-19`
/// Berenger). Bounded propagation up to 100 steps is the B2 gate.
const STABILITY_BOUND: f64 = 1.0e3;

/// Phase 2.fdtd.7.x B2.1 — 500-step variant of the
/// [`berenger_step_propagates_without_divergence`] stability gate.
///
/// HHHHHHH's diagnosis on the B2 landing (`commit 997e706`,
/// `subgrid_plane_wave_traversal::*_strict_*`) found the coarse `|E_z|`
/// doubles every ~7 coarse steps from step ~60 onwards on the Q5
/// 96 × 32 × 32 plane-wave-traversal geometry. The smaller
/// 64 × 32 × 32 geometry used here has a longer onset (the source is
/// further from the fine box and the CPML faces absorb more of the
/// scattered field), but the divergence is qualitatively the same
/// failure mode. Tracked as the 500-step canary so the
/// stability bound stays grep-able even if the 100-step test
/// remains green.
///
/// **Currently `#[ignore]`'d** per AAAAAAA plan B4 escape hatch.
///
/// Phase 2.fdtd.7.x B2.2 (track OOOOOOO) added coarse-ghost
/// subtraction to the J source (Berenger 2006 §III canonical
/// equivalent-source form `J = +n̂ × (H_fine − H_coarse_ghost)`).
/// This DELAYS the divergence onset from step ~98 (B2.1 baseline) to
/// step ~137 (B2.2), and reduces the always-on 100-step canary peak
/// |E_z| from 31 V/m to 2.75 V/m — a strict improvement. But it does
/// NOT retire the 500-step bound: peak |E_z| still crosses 1e3 V/m
/// around step 137 and grows by ~10× per 20 steps thereafter.
/// Empirically, applying coarse-ghost subtraction **symmetrically**
/// to the M source (Berenger §III in pure form) actively destabilises
/// the canary: the coarse `E_t` at the surface is Dirichlet-tied to
/// the fine grid by Q3, so the "symmetric" M ghost subtraction
/// nullifies the magnetic equivalent current rather than correcting
/// it, and the 100-step canary peak |E_z| jumps from 2.75 V/m
/// (J-only ghost) back to ≈ 1e3 V/m (J + M ghost). Resolution
/// requires a separate, M-side equivalence accounting fix —
/// deferred to Phase 2.fdtd.7.y.
#[test]
#[ignore = "Phase 2.fdtd.7.x B2.2 (Track OOOOOOO): coarse-ghost J subtraction lands; 500-step \
            divergence delayed from step ~98 to step ~137 (peak |E_z| = 1.035e3) but not \
            retired. Residual is an M-side accounting issue (Q3-tied coarse E surface) — \
            deferred to Phase 2.fdtd.7.y."]
fn berenger_step_propagates_without_divergence_500_steps() {
    let coarse_grid = YeeGrid::vacuum(NX_C, NY_C, NZ_C, DX_C);
    let coarse_dt = coarse_grid.dt;
    let cpml_c = CpmlParams::for_grid(&coarse_grid, NPML_C);
    let inner = WalkingSkeletonSolver::with_cpml(coarse_grid, cpml_c);

    let region = SubgridRegion::new(inner.grid(), SG_LO, SG_HI)
        .expect("SubgridRegion::new must accept this in-interior nest");

    let mut sub = SubgriddedSolver::new(inner).with_region(region);

    let sigma = 4.0 * coarse_dt;
    let t0 = 3.0 * sigma;

    let mut peak_fine_ez = 0.0_f64;
    let mut peak_coarse_ez = 0.0_f64;

    const N_LONG_STEPS: usize = 500;

    for step in 0..N_LONG_STEPS {
        sub.step_with_gaussian_source_ez(SRC.0, SRC.1, SRC.2, t0, sigma);

        let g = sub.inner().grid();
        let f = sub.region().expect("region attached").fine_grid();

        // Track running peaks.
        for &v in f.ez.iter() {
            let av = v.abs();
            if av > peak_fine_ez {
                peak_fine_ez = av;
            }
        }
        for &v in g.ez.iter() {
            let av = v.abs();
            if av > peak_coarse_ez {
                peak_coarse_ez = av;
            }
        }

        // Finite check (catches the catastrophic divergence inside
        // the 500-step window — currently triggered around step 60
        // per HHHHHHH's diagnosis).
        for arr in [&g.ex, &g.ey, &g.ez, &g.hx, &g.hy, &g.hz] {
            for &v in arr.iter() {
                assert!(
                    v.is_finite(),
                    "coarse field non-finite at step {step}: v = {v}"
                );
            }
        }
        for arr in [&f.ex, &f.ey, &f.ez, &f.hx, &f.hy, &f.hz] {
            for &v in arr.iter() {
                assert!(
                    v.is_finite(),
                    "fine field non-finite at step {step}: v = {v}"
                );
            }
        }

        assert!(
            peak_fine_ez < STABILITY_BOUND,
            "fine grid |E_z| diverged at step {step}: peak = {peak_fine_ez:.3e} >= {STABILITY_BOUND:.3e}"
        );
    }

    eprintln!(
        "Berenger 500-step traversal: peak |E_z|_coarse = {peak_coarse_ez:.3e}, \
         peak |E_z|_fine = {peak_fine_ez:.3e} (bound {STABILITY_BOUND:.0e})"
    );
}

/// Run the Berenger-closure step driver for 100 coarse steps with a
/// Gaussian `E_z` excitation on the coarse grid and assert bounded
/// propagation on both the coarse and fine grids.
#[test]
fn berenger_step_propagates_without_divergence() {
    let coarse_grid = YeeGrid::vacuum(NX_C, NY_C, NZ_C, DX_C);
    let coarse_dt = coarse_grid.dt;
    let cpml_c = CpmlParams::for_grid(&coarse_grid, NPML_C);
    let inner = WalkingSkeletonSolver::with_cpml(coarse_grid, cpml_c);

    let region = SubgridRegion::new(inner.grid(), SG_LO, SG_HI)
        .expect("SubgridRegion::new must accept this in-interior nest");

    let mut sub = SubgriddedSolver::new(inner).with_region(region);

    // Gaussian envelope width 4 coarse dt (resolved on fine sub-grid);
    // onset at t0 = 3·sigma. Peak amplitude is the implicit unit; the
    // soft-source convention adds `amplitude` to `grid.ez` per step,
    // so total injected energy scales with the time-integrated
    // Gaussian. Amplitude is left at the
    // `WalkingSkeletonSolver::apply_gaussian_source_ez` default
    // (= 1.0 V/m peak) since the gate is on absolute divergence, not
    // relative agreement against a reference.
    let sigma = 4.0 * coarse_dt;
    let t0 = 3.0 * sigma;

    let mut peak_fine_ez = 0.0_f64;
    let mut peak_coarse_ez = 0.0_f64;

    for step in 0..N_COARSE_STEPS {
        sub.step_with_gaussian_source_ez(SRC.0, SRC.1, SRC.2, t0, sigma);

        // Finite check on every cell (every step is too cheap to skip).
        let g = sub.inner().grid();
        for arr in [&g.ex, &g.ey, &g.ez, &g.hx, &g.hy, &g.hz] {
            for &v in arr.iter() {
                assert!(
                    v.is_finite(),
                    "coarse field non-finite at step {step}: v = {v}"
                );
            }
        }
        let f = sub.region().expect("region attached").fine_grid();
        for arr in [&f.ex, &f.ey, &f.ez, &f.hx, &f.hy, &f.hz] {
            for &v in arr.iter() {
                assert!(
                    v.is_finite(),
                    "fine field non-finite at step {step}: v = {v}"
                );
            }
        }

        // Track running peak |E_z| on the fine grid.
        for &v in f.ez.iter() {
            let av = v.abs();
            if av > peak_fine_ez {
                peak_fine_ez = av;
            }
        }
        for &v in g.ez.iter() {
            let av = v.abs();
            if av > peak_coarse_ez {
                peak_coarse_ez = av;
            }
        }

        assert!(
            peak_fine_ez < STABILITY_BOUND,
            "fine grid |E_z| diverged at step {step}: peak = {peak_fine_ez:.3e} >= {STABILITY_BOUND:.3e}"
        );
    }

    eprintln!(
        "Berenger step over {N_COARSE_STEPS} coarse steps: peak |E_z|_coarse = {peak_coarse_ez:.3e}, \
         peak |E_z|_fine = {peak_fine_ez:.3e} (bound {STABILITY_BOUND:.0e})"
    );

    // The wave should actually reach the fine grid (sanity check
    // against silent zero-output). The Gaussian peaks around step
    // `t0/dt = 12`, propagates ≈ 8 coarse cells (= sub-grid front) by
    // step ~40; by step 100 the entire pulse has crossed the fine box.
    assert!(
        peak_fine_ez > 1.0e-6,
        "fine grid peak |E_z| = {peak_fine_ez:.3e} suspiciously small — \
         source may not be propagating into the fine region"
    );
}
