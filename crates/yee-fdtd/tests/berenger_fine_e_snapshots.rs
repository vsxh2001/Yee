//! Phase 2.fdtd.7.y Step C1 — pre/post fine-`E` snapshot capture tests.
//!
//! Pins the two new snapshot helpers added by Step C1 of the Phase
//! 2.fdtd.7.y Berenger M-side compensating-source amendment
//! (ADR-0038 / `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md`):
//!
//! - [`SubgridRegion::snapshot_fine_e_pre_update`] — captures fine
//!   `E_t` immediately after the Q3 Dirichlet
//!   [`SubgridRegion::interpolate_coarse_e_to_fine`] writes the outer
//!   layer and **before** [`SubgridRegion::update_fine_e`] runs.
//! - [`SubgridRegion::snapshot_fine_e_post_update`] — captures fine
//!   `E_t` immediately **after** [`SubgridRegion::update_fine_e`]
//!   advances the fine grid by one fine sub-step.
//!
//! At Step C1 the snapshots are populated but not yet consumed; Step C2
//! wires the difference `E_post − E_pre` into the compensating M source
//! `M = -n̂ × (E_post − E_pre)`. These tests guard the time-level
//! semantics so the C2 wiring inherits a known-good snapshot pair.
//!
//! Spec: `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
//! Plan: `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-y-m-coupling.md` Step C1.
//! ADR:  `docs/src/decisions/0038-berenger-m-coupling-spec-amendment.md`.

use yee_fdtd::{SubgridRegion, SubgriddedSolver, WalkingSkeletonSolver, YeeGrid};

const N: usize = 16;
const DX: f64 = 1.0e-3;
const SG_LO: (usize, usize, usize) = (4, 4, 4);
const SG_HI: (usize, usize, usize) = (10, 10, 10);

/// Build a fresh `N³` vacuum parent grid.
fn parent() -> YeeGrid {
    YeeGrid::vacuum(N, N, N, DX)
}

/// Drive every coarse `E` array to a non-trivial position-dependent
/// pattern so the Q3 Dirichlet interpolation produces a non-zero fine
/// `E_t`, and the subsequent `update_fine_e` has something curl-of-`H`
/// to advance against (we seed a matching `H` perturbation on the fine
/// grid below).
fn seed_parent_e(grid: &mut YeeGrid) {
    let nx = grid.ex.shape()[0];
    let ny = grid.ex.shape()[1];
    let nz = grid.ex.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                grid.ex[[i, j, k]] = 0.10 + 0.01 * (i as f64) + 0.003 * (j as f64);
            }
        }
    }
    let nx = grid.ey.shape()[0];
    let ny = grid.ey.shape()[1];
    let nz = grid.ey.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                grid.ey[[i, j, k]] = 0.20 + 0.007 * (j as f64);
            }
        }
    }
    let nx = grid.ez.shape()[0];
    let ny = grid.ez.shape()[1];
    let nz = grid.ez.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                grid.ez[[i, j, k]] = 0.30 + 0.005 * (k as f64);
            }
        }
    }
}

/// Seed the fine grid's `H` field with a non-zero pattern so that
/// `update_fine_e` (whose update is curl(H) ∝ ΔH) actually moves the
/// fine `E` between the pre- and post-update snapshots. Without this,
/// the fine grid starts entirely zero except on the Dirichlet outer
/// layer and `update_fine_e`'s contribution is dominated by the
/// boundary rather than the interior — still non-zero in general, but
/// `seed_fine_h` makes the differential signal robust.
fn seed_fine_h(region: &mut SubgridRegion) {
    let g = region.fine_grid_mut();
    let nx = g.hx.shape()[0];
    let ny = g.hx.shape()[1];
    let nz = g.hx.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                g.hx[[i, j, k]] = 1.0e-4 * ((i + 1) as f64) * 0.5;
            }
        }
    }
    let nx = g.hy.shape()[0];
    let ny = g.hy.shape()[1];
    let nz = g.hy.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                g.hy[[i, j, k]] = 1.0e-4 * ((j + 1) as f64) * 0.5;
            }
        }
    }
    let nx = g.hz.shape()[0];
    let ny = g.hz.shape()[1];
    let nz = g.hz.shape()[2];
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                g.hz[[i, j, k]] = 1.0e-4 * ((k + 1) as f64) * 0.5;
            }
        }
    }
}

// -----------------------------------------------------------------------
// Test 1 — `snapshot_fine_e_pre_update` captures the post-Q3 Dirichlet
// state bit-for-bit.
// -----------------------------------------------------------------------

#[test]
fn pre_snapshot_captures_post_q3_dirichlet() {
    let mut p = parent();
    seed_parent_e(&mut p);

    let mut region = SubgridRegion::new(&p, SG_LO, SG_HI).expect("valid subgrid bounds");

    // Q3 surface pair: start- and end-of-coarse-step snapshots so the
    // temporal blend in `interpolate_coarse_e_to_fine` is well-defined.
    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.75);

    region.snapshot_fine_e_pre_update();

    let (pre_ex, pre_ey, pre_ez) = region
        .fine_e_pre_snapshot()
        .expect("pre snapshot populated by call above");
    let f = region.fine_grid();

    assert_eq!(pre_ex, &f.ex, "pre-snapshot ex must clone the fine ex");
    assert_eq!(pre_ey, &f.ey, "pre-snapshot ey must clone the fine ey");
    assert_eq!(pre_ez, &f.ez, "pre-snapshot ez must clone the fine ez");

    // The Q3 Dirichlet must have written *something* non-trivial onto
    // the outer fine `E_t` — guards against the trivial all-zero
    // failure mode in which the snapshot match is meaningless.
    let any_nonzero = pre_ex.iter().any(|&v| v != 0.0)
        || pre_ey.iter().any(|&v| v != 0.0)
        || pre_ez.iter().any(|&v| v != 0.0);
    assert!(
        any_nonzero,
        "Q3 Dirichlet should have left the fine E_t non-zero"
    );
}

// -----------------------------------------------------------------------
// Test 2 — `snapshot_fine_e_post_update` captures the post-update fine
// `E` state bit-for-bit.
// -----------------------------------------------------------------------

#[test]
fn post_snapshot_captures_post_update_fine_e() {
    let mut p = parent();
    seed_parent_e(&mut p);

    let mut region = SubgridRegion::new(&p, SG_LO, SG_HI).expect("valid subgrid bounds");
    seed_fine_h(&mut region);

    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.75);
    region.update_fine_h();
    region.update_fine_e();

    region.snapshot_fine_e_post_update();

    let (post_ex, post_ey, post_ez) = region
        .fine_e_post_snapshot()
        .expect("post snapshot populated by call above");
    let f = region.fine_grid();

    assert_eq!(post_ex, &f.ex, "post-snapshot ex must clone the fine ex");
    assert_eq!(post_ey, &f.ey, "post-snapshot ey must clone the fine ey");
    assert_eq!(post_ez, &f.ez, "post-snapshot ez must clone the fine ez");
}

// -----------------------------------------------------------------------
// Test 3 — pre and post snapshots differ by exactly the `update_fine_e`
// delta, bounded but non-zero somewhere.
// -----------------------------------------------------------------------

#[test]
fn pre_post_differ_by_update_e_delta() {
    let mut p = parent();
    seed_parent_e(&mut p);

    let mut region = SubgridRegion::new(&p, SG_LO, SG_HI).expect("valid subgrid bounds");
    seed_fine_h(&mut region);

    region.snapshot_coarse_e_t(&p);
    region.snapshot_coarse_e_t_end(&p);
    region.interpolate_coarse_e_to_fine(0.75);

    // Capture pre-update snapshot, then run the H-then-E half-step.
    region.snapshot_fine_e_pre_update();
    region.update_fine_h();
    region.update_fine_e();
    region.snapshot_fine_e_post_update();

    let (pre_ex, pre_ey, pre_ez) = region
        .fine_e_pre_snapshot()
        .expect("pre snapshot populated");
    let (post_ex, post_ey, post_ez) = region
        .fine_e_post_snapshot()
        .expect("post snapshot populated");

    // Component-wise max-abs delta.
    let max_d_ex = pre_ex
        .iter()
        .zip(post_ex.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let max_d_ey = pre_ey
        .iter()
        .zip(post_ey.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    let max_d_ez = pre_ez
        .iter()
        .zip(post_ez.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);

    let max_d = max_d_ex.max(max_d_ey).max(max_d_ez);

    // Non-zero somewhere: `update_fine_e` (curl(H) ∝ ΔH) must have
    // moved at least one fine `E` cell given the seeded fine `H`
    // perturbation.
    assert!(
        max_d > f64::EPSILON,
        "expected non-zero |E_post − E_pre|; saw max_d = {max_d}"
    );

    // Bounded above: the fine sub-step's E update is `dt_fine / ε_0 ·
    // curl(H)`; for the seeded `|H| ~ 1e-3` and `dt_fine ~ DX/(2·c)
    // ~ 1.7e-12 s`, the per-step `|ΔE|` is bounded by `~ dt_fine / ε_0
    // · |H| / dx_fine ~ 4e3` worst-case. A generous `1e6` bound
    // captures any reasonable seeded amplitude and still rejects the
    // catastrophic all-cells-blown-up failure mode.
    assert!(
        max_d < 1.0e6,
        "expected bounded |E_post − E_pre|; saw max_d = {max_d}"
    );
}

// -----------------------------------------------------------------------
// Test 4 — `SubgriddedSolver::step` populates both snapshots in one
// pipeline pass.
// -----------------------------------------------------------------------

#[test]
fn step_pipeline_calls_both_snapshots() {
    let mut p = parent();
    seed_parent_e(&mut p);

    let region = SubgridRegion::new(&p, SG_LO, SG_HI).expect("valid subgrid bounds");
    let inner = WalkingSkeletonSolver::new(p);
    let mut sub = SubgriddedSolver::new(inner).with_region(region);

    // Before the first step neither snapshot should exist.
    {
        let r = sub.region().expect("region attached");
        assert!(
            r.fine_e_pre_snapshot().is_none(),
            "pre snapshot must be None before any step"
        );
        assert!(
            r.fine_e_post_snapshot().is_none(),
            "post snapshot must be None before any step"
        );
    }

    sub.step();

    let r = sub.region().expect("region attached");
    assert!(
        r.fine_e_pre_snapshot().is_some(),
        "step() must have called snapshot_fine_e_pre_update"
    );
    assert!(
        r.fine_e_post_snapshot().is_some(),
        "step() must have called snapshot_fine_e_post_update"
    );
}
