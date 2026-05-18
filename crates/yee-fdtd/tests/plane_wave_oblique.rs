//! Oblique-incidence TF/SF validation (Phase 2.fdtd.5.3 / 2.fdtd.5.3.1).
//!
//! Phase 2.fdtd.5.0/5.1/5.2 shipped the `+x` `E_z` polarized TF/SF
//! source (a normal-incidence special case). Phase 2.fdtd.5.3 lifted
//! that restriction to arbitrary `(θ, φ, ψ)` via a 1-D auxiliary
//! incident-field grid along `k_hat` plus per-face vector projection
//! onto the box-stencil Yee nodes (see the in-source design notes at
//! the top of `crates/yee-fdtd/src/sources.rs`). Phase 2.fdtd.5.3.1
//! added Taflove §5.10.5 dispersion matching of the 1-D aux step so
//! the 1-D Yee phase velocity tracks the 3-D Yee phase velocity along
//! `k_hat` — this raises 30°/45° finite-box contrast from ~14.5× to
//! >1000×, hitting the original 5.3 DoD.
//!
//! Four cases:
//!
//! 1. **Normal-incidence regression** — `(θ=π/2, φ=0, ψ=π)` (the
//!    `(k̂=+x̂, E_inc_hat=+ẑ)` mapping of the 5.2 setup) constructed
//!    via `with_oblique_incidence` (NOT the legacy `new` path) must
//!    reach within 1% of the 5.2 finite-box contrast floor
//!    (≥ 1e10× per the brief; the empirical value is ~9.8e14×, only
//!    1.5× below the legacy 5.2 floor). For on-axis propagation the
//!    1-D / 3-D Yee dispersion relations coincide exactly, so the
//!    dispersion-matched `ds_aux` collapses to `dx` and this test is
//!    indifferent to the 5.3.1 change.
//! 2. **Oblique sanity (5.3 reproduction)** — `θ = 30°, φ = 45°,
//!    ψ = π/2` with `dispersion_match = false`. Must reproduce the
//!    Phase 2.fdtd.5.3 ~14.5× contrast within 5%, verifying the new
//!    code path doesn't change the old (no-match) result.
//! 3. **Oblique DoD (5.3.1)** — same angles with the default
//!    `with_oblique_incidence` (dispersion match enabled). Must
//!    clear the brief's 1000× DoD.
//! 4. **Grazing rejection** — `θ = 85°` must complete without
//!    NaN / panic. Contrast is expected to degrade.
//!
//! Test format mirrors `plane_wave_finite_box.rs`.

use std::ops::Range;

use yee_fdtd::{CpmlParams, PlaneWaveSource, WalkingSkeletonSolver, YeeGrid};

fn max_abs_field_in_region(
    grid: &YeeGrid,
    is: Range<usize>,
    js: Range<usize>,
    ks: Range<usize>,
) -> f64 {
    let mut peak: f64 = 0.0;
    for i in is.clone() {
        for j in js.clone() {
            for k in ks.clone() {
                // Examine all three E components at this cell index.
                // The arrays have different shapes; clamp via bounds
                // checks. We only care about the maximum magnitude, so
                // any out-of-bounds index is silently skipped.
                if i < grid.ex.shape()[0] && j < grid.ex.shape()[1] && k < grid.ex.shape()[2] {
                    let v = grid.ex[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
                if i < grid.ey.shape()[0] && j < grid.ey.shape()[1] && k < grid.ey.shape()[2] {
                    let v = grid.ey[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
                if i < grid.ez.shape()[0] && j < grid.ez.shape()[1] && k < grid.ez.shape()[2] {
                    let v = grid.ez[(i, j, k)].abs();
                    if v > peak {
                        peak = v;
                    }
                }
            }
        }
    }
    peak
}

fn assert_all_finite(grid: &YeeGrid, label: &str) {
    for (name, arr) in [
        ("ex", &grid.ex),
        ("ey", &grid.ey),
        ("ez", &grid.ez),
        ("hx", &grid.hx),
        ("hy", &grid.hy),
        ("hz", &grid.hz),
    ] {
        for v in arr.iter() {
            assert!(v.is_finite(), "{label}: {name} contained non-finite value");
        }
    }
}

#[test]
#[ignore = "slow: ~10s for 80^3 x 600 steps"]
fn oblique_normal_incidence_regression() {
    // 5.2-equivalent normal-incidence angles: (θ=π/2, φ=0, ψ=π) gives
    // k_hat = (1, 0, 0) (i.e. +x̂) and E_inc_hat = +ẑ, matching the
    // 5.2 +x / E_z polarization bit-for-bit. The oblique kernel
    // dispatches through `correct_h_oblique` / `correct_e_oblique`
    // (not the legacy four-face shortcut), so this test verifies the
    // general 12-face stencil collapses to the same physics as the
    // four-face one at exact-on-axis incidence.
    use std::f64::consts::PI;

    const N: usize = 80;
    const DX: f64 = 5.0e-3;
    const N_STEPS: usize = 600;
    const FREQ_HZ: f64 = 3.0e9;
    const RAMP: usize = 40;
    const PAD: usize = 8;
    const NPML: usize = 8;

    const I0: usize = 25;
    const I1: usize = 55;
    const J0: usize = 25;
    const J1: usize = 55;
    const K0: usize = 25;
    const K1: usize = 55;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;
    let cpml_params = CpmlParams::for_grid(&grid, NPML);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, cpml_params);

    let mut pw = PlaneWaveSource::with_oblique_incidence(
        I0,
        I1,
        J0,
        J1,
        K0,
        K1,
        PI / 2.0,
        0.0,
        PI, // θ=π/2, φ=0, ψ=π → k=+x̂, E=+ẑ (5.2 setup)
        FREQ_HZ,
        RAMP,
        DX,
        dt,
        PAD,
    );
    assert!(
        !pw.is_legacy_normal_incidence(),
        "oblique constructor must NOT mark source as legacy"
    );

    for _ in 0..N_STEPS {
        solver.step_with_plane_wave(&mut pw);
    }
    assert_all_finite(solver.grid(), "oblique-normal");

    // Inside TF: well clear of the box-face stencil layers. The wave
    // here is +x propagating with E_z polarization, so |E_z| should
    // dominate.
    let inside_amp = max_abs_field_in_region(
        solver.grid(),
        (I0 + 5)..(I1 - 4),
        (J0 + 5)..(J1 - 4),
        (K0 + 5)..(K1 - 4),
    );
    assert!(
        inside_amp > 0.5,
        "expected TF region to carry incident wave, got {inside_amp}"
    );

    // SF "quiet zone" — between the low-x CPML inner edge and the TF
    // front x-face. Cross-section clipped to the TF box y, z extent so
    // we measure direct-through (not side-face) leakage.
    const CPML_INTERIOR: usize = NPML + 2;
    let outside_amp = max_abs_field_in_region(
        solver.grid(),
        CPML_INTERIOR..(I0 - 1),
        J0..(J1 + 1),
        K0..(K1 + 1),
    );
    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!(
        "oblique-normal inside  max |E| = {inside_amp:.6e}\n\
         oblique-normal outside max |E| = {outside_amp:.6e}\n\
         oblique-normal contrast        = {contrast:.6e} ({:.2} dB)",
        20.0 * contrast.log10().max(-1000.0)
    );

    // 1% of the 5.2 floor (~7e14× → 7e12×). The brief calls for "still
    // > 1e10×" as the explicit gate; we hold to that.
    assert!(
        contrast > 1.0e10,
        "oblique-normal-incidence contrast {contrast:.2e} too low \
         (expected > 1e10× to satisfy the 1% Phase 2.fdtd.5.3 DoD)"
    );
}

/// Common 30°/45° geometry used by both the sanity (no-match) and the
/// DoD (dispersion-match) tests. Returns `(inside_amp, outside_amp)`
/// for the configured `dispersion_match` flag.
fn run_oblique_30_45_case(dispersion_match: bool) -> (f64, f64) {
    use std::f64::consts::PI;

    const N: usize = 60;
    const DX: f64 = 5.0e-3;
    const N_STEPS: usize = 400;
    const FREQ_HZ: f64 = 3.0e9;
    const RAMP: usize = 40;
    const PAD: usize = 8;
    const NPML: usize = 8;

    const I0: usize = 20;
    const I1: usize = 40;
    const J0: usize = 20;
    const J1: usize = 40;
    const K0: usize = 20;
    const K1: usize = 40;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;
    let cpml_params = CpmlParams::for_grid(&grid, NPML);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, cpml_params);

    let theta = 30.0_f64.to_radians();
    let phi = 45.0_f64.to_radians();
    let psi = PI / 2.0; // E along e_phi
    let mut pw = PlaneWaveSource::with_oblique_incidence_match(
        I0,
        I1,
        J0,
        J1,
        K0,
        K1,
        theta,
        phi,
        psi,
        FREQ_HZ,
        RAMP,
        DX,
        dt,
        PAD,
        dispersion_match,
    );

    for _ in 0..N_STEPS {
        solver.step_with_plane_wave(&mut pw);
    }
    assert_all_finite(
        solver.grid(),
        if dispersion_match {
            "oblique-30-45-match"
        } else {
            "oblique-30-45-nomatch"
        },
    );

    let inside_amp = max_abs_field_in_region(
        solver.grid(),
        (I0 + 4)..(I1 - 3),
        (J0 + 4)..(J1 - 3),
        (K0 + 4)..(K1 - 3),
    );
    assert!(
        inside_amp > 0.1,
        "expected oblique TF region to carry incident wave, got {inside_amp}"
    );

    const CPML_INTERIOR: usize = NPML + 2;
    let sf_lo_x = max_abs_field_in_region(
        solver.grid(),
        CPML_INTERIOR..(I0 - 1),
        J0..(J1 + 1),
        K0..(K1 + 1),
    );
    let sf_lo_y = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        CPML_INTERIOR..(J0 - 1),
        K0..(K1 + 1),
    );
    let sf_lo_z = max_abs_field_in_region(
        solver.grid(),
        I0..(I1 + 1),
        J0..(J1 + 1),
        CPML_INTERIOR..(K0 - 1),
    );
    let outside_amp = sf_lo_x.max(sf_lo_y).max(sf_lo_z);
    eprintln!(
        "oblique-30°/45° (match={dispersion_match}) inside  max |E|       = {inside_amp:.6e}\n\
         oblique-30°/45° (match={dispersion_match}) SF max (lo-x, lo-y, lo-z) = ({sf_lo_x:.3e}, {sf_lo_y:.3e}, {sf_lo_z:.3e})\n\
         oblique-30°/45° (match={dispersion_match}) outside_amp           = {outside_amp:.6e}"
    );

    (inside_amp, outside_amp)
}

#[test]
#[ignore = "slow: ~10s for 60^3 x 400 steps"]
fn oblique_30deg_45deg_no_match_reproduces_phase_2_fdtd_5_3_baseline() {
    // Originally asserted ~14.5× (the empirical Phase 2.fdtd.5.3
    // ship value) with ±50% tolerance. Phase 2.fdtd.5.3.2 upgraded
    // the 1-D auxiliary incident-field interpolation from linear to
    // 4-point cubic Lagrange (commit prior). That upgrade improves
    // *both* match=true and match=false paths, so the no-match
    // contrast now lands at ~340× rather than ~14.5×. The original
    // 14.5× baseline is therefore obsolete; widen the window to a
    // floor-only regression guard that catches a true regression
    // without re-asserting the pre-cubic figure.
    let (inside_amp, outside_amp) = run_oblique_30_45_case(false);
    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!(
        "oblique-30°/45° (no dispersion match) contrast = {contrast:.6e} \
         ({:.2} dB) — post-Phase 2.fdtd.5.3.2 cubic interpolation",
        20.0 * contrast.log10().max(-1000.0)
    );

    // Floor-only guard: contrast must exceed the pre-cubic 14.5×
    // figure by a comfortable margin. The old [7×, 22×] window is
    // gone because cubic interpolation benefits both paths.
    assert!(
        contrast > 50.0,
        "oblique 30°/45° (no match) contrast {contrast:.3} regressed \
         below the post-Phase 2.fdtd.5.3.2 floor of 50× — the cubic \
         interpolation upgrade may have been reverted"
    );
}

#[test]
#[ignore = "slow: ~10s for 60^3 x 400 steps"]
fn oblique_30deg_45deg_ephi_polarization() {
    // θ = 30°, φ = 45°, ψ = π/2 (E along e_phi). k_hat in the
    // first octant (all positive components), so the (i0, j0, k0)
    // corner is upstream and every face-projected distance is ≥ 0.
    //
    // Phase 2.fdtd.5.3.1 enables Taflove §5.10.5 dispersion matching
    // of the 1-D auxiliary-grid step in `with_oblique_incidence` by
    // default. With the matched `ds_aux` the 1-D and 3-D numerical
    // phase velocities agree at the source carrier, and the residual
    // TF/SF leakage drops by ~2 orders of magnitude versus the
    // pre-5.3.1 baseline (which sat at ~14.5×).
    //
    // Phase 2.fdtd.5.3 DoD: >1000×. Phase 2.fdtd.5.3.2 face-stencil
    // audit landed cubic Lagrange interpolation of the 1-D auxiliary
    // incident-field grid (the dominant SF leakage source turned out
    // to be the linear-interpolation residual, not the face stencils
    // themselves). Measured contrast with dispersion match is now
    // ~1027× (60 dB), clearing the original DoD.

    let (inside_amp, outside_amp) = run_oblique_30_45_case(true);
    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!(
        "oblique-30°/45° (dispersion match)    contrast = {contrast:.6e} \
         ({:.2} dB)",
        20.0 * contrast.log10().max(-1000.0)
    );

    assert!(
        contrast > 1000.0,
        "oblique 30°/45° contrast {contrast:.3e} below the Phase 2.fdtd.5.3 \
         DoD of >1000× — Phase 2.fdtd.5.3.2 cubic-interpolation upgrade may \
         have been reverted."
    );
}

#[test]
#[ignore = "slow: ~10s for 60^3 x 300 steps"]
fn oblique_grazing_85deg_runs_without_panic() {
    // θ = 85° (nearly tangential to z = const). Just verify no NaN /
    // panic — contrast is expected to degrade for grazing angles.
    const N: usize = 60;
    const DX: f64 = 5.0e-3;
    const N_STEPS: usize = 300;
    const FREQ_HZ: f64 = 3.0e9;
    const RAMP: usize = 40;
    const PAD: usize = 8;
    const NPML: usize = 8;

    const I0: usize = 20;
    const I1: usize = 40;
    const J0: usize = 20;
    const J1: usize = 40;
    const K0: usize = 20;
    const K1: usize = 40;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;
    let cpml_params = CpmlParams::for_grid(&grid, NPML);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, cpml_params);

    let theta = 85.0_f64.to_radians();
    let phi = 30.0_f64.to_radians();
    let psi = 0.0;
    let mut pw = PlaneWaveSource::with_oblique_incidence(
        I0, I1, J0, J1, K0, K1, theta, phi, psi, FREQ_HZ, RAMP, DX, dt, PAD,
    );

    for _ in 0..N_STEPS {
        solver.step_with_plane_wave(&mut pw);
    }
    assert_all_finite(solver.grid(), "oblique-85");
    eprintln!("oblique-85° survived 300 steps without NaN");
}
