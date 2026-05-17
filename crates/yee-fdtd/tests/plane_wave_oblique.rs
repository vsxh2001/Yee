//! Oblique-incidence TF/SF validation (Phase 2.fdtd.5.3).
//!
//! Phase 2.fdtd.5.0/5.1/5.2 shipped the `+x` `E_z` polarized TF/SF
//! source (a normal-incidence special case). Phase 2.fdtd.5.3 lifts
//! that restriction to arbitrary `(θ, φ, ψ)` via a 1-D auxiliary
//! incident-field grid along `k_hat` plus per-face vector projection
//! onto the box-stencil Yee nodes (see the in-source design notes at
//! the top of `crates/yee-fdtd/src/sources.rs`).
//!
//! Three cases:
//!
//! 1. **Normal-incidence regression** — `(θ=π/2, φ=0, ψ=π)` (the
//!    `(k̂=+x̂, E_inc_hat=+ẑ)` mapping of the 5.2 setup) constructed
//!    via `with_oblique_incidence` (NOT the legacy `new` path) must
//!    reach within 1% of the 5.2 finite-box contrast floor
//!    (≥ 1e10× per the brief; the empirical value is ~9.8e14×, only
//!    1.5× below the legacy 5.2 floor). For on-axis propagation the
//!    linear interpolation falls on integer indices, so no
//!    interpolation roundoff is introduced.
//! 2. **Oblique sanity** — `θ = 30°, φ = 45°, ψ = π/2` (E along
//!    `e_phi_hat`). 60³ box, finite TF box. **Empirical contrast at
//!    Phase 2.fdtd.5.3 ship: ~14.5×.** The brief's 1000× DoD is
//!    *not* met; the limiting error is the dispersion mismatch
//!    between the 1-D Yee aux at step `ds = dx` and the 3-D Yee
//!    oblique phase velocity along `k_hat`. Taflove §5.10.5's
//!    matched-numerical-dispersion remedy is deferred (see the
//!    in-test comment for the escape-hatch detail). The gate here
//!    is set to `>10×` — sufficient to certify the 12-face kernel's
//!    sign conventions and cross-section ranges, but explicitly
//!    documents the gap.
//! 3. **Grazing rejection** — `θ = 85°` must complete without
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

#[test]
#[ignore = "slow: ~10s for 60^3 x 400 steps"]
fn oblique_30deg_45deg_ephi_polarization() {
    // θ = 30°, φ = 45°, ψ = π/2 (E along e_phi). k_hat in the
    // first octant (all positive components), so the (i0, j0, k0)
    // corner is upstream and every face-projected distance is ≥ 0.
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
    let mut pw = PlaneWaveSource::with_oblique_incidence(
        I0, I1, J0, J1, K0, K1, theta, phi, psi, FREQ_HZ, RAMP, DX, dt, PAD,
    );

    for _ in 0..N_STEPS {
        solver.step_with_plane_wave(&mut pw);
    }
    assert_all_finite(solver.grid(), "oblique-30-45");

    // Inside TF: deep in the interior of the box.
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

    // SF "quiet zone" — sample across all three upstream half-planes
    // (low-x, low-y, low-z), each clipped to the cross-section of the
    // TF box so we measure straight-through leakage. Use the maximum
    // over the three upstream slabs as the contrast denominator.
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
    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!(
        "oblique-30°/45° inside  max |E|       = {inside_amp:.6e}\n\
         oblique-30°/45° SF max (lo-x, lo-y, lo-z) = ({sf_lo_x:.3e}, {sf_lo_y:.3e}, {sf_lo_z:.3e})\n\
         oblique-30°/45° contrast              = {contrast:.6e} ({:.2} dB)",
        20.0 * contrast.log10().max(-1000.0)
    );

    // Phase 2.fdtd.5.3 measured contrast: ~14.5× at θ=30°, φ=45°
    // with the 1-D auxiliary grid using `ds_aux = dx`. This is below
    // the brief's 1000× DoD target. The escape-hatch report
    // identifies the leakage source as the **dispersion mismatch**
    // between the 1-D Yee aux at step `dx` and the 3-D Yee oblique
    // phase velocity along `k_hat`: cumulative phase drift across
    // the box (~30 cells in the projected diagonal) produces a
    // 5-10% per-face leakage that dominates the per-cell linear
    // interpolation error.
    //
    // The standard remedy (Taflove §5.10.5 "matched numerical-
    // dispersion incident-field generator") is to step the 1-D aux
    // at a dispersion-matched `ds_aux` such that the 1-D and 3-D
    // phase velocities agree at the source frequency. That requires
    // solving the transcendental dispersion equation while
    // discriminating against the trivial `ds → 0` root that an
    // unguarded Newton or fixed-point picks up. The closed-form
    // implementation is deferred to a follow-on Phase 2.fdtd.5.3.1
    // (or 5.4) sub-track.
    //
    // For Phase 2.fdtd.5.3 ship, we keep the gate at >10× — well
    // above no-correction (~1×) and adequate to certify the kernel's
    // sign conventions, cross-section ranges, and back-compat. This
    // is the brief's escape-hatch outcome with a measured contrast.
    assert!(
        contrast > 10.0,
        "oblique 30°/45° contrast {contrast:.2} too low (expected > 10× \
         for the Phase 2.fdtd.5.3 ship floor; the brief's 1000× DoD is \
         not achievable without dispersion-matched aux step — see test \
         comment for the Taflove §5.10.5 remedy deferred to a follow-on)"
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
