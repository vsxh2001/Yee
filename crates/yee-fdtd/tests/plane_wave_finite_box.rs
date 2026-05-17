//! TF/SF finite-box (not slab) contrast validation (Phase 2.fdtd.5.2).
//!
//! Phase 2.fdtd.5 (Track OO) shipped the `+x` `E_z` polarized TF/SF
//! plane-wave source with the TF region spanning the full `j`, `k`
//! extent of the grid (a slab), and the slab case measures
//! ~2676× (~68 dB) inside/outside contrast.
//!
//! This test exercises the same kernel on a **finite** TF box
//! (bounded on all six faces). Naïve polarization analysis says the
//! `j0`/`j1`/`k0`/`k1` faces should need no corrections at normal
//! `+x` `E_z` incidence (since `H_inc_x = H_inc_z = E_inc_x = E_inc_y = 0`),
//! but in fact two more pairs of corrections are needed:
//!
//! - `H_x` at the j-faces (`j = j0 - 1` and `j = j1`) needs the
//!   `−∂E_z/∂y` term corrected for the TF/SF discontinuity in `E_z`
//!   across the j-face. This is the same physics as the existing
//!   i-face `H_y` correction, applied to a different curl component.
//! - `E_x` at the k-faces (`k = k0` and `k = k1 + 1`) needs the
//!   `−∂H_y/∂z` term corrected for the TF/SF discontinuity in `H_y`
//!   across the k-face.
//!
//! With those four extra face corrections in place (Phase 2.fdtd.5.2),
//! the finite-box configuration measures contrast at machine
//! precision — the SF "quiet zone" outside the TF box drops to
//! ~1e-15 (the same `f64` roundoff floor the slab case bottoms out
//! at if CPML weren't the dominant residual).
//!
//! This test pins the finite-box contrast at ≥ 100× — that's the
//! Phase 2.fdtd.5.2 deliverable gate. The empirical value with all
//! six-face corrections in place is ~7×10¹⁴ (~298 dB; effectively
//! roundoff-limited), so the 100× gate is a *very* loose guardrail
//! whose only purpose is to catch a regression that re-introduces
//! incident-leakage at the j/k faces.
//!
//! Setup mirrors `plane_wave_propagation.rs` but with a tightly
//! bounded TF box in `y` and `z`.

use std::ops::Range;

use yee_fdtd::{CpmlParams, PlaneWaveDirection, PlaneWaveSource, WalkingSkeletonSolver, YeeGrid};

fn max_abs_ez_in_region(
    grid: &YeeGrid,
    is: Range<usize>,
    js: Range<usize>,
    ks: Range<usize>,
) -> f64 {
    let mut peak: f64 = 0.0;
    for i in is.clone() {
        for j in js.clone() {
            for k in ks.clone() {
                let v = grid.ez[(i, j, k)].abs();
                if v > peak {
                    peak = v;
                }
            }
        }
    }
    peak
}

#[test]
#[ignore = "slow: ~10s for 80^3 x 600 steps"]
fn finite_box_tfsf_preserves_contrast() {
    // 80³ vacuum grid; finite TF box that is bounded in *all three*
    // axes — i.e. not a j/k-spanning slab. The box is centred so that
    // CPML, TF-front, TF-interior, TF-back, and SF "outside-the-box"
    // regions are all comfortably separated.
    const N: usize = 80;
    const DX: f64 = 5.0e-3;
    const N_STEPS: usize = 600;
    const FREQ_HZ: f64 = 3.0e9; // λ ≈ 100 mm = 20 cells
    const RAMP: usize = 40;
    const PAD: usize = 8;
    const NPML: usize = 8;

    // Finite TF box: bounded on all six faces.
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
    let mut pw = PlaneWaveSource::new(
        I0,
        I1,
        J0,
        J1,
        K0,
        K1,
        PlaneWaveDirection::PlusX,
        FREQ_HZ,
        RAMP,
        DX,
        dt,
        PAD,
    );

    for _ in 0..N_STEPS {
        solver.step_with_plane_wave(&mut pw);
    }

    // Sanity: nothing exploded.
    for k in 0..N {
        for j in 0..=N {
            for i in 0..=N {
                let v = solver.grid().ez[(i, j, k)];
                assert!(v.is_finite(), "E_z went non-finite at ({i},{j},{k})");
            }
        }
    }

    // Inside TF region: well clear of the i0/i1/j0/j1/k0/k1 stencil
    // layers so we measure the developed incident wave rather than the
    // immediate-face correction artefacts.
    let inside_amp = max_abs_ez_in_region(
        solver.grid(),
        (I0 + 5)..(I1 - 4),
        (J0 + 5)..(J1 - 4),
        (K0 + 5)..(K1 - 4),
    );

    // Outside TF region (SF "quiet zone"): sampled in the source-side
    // SF region between the low-x CPML inner edge and the TF front
    // face, with j, k slabs constrained to lie inside the j, k extent
    // of the TF box. This is the cleanest SF zone — it sits in the
    // wave's direct path (so any spurious through-leakage would show
    // up here) without crossing the j/k side faces of the TF box.
    const CPML_INTERIOR: usize = NPML + 2;
    let outside_amp =
        max_abs_ez_in_region(solver.grid(), CPML_INTERIOR..(I0 - 1), J0..(J1 + 1), K0..K1);

    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!("finite-box inside  max |E_z| = {inside_amp:.6e}");
    eprintln!("finite-box outside max |E_z| = {outside_amp:.6e}");
    eprintln!("finite-box contrast (inside/outside) = {contrast:.2}");
    eprintln!(
        "finite-box contrast (dB)             = {:.2}",
        20.0 * contrast.log10().max(-1000.0)
    );

    assert!(
        inside_amp > 0.5,
        "expected TF region to carry incident wave, got {inside_amp}"
    );
    // Finite-box outside-amplitude bound. With the Phase 2.fdtd.5.2
    // j/k-face corrections in place the empirical value is at
    // roundoff (~1e-15); we leave the 0.01 guardrail so a regression
    // that takes the outside amplitude back into the visible range
    // is caught loudly.
    assert!(
        outside_amp < 0.01,
        "expected SF region to be quiet, got {outside_amp}"
    );
    // 100× gate from the Phase 2.fdtd.5.2 brief. Empirical value with
    // all six-face corrections in place is ~7×10¹⁴ (effectively
    // floating-point roundoff-limited). The 100× threshold is a
    // *very* loose guardrail whose only job is to detect a regression
    // that re-introduces incident-leakage at the j- or k-faces.
    assert!(
        contrast > 100.0,
        "finite-box TF/SF contrast {contrast:.2} too low (expected > 100)"
    );
}
