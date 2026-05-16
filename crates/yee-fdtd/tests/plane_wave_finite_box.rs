//! TF/SF finite-box (not slab) contrast validation (Phase 2.fdtd.5.1).
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
//! but the discrete `H_x` stencil at `j = j0 ± 1/2` reads `E_z`
//! across the TF/SF boundary, and `E_inc_z ≠ 0`. The standard
//! Yee update there does not pick up an "incident-H_x" curl term
//! to cancel — but the *scattered* `H_x` it produces from the
//! `E_z` jump is a real (and spurious) artefact of the finite-box
//! geometry that the slab case avoids by terminating those faces
//! in CPML.
//!
//! Empirically (this test), the finite-box configuration achieves
//! ~6× contrast — well above the 5× **prose** gate documented in
//! the Phase 2.fdtd.5.1 brief, but substantially below the slab
//! configuration's 2676×. This confirms the brief's expectation
//! that the corners are "imperfect" without full side-face
//! corrections; the missing physics is the j/k-face `E_z`
//! discontinuity correction described above.
//!
//! Side-face physics corrections — the proper Phase 2.fdtd.5.2
//! deliverable — are deferred. This test pins the present
//! behaviour with the loose 5× bound so any regression in the
//! existing i0/i1 kernel is caught.
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
    // Finite-box outside-amplitude bound is *loose* (0.5) — see module
    // docstring. The slab geometry would pin this at well under 1e-3
    // because all four transverse faces sit in CPML; here, the j/k
    // side faces emit scattered fields the i-only correction kernel
    // cannot cancel. Tightening this requires Phase 2.fdtd.5.2
    // side-face corrections.
    assert!(
        outside_amp < 0.5,
        "expected SF region to be relatively quiet, got {outside_amp}"
    );
    // 5× prose gate from the Phase 2.fdtd.5.1 brief. Empirical value
    // is ~6× (~15.6 dB). This is the pinning test for the existing
    // i0/i1-only correction kernel applied to a finite TF box; any
    // regression in the i-face corrections (or in CPML, which has to
    // hold up the outside-box ambient) will break it.
    assert!(
        contrast > 5.0,
        "finite-box TF/SF contrast {contrast:.2} too low (expected > 5)"
    );
}
