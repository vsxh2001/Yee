//! TF/SF plane-wave-source integration test (Phase 2.fdtd.5).
//!
//! Verifies that a normally-incident plane wave injected by
//! [`yee_fdtd::PlaneWaveSource`] propagates inside the total-field
//! region and is suppressed in the scattered-field region (the TF/SF
//! "quiet zone" property).
//!
//! Setup:
//!
//! - 60³ vacuum grid with `dx = 5 mm` (`dt ≈ 0.9 × Courant`).
//! - TF region: `[15..=45, 15..=45, 15..=45]`.
//! - Source: `+x` propagating, `E_z` polarized, sinusoid at `f = 3 GHz`
//!   (free-space `λ ≈ 100 mm = 20 cells`) with a Hann ramp over the
//!   first 40 time steps.
//! - CPML absorbing outer boundaries with `npml = 8` — without this,
//!   any small TF/SF leakage bounces off the PEC walls and accumulates
//!   into the SF region's amplitude, swamping the "quiet zone" check.
//!
//! After 400 steps we sample `|E_z|` inside the TF box (interior) and
//! outside it on the source side. The TF region carries the incident
//! wave; the SF region should be near-quiet. A contrast of ≥ 10× is
//! required by Phase 2.fdtd.5; well-tuned TF/SF implementations hit
//! 100–1000× (40–60 dB).
//!
//! Marked `#[ignore]` because the run takes ~10 s in release.

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
#[ignore = "slow: ~10s for 60^3 x 400 steps"]
fn plane_wave_propagates_with_tfsf_quiet_outside_box() {
    const N: usize = 60;
    const DX: f64 = 5.0e-3;
    const N_STEPS: usize = 400;
    const FREQ_HZ: f64 = 3.0e9; // λ ≈ 100 mm = 20 cells in vacuum
    const RAMP: usize = 40;
    const PAD: usize = 8;

    // TF region: a slab spanning the full grid in y and z (so the only
    // TF/SF interfaces are the i0 / i1 faces), with x bounds well inside
    // the CPML interior. Phase 2.fdtd.5 only implements i0 / i1 corrections
    // (see [`PlaneWaveSource`] docs) — making the TF region full-width in
    // y and z removes the finite-extent side-face artefacts that a 3D
    // TF/SF box would otherwise emit at the j0 / j1 / k0 / k1 jumps.
    const I0: usize = 15;
    const I1: usize = 45;
    const J0: usize = 0;
    const J1: usize = N;
    const K0: usize = 0;
    const K1: usize = N;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;
    let cpml_params = CpmlParams::for_grid(&grid, 8);
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

    // Inside TF region (clear of the i0/i1 stencil layers and the CPML
    // depth in y,z so we see a clean plane wave).
    const CPML_INTERIOR: usize = 10; // npml=8 + 2 cells of margin
    let inside_amp = max_abs_ez_in_region(
        solver.grid(),
        (I0 + 5)..(I1 - 4),
        CPML_INTERIOR..(N - CPML_INTERIOR),
        CPML_INTERIOR..(N - CPML_INTERIOR),
    );

    // Outside TF region: between the low-x CPML inner edge and the TF
    // front face. This is the SF "quiet zone" the test verifies.
    let outside_amp = max_abs_ez_in_region(
        solver.grid(),
        CPML_INTERIOR..(I0 - 1),
        CPML_INTERIOR..(N - CPML_INTERIOR),
        CPML_INTERIOR..(N - CPML_INTERIOR),
    );

    let contrast = inside_amp / outside_amp.max(1e-30);
    eprintln!("inside (TF)  max |E_z| = {inside_amp:.6e}");
    eprintln!("outside (SF) max |E_z| = {outside_amp:.6e}");
    eprintln!("contrast (inside/outside) = {contrast:.2}");
    eprintln!(
        "contrast (dB) = {:.2}",
        20.0 * contrast.log10().max(-1000.0)
    );

    assert!(
        inside_amp > 0.5,
        "expected TF region to carry incident wave, got {inside_amp}"
    );
    assert!(
        outside_amp < 0.01,
        "expected SF region to be quiet, got {outside_amp}"
    );
    // 1000× gate (Phase 2.fdtd.5.2 DoD). The slab configuration is
    // CPML-bounded on the j/k faces; the only sources of SF amplitude
    // are i-face round-off and CPML reflection. Empirical value
    // ~2676× (~68 dB). The 1000× threshold guards against a
    // regression in either the i-face TF/SF kernel or the CPML.
    assert!(
        contrast > 1000.0,
        "slab TF/SF contrast {contrast:.2} too low (expected > 1000)"
    );
}
