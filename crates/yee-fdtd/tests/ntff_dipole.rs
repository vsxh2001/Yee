//! Phase 2.fdtd.2 integration test for the near-to-far-field
//! transformation.
//!
//! Drives a centred z-polarised Gaussian source on `E_z` of a 50³
//! vacuum grid with CPML on every face. The [`NtffState`] DFT
//! accumulator records the equivalent surface currents at
//! `f = 15 GHz` (λ ≈ 20 mm) over a 2000-step run.
//!
//! Once the FDTD loop finishes we project the accumulated currents to
//! two observation directions:
//!
//! - **broadside** at `θ = π/2, φ = 0` (in the dipole's equatorial
//!   plane);
//! - **endfire** at `θ = 0` (along the dipole axis).
//!
//! For an `E_z`-polarised source the far-field amplitude varies as
//! `sin θ`, so broadside is the maximum and endfire is the null. We
//! expect the broadside/endfire ratio to be ≥ 20 dB.
//!
//! Wall-time budget on a single core (release): < 60 s.

use std::f64::consts::FRAC_PI_2;

use yee_fdtd::{CpmlParams, NtffParams, NtffState, WalkingSkeletonSolver, YeeGrid};

#[test]
#[ignore = "heavy solver test, release-gated in CI (fdtd-heavy-validation-gate); \
            skipped in the default debug `cargo test --workspace` which would time out \
            (CLAUDE.md §10). Run via `cargo test -p yee-fdtd --release --test ntff_dipole \
            -- --ignored ntff_recovers_dipole_pattern_broadside`."]
fn ntff_recovers_dipole_pattern_broadside() {
    const N: usize = 50;
    const DX: f64 = 1.0e-3; // 1 mm
    const NPML: usize = 10;
    const N_STEPS: usize = 2000;
    const F_PROBE: f64 = 15.0e9; // 15 GHz → λ ≈ 20 mm
    // Source at the centre of the grid, on E_z.
    const SRC: (usize, usize, usize) = (25, 25, 25);
    // Integration surface 5 cells inside the inner PML edge → margin
    // npml + 5 from the outer face. Surface bounds: [15, 35] on every
    // axis → 21×21×21 box (well outside the centred source).
    const BOX_MARGIN_CELLS: usize = NPML + 5;

    let grid = YeeGrid::vacuum(N, N, N, DX);
    let dt = grid.dt;

    // Soft Gaussian-in-time source on E_z. The Gaussian's spectrum is
    // also a Gaussian; we want F_PROBE = 15 GHz to be inside its main
    // lobe but not at the spectral peak (which would be DC for an
    // unmodulated Gaussian). Choosing sigma so that the spectral
    // amplitude at 15 GHz is ≈ exp(-(ω σ)²/2) of the DC peak; for
    // sigma = 4·dt ≈ 6.9 ps, ω σ = 2π·15 GHz · 6.9 ps ≈ 0.65 → spectral
    // attenuation ≈ exp(-0.21) ≈ 0.81 at 15 GHz. Plenty of signal.
    let t0 = 12.0 * dt;
    let sigma = 4.0 * dt;

    let params = CpmlParams::for_grid(&grid, NPML);
    let mut solver = WalkingSkeletonSolver::with_cpml(grid, params);

    let ntff_params = NtffParams {
        f_probe: F_PROBE,
        box_margin_cells: BOX_MARGIN_CELLS,
        // Placeholder; we sweep below via far_field_at.
        theta_rad: FRAC_PI_2,
        phi_rad: 0.0,
    };
    let mut ntff = NtffState::new(solver.grid(), ntff_params);

    for _ in 0..N_STEPS {
        solver.step_with_source_and_ntff(SRC.0, SRC.1, SRC.2, t0, sigma, &mut ntff);
    }

    // Sanity: NTFF saw every step.
    assert_eq!(ntff.n_samples(), N_STEPS as u64);

    // Project at two directions for the E_z-polarised "dipole":
    //   broadside: θ = π/2, φ = 0   (xy plane)
    //   endfire : θ = 0             (+z axis)
    let e_broadside = ntff.far_field_at(FRAC_PI_2, 0.0);
    let e_endfire = ntff.far_field_at(0.0, 0.0);

    let mag_broad = e_broadside.norm();
    let mag_end = e_endfire.norm();

    eprintln!("|E_far(broadside, θ=π/2, φ=0)| = {mag_broad:.3e}");
    eprintln!("|E_far(endfire,  θ=0)|         = {mag_end:.3e}");
    assert!(mag_broad.is_finite(), "broadside non-finite");
    assert!(mag_end.is_finite(), "endfire non-finite");
    assert!(
        mag_broad > 0.0,
        "broadside is zero — no radiation captured?"
    );

    let ratio = if mag_end > 0.0 {
        mag_broad / mag_end
    } else {
        f64::INFINITY
    };
    let db = 20.0 * ratio.log10();
    eprintln!("broadside / endfire = {ratio:.3e}  ({db:.2} dB)");

    // Theoretical short-dipole null at θ = 0 is infinite. With a
    // coarse grid, finite CPML, and a single-frequency probe of a
    // broadband Gaussian we expect ≥ 20 dB.
    assert!(
        db >= 20.0,
        "dipole broadside/endfire = {db:.2} dB is below the 20 dB target",
    );
}
