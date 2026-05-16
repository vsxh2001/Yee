//! Phase 2.fdtd.4 end-to-end driver validation: short-dipole pattern.
//!
//! Drives the [`FdtdDriver`] with a z-polarised Hann-windowed sinusoidal
//! `J_z` source distributed over a few cells along z, runs the FDTD
//! loop with CPML on all six faces, and validates the resulting NTFF
//! far-field pattern against the analytic short-dipole form
//! `|E_θ| ∝ sin θ`.
//!
//! Tolerances are loose because:
//!
//! - The grid is coarse (60³, dx = 5 mm at 1 GHz → λ/dx ≈ 60).
//! - The NTFF integrates over a finite box so geometric symmetry is only
//!   approximate.
//! - The Hann ramp produces a transient at startup that the DFT bin
//!   still averages over.
//!
//! Wall-time budget: < 60 s single-threaded release build. Marked
//! `#[ignore]` so the default `cargo test` stays fast — run with
//! `--include-ignored` to exercise.
//!
//! References:
//!
//! > C. A. Balanis, *Antenna Theory: Analysis and Design*, 4th ed.,
//! > Wiley, 2016, §4.2 (infinitesimal dipole far-field
//! > `E_θ = j η₀ k I₀ ℓ e^{−jkr}/(4πr) · sin θ`).

use yee_fdtd::{FdtdDriver, FdtdDriverConfig, YeeGrid};

#[test]
#[ignore = "slow: ~30s for 60³ grid × 800 steps"]
fn short_dipole_radiates_sin_theta_pattern() {
    let grid = YeeGrid::vacuum(60, 60, 60, 5.0e-3);
    let cfg = FdtdDriverConfig {
        n_steps: 800,
        dipole_center_cells: (30, 30, 30),
        dipole_length_cells: 5,
        source_freq_hz: 1.0e9,
        ntff_surface_pad_cells: 4,
        cpml_thickness_cells: 10,
    };
    let pattern = FdtdDriver::new(grid, cfg).run();

    // For a short dipole along z, |E_θ| ∝ sin θ.
    // Sample at 0°, 45°, 90°, 135°, 180°.
    let find = |deg: f64| {
        *pattern
            .theta_deg
            .iter()
            .zip(&pattern.e_theta_phi0)
            .min_by(|a, b| (a.0 - deg).abs().partial_cmp(&(b.0 - deg).abs()).unwrap())
            .unwrap()
            .1
    };

    let e_0 = find(0.0);
    let e_45 = find(45.0);
    let e_90 = find(90.0);
    let e_135 = find(135.0);
    let e_180 = find(180.0);

    eprintln!("|E_θ|(  0°) = {e_0:.4}");
    eprintln!("|E_θ|( 45°) = {e_45:.4}");
    eprintln!("|E_θ|( 90°) = {e_90:.4}");
    eprintln!("|E_θ|(135°) = {e_135:.4}");
    eprintln!("|E_θ|(180°) = {e_180:.4}");

    // sin(0) = 0, sin(45) = sin(135) = 0.707, sin(90) = 1, sin(180) = 0.
    // Loose tolerances: FDTD grid is coarse, NTFF integrates over a
    // finite box.
    assert!(e_0 < 0.15, "expected null at θ=0°, got {e_0}");
    assert!(e_180 < 0.15, "expected null at θ=180°, got {e_180}");
    assert!(
        (e_90 - 1.0).abs() < 0.05,
        "expected peak at θ=90°, got {e_90}"
    );
    assert!(
        (e_45 - 0.707).abs() < 0.15,
        "expected ~0.707 at θ=45°, got {e_45}"
    );
    assert!(
        (e_135 - 0.707).abs() < 0.15,
        "expected ~0.707 at θ=135°, got {e_135}"
    );
}
