//! Validation gate **fdtd-205** — Ohmic skin-depth penetration profile.
//!
//! # Physics
//!
//! A sinusoidal E_y source drives a wave that propagates in the +x direction
//! through vacuum and into a conducting half-space (x ≥ X_SURFACE) with
//! conductivity σ. In the good-conductor approximation the field amplitude
//! decays as (Griffiths §9.4.1):
//!
//! ```text
//! |E(x)| = |E₀| exp(-(x - x_surface) / δ),   δ = √(2 / (ω μ₀ σ))
//! ```
//!
//! This gate validates the CA/CB Ohmic-loss E-update (Taflove §3.7).
//!
//! # Geometry
//!
//! 80×1×200 grid, dx = 1 mm, f = 1 GHz, σ = 2.533 S/m → δ = 10 mm = 10 cells.
//! Vacuum: x = 0..49; conductor: x = 50..79 (30 cells = 3δ).
//!
//! E_y propagates in the +x direction via the H_z ↔ E_y leapfrog stencil.
//! The z-axis (NZ = 200 cells, 200 mm) has PEC walls at z = 0 and z = NZ.
//! The dominant propagating mode has kz = π/(200 mm) = 15.7 /m, which satisfies
//! kz << 1/δ = 100 /m → the skin-depth in x matches e^{-z/δ} to < 1 %.
//!
//! NY = 1 makes the grid effectively 2-D in x-z. PEC at the y-faces does not
//! zero E_y (E_y is normal to y-faces), so the y-dimension is irrelevant.
//!
//! # Gate criteria
//!
//! Over steps 6 000–7 999 (steady-state measurement window):
//!
//! ```text
//! ratio_1δ = max|E_y(60,0,100)| / max|E_y(50,0,100)|  →  target e⁻¹ ≈ 0.3679
//! ratio_2δ = max|E_y(70,0,100)| / max|E_y(50,0,100)|  →  target e⁻² ≈ 0.1353
//!
//! Gate A: |ratio_1δ - e⁻¹| / e⁻¹ < 10 %
//! Gate B: |ratio_2δ - e⁻²| / e⁻² < 15 %
//! ```
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test ohmic_skin_depth --release -- --nocapture
//! ```

use std::f64::consts::PI;

use yee_fdtd::{WalkingSkeletonSolver, YeeGrid};

// ---------------------------------------------------------------------------
// Grid parameters
// ---------------------------------------------------------------------------

const NX: usize = 80;
const NY: usize = 1;
const NZ: usize = 200;
const DX: f64 = 1.0e-3; // 1 mm cells

// ---------------------------------------------------------------------------
// Material parameters
// ---------------------------------------------------------------------------

/// Target frequency (Hz).
const FREQ: f64 = 1.0e9; // 1 GHz
/// μ₀ (H/m).
const MU0: f64 = 1.256_637_061_4e-6;
/// Conductivity chosen so δ = 10 mm = 10 cells at 1 GHz.
///
/// σ = 2 / (ω μ₀ δ²) = 2 / (2π·10⁹ · 4π·10⁻⁷ · (10⁻²)²) ≈ 2.533 S/m.
const SIGMA: f64 = 2.5331;

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

/// First conductor cell (x index). Cells 0..49 are vacuum; 50..79 are conductor.
const X_SURFACE: usize = 50;
/// Source cell — E_y at this (x, y, z) position in the vacuum region.
const SRC_X: usize = 25;
const SRC_Y: usize = 0;
/// Place the source at z = NZ/2 so both PEC walls are equidistant; this
/// maximises the amplitude of the dominant kz = π/(NZ·dx) standing-wave mode.
const SRC_Z: usize = NZ / 2; // 100

// ---------------------------------------------------------------------------
// Time-stepping parameters
// ---------------------------------------------------------------------------

/// Steps to let the field couple into the conductor and reach sinusoidal
/// quasi-steady state.
const N_TRANSIENT: usize = 6_000;
/// Steps over which the peak amplitude is recorded (≈ 3.9 periods at 1 GHz).
const N_MEASURE: usize = 2_000;

// ---------------------------------------------------------------------------
// Physics helpers
// ---------------------------------------------------------------------------

/// Analytic skin depth: `δ = √(2 / (ω μ₀ σ))` (Griffiths §9.4.1).
fn analytic_skin_depth(sigma: f64, freq: f64) -> f64 {
    let omega = 2.0 * PI * freq;
    (2.0 / (omega * MU0 * sigma)).sqrt()
}

// ---------------------------------------------------------------------------
// Simulation runner
// ---------------------------------------------------------------------------

/// Run the skin-depth simulation and return peak E_y amplitudes at:
/// - `(X_SURFACE, 0, SRC_Z)` — conductor surface
/// - `(X_SURFACE + 10, 0, SRC_Z)` — 1δ depth
/// - `(X_SURFACE + 20, 0, SRC_Z)` — 2δ depth
///
/// E_y propagates in x via the H_z ↔ E_y leapfrog stencil. The grid is
/// 80×1×200: NY=1 makes the problem effectively 2-D in x-z, and NZ=200 mm
/// is wide enough that the dominant standing-wave mode (kz=π/200 mm = 15.7/m)
/// satisfies kz << 1/δ = 100/m, so the skin-depth decay dominates in x.
fn run_skin_depth_sim() -> (f64, f64, f64) {
    let mut grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    // set_sigma_box uses exclusive upper bounds; cover x = X_SURFACE..NX.
    grid.set_sigma_box(X_SURFACE, NX + 1, 0, NY + 1, 0, NZ + 1, SIGMA);

    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    let mut amp_surface: f64 = 0.0;
    let mut amp_1delta: f64 = 0.0;
    let mut amp_2delta: f64 = 0.0;

    for n in 0..N_TRANSIENT + N_MEASURE {
        let t = n as f64 * dt;

        solver.update_h_only();

        // Soft sinusoidal E_y source in the vacuum region.
        // E_y propagates in x via the H_z ↔ E_y leapfrog, driving the
        // conductor's CA/CB update and establishing a skin-depth profile.
        solver.grid_mut().ey[(SRC_X, SRC_Y, SRC_Z)] += (2.0 * PI * FREQ * t).sin();

        solver.update_e_only();
        // apply_cpml_e applies PEC boundary when no CPML is configured.
        solver.apply_cpml_e();
        solver.advance_clock();

        // Record peak during the measurement window.
        if n >= N_TRANSIENT {
            amp_surface = amp_surface.max(solver.grid().ey[(X_SURFACE, SRC_Y, SRC_Z)].abs());
            amp_1delta = amp_1delta.max(solver.grid().ey[(X_SURFACE + 10, SRC_Y, SRC_Z)].abs());
            amp_2delta = amp_2delta.max(solver.grid().ey[(X_SURFACE + 20, SRC_Y, SRC_Z)].abs());
        }
    }

    (amp_surface, amp_1delta, amp_2delta)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Sanity-check for the analytic skin-depth formula.
#[test]
fn analytic_skin_depth_is_ten_cells() {
    let delta = analytic_skin_depth(SIGMA, FREQ);
    let cells = delta / DX;
    assert!(
        (cells - 10.0).abs() < 0.01,
        "expected δ = 10.00 cells, got {cells:.4}"
    );
}

/// **fdtd-205** gate: exponential skin-depth penetration profile.
///
/// Verifies that the CA/CB Ohmic-loss E-update (Taflove §3.7) reproduces
/// `|E(x)| = |E₀| exp(-x/δ)` inside a conducting half-space.
///
/// Gate A (1δ): `|ratio_1δ - e⁻¹| / e⁻¹ < 10 %`
/// Gate B (2δ): `|ratio_2δ - e⁻²| / e⁻² < 15 %`
#[test]
fn skin_depth_ratios_match_analytic() {
    let delta = analytic_skin_depth(SIGMA, FREQ);
    // Confirm the geometry before running.
    assert!(
        (delta / DX - 10.0).abs() < 0.01,
        "δ should be 10 cells, got {:.4} cells",
        delta / DX
    );

    let (amp_s, amp_1, amp_2) = run_skin_depth_sim();

    assert!(
        amp_s > 1e-30,
        "amp_surface is negligible ({amp_s:.2e}); E_y did not couple to conductor surface"
    );

    let ratio_1 = amp_1 / amp_s;
    let ratio_2 = amp_2 / amp_s;
    let target_1 = (-1.0_f64).exp(); // e⁻¹ ≈ 0.3679
    let target_2 = (-2.0_f64).exp(); // e⁻² ≈ 0.1353

    let err_1 = (ratio_1 - target_1).abs() / target_1;
    let err_2 = (ratio_2 - target_2).abs() / target_2;

    println!(
        "\nfdtd-205 skin-depth penetration gate
  δ_analytic   = {:.2} mm = {:.1} cells
  amp_surface  = {amp_s:.4e}
  ratio_1δ     = {ratio_1:.4}  (target e⁻¹ = {target_1:.4},  rel_err = {:.2} %)
  ratio_2δ     = {ratio_2:.4}  (target e⁻² = {target_2:.4},  rel_err = {:.2} %)
  Gate A threshold: 10 %
  Gate B threshold: 15 %",
        delta * 1e3,
        delta / DX,
        err_1 * 100.0,
        err_2 * 100.0,
    );

    assert!(
        err_1 < 0.10,
        "Gate A FAILED: ratio_1δ = {ratio_1:.4}, target = {target_1:.4}, \
         rel_err = {:.2} % (threshold 10 %)",
        err_1 * 100.0,
    );
    assert!(
        err_2 < 0.15,
        "Gate B FAILED: ratio_2δ = {ratio_2:.4}, target = {target_2:.4}, \
         rel_err = {:.2} % (threshold 15 %)",
        err_2 * 100.0,
    );
}
