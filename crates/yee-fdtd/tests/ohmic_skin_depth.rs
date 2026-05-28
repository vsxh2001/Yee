//! fdtd-205: Ohmic skin-depth penetration gate.
//!
//! # Physics
//!
//! Validates that the FDTD CA/CB Ohmic-loss E-update reproduces the
//! exponential spatial decay of a plane wave inside a conducting half-space:
//!
//! ```text
//! |E(z_surface + n·δ)| / |E(z_surface)| = e^{-n}
//! ```
//!
//! where `δ = √(2 / (ω μ₀ σ))` is the skin depth (Griffiths §9.4.1,
//! Taflove §3.7).
//!
//! # Geometry
//!
//! 1D-like slab: `NX × NY × NZ = 5 × 5 × 130`, `DX = 1 mm`.
//! Vacuum region: `z ∈ [0, 50)` cells.
//! Conductor: `z ∈ [50, 130)` cells, `σ = 2.5331 S/m → δ = 10 mm = 10 cells`.
//!
//! Source: soft sinusoidal `E_x` injection spanning the full transverse
//! cross-section at k = `SRC_Z` = 25. Injecting across all interior `(i, j)`
//! nodes (`j ∈ [1, NY)`) ensures a uniform transverse profile — this is
//! required to launch a 1D TEM-like
//! plane wave in the +z direction.
//!
//! # Boundary conditions
//!
//! The transverse (x,y) boundaries use PMC (Perfect Magnetic Conductor)
//! at the y-faces: `H_z` is zeroed at `j = 0` and `j = NY − 1` after each
//! H-field update. This prevents the PEC-wall-induced `H_z` cascade that
//! would otherwise distort the (E_x, H_y) TEM wave: in a PEC box the
//! sinusoidal transverse mode has cutoff ≈ 42 GHz >> 1 GHz and is
//! evanescent, not propagating. The PMC-y condition removes this evanescent
//! transverse mode contamination and recovers the true 1D skin-depth physics.
//!
//! The (E_x, H_y) pair propagates in +z as a quasi-TEM wave: E_x drives
//! ∂H_y/∂z via Faraday's law and H_y drives ∂E_x/∂z via Ampère's law,
//! identical to the textbook 1D FDTD TEM polarisation. E_z, H_x and E_y
//! remain zero by symmetry; H_z is enforced zero at y-boundaries by PMC.
//!
//! Probes: average over `i ∈ [0, NX)` at `j = NY/2`: k=50 surface,
//!         k=60 (1δ), k=70 (2δ).
//!
//! # Gate tolerances
//!
//! * Gate A: `|ratio_1δ − e⁻¹| / e⁻¹ < 10 %`
//! * Gate B: `|ratio_2δ − e⁻²| / e⁻² < 15 %`
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test ohmic_skin_depth -- --nocapture
//! ```

use std::f64::consts::{E, PI};

use yee_fdtd::{WalkingSkeletonSolver, YeeGrid};

// ---------------------------------------------------------------------------
// Grid parameters
// ---------------------------------------------------------------------------

const NX: usize = 5;
const NY: usize = 5;
const NZ: usize = 130;
const DX: f64 = 1.0e-3; // 1 mm cells

// ---------------------------------------------------------------------------
// Material parameters
// ---------------------------------------------------------------------------

/// Electric conductivity. Chosen so that δ = 10 mm = 10 cells exactly:
///   δ = √(2 / (ω μ₀ σ)), so σ = 2 / (ω μ₀ δ²)
///   with ω = 2π·1e9, μ₀ = 4π×10⁻⁷, δ = 10e-3.
const SIGMA: f64 = 2.5331; // S/m

// ---------------------------------------------------------------------------
// Source parameters
// ---------------------------------------------------------------------------

/// Excitation frequency.
const FREQ: f64 = 1.0e9; // 1 GHz
/// Source z-index (k) for the Ex plane-wave injection.
const SRC_Z: usize = 25;

// ---------------------------------------------------------------------------
// Conductor interface and measurement points
// ---------------------------------------------------------------------------

/// First conductor cell index (inclusive, z direction).
const Z_SURFACE: usize = 50;

// ---------------------------------------------------------------------------
// Run parameters
// ---------------------------------------------------------------------------

/// Transient burn-in: allow the CW to build up and reach steady state at the
/// surface before recording amplitudes.
const N_TRANSIENT: usize = 6_000;
/// Measurement window: track peak |E_x| over this many steps after transient.
const N_MEASURE: usize = 2_000;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Free-space permeability (H/m).
const MU0: f64 = 1.256_637_061_4e-6;

// ---------------------------------------------------------------------------
// Physics helpers
// ---------------------------------------------------------------------------

/// Analytic skin depth (m) for a uniform conductor at frequency `freq`.
///
/// ```text
/// δ = √(2 / (ω μ₀ σ))
/// ```
fn analytic_skin_depth(sigma: f64, freq: f64) -> f64 {
    let omega = 2.0 * PI * freq;
    (2.0 / (omega * MU0 * sigma)).sqrt()
}

// ---------------------------------------------------------------------------
// Core simulation runner
// ---------------------------------------------------------------------------

/// Run the Ohmic-skin-depth simulation and return
/// `(amp_surface, amp_1delta, amp_2delta)`.
///
/// Injects a CW sinusoidal `E_x` source spanning the full transverse
/// cross-section (all `i ∈ [0, NX)`, `j ∈ [1, NY)`) at `k = SRC_Z` to
/// create a quasi-1D TEM plane wave propagating in +z. Probes the `E_x`
/// amplitude at three z-depths inside the conductor (surface, 1δ, 2δ).
///
/// # Boundary conditions
///
/// PMC (Perfect Magnetic Conductor) is enforced at the y-faces by zeroing
/// `H_z` at `j = 0` and `j = NY − 1` after each H-update. This prevents
/// the PEC-box evanescent-mode cascade and recovers the 1D TEM skin-depth
/// physics. PEC at z-faces is enforced by zeroing `E_x` at `k = 0` and
/// `k = NZ`. The x-faces require no special treatment because `E_y = 0`
/// and `E_z = 0` throughout (1D TEM symmetry).
///
/// * `amp_surface` — peak |E_x| at `(i, NY/2, Z_SURFACE)` (conductor face)
/// * `amp_1delta`  — peak |E_x| at `(i, NY/2, Z_SURFACE + 10)` (1δ inside)
/// * `amp_2delta`  — peak |E_x| at `(i, NY/2, Z_SURFACE + 20)` (2δ inside)
fn run_skin_depth_sim() -> (f64, f64, f64) {
    // Build vacuum grid.
    let mut grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    // Fill conductor region: z ∈ [Z_SURFACE, NZ).
    grid.set_sigma_box(0, NX + 1, 0, NY + 1, Z_SURFACE, NZ + 1, SIGMA);

    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    let mut amp_surface = 0.0_f64;
    let mut amp_1delta = 0.0_f64;
    let mut amp_2delta = 0.0_f64;

    // Probe y-index: centre of the transverse cross-section.
    let j_probe = NY / 2; // = 2

    for n in 0..(N_TRANSIENT + N_MEASURE) {
        let t = n as f64 * dt;

        // 1. H update (curl of E).
        solver.update_h_only();

        // 2. Enforce PMC at y-faces: zero H_z at j=0 and j=NY-1.
        //    H_z shape is (NX, NY, NZ+1), valid j ∈ [0, NY).
        //    This removes the evanescent transverse-mode H_z that would
        //    otherwise contaminate the 1D TEM skin-depth physics.
        for i in 0..NX {
            for k in 0..=NZ {
                solver.grid_mut().hz[(i, 0, k)] = 0.0;
                solver.grid_mut().hz[(i, NY - 1, k)] = 0.0;
            }
        }

        // 3. Soft sinusoidal source: inject E_x over the interior transverse
        //    slice at k = SRC_Z. j runs 1..NY (interior y).
        let src_amp = (2.0 * PI * FREQ * t).sin();
        for i in 0..NX {
            for j in 1..NY {
                solver.grid_mut().ex[(i, j, SRC_Z)] += src_amp;
            }
        }

        // 4. E update (curl of H, with CA/CB in conductor).
        solver.update_e_only();

        // 5. PEC at z-faces: zero E_x at k=0 and k=NZ.
        //    (update_e_only skips j=0,NY and k=0,NZ so E_x at y-faces
        //    stays zero from initialization — no explicit zeroing needed there.)
        for i in 0..NX {
            for j in 0..=NY {
                solver.grid_mut().ex[(i, j, 0)] = 0.0;
                solver.grid_mut().ex[(i, j, NZ)] = 0.0;
            }
        }

        // 6. Advance clock.
        solver.advance_clock();

        if n >= N_TRANSIENT {
            // Average over i to cancel any residual transverse asymmetry.
            let ex_surf: f64 = (0..NX)
                .map(|i| solver.grid().ex[(i, j_probe, Z_SURFACE)].abs())
                .sum::<f64>()
                / NX as f64;
            let ex_1d: f64 = (0..NX)
                .map(|i| solver.grid().ex[(i, j_probe, Z_SURFACE + 10)].abs())
                .sum::<f64>()
                / NX as f64;
            let ex_2d: f64 = (0..NX)
                .map(|i| solver.grid().ex[(i, j_probe, Z_SURFACE + 20)].abs())
                .sum::<f64>()
                / NX as f64;
            amp_surface = amp_surface.max(ex_surf);
            amp_1delta = amp_1delta.max(ex_1d);
            amp_2delta = amp_2delta.max(ex_2d);
        }
    }

    (amp_surface, amp_1delta, amp_2delta)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// **fdtd-205** gate: exponential skin-depth penetration in a conducting
/// half-space.
///
/// σ = 2.5331 S/m at f = 1 GHz → δ_analytic = 10 mm = 10 cells.
/// Source: E_x plane-wave injection spanning the full transverse cross-section.
/// Boundary: PMC at y-faces (H_z = 0) to suppress evanescent transverse modes.
///
/// Gate A: `|ratio_1δ − e⁻¹| / e⁻¹ < 10 %`
/// Gate B: `|ratio_2δ − e⁻²| / e⁻² < 15 %`
#[test]
fn skin_depth_ratios_match_analytic() {
    let delta = analytic_skin_depth(SIGMA, FREQ);

    // Guard: δ ≈ 10 cells (sanity check on constants).
    assert!(
        (delta / DX - 10.0).abs() < 0.01,
        "δ should be 10 cells, got {:.4}",
        delta / DX
    );

    let (amp_s, amp_1, amp_2) = run_skin_depth_sim();

    assert!(
        amp_s > 1e-10,
        "amp_surface = {amp_s:.4e} — field is not reaching the conductor \
         surface; check sigma region and N_TRANSIENT"
    );

    let ratio_1 = amp_1 / amp_s;
    let ratio_2 = amp_2 / amp_s;
    let target_1 = 1.0_f64 / E; // e^{-1}
    let target_2 = (-2.0_f64).exp(); // e^{-2}
    let err_1 = (ratio_1 - target_1).abs() / target_1;
    let err_2 = (ratio_2 - target_2).abs() / target_2;

    eprintln!(
        "\nfdtd-205: Ohmic skin-depth penetration gate
  δ_analytic  = {:.2} mm  ({:.4} cells)
  amp_surface = {amp_s:.4e}
  amp_1δ      = {amp_1:.4e}
  amp_2δ      = {amp_2:.4e}
  ratio_1δ    = {ratio_1:.6}  (target e⁻¹ = {target_1:.6})
  ratio_2δ    = {ratio_2:.6}  (target e⁻² = {target_2:.6})
  err_1δ      = {:.2} %  (gate < 10 %)
  err_2δ      = {:.2} %  (gate < 15 %)
",
        delta * 1e3,
        delta / DX,
        err_1 * 100.0,
        err_2 * 100.0,
    );

    // Gate A: 10% tolerance at 1δ.
    assert!(
        err_1 < 0.10,
        "Gate A FAILED: ratio_1δ = {ratio_1:.6} (target {target_1:.6}), \
         rel_err = {:.2} % (threshold 10 %)",
        err_1 * 100.0
    );
    // Gate B: 15% tolerance at 2δ.
    assert!(
        err_2 < 0.15,
        "Gate B FAILED: ratio_2δ = {ratio_2:.6} (target {target_2:.6}), \
         rel_err = {:.2} % (threshold 15 %)",
        err_2 * 100.0
    );
}
