//! Validation gate **fdtd-201** — rectangular PEC-cavity TE₁₀₁ resonant
//! frequency extracted from a time-domain FDTD run.
//!
//! # Physics
//!
//! An air-filled rectangular cavity with PEC walls at dimensions
//! `a × b × d` (x × y × z) supports resonant modes at frequencies
//! (Pozar §6.3):
//!
//! ```text
//! f_{mnp} = (c / 2) * sqrt((m/a)² + (n/b)² + (p/d)²)
//! ```
//!
//! The dominant TE₁₀₁ mode (`m=1, n=0, p=1`) has:
//!
//! ```text
//! f₁₀₁ = (c / 2) * sqrt((1/a)² + (1/d)²)
//! ```
//!
//! With `a = d > b` this is the lowest-frequency mode and is well-separated
//! from the next (TE₂₀₁ / TE₁₀₃ at ~2× f₁₀₁, and TE₁₁₁ which requires
//! a non-zero n/b term).
//!
//! # Mode structure and source alignment
//!
//! For TE₁₀₁ the dominant E-field component is **E_y**
//! (tangential to the broad faces at y = 0 and y = b, which are the driven
//! wall pair, matching the TE₁₀ waveguide derivation):
//!
//! ```text
//! E_y  ∝  sin(π·x/a) · sin(π·z/d)
//! H_x  ∝  cos(π·x/a) · sin(π·z/d)  (roughly)
//! H_z  ∝  sin(π·x/a) · cos(π·z/d)  (roughly)
//! ```
//!
//! The maximum of E_y sits at x = a/2, z = d/2 — the geometric centre of
//! the broad face.  The source and probe are placed there (or nearby) to
//! maximise coupling to this mode.
//!
//! We inject a soft Gaussian pulse directly into `grid.ey` (via the public
//! [`yee_fdtd::YeeGrid::ey`] field accessed through `solver.grid_mut()`),
//! using the custom-step-body pattern demonstrated in
//! `tests/lumped_resistor.rs`.  This stays inside the lane — no `src/`
//! changes, pure consumer of the public grid API.
//!
//! # Method
//!
//! 1. Build the cavity via `YeeGrid::vacuum(nx, ny, nz, dx)` with hard PEC
//!    outer walls (the `WalkingSkeletonSolver::new` default).
//! 2. Inject an off-centre Gaussian pulse into E_y via the public grid field.
//! 3. Step N times using the custom-body pattern, recording an E_y probe
//!    time series at an interior point.
//! 4. Extract the resonant frequency by scanning a dense candidate grid of
//!    single-bin DFTs over `[f₁₀₁·0.65, f₁₀₁·1.50]`, then peak-finding.
//!    No FFT library is needed; this mirrors the `ntff.rs:253` single-
//!    frequency DFT accumulator idiom.
//!
//! # Tolerance
//!
//! The gate asserts the extracted resonance matches the analytic TE₁₀₁
//! frequency within **±2.5 %**.  Grid dispersion on a ≈28-cells-per-
//! wavelength Yee mesh in vacuum typically contributes ~0.5–1 % numerical
//! phase velocity error; the ±2.5 % band gives comfortable margin while
//! being clearly non-trivial (the next mode is >40 % above f₁₀₁).
//!
//! ## Strict ±0.5 % refinement path
//!
//! To tighten to ±0.5 % (the `validation/README.md` target), halve `dx`
//! (doubling `nx/ny/nz`) and double `N_STEPS`.  With `dx = 5 mm` /
//! ≈56 cells·λ⁻¹ / `N = 60 000` the grid-dispersion error falls below 0.3 %
//! and the DFT frequency resolution is ~0.08 % of f₁₀₁.  The Q-factor
//! extraction (the other column in the README row) is deferred to a
//! follow-on slice; it requires fitting the damped-exponential decay of the
//! free-running cavity resonance and is out of scope for fdtd-201.
//!
//! # Wall-time budget
//!
//! Grid: `20 × 10 × 20 = 4 000` cells; `30 000` steps; 400-candidate DFT
//! scan.  Observed wall-time: ~5–15 s release; gated with `#[ignore]` like
//! the sibling slow integration tests.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test cavity_resonance --release -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_fdtd::boundary;
use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver, YeeGrid};

// --------------------------------------------------------------------------
// Cavity geometry
// --------------------------------------------------------------------------
//
// Choose a = d > b so TE₁₀₁ is cleanly dominant and well-separated.
// Physical: a = d = 0.20 m, b = 0.10 m.
//
// Grid: dx = 0.010 m → nx = 20, ny = 10, nz = 20.
// Physical cavity dimensions = N_cells · dx for each axis.
//
// Note: YeeGrid::vacuum(nx, ny, nz, dx) creates `nx×ny×nz` primary cells.
// The PEC wall clamps tangential E on faces at [0] and [nx], so the
// physical interior dimension is a = nx·dx.
//
const NX: usize = 20;
const NY: usize = 10;
const NZ: usize = 20;
const DX: f64 = 0.010; // metres

// Physical interior dimensions (= N_cells · dx):
const PHYS_A: f64 = NX as f64 * DX; // 0.20 m  (x)
const PHYS_B: f64 = NY as f64 * DX; // 0.10 m  (y)
const PHYS_D: f64 = NZ as f64 * DX; // 0.20 m  (z)
// (PHYS_B kept for documentation clarity; not used in the frequency formula)
const _PHYS_B: f64 = PHYS_B;

// --------------------------------------------------------------------------
// Run parameters
// --------------------------------------------------------------------------

/// Total number of FDTD time steps.
///
/// At dt ≈ 1.732 × 10⁻¹¹ s (0.9 × CFL on 10 mm cubic cells), 30 000 steps
/// gives a total simulation time of ~0.520 µs.
/// Frequency resolution: Δf = 1 / (N · dt) ≈ 1.92 MHz, which is ~0.18 %
/// of f₁₀₁ ≈ 1.061 GHz — well below the ±2.5 % assertion window.
const N_STEPS: usize = 30_000;

/// Source cell index into `ey` (shape `[nx+1, ny, nz+1]`).
///
/// TE₁₀₁ has E_y ∝ sin(π·x/a)·sin(π·z/d).  We place the source at roughly
/// (nx/4, ny/2, nz/4) for strong coupling (sin²(π/4) ≈ 0.5 of peak) while
/// being off-centre so higher modes also see the source.
///
/// E_y staggering: ey[i, j, k] lives at (i+0.5)·dx, j·dy, (k+0.5)·dz.
/// Index (NX/4, NY/2, NZ/4) = (5, 5, 5) maps to x = 5.5·dx = 0.055 m,
/// z = 5.5·dx = 0.055 m.  sin(π·0.055/0.2) ≈ sin(0.275π) ≈ 0.757.
const SRC_I: usize = NX / 4; // 5
const SRC_J: usize = NY / 2; // 5
const SRC_K: usize = NZ / 4; // 5

/// Probe cell for E_y, away from the source and at a strong antinode.
///
/// Index (NX*3/4, NY/2, NZ*3/4) = (15, 5, 15).
/// x = 15.5·dx = 0.155 m, z = 15.5·dx = 0.155 m.
/// sin(π·0.155/0.2) ≈ sin(0.775π) ≈ 0.757.
/// Source and probe are symmetric about the centre; TE₁₀₁ drives both in phase.
const PRB_I: usize = NX * 3 / 4; // 15
const PRB_J: usize = NY / 2; // 5
const PRB_K: usize = NZ * 3 / 4; // 15

/// Number of candidate frequencies in the DFT scan.
const N_FREQ_BINS: usize = 400;

// --------------------------------------------------------------------------
// Speed of light (matches yee-core::units::C0)
// --------------------------------------------------------------------------
const C0: f64 = 299_792_458.0; // m/s

/// Analytic TE₁₀₁ resonant frequency:
///
/// ```text
/// f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)
/// ```
fn analytic_f101() -> f64 {
    0.5 * C0 * ((1.0 / (PHYS_A * PHYS_A)) + (1.0 / (PHYS_D * PHYS_D))).sqrt()
}

/// Gaussian pulse amplitude at time `t` centred at `t0` with width `sigma`.
#[inline]
fn gaussian(t: f64, t0: f64, sigma: f64) -> f64 {
    let arg = (t - t0) / sigma;
    (-arg * arg).exp()
}

/// `fdtd-201` gate: extract TE₁₀₁ resonant frequency via single-bin DFT scan.
///
/// Gate tolerance: ±2.5 %.
/// See module docstring for the ±0.5 % refinement path.
///
/// Run with:
/// ```bash
/// cargo test -p yee-fdtd --test cavity_resonance --release -- --ignored --nocapture
/// ```
#[test]
#[ignore = "slow: ~5-15 s release; fdtd-201 TE101 cavity resonance gate (Phase 2.fdtd)"]
fn te101_resonance_matches_analytic_within_two_point_five_percent() {
    // ----------------------------------------------------------------
    // Build cavity: vacuum grid + hard PEC outer walls.
    // WalkingSkeletonSolver::new uses the deprecated apply_pec boundary
    // (reflecting), which is exactly right for a closed cavity.
    // ----------------------------------------------------------------
    let grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    // ----------------------------------------------------------------
    // Source parameters: broadband Gaussian centred at t0 with σ small
    // enough to cover f₁₀₁ and several cavity modes.
    //
    // σ = 4·dt → bandwidth ~1/(2πσ) ≈ 2.3 GHz, easily covering f₁₀₁.
    // t0 = 12·dt → Gaussian tail at t=0 is e^{-(12/4)²} ≈ 10⁻³·e⁻⁹ ≈ 10⁻⁷.
    // ----------------------------------------------------------------
    let t0 = 12.0 * dt;
    let sigma = 4.0 * dt;

    // ----------------------------------------------------------------
    // Run the cavity using a custom step body that injects E_y.
    //
    // This mirrors the pattern in tests/lumped_resistor.rs: we call
    // the individual sub-step helpers (update_h, apply_cpml_h,
    // update_e, apply_cpml_e, advance_clock) and insert the source
    // injection between the H and E updates.
    //
    // The E_y field at ey[(i, j, k)] is public on YeeGrid, so writing
    // to it is a pure consumer call — no src/ change required.
    // ----------------------------------------------------------------
    let mut probe_series: Vec<f64> = Vec::with_capacity(N_STEPS);

    for _n in 0..N_STEPS {
        let t = solver.current_time();

        // H update.
        solver.update_h_only();
        // PEC outer-face clamp on H (no CPML in the vanilla solver).
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());

        // Inject broadband Gaussian into E_y at the source cell.
        // ey shape: [nx+1, ny, nz+1]; SRC indices are in-bounds.
        {
            let amp = gaussian(t, t0, sigma);
            solver.grid_mut().ey[(SRC_I, SRC_J, SRC_K)] += amp;
        }

        // E update.
        solver.update_e_only();
        // PEC outer-face clamp on E (also applies any interior PEC masks).
        solver.apply_cpml_e();

        // Advance the step counter.
        solver.advance_clock();

        // Record probe.
        probe_series.push(solver.grid().ey[(PRB_I, PRB_J, PRB_K)]);
    }

    // ----------------------------------------------------------------
    // Frequency scan via single-bin DFT (Goertzel accumulation).
    //
    // For each candidate frequency f we compute:
    //
    //   |DFT(f)|² = (Σ_n x[n]·cos(ω·n·dt))² + (Σ_n x[n]·sin(ω·n·dt))²
    //
    // This is the same pattern as ntff.rs:253 (phase = exp(-jωt)·dt).
    // Scan over N_FREQ_BINS linearly-spaced candidates in the band
    // [0.65·f₁₀₁, 1.50·f₁₀₁] — wide enough to capture f₁₀₁ and the
    // first few modes above it.
    // ----------------------------------------------------------------
    let f_ref = analytic_f101();
    let f_lo = 0.65 * f_ref;
    let f_hi = 1.50 * f_ref;
    let df_scan = (f_hi - f_lo) / (N_FREQ_BINS - 1) as f64;

    let mut peak_power = 0.0_f64;
    let mut peak_freq = f_lo;

    for bin in 0..N_FREQ_BINS {
        let f_candidate = f_lo + bin as f64 * df_scan;
        let omega = 2.0 * PI * f_candidate;

        let mut re_acc = 0.0_f64;
        let mut im_acc = 0.0_f64;
        for (n, &x) in probe_series.iter().enumerate() {
            let phase = omega * n as f64 * dt;
            re_acc += x * phase.cos();
            im_acc -= x * phase.sin();
        }

        let power = re_acc * re_acc + im_acc * im_acc;
        if power > peak_power {
            peak_power = power;
            peak_freq = f_candidate;
        }
    }

    // ----------------------------------------------------------------
    // Diagnostics — printed when run with --nocapture.
    // ----------------------------------------------------------------
    let rel_error = (peak_freq - f_ref) / f_ref;
    eprintln!(
        "\nfdtd-201 TE₁₀₁ rectangular-cavity resonance gate
  cavity:         a = {:.4} m, b = {:.4} m, d = {:.4} m
  grid:           {}×{}×{}, dx = {:.1} mm, dt = {:.4e} s
  steps:          {} (T_total = {:.4e} s)
  DFT scan:       {N_FREQ_BINS} bins in [{:.4} GHz, {:.4} GHz]
  analytic f₁₀₁:  {:.6} GHz
  extracted f:    {:.6} GHz  (|DFT|² = {:.3e})
  relative error: {:.4} %
",
        PHYS_A,
        PHYS_B,
        PHYS_D,
        NX,
        NY,
        NZ,
        DX * 1e3,
        dt,
        N_STEPS,
        N_STEPS as f64 * dt,
        f_lo * 1e-9,
        f_hi * 1e-9,
        f_ref * 1e-9,
        peak_freq * 1e-9,
        peak_power,
        rel_error * 100.0,
    );

    // ----------------------------------------------------------------
    // Gate: |rel_error| ≤ 2.5 %.
    //
    // Grid dispersion on a ~28-cells/λ Yee mesh in vacuum shifts the
    // numerical phase velocity by ~0.5–1 %, so the true tight ±0.5 %
    // bound is reached only on a refined grid (see module docstring).
    // ----------------------------------------------------------------
    assert!(
        rel_error.abs() < 0.025,
        "fdtd-201 FAILED: extracted f₁₀₁ = {:.6} GHz, analytic = {:.6} GHz, \
         rel_error = {:.4} % (threshold ±2.5 %)",
        peak_freq * 1e-9,
        f_ref * 1e-9,
        rel_error * 100.0,
    );
}
