//! Validation gate **fdtd-202** — Q-factor extraction from a lossy
//! rectangular PEC cavity driven by the CA/CB E-update.
//!
//! # Physics
//!
//! An air-filled rectangular cavity with PEC walls and a uniform electric
//! conductivity σ throughout the volume. The TE₁₀₁ mode excited at frequency
//! f₁₀₁ decays exponentially with time constant τ = 2ε₀/σ, giving a Q-factor
//! (Pozar §6.7, Taflove §3.7):
//!
//! ```text
//! Q = ω₁₀₁ · ε₀ / σ  =  2π · f₁₀₁ · ε₀ / σ
//! ```
//!
//! This test validates the lossy Yee E-update (Taflove §3.7 CA/CB form) by:
//!
//! 1. Exciting the cavity with a broadband Gaussian pulse.
//! 2. Letting the TE₁₀₁ mode ring down freely after the source is off.
//! 3. Fitting an exponential decay to the ring-down and extracting τ.
//! 4. Computing Q = π · f₁₀₁ · τ and comparing to the analytic value.
//!
//! # Cavity geometry
//!
//! `a × b × d` (x × y × z) with `a = d = 0.20 m`, `b = 0.10 m`.
//! Grid: `nx=20, ny=10, nz=20`, `dx=0.01 m`.
//!
//! Analytic TE₁₀₁ frequency:
//!
//! ```text
//! f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)  ≈  1.0607 GHz
//! ```
//!
//! Analytic Q (σ₀ = 2.96e-3 S/m):
//!
//! ```text
//! Q = 2π · f₁₀₁ · ε₀ / σ₀  ≈  20
//! ```
//!
//! # Stability check
//!
//! For σ₀ = 2.96e-3 S/m at dt ≈ 17.3 ps:
//!   CA = (2ε₀ − σΔt) / (2ε₀ + σΔt) ≈ 0.9971  (stable, |CA| < 1)
//!
//! # Gate tolerance
//!
//! `|Q_measured − Q_analytic| / Q_analytic < 5 %`
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test cavity_q --release -- --nocapture
//! ```

use std::f64::consts::PI;

use yee_fdtd::boundary;
use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver, YeeGrid};

// ---------------------------------------------------------------------------
// Cavity geometry — same as fdtd-201 (cavity_resonance.rs)
// ---------------------------------------------------------------------------

const NX: usize = 20;
const NY: usize = 10;
const NZ: usize = 20;
const DX: f64 = 0.010; // metres

const PHYS_A: f64 = NX as f64 * DX; // 0.20 m (x)
const PHYS_D: f64 = NZ as f64 * DX; // 0.20 m (z)

// Speed of light (matches yee-core::units::C0)
const C0: f64 = 299_792_458.0; // m/s

// ε₀ in SI units
const EPS0: f64 = 8.854_187_817e-12; // F/m

// ---------------------------------------------------------------------------
// Material parameters
// ---------------------------------------------------------------------------

/// Electric conductivity used in fdtd-202 main gate.
///
/// Analytic Q = 2π · f₁₀₁ · ε₀ / σ₀ ≈ 20.
const SIGMA0: f64 = 2.96e-3; // S/m

// ---------------------------------------------------------------------------
// Run parameters
// ---------------------------------------------------------------------------

/// Number of source-injection steps (broadband pulse).
const N_SRC: usize = 200;

/// Number of ring-down steps recorded for Q extraction.
///
/// At dt ≈ 17.3 ps, 6000 steps ≈ 104 ns ≈ 17.4 τ (where τ = 2ε₀/σ ≈ 6 ns).
/// Using the last 2/3 of the window places the fit at 5.8–17.4 τ, well into
/// the single-mode ring-down where higher modes (TE₁₀₃ etc.) are negligible.
const N_RING: usize = 6_000;

// ---------------------------------------------------------------------------
// Physics helpers
// ---------------------------------------------------------------------------

/// Analytic TE₁₀₁ resonant frequency.
///
/// ```text
/// f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)
/// ```
fn analytic_f101() -> f64 {
    0.5 * C0 * ((1.0 / (PHYS_A * PHYS_A)) + (1.0 / (PHYS_D * PHYS_D))).sqrt()
}

/// Analytic Q-factor for a uniformly lossy cavity (Taflove §3.7):
///
/// ```text
/// Q = 2π · f₁₀₁ · ε₀ / σ
/// ```
fn analytic_q(sigma: f64) -> f64 {
    2.0 * PI * analytic_f101() * EPS0 / sigma
}

/// Gaussian pulse amplitude at time `t` centred at `t0` with width `sigma_t`.
#[inline]
fn gaussian(t: f64, t0: f64, sigma_t: f64) -> f64 {
    let arg = (t - t0) / sigma_t;
    (-arg * arg).exp()
}

/// Fit an exponential decay to `time_series` and return the decay time
/// constant τ (seconds).
///
/// Uses log-linear (least-squares) regression:
///
/// ```text
/// log|y[n]| = A + slope · t[n],   slope = −1/τ
/// ```
///
/// `t_start` is the absolute simulation time of `time_series[0]`.
/// Samples with `|y| < 1e-30` are excluded to avoid `ln(0)`.
fn fit_log_decay(time_series: &[f64], dt: f64, t_start: f64) -> f64 {
    let n = time_series.len() as f64;
    let ts: Vec<f64> = (0..time_series.len())
        .map(|i| t_start + i as f64 * dt)
        .collect();
    let ys: Vec<f64> = time_series
        .iter()
        .map(|&v| v.abs().max(1e-30).ln())
        .collect();

    let t_mean = ts.iter().sum::<f64>() / n;
    let y_mean = ys.iter().sum::<f64>() / n;

    let num = ts
        .iter()
        .zip(ys.iter())
        .map(|(&t, &y)| (t - t_mean) * (y - y_mean))
        .sum::<f64>();
    let den = ts.iter().map(|&t| (t - t_mean).powi(2)).sum::<f64>();

    // slope should be negative (decay); -1/slope = τ
    let slope = num / den;
    -1.0 / slope
}

// ---------------------------------------------------------------------------
// Core simulation runner
// ---------------------------------------------------------------------------

/// Run the lossy-cavity simulation and return the probe ring-down series.
///
/// 1. Builds a vacuum grid and fills the entire domain with conductivity
///    `sigma`.
/// 2. Injects a Gaussian pulse into E_y for `n_src` steps.
/// 3. Lets the field ring down for `n_ring` steps, recording E_y at the
///    probe location.
///
/// Returns `(probe, dt)` — the ring-down time series and the grid time step.
fn run_lossy_cavity(sigma: f64, n_src: usize, n_ring: usize) -> (Vec<f64>, f64) {
    // Build grid and attach conductivity over the full domain.
    let mut grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    // set_sigma_box uses inclusive-exclusive indexing; pass NX+1 etc. to
    // cover the full [nx+1, ny+1, nz+1] extent.
    grid.set_sigma_box(0, NX + 1, 0, NY + 1, 0, NZ + 1, sigma);

    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    // Source placed off-centre to couple to TE₁₀₁ (same as cavity_resonance.rs).
    // E_y shape is [NX+1, NY, NZ+1] = [21, 10, 21].
    let src_i = NX / 4; // 5  — sin(π/4) ≈ 0.707 coupling to TE₁₀₁
    let src_j = NY / 2; // 5
    let src_k = NZ / 4; // 5

    // Probe at a strong TE₁₀₁ antinode (opposite quarter to reduce source
    // near-field bias in the ring-down).
    let prb_i = NX * 3 / 4; // 15
    let prb_j = NY / 2; // 5
    let prb_k = NZ * 3 / 4; // 15

    // Source parameters: Gaussian centred at t0 = 12·dt, width 4·dt.
    let t0 = 12.0 * dt;
    let sigma_t = 4.0 * dt;
    let amp = 1.0_f64;

    // Phase 1: inject source for n_src steps.
    for _n in 0..n_src {
        let t = solver.current_time();

        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());

        // Soft Gaussian source on E_y.
        solver.grid_mut().ey[(src_i, src_j, src_k)] += amp * gaussian(t, t0, sigma_t);

        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
    }

    // Phase 2: ring-down, record probe.
    let mut probe = Vec::with_capacity(n_ring);
    for _n in 0..n_ring {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());

        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();

        probe.push(solver.grid().ey[(prb_i, prb_j, prb_k)]);
    }

    (probe, dt)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// **fdtd-202** gate: Q-factor extracted from lossy cavity ring-down.
///
/// σ₀ = 2.96e-3 S/m → Q_analytic ≈ 20 at f₁₀₁ ≈ 1.0607 GHz.
/// Gate tolerance: ±5 %.
#[test]
fn fdtd_202_q_factor_lossy_cavity() {
    let (probe, dt) = run_lossy_cavity(SIGMA0, N_SRC, N_RING);

    // Use last 2/3 of ring-down (skip early transient where higher modes
    // are still present).
    let skip = N_RING / 3;
    let window = &probe[skip..];

    // Absolute simulation time of window[0].
    // During source phase: N_SRC steps.  During ring-down skip: skip steps.
    let t_start = (N_SRC + skip) as f64 * dt;

    let tau = fit_log_decay(window, dt, t_start);
    let f101 = analytic_f101();
    let q_measured = PI * f101 * tau;
    let q_analytic = analytic_q(SIGMA0);
    let rel_err = (q_measured - q_analytic).abs() / q_analytic;

    eprintln!(
        "\nfdtd-202 Q-factor lossy-cavity gate
  cavity:       a = {:.4} m, d = {:.4} m
  grid:         {}×{}×{}, dx = {:.1} mm, dt = {:.4e} s
  sigma:        {:.4e} S/m
  f₁₀₁:         {:.6} GHz
  Q_analytic:   {:.4}
  tau_fit:      {:.6e} s
  Q_measured:   {:.4}
  rel_error:    {:.4} %
  gate:         ±5 %
",
        PHYS_A,
        PHYS_D,
        NX,
        NY,
        NZ,
        DX * 1e3,
        dt,
        SIGMA0,
        f101 * 1e-9,
        q_analytic,
        tau,
        q_measured,
        rel_err * 100.0,
    );

    assert!(
        rel_err < 0.05,
        "fdtd-202 FAILED: Q_measured = {:.4}, Q_analytic = {:.4}, \
         rel_error = {:.4} % (threshold ±5 %)",
        q_measured,
        q_analytic,
        rel_err * 100.0,
    );
}

/// Regression: a grid with σ = 0 via `set_sigma_box` must produce field
/// values bit-identical to a grid with no `sigma_cells` at all.
///
/// When CA = 1 and CB = dt/(ε₀ε_r), the lossy form reduces to the lossless
/// `E += coeff * curl_H` — this test enforces that identity.
#[test]
fn sigma_zero_matches_lossless_update() {
    let n_steps = 5;
    let src_i = NX / 4;
    let src_j = NY / 2;
    let src_k = NZ / 4;

    // Grid A: no sigma_cells at all (pure lossless path).
    let grid_a = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid_a.dt;
    let mut solver_a = WalkingSkeletonSolver::new(grid_a);

    // Grid B: sigma_cells set to 0.0 everywhere (lossy path, but CA=1 CB=dt/ε₀).
    let mut grid_b = YeeGrid::vacuum(NX, NY, NZ, DX);
    grid_b.set_sigma_box(0, NX + 1, 0, NY + 1, 0, NZ + 1, 0.0);
    let mut solver_b = WalkingSkeletonSolver::new(grid_b);

    // Source parameters (same in both).
    let t0 = 12.0 * dt;
    let sigma_t = 4.0 * dt;
    let amp = 1.0_f64;

    for _n in 0..n_steps {
        let t_a = solver_a.current_time();
        let t_b = solver_b.current_time();

        solver_a.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver_a.grid_mut());
        solver_a.grid_mut().ey[(src_i, src_j, src_k)] += amp * gaussian(t_a, t0, sigma_t);
        solver_a.update_e_only();
        solver_a.apply_cpml_e();
        solver_a.advance_clock();

        solver_b.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver_b.grid_mut());
        solver_b.grid_mut().ey[(src_i, src_j, src_k)] += amp * gaussian(t_b, t0, sigma_t);
        solver_b.update_e_only();
        solver_b.apply_cpml_e();
        solver_b.advance_clock();
    }

    // Compare all E field components for bit-exact equality.
    let ga = solver_a.grid();
    let gb = solver_b.grid();

    for ((va, vb), idx) in ga.ex.iter().zip(gb.ex.iter()).zip(0_usize..) {
        assert_eq!(
            va, vb,
            "ex[{}]: lossless={} vs sigma=0 lossy={}",
            idx, va, vb
        );
    }
    for ((va, vb), idx) in ga.ey.iter().zip(gb.ey.iter()).zip(0_usize..) {
        assert_eq!(
            va, vb,
            "ey[{}]: lossless={} vs sigma=0 lossy={}",
            idx, va, vb
        );
    }
    for ((va, vb), idx) in ga.ez.iter().zip(gb.ez.iter()).zip(0_usize..) {
        assert_eq!(
            va, vb,
            "ez[{}]: lossless={} vs sigma=0 lossy={}",
            idx, va, vb
        );
    }
}

/// Higher-Q variant: σ = σ₀/10 → Q_analytic ≈ 200. N_RING = 30 000 steps.
///
/// Marked `#[ignore]` because of the longer ring-down time required;
/// run explicitly with `-- --ignored` or `-- --include-ignored`.
#[test]
#[ignore = "slow: ~5 s release; fdtd-202 high-Q (Q≈200) cavity ring-down gate"]
fn fdtd_202_q_factor_hi_q_ignored() {
    const SIGMA_LO: f64 = SIGMA0 / 10.0; // 2.96e-4 S/m → Q ≈ 200
    const N_RING_HI: usize = 30_000;

    let (probe, dt) = run_lossy_cavity(SIGMA_LO, N_SRC, N_RING_HI);

    let skip = N_RING_HI / 3;
    let window = &probe[skip..];
    let t_start = (N_SRC + skip) as f64 * dt;

    let tau = fit_log_decay(window, dt, t_start);
    let f101 = analytic_f101();
    let q_measured = PI * f101 * tau;
    let q_analytic = analytic_q(SIGMA_LO);
    let rel_err = (q_measured - q_analytic).abs() / q_analytic;

    eprintln!(
        "\nfdtd-202 high-Q variant
  sigma:        {:.4e} S/m
  Q_analytic:   {:.4}
  tau_fit:      {:.6e} s
  Q_measured:   {:.4}
  rel_error:    {:.4} %
  gate:         ±5 %
",
        SIGMA_LO,
        q_analytic,
        tau,
        q_measured,
        rel_err * 100.0,
    );

    assert!(
        rel_err < 0.05,
        "fdtd-202 hi-Q FAILED: Q_measured = {:.4}, Q_analytic = {:.4}, \
         rel_error = {:.4} % (threshold ±5 %)",
        q_measured,
        q_analytic,
        rel_err * 100.0,
    );
}
