//! Validation gate **fdtd-201.x** — rectangular PEC-cavity TE₂₀₁ resonant
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
//! The **higher-order** mode under test is **TE₂₀₁** (`m=2, n=0, p=1`):
//!
//! ```text
//! f₂₀₁ = (c / 2) * sqrt((2/a)² + (1/d)²)
//! ```
//!
//! # Why a ≠ d (non-degeneracy requirement)
//!
//! For a box with `a = d` (i.e. equal x and z extents), TE₂₀₁ and TE₁₀₂
//! are **degenerate**: both satisfy `sqrt((2/a)² + (1/a)²) = sqrt(5)/a`, so
//! a single peak in the DFT cannot be attributed to a *named* mode.  We
//! break the degeneracy by choosing `a ≠ d`: with `a = 0.24 m, d = 0.16 m`
//! the two frequencies are:
//!
//! ```text
//! f₂₀₁ = 1.5614 GHz   (a=0.24, d=0.16)
//! f₁₀₂ = 1.9751 GHz   (a=0.24, d=0.16)
//! ```
//!
//! — separated by 26.5 %, making TE₂₀₁ cleanly identifiable.
//!
//! # Cavity geometry and mode isolation
//!
//! The box is chosen so TE₂₀₁ is the **sole** mode in the DFT scan band
//! `[1.40, 1.72] GHz`.  The nearby-mode table (all TE₋ₘ₀ₚ resonances,
//! TE₋ₘₙₚ with n=0 being the relevant family; n≥1 modes are pushed well
//! above 3.7 GHz by the narrow-y dimension b=0.04 m):
//!
//! ```text
//! TE₁₀₁:  1.1260 GHz  (dominant mode, well below band floor 1.40 GHz)
//! TE₂₀₁:  1.5614 GHz  *** TARGET — sole peak in [1.40, 1.72] GHz ***
//! TE₁₀₂:  1.9751 GHz  (well above band ceiling 1.72 GHz)
//! TE₃₀₁:  2.0949 GHz
//! TM₁₁₀:  3.7991 GHz  (n=1 family, pushed high by small b; p=0 ⇒ TM, not TE)
//! ```
//!
//! The scan band `[1.40, 1.72] GHz` brackets TE₂₀₁ by ±(10–18)% and
//! excludes every other resonance.
//!
//! # Mode structure and source alignment
//!
//! For TE₂₀₁ the dominant E-field component is **E_y**:
//!
//! ```text
//! E_y  ∝  sin(2π·x/a) · sin(π·z/d)
//! ```
//!
//! The antinodes sit at `x = a/4` and `x = 3a/4`, `z = d/2`.  We place
//! the Gaussian source at `(a/4, b/2, d/2)` (strong coupling, sin ≈ 0.99)
//! and the E_y probe at the symmetric antinode `(3a/4, b/2, d/2)` so that
//! TE₂₀₁ is the dominant resonance in the probe trace.
//!
//! Importantly, the source x-position `a/4 = 0.06 m` is off the TE₁₀₁
//! antinode at `x = a/2 = 0.12 m` (sin(π·a/4·/a) = sin(π/4) ≈ 0.71 ≠ 1).
//! TE₁₀₁ is still excited (broadband injection), but the scan band is
//! chosen to exclude it entirely.
//!
//! # Grid dispersion note
//!
//! At `dx = 10 mm` and `f₂₀₁ = 1.5614 GHz`, the wavelength is
//! `λ = c/f₂₀₁ ≈ 0.192 m` → ~19.2 cells/λ.  The 3D Yee-FDTD numerical
//! phase velocity error scales as `−(π/N)²/3`, giving an estimated
//! dispersion of **~−0.9 %** at this resolution — comfortably within the
//! ±2.5 % gate.  The TE₁₀₁ fdtd-201 gate operates at ~28 cells/λ; the
//! higher frequency here reduces the effective grid sampling, which is
//! precisely the dispersion sensitivity being validated.
//!
//! # Method
//!
//! 1. Build the cavity via `YeeGrid::vacuum(nx, ny, nz, dx)` with hard PEC
//!    outer walls (the `WalkingSkeletonSolver::new` default).
//! 2. Inject an off-centre Gaussian pulse into E_y via the public grid field.
//! 3. Step N times using the custom-body pattern, recording an E_y probe
//!    time series at an interior point.
//! 4. Extract the resonant frequency by scanning a dense candidate grid of
//!    single-bin DFTs over `[1.40, 1.72] GHz`, then peak-finding.
//!    No FFT library is needed; mirrors the `ntff.rs:253` idiom and the
//!    fdtd-201 harness in `cavity_resonance.rs`.
//!
//! # Tolerance
//!
//! The gate asserts the extracted resonance matches the analytic TE₂₀₁
//! frequency within **±2.5 %**.  Grid dispersion at ~19 cells/λ contributes
//! roughly −0.9 % numerical phase error; the ±2.5 % band gives comfortable
//! margin while being clearly non-trivial (the next mode in the band is >18 %
//! away).
//!
//! ## Strict ±0.5 % refinement path
//!
//! To tighten to ±0.5 %, halve `dx` to 5 mm (doubling `nx/ny/nz`) and
//! double `N_STEPS` to 60 000.  With ~38 cells/λ at f₂₀₁ the dispersion
//! error falls below 0.3 % and the DFT frequency resolution is ≈0.06 % of
//! f₂₀₁.  Q-factor extraction (damped-exponential decay fitting) is deferred
//! as a follow-on, as in fdtd-201.
//!
//! # Wall-time budget
//!
//! Grid: `24 × 4 × 16 = 1 536` cells; `30 000` steps; 400-candidate DFT
//! scan.  Observed wall-time: similar to fdtd-201 (~5–15 s release).
//! Gated with `#[ignore]` like sibling slow integration tests.
//!
//! # Running
//!
//! ```bash
//! cargo test -p yee-fdtd --test cavity_higher_mode --release -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_fdtd::boundary;
use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver, YeeGrid};

// --------------------------------------------------------------------------
// Cavity geometry — a ≠ d to break TE₂₀₁ / TE₁₀₂ degeneracy.
// --------------------------------------------------------------------------
//
// Physical: a = 0.24 m (x), b = 0.04 m (y), d = 0.16 m (z).
//
// Grid: dx = 0.010 m → nx = 24, ny = 4, nz = 16.  Cells stay cubic;
// only the counts differ from the fdtd-201 geometry.
//
// b is kept narrow (ny = 4 → b = 0.04 m) so that all n≥1 modes
// (TM₁₁₀, TE₁₁₁, TE₀₁₁, …) are pushed above 3.7 GHz, leaving the
// scan band [1.40, 1.72] GHz occupied by TE₂₀₁ alone.
//
// Note: YeeGrid::vacuum(nx, ny, nz, dx) creates `nx×ny×nz` primary cells.
// The PEC wall clamps tangential E on faces at [0] and [nx], so the
// physical interior dimension is a = nx·dx.
//
const NX: usize = 24;
const NY: usize = 4;
const NZ: usize = 16;
const DX: f64 = 0.010; // metres

// Physical interior dimensions (= N_cells · dx):
const PHYS_A: f64 = NX as f64 * DX; // 0.24 m  (x)
const PHYS_B: f64 = NY as f64 * DX; // 0.04 m  (y)
const PHYS_D: f64 = NZ as f64 * DX; // 0.16 m  (z)

// --------------------------------------------------------------------------
// Run parameters
// --------------------------------------------------------------------------

/// Total number of FDTD time steps.
///
/// At dt ≈ 1.733 × 10⁻¹¹ s (0.9 × CFL on 10 mm cubic cells), 30 000 steps
/// gives a total simulation time of ~0.520 µs.
/// Frequency resolution: Δf = 1 / (N · dt) ≈ 1.92 MHz ≈ 0.123 % of f₂₀₁.
const N_STEPS: usize = 30_000;

/// Source cell index into `ey` (shape `[nx+1, ny, nz+1]`).
///
/// TE₂₀₁ has E_y ∝ sin(2π·x/a) · sin(π·z/d).
///
/// The antinode nearest x = a/4 = 0.06 m:
/// - E_y staggering: ey[i,j,k] sits at (i·dx, (j+0.5)·dy, k·dz)
///   (integer i,k; half-integer j — E_y lives on y-parallel edges).
/// - Index i=6: x = 6·dx = 0.060 m (= a/4);  sin(2π·0.060/0.24) = 1.000.
///
/// The antinode at z = d/2 = 0.08 m:
/// - Index k=7: z = 7·dx = 0.070 m;  sin(π·0.070/0.16) ≈ 0.981.
///
/// Index j=2 places the source in the interior. E_y is NORMAL to the
/// y = 0 / y = b walls, so `apply_pec` does not clamp it there (any j is
/// valid); TE₂₀₁ is uniform in y anyway since n = 0.
const SRC_I: usize = NX / 4; // 6  ← x = 6·dx = 0.060 m = a/4; sin(2πx/a) = 1.000
const SRC_J: usize = NY / 2; // 2
const SRC_K: usize = NZ / 2 - 1; // 7  ← z = 7·dx = 0.070 m ≈ d/2 (0.080 m); sin(πz/d) ≈ 0.981

/// Probe cell for E_y at the symmetric antinode (3a/4, b/2, d/2).
///
/// - Index i=17: x = 17·dx = 0.170 m;  sin(2π·0.170/0.24) ≈ −0.966.
///   (Opposite sign to source; TE₂₀₁ drives both in phase for the standing wave.)
/// - Index k=7: same z as source (z ≈ d/2 antinode).
///
/// Source and probe are at symmetric antinodes of TE₂₀₁; TE₁₀₁ (antinode
/// at x = a/2) still couples to both but is excluded by the scan band.
const PRB_I: usize = NX * 3 / 4 - 1; // 17  ← x = 17·dx = 0.170 m ≈ 3a/4 (0.180 m); sin(2πx/a) ≈ -0.966
const PRB_J: usize = NY / 2; // 2
const PRB_K: usize = NZ / 2 - 1; // 7

/// Number of candidate frequencies in the DFT scan.
const N_FREQ_BINS: usize = 400;

// --------------------------------------------------------------------------
// Speed of light (matches yee-core::units::C0)
// --------------------------------------------------------------------------
const C0: f64 = 299_792_458.0; // m/s

/// Analytic TE₂₀₁ resonant frequency (Pozar §6.3):
///
/// ```text
/// f₂₀₁ = (c/2) · √((2/a)² + (1/d)²)
/// ```
fn analytic_f201() -> f64 {
    0.5 * C0 * ((4.0 / (PHYS_A * PHYS_A)) + (1.0 / (PHYS_D * PHYS_D))).sqrt()
}

/// Analytic TE₁₀₁ resonant frequency (reference for nearby-mode table):
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

/// `fdtd-201.x` gate: extract TE₂₀₁ resonant frequency via single-bin DFT scan.
///
/// Validates the higher-order-mode selectivity and higher-frequency grid-dispersion
/// behaviour of the FDTD solver — a distinct claim from the dominant-mode fdtd-201
/// gate (`cavity_resonance.rs`).
///
/// Gate tolerance: ±2.5 %.
/// See module docstring for the ±0.5 % refinement path.
///
/// Nearby modes (geometry a=0.24 m, b=0.04 m, d=0.16 m; all TE₋ₘ₀ₚ):
/// - TE₁₀₁: 1.1260 GHz (dominant mode, 28 % below band floor)
/// - TE₂₀₁: 1.5614 GHz (TARGET — sole peak in [1.40, 1.72] GHz)
/// - TE₁₀₂: 1.9751 GHz (15 % above band ceiling)
/// - TE₃₀₁: 2.0949 GHz
///
/// Run with:
/// ```bash
/// cargo test -p yee-fdtd --test cavity_higher_mode --release -- --ignored --nocapture
/// ```
#[test]
#[ignore = "slow: ~5-15 s release; fdtd-201.x TE201 higher-cavity-mode gate (Phase 2.fdtd)"]
fn te201_resonance_matches_analytic_within_two_point_five_percent() {
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
    // enough to cover f₂₀₁ and several cavity modes.
    //
    // σ = 4·dt → bandwidth ~1/(2πσ) ≈ 2.3 GHz, easily covering f₂₀₁.
    // t0 = 12·dt → Gaussian tail at t=0 is e^{-(12/4)²} ≈ 10⁻⁷.
    // ----------------------------------------------------------------
    let t0 = 12.0 * dt;
    let sigma = 4.0 * dt;

    // ----------------------------------------------------------------
    // Run the cavity using a custom step body that injects E_y.
    //
    // This mirrors the pattern in `cavity_resonance.rs` (fdtd-201) and
    // `tests/lumped_resistor.rs`: we call the individual sub-step
    // helpers (update_h_only, apply_pec, update_e_only, apply_cpml_e,
    // advance_clock) and insert the source injection between H and E.
    //
    // The E_y field at ey[(i, j, k)] is public on YeeGrid, so writing
    // to it is a pure consumer call — no src/ change required.
    // ----------------------------------------------------------------
    let mut probe_series: Vec<f64> = Vec::with_capacity(N_STEPS);

    for _n in 0..N_STEPS {
        let t = solver.current_time();

        // H update.
        solver.update_h_only();
        // Outer-wall boundary: `apply_pec` zeroes tangential E on the six
        // outer faces — the H-half-step boundary.
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());

        // Inject broadband Gaussian into E_y at the TE₂₀₁ antinode.
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

        // Record probe at the symmetric TE₂₀₁ antinode.
        probe_series.push(solver.grid().ey[(PRB_I, PRB_J, PRB_K)]);
    }

    // ----------------------------------------------------------------
    // Frequency scan via single-bin DFT (Goertzel accumulation).
    //
    // For each candidate frequency f we compute:
    //
    //   |DFT(f)|² = (Σ_n x[n]·cos(ω·n·dt))² + (Σ_n x[n]·sin(ω·n·dt))²
    //
    // Scan band [1.40, 1.72] GHz is chosen to bracket TE₂₀₁ (1.5614 GHz)
    // and exclude all other resonances:
    //  - TE₁₀₁ at 1.126 GHz is 24 % below the band floor.
    //  - TE₁₀₂ at 1.975 GHz is 15 % above the band ceiling.
    // ----------------------------------------------------------------
    let f_lo = 1.40e9_f64; // 1.40 GHz — 24 % above TE₁₀₁, 10 % below TE₂₀₁
    let f_hi = 1.72e9_f64; // 1.72 GHz — 10 % above TE₂₀₁, 15 % below TE₁₀₂
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
    let f_ref = analytic_f201();
    let f_ref_101 = analytic_f101();
    let f_ref_102 = 0.5 * C0 * ((1.0 / (PHYS_A * PHYS_A)) + (4.0 / (PHYS_D * PHYS_D))).sqrt();
    let rel_error = (peak_freq - f_ref) / f_ref;
    eprintln!(
        "\nfdtd-201.x TE₂₀₁ higher-cavity-mode resonance gate
  cavity:         a = {:.4} m, b = {:.4} m, d = {:.4} m
  grid:           {}×{}×{}, dx = {:.1} mm, dt = {:.4e} s
  steps:          {} (T_total = {:.4e} s)
  DFT scan:       {N_FREQ_BINS} bins in [{:.4} GHz, {:.4} GHz]

  nearby modes (TE_m0p, a=0.24 m, b=0.04 m, d=0.16 m):
    TE₁₀₁ = {:.6} GHz  (below band, dominant)
    TE₂₀₁ = {:.6} GHz  (TARGET — sole peak in scan band)
    TE₁₀₂ = {:.6} GHz  (above band)

  analytic f₂₀₁:  {:.6} GHz
  extracted f:    {:.6} GHz  (|DFT|² = {:.3e})
  relative error: {:+.4} %
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
        f_ref_101 * 1e-9,
        f_ref * 1e-9,
        f_ref_102 * 1e-9,
        f_ref * 1e-9,
        peak_freq * 1e-9,
        peak_power,
        rel_error * 100.0,
    );

    // ----------------------------------------------------------------
    // Gate: |rel_error| ≤ 2.5 %.
    //
    // Grid dispersion on a ~19-cells/λ Yee mesh at f₂₀₁ contributes an
    // estimated −0.9 % numerical phase velocity error.  The ±2.5 % band
    // is the appropriate loose gate for this resolution; the ±0.5 %
    // refinement path requires halving dx (see module docstring).
    // ----------------------------------------------------------------
    assert!(
        rel_error.abs() < 0.025,
        "fdtd-201.x FAILED: extracted f₂₀₁ = {:.6} GHz, analytic = {:.6} GHz, \
         rel_error = {:+.4} % (threshold ±2.5 %)",
        peak_freq * 1e-9,
        f_ref * 1e-9,
        rel_error * 100.0,
    );
}
