//! FDTD driver Python bindings.
//!
//! Exposes [`yee_fdtd::FdtdDriver`], [`yee_fdtd::FdtdDriverConfig`], and
//! [`yee_fdtd::RadiationPattern`] to Python so notebook / scripting callers
//! can run a short-dipole radiation-pattern simulation end-to-end and
//! receive numpy arrays for the θ-cut.
//!
//! Also exposes [`run_cavity_q`] / [`PyCavityQResult`] for the fdtd-202
//! lossy-cavity Q-factor ring-down gate (Phase 2.fdtd.py.0),
//! [`run_cavity_resonance`] / [`PyCavityResonanceResult`] for the fdtd-201
//! TE₁₀₁ resonance gate (Phase 2.fdtd.py.1),
//! [`run_dipole_pattern`] / [`PyDipolePatternResult`] for the fdtd-203
//! short-dipole sin-θ NTFF radiation-pattern gate (Phase 2.fdtd.py.2), and
//! [`run_fresnel_tfsf`] / [`PyFresnelTfsfResult`] for the fdtd-204
//! TF/SF Fresnel-transmission gate (Phase 2.fdtd.py.3).
//!
//! The Python wrapper holds only the configuration plus the grid
//! parameters (`nx`, `ny`, `nz`, `dx`); the underlying Rust
//! [`yee_fdtd::FdtdDriver`] is constructed fresh inside every
//! [`PyFdtdDriver::run`] call because the Rust `run` method consumes
//! `self`. As a result a `PyFdtdDriver` can be `.run()` multiple times,
//! each call returning an independent [`PyRadiationPattern`].

use std::f64::consts::PI;

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use yee_fdtd::{FdtdDriver as RustDriver, FdtdDriverConfig as RustCfg, YeeGrid};

/// Python wrapper for [`yee_fdtd::FdtdDriverConfig`].
///
/// All cell indices are 0-based and refer to the integer-E node lattice
/// of the underlying [`yee_fdtd::YeeGrid`]. Frequencies are in Hz.
#[pyclass(name = "FdtdDriverConfig", module = "yee._yee", from_py_object)]
#[derive(Clone)]
pub struct PyFdtdDriverConfig {
    /// Number of FDTD time steps to run.
    #[pyo3(get)]
    pub n_steps: usize,
    /// `(i, j, k)` of the dipole centre cell.
    #[pyo3(get)]
    pub dipole_center_cells: (usize, usize, usize),
    /// Length of the dipole, in cells along z.
    #[pyo3(get)]
    pub dipole_length_cells: usize,
    /// Sinusoid drive frequency in Hz. Also the NTFF DFT bin.
    #[pyo3(get)]
    pub source_freq_hz: f64,
    /// Cells between the inner CPML edge and the NTFF integration
    /// surface.
    #[pyo3(get)]
    pub ntff_surface_pad_cells: usize,
    /// CPML thickness on every face, in cells. Typical value `10`.
    #[pyo3(get)]
    pub cpml_thickness_cells: usize,
}

#[pymethods]
impl PyFdtdDriverConfig {
    /// Construct a config.
    ///
    /// Args:
    ///     n_steps: number of FDTD time steps to run.
    ///     dipole_center_cells: `(i, j, k)` integer index of the dipole
    ///         centre cell.
    ///     dipole_length_cells: dipole length along z, in cells (>= 1).
    ///     source_freq_hz: sinusoid drive frequency (Hz, > 0).
    ///     ntff_surface_pad_cells: cells between the inner CPML edge and
    ///         the NTFF integration surface. Defaults to `4`.
    ///     cpml_thickness_cells: CPML thickness on every face, in cells.
    ///         Defaults to `10`.
    #[new]
    #[pyo3(signature = (
        n_steps,
        dipole_center_cells,
        dipole_length_cells,
        source_freq_hz,
        ntff_surface_pad_cells = 4,
        cpml_thickness_cells = 10,
    ))]
    fn new(
        n_steps: usize,
        dipole_center_cells: (usize, usize, usize),
        dipole_length_cells: usize,
        source_freq_hz: f64,
        ntff_surface_pad_cells: usize,
        cpml_thickness_cells: usize,
    ) -> Self {
        Self {
            n_steps,
            dipole_center_cells,
            dipole_length_cells,
            source_freq_hz,
            ntff_surface_pad_cells,
            cpml_thickness_cells,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "FdtdDriverConfig(n_steps={}, dipole_center_cells={:?}, \
             dipole_length_cells={}, source_freq_hz={}, \
             ntff_surface_pad_cells={}, cpml_thickness_cells={})",
            self.n_steps,
            self.dipole_center_cells,
            self.dipole_length_cells,
            self.source_freq_hz,
            self.ntff_surface_pad_cells,
            self.cpml_thickness_cells,
        )
    }
}

impl PyFdtdDriverConfig {
    /// Internal: convert this config into the Rust `FdtdDriverConfig`.
    pub(crate) fn to_rust(&self) -> RustCfg {
        RustCfg {
            n_steps: self.n_steps,
            dipole_center_cells: self.dipole_center_cells,
            dipole_length_cells: self.dipole_length_cells,
            source_freq_hz: self.source_freq_hz,
            ntff_surface_pad_cells: self.ntff_surface_pad_cells,
            cpml_thickness_cells: self.cpml_thickness_cells,
        }
    }
}

/// Python wrapper for [`yee_fdtd::RadiationPattern`].
///
/// `theta_deg[i]` corresponds to `e_theta_phi0[i]`. The magnitudes are
/// normalized so that `max e_theta_phi0 == 1.0`.
#[pyclass(name = "FdtdRadiationPattern", module = "yee._yee")]
pub struct PyRadiationPattern {
    pub(crate) theta_deg: Vec<f64>,
    pub(crate) e_theta_phi0: Vec<f64>,
}

#[pymethods]
impl PyRadiationPattern {
    /// Polar angles in **degrees**, sampled from 0 to 180 inclusive in
    /// 5° steps (37 points total).
    #[getter]
    fn theta_deg<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.theta_deg.clone().into_pyarray(py)
    }

    /// `|E_θ|` at each `theta_deg[i]` at `φ = 0`, normalized so the
    /// maximum across the cut equals `1.0`.
    #[getter]
    fn e_theta_phi0<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.e_theta_phi0.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!("FdtdRadiationPattern(n_points={})", self.theta_deg.len())
    }
}

/// Python wrapper for [`yee_fdtd::FdtdDriver`].
///
/// Holds the grid dimensions plus the configuration; the underlying
/// Rust driver is built fresh on every [`PyFdtdDriver::run`] call (the
/// Rust `run` consumes `self`). This means a single `PyFdtdDriver`
/// instance can be `.run()` repeatedly — each call performs an
/// independent simulation.
#[pyclass(name = "FdtdDriver", module = "yee._yee")]
pub struct PyFdtdDriver {
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    pub(crate) nz: usize,
    pub(crate) dx: f64,
    pub(crate) cfg: PyFdtdDriverConfig,
}

#[pymethods]
impl PyFdtdDriver {
    /// Build a vacuum grid + configure the driver.
    ///
    /// Args:
    ///     nx, ny, nz: cell counts along x, y, z (all >= 1).
    ///     dx: cell size in metres (must be > 0).
    ///     cfg: `FdtdDriverConfig` describing the source, NTFF, and CPML.
    ///
    /// The grid is allocated lazily inside [`PyFdtdDriver.run`]; this
    /// constructor only stores the parameters.
    #[new]
    fn new(nx: usize, ny: usize, nz: usize, dx: f64, cfg: PyFdtdDriverConfig) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            cfg,
        }
    }

    /// Run the FDTD simulation to completion and return the far-field
    /// radiation pattern.
    ///
    /// This method does **not** consume the Python driver: a single
    /// `FdtdDriver` instance may be `.run()` multiple times, each call
    /// building a fresh underlying [`yee_fdtd::FdtdDriver`] from the
    /// stored grid and configuration parameters.
    ///
    /// Returns:
    ///     `FdtdRadiationPattern` — `.theta_deg` and `.e_theta_phi0` are
    ///     1-D `float64` numpy arrays of length 37 (5° steps from 0 to
    ///     180 inclusive).
    fn run(&self, py: Python<'_>) -> PyRadiationPattern {
        // Release the GIL while running the FDTD time loop: the inner
        // call is pure-Rust CPU work and can take several seconds for
        // realistic grids.
        py.detach(|| {
            let grid = YeeGrid::vacuum(self.nx, self.ny, self.nz, self.dx);
            let driver = RustDriver::new(grid, self.cfg.to_rust());
            let pat = driver.run();
            PyRadiationPattern {
                theta_deg: pat.theta_deg,
                e_theta_phi0: pat.e_theta_phi0,
            }
        })
    }
}

// ---------------------------------------------------------------------------
// Phase 2.fdtd.py.0 — lossy-cavity Q-factor driver
// ---------------------------------------------------------------------------

/// ε₀ in SI units (F/m). Mirrors the constant used by the fdtd-202 gate in
/// `crates/yee-fdtd/tests/cavity_q.rs`.
const CQ_EPS0: f64 = 8.854_187_817e-12;

/// Speed of light in vacuum (m/s).
const CQ_C0: f64 = 299_792_458.0;

/// Analytic TE₁₀₁ resonant frequency for a square `nx × nz` cavity at cell
/// size `dx` (Pozar §6.3):
///
/// ```text
/// f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)
/// ```
fn cq_f101(nx: usize, nz: usize, dx: f64) -> f64 {
    let a = nx as f64 * dx;
    let d = nz as f64 * dx;
    0.5 * CQ_C0 * ((1.0 / (a * a)) + (1.0 / (d * d))).sqrt()
}

/// Analytic Q-factor for a uniformly lossy PEC cavity (Taflove §3.7):
///
/// ```text
/// Q = 2π · f₁₀₁ · ε₀ / σ
/// ```
fn cq_q_analytic(nx: usize, nz: usize, dx: f64, sigma: f64) -> f64 {
    2.0 * PI * cq_f101(nx, nz, dx) * CQ_EPS0 / sigma
}

/// Gaussian pulse amplitude at time `t` centred at `t0` with width `sigma_t`.
#[inline]
fn cq_gaussian(t: f64, t0: f64, sigma_t: f64) -> f64 {
    let arg = (t - t0) / sigma_t;
    (-arg * arg).exp()
}

/// Fit an exponential decay to `series` and return the decay time constant τ.
///
/// Uses log-linear least-squares regression:
///
/// ```text
/// log|y[n]| = A + slope · t[n],   slope = −1/τ
/// ```
fn cq_fit_log_decay(series: &[f64], dt: f64, t_start: f64) -> f64 {
    let n = series.len() as f64;
    let ts: Vec<f64> = (0..series.len()).map(|i| t_start + i as f64 * dt).collect();
    let ys: Vec<f64> = series.iter().map(|&v| v.abs().max(1e-30).ln()).collect();

    let t_mean = ts.iter().sum::<f64>() / n;
    let y_mean = ys.iter().sum::<f64>() / n;

    let num: f64 = ts
        .iter()
        .zip(ys.iter())
        .map(|(&t, &y)| (t - t_mean) * (y - y_mean))
        .sum();
    let den: f64 = ts.iter().map(|&t| (t - t_mean).powi(2)).sum();

    // slope is negative (decay); τ = −1/slope.
    -1.0 / (num / den)
}

/// Result of a lossy-cavity Q-factor simulation (fdtd-202 gate).
///
/// Returned by [`run_cavity_q`]. All scalar fields are exposed as Python
/// read-only properties via `#[pyo3(get)]`. The ring-down time series is
/// available through [`PyCavityQResult::probe_array`].
#[pyclass(name = "CavityQResult", module = "yee._yee")]
pub struct PyCavityQResult {
    /// Q extracted from the log-linear ring-down fit.
    #[pyo3(get)]
    pub q_measured: f64,
    /// Analytic Q = 2π · f₁₀₁ · ε₀ / σ.
    #[pyo3(get)]
    pub q_analytic: f64,
    /// Analytic TE₁₀₁ resonant frequency (Hz).
    #[pyo3(get)]
    pub f101_hz: f64,
    /// |q_measured − q_analytic| / q_analytic.
    #[pyo3(get)]
    pub rel_err: f64,
    /// `true` iff `rel_err < 0.05` (fdtd-202 gate: ±5 %).
    #[pyo3(get)]
    pub passed: bool,
    /// Internal ring-down probe time series (n_ring samples). Not directly
    /// exposed as a Python attribute; access via [`Self::probe_array`].
    pub probe_vec: Vec<f64>,
}

#[pymethods]
impl PyCavityQResult {
    /// Return the ring-down probe time series as a 1-D numpy `float64` array.
    ///
    /// The array has length `n_ring` (the value passed to [`run_cavity_q`]).
    /// Sample `i` corresponds to the E_y probe value recorded at time step
    /// `n_src + i` during the free ring-down phase.
    pub fn probe_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.probe_vec.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "CavityQResult(q_measured={:.4}, q_analytic={:.4}, \
             rel_err={:.4e}, passed={})",
            self.q_measured, self.q_analytic, self.rel_err, self.passed,
        )
    }
}

// ---------------------------------------------------------------------------
// Phase 2.fdtd.py.1 — cavity resonance frequency driver
// ---------------------------------------------------------------------------

/// Analytic TE₁₀₁ resonant frequency for a rectangular `nx × nz` cavity at
/// cell size `dx` (Pozar §6.3):
///
/// ```text
/// f₁₀₁ = (c/2) · √((1/a)² + (1/d)²)
/// ```
fn cr_analytic_f101(nx: usize, nz: usize, dx: f64) -> f64 {
    const C0: f64 = 299_792_458.0;
    let a = nx as f64 * dx;
    let d = nz as f64 * dx;
    0.5 * C0 * ((1.0 / (a * a)) + (1.0 / (d * d))).sqrt()
}

/// Result of a rectangular PEC cavity resonance simulation (fdtd-201 gate).
///
/// Returned by [`run_cavity_resonance`]. All scalar fields are exposed as
/// Python read-only properties via `#[pyo3(get)]`. The full probe time series
/// is available through [`PyCavityResonanceResult::probe_array`].
#[pyclass(name = "CavityResonanceResult", module = "yee._yee")]
pub struct PyCavityResonanceResult {
    /// Peak frequency extracted from the DFT magnitude scan (Hz).
    #[pyo3(get)]
    pub f_extracted_hz: f64,
    /// Analytic TE₁₀₁ resonant frequency from Pozar §6.3 (Hz).
    #[pyo3(get)]
    pub f_analytic_hz: f64,
    /// `|f_extracted_hz − f_analytic_hz| / f_analytic_hz`.
    #[pyo3(get)]
    pub rel_err: f64,
    /// `true` iff `rel_err < 0.025` (fdtd-201 gate: ±2.5 %).
    #[pyo3(get)]
    pub passed: bool,
    /// Full probe time series (n_steps samples). Not directly exposed as a
    /// Python attribute; access via [`Self::probe_array`].
    pub probe_vec: Vec<f64>,
}

#[pymethods]
impl PyCavityResonanceResult {
    /// Return the full probe time series as a 1-D numpy `float64` array.
    ///
    /// The array has length `n_steps` (the value passed to
    /// [`run_cavity_resonance`]). Sample `i` corresponds to the E_y probe
    /// value recorded at time step `i` of the simulation.
    pub fn probe_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.probe_vec.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "CavityResonanceResult(f_extracted_hz={:.6e}, f_analytic_hz={:.6e}, \
             rel_err={:.4e}, passed={})",
            self.f_extracted_hz, self.f_analytic_hz, self.rel_err, self.passed,
        )
    }
}

/// Run a rectangular PEC cavity simulation and extract the TE₁₀₁ resonant
/// frequency via a broadband DFT scan.
///
/// Builds a vacuum grid of `nx × ny × nz` cells at cell size `dx` metres,
/// injects a broadband Gaussian pulse into E_y, runs for `n_steps` time steps
/// with PEC boundary conditions, and scans `n_freq_bins` candidate frequencies
/// in `[0.65·f_ref, 1.50·f_ref]` (where `f_ref` is the analytic TE₁₀₁
/// frequency) to find the spectral peak. Returns the frequency of maximum DFT
/// magnitude alongside the analytic reference and the relative error.
///
/// # Arguments
///
/// * `nx` — cells in x (default 20)
/// * `ny` — cells in y (default 10)
/// * `nz` — cells in z (default 20)
/// * `dx` — cell size in metres (default 0.01 = 10 mm)
/// * `n_steps` — number of FDTD time steps (default 30 000)
/// * `n_freq_bins` — number of DFT candidate frequencies (default 400)
///
/// # Returns
///
/// A [`CavityResonanceResult`] with `.f_extracted_hz`, `.f_analytic_hz`,
/// `.rel_err`, `.passed`, and `.probe_array()`.
#[pyfunction]
#[pyo3(signature = (
    nx = 20,
    ny = 10,
    nz = 20,
    dx = 0.01,
    n_steps = 30_000usize,
    n_freq_bins = 400usize,
))]
pub fn run_cavity_resonance(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    n_steps: usize,
    n_freq_bins: usize,
) -> PyCavityResonanceResult {
    use yee_fdtd::boundary;
    use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver};

    let grid = YeeGrid::vacuum(nx, ny, nz, dx);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    // Source: (nx/4, ny/2, nz/4) in ey — off-centre for TE₁₀₁ coupling.
    let src_i = nx / 4;
    let src_j = ny / 2;
    let src_k = nz / 4;
    // Probe: (3*nx/4, ny/2, 3*nz/4) — opposite quarter, reduces source
    // near-field bias.
    let prb_i = nx * 3 / 4;
    let prb_j = ny / 2;
    let prb_k = nz * 3 / 4;

    // Gaussian parameters: t0 = 12·dt, σ_t = 4·dt.
    let t0 = 12.0 * dt;
    let sigma_t = 4.0 * dt;

    let mut probe_series: Vec<f64> = Vec::with_capacity(n_steps);

    for _ in 0..n_steps {
        let t = solver.current_time();
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        {
            let arg = (t - t0) / sigma_t;
            solver.grid_mut().ey[(src_i, src_j, src_k)] += (-arg * arg).exp();
        }
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
        probe_series.push(solver.grid().ey[(prb_i, prb_j, prb_k)]);
    }

    // DFT scan: n_freq_bins candidates in [0.65·f_ref, 1.50·f_ref].
    let f_ref = cr_analytic_f101(nx, nz, dx);
    let f_lo = 0.65 * f_ref;
    let f_hi = 1.50 * f_ref;
    let df_scan = (f_hi - f_lo) / (n_freq_bins - 1) as f64;

    let mut peak_power = 0.0_f64;
    let mut peak_freq = f_lo;

    for bin in 0..n_freq_bins {
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

    let rel_err = (peak_freq - f_ref).abs() / f_ref;
    PyCavityResonanceResult {
        f_extracted_hz: peak_freq,
        f_analytic_hz: f_ref,
        rel_err,
        passed: rel_err < 0.025,
        probe_vec: probe_series,
    }
}

// ---------------------------------------------------------------------------
// Phase 2.fdtd.py.2 — short-dipole NTFF radiation-pattern driver
// ---------------------------------------------------------------------------

/// Result of a short-dipole NTFF radiation-pattern simulation (fdtd-203 gate).
///
/// Returned by [`run_dipole_pattern`]. The five scalar fields are exposed as Python
/// read-only properties. The full θ-cut is available through `.theta_deg_array()` and
/// `.e_theta_array()`.
#[pyclass(name = "DipolePatternResult", module = "yee._yee")]
pub struct PyDipolePatternResult {
    /// `|E_θ|` at θ = 0° (endfire null; gate requires < 0.15).
    #[pyo3(get)]
    pub e_theta_0: f64,
    /// `|E_θ|` at θ = 45° (gate: |e − 0.707| < 0.15).
    #[pyo3(get)]
    pub e_theta_45: f64,
    /// `|E_θ|` at θ = 90° (broadside peak; normalized to 1.0; gate: |e − 1.0| < 0.05).
    #[pyo3(get)]
    pub e_theta_90: f64,
    /// `|E_θ|` at θ = 135° (gate: |e − 0.707| < 0.15).
    #[pyo3(get)]
    pub e_theta_135: f64,
    /// `|E_θ|` at θ = 180° (endfire null; gate requires < 0.15).
    #[pyo3(get)]
    pub e_theta_180: f64,
    /// `true` iff all five gate criteria pass (fdtd-203, Balanis §4.2 sin θ).
    #[pyo3(get)]
    pub passed: bool,
    /// Internal: full θ array (private; access via `theta_deg_array()`).
    pub theta_deg_vec: Vec<f64>,
    /// Internal: full `|E_θ|` array (private; access via `e_theta_array()`).
    pub e_theta_vec: Vec<f64>,
}

#[pymethods]
impl PyDipolePatternResult {
    /// Return the full θ array as a 1-D numpy `float64` array (37 points, 0°–180° in 5° steps).
    pub fn theta_deg_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.theta_deg_vec.clone().into_pyarray(py)
    }

    /// Return the full `|E_θ|` array as a 1-D numpy `float64` array (normalized to max=1.0).
    pub fn e_theta_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.e_theta_vec.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "DipolePatternResult(e_theta_0={:.4}, e_theta_90={:.4}, \
             e_theta_180={:.4}, passed={})",
            self.e_theta_0, self.e_theta_90, self.e_theta_180, self.passed,
        )
    }
}

/// Run a short-dipole FDTD simulation and return the far-field radiation pattern.
///
/// Builds a vacuum grid of `nx × ny × nz` cells at cell size `dx` metres, drives a
/// z-polarised sinusoidal current source over a 5-cell dipole at the grid centre,
/// runs for `n_steps` time steps with CPML on all six faces, and returns the NTFF
/// θ-cut of `|E_θ|` at `φ = 0` normalized to its maximum.
///
/// Default scenario matches the fdtd-203 gate (60³, dx=5mm, 800 steps, 1 GHz).
///
/// # Reference
///
/// C. A. Balanis, *Antenna Theory*, 4th ed., §4.2 — short-dipole far-field
/// `E_θ ∝ sin θ` (maximum at θ=90°, nulls at θ=0°/180°).
///
/// # Gate criteria (fdtd-203)
///
/// | θ    | Expected | Gate               |
/// |------|----------|--------------------|
/// |  0°  | 0        | e < 0.15           |
/// | 45°  | 0.707    | |e − 0.707| < 0.15 |
/// | 90°  | 1.0      | |e − 1.0|  < 0.05  |
/// | 135° | 0.707    | |e − 0.707| < 0.15 |
/// | 180° | 0        | e < 0.15           |
#[pyfunction]
#[pyo3(signature = (
    nx = 60usize,
    ny = 60usize,
    nz = 60usize,
    dx = 5.0e-3_f64,
    n_steps = 800usize,
    source_freq_hz = 1.0e9_f64,
))]
pub fn run_dipole_pattern(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    n_steps: usize,
    source_freq_hz: f64,
    py: Python<'_>,
) -> PyDipolePatternResult {
    let (theta_deg, e_theta_phi0) = py.detach(|| {
        let grid = YeeGrid::vacuum(nx, ny, nz, dx);
        let ci = nx / 2;
        let cj = ny / 2;
        let ck = nz / 2;
        let cfg = RustCfg {
            n_steps,
            dipole_center_cells: (ci, cj, ck),
            dipole_length_cells: 5,
            source_freq_hz,
            ntff_surface_pad_cells: 4,
            cpml_thickness_cells: 10,
        };
        let pat = RustDriver::new(grid, cfg).run();
        (pat.theta_deg, pat.e_theta_phi0)
    });

    // Helper closure: find the nearest θ in the discrete array to `target_deg`.
    let sample = |target_deg: f64| -> f64 {
        theta_deg
            .iter()
            .zip(e_theta_phi0.iter())
            .min_by(|a, b| {
                (a.0 - target_deg)
                    .abs()
                    .partial_cmp(&(b.0 - target_deg).abs())
                    .unwrap()
            })
            .map(|(_, &v)| v)
            .unwrap_or(0.0)
    };

    let e0 = sample(0.0);
    let e45 = sample(45.0);
    let e90 = sample(90.0);
    let e135 = sample(135.0);
    let e180 = sample(180.0);

    let passed = e0 < 0.15
        && (e45 - 0.707_f64).abs() < 0.15
        && (e90 - 1.0_f64).abs() < 0.05
        && (e135 - 0.707_f64).abs() < 0.15
        && e180 < 0.15;

    PyDipolePatternResult {
        e_theta_0: e0,
        e_theta_45: e45,
        e_theta_90: e90,
        e_theta_135: e135,
        e_theta_180: e180,
        passed,
        theta_deg_vec: theta_deg,
        e_theta_vec: e_theta_phi0,
    }
}

/// Run a lossy rectangular PEC cavity simulation and return the Q-factor.
///
/// Builds a vacuum grid of `nx × ny × nz` cells at cell size `dx` metres,
/// fills it uniformly with electric conductivity `sigma` S/m (via the
/// CA/CB Yee E-update, Taflove §3.7), injects a broadband Gaussian pulse
/// into E_y for `n_src` steps, records the ring-down for `n_ring` steps,
/// fits an exponential decay to the last 2/3 of the ring-down, and
/// extracts Q = π · f₁₀₁ · τ.
///
/// # Arguments
///
/// * `nx` — cells in x (default 20)
/// * `ny` — cells in y (default 10)
/// * `nz` — cells in z (default 20)
/// * `dx` — cell size in metres (default 0.01 = 10 mm)
/// * `sigma` — uniform electric conductivity in S/m (default 2.96e-3,
///   giving Q ≈ 20)
/// * `n_src` — number of source-injection steps (default 200)
/// * `n_ring` — number of ring-down steps to record (default 6000)
///
/// # Returns
///
/// A [`CavityQResult`] with `.q_measured`, `.q_analytic`, `.f101_hz`,
/// `.rel_err`, `.passed`, and `.probe_array()`.
#[pyfunction]
#[pyo3(signature = (
    nx = 20,
    ny = 10,
    nz = 20,
    dx = 0.01,
    sigma = 2.96e-3_f64,
    n_src = 200,
    n_ring = 6000,
))]
pub fn run_cavity_q(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    sigma: f64,
    n_src: usize,
    n_ring: usize,
) -> PyCavityQResult {
    use yee_fdtd::boundary;
    use yee_fdtd::{FdtdSolver, WalkingSkeletonSolver};

    // Build a vacuum grid and attach uniform conductivity over the full domain.
    // set_sigma_box uses inclusive-exclusive indexing; NX+1 etc. covers the
    // full [nx+1, ny+1, nz+1] staggered extent.
    let mut grid = YeeGrid::vacuum(nx, ny, nz, dx);
    grid.set_sigma_box(0, nx + 1, 0, ny + 1, 0, nz + 1, sigma);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    // Source: off-centre for TE₁₀₁ coupling (sin(π/4) ≈ 0.707).
    let src_i = nx / 4;
    let src_j = ny / 2;
    let src_k = nz / 4;

    // Probe: opposite quarter to reduce source near-field bias.
    let prb_i = nx * 3 / 4;
    let prb_j = ny / 2;
    let prb_k = nz * 3 / 4;

    // Gaussian parameters: t0 = 12·dt, σ_t = 4·dt.
    let t0 = 12.0 * dt;
    let sigma_t = 4.0 * dt;

    // Phase 1: inject Gaussian pulse into E_y for n_src steps.
    for _ in 0..n_src {
        let t = solver.current_time();
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.grid_mut().ey[(src_i, src_j, src_k)] += cq_gaussian(t, t0, sigma_t);
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
    }

    // Phase 2: ring-down, record E_y probe.
    let mut probe = Vec::with_capacity(n_ring);
    for _ in 0..n_ring {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
        probe.push(solver.grid().ey[(prb_i, prb_j, prb_k)]);
    }

    // Fit log-decay to the last 2/3 of the ring-down window.
    let skip = n_ring / 3;
    let window = &probe[skip..];
    let t_start = (n_src + skip) as f64 * dt;
    let tau = cq_fit_log_decay(window, dt, t_start);

    let f101 = cq_f101(nx, nz, dx);
    let q_measured = PI * f101 * tau;
    let q_analytic = cq_q_analytic(nx, nz, dx, sigma);
    let rel_err = (q_measured - q_analytic).abs() / q_analytic;

    PyCavityQResult {
        q_measured,
        q_analytic,
        f101_hz: f101,
        rel_err,
        passed: rel_err < 0.05,
        probe_vec: probe,
    }
}

// ---------------------------------------------------------------------------
// Phase 2.fdtd.py.3 — TF/SF Fresnel-transmission driver (fdtd-204 gate)
// ---------------------------------------------------------------------------

/// Result of a TF/SF Fresnel-transmission simulation (fdtd-204 gate).
///
/// The simulation runs a normal-incidence plane wave through a dielectric
/// slab (ε_r = 2.2, thickness 5 mm) at 10 GHz and measures the amplitude
/// transmission coefficient.
#[pyclass(name = "FresnelTfsfResult", module = "yee._yee")]
pub struct PyFresnelTfsfResult {
    /// Measured amplitude transmission coefficient (A_trans / A_inc).
    #[pyo3(get)]
    pub t_measured: f64,
    /// Analytic amplitude transmission coefficient from the Born-Wolf
    /// transfer-matrix formula for the slab geometry.
    #[pyo3(get)]
    pub t_analytic: f64,
    /// Relative error |t_measured − t_analytic| / t_analytic.
    #[pyo3(get)]
    pub rel_err: f64,
    /// `true` iff `rel_err < 0.05` (fdtd-204 gate: ±5 %).
    #[pyo3(get)]
    pub passed: bool,
}

#[pymethods]
impl PyFresnelTfsfResult {
    fn __repr__(&self) -> String {
        format!(
            "FresnelTfsfResult(t_measured={:.4}, t_analytic={:.4}, \
             rel_err={:.4}, passed={})",
            self.t_measured, self.t_analytic, self.rel_err, self.passed
        )
    }
}

/// Run the fdtd-204 TF/SF Fresnel-transmission gate from Python.
///
/// Default arguments reproduce the published fdtd-204 scenario:
/// 80³ grid, ε_r=2.2 slab (5 cells × 1 mm), f=10 GHz, 600 steps.
///
/// # Gate criteria (fdtd-204, Born & Wolf §1.6.2)
/// `|t_measured / t_analytic − 1| < 0.05`  (5 %)
///
/// # Example
/// ```python
/// from yee import run_fresnel_tfsf
/// result = run_fresnel_tfsf()
/// assert result.passed, f"fdtd-204 gate failed: rel_err={result.rel_err:.3f}"
/// ```
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (
    *,
    nx = 80_usize,
    ny = 80_usize,
    nz = 80_usize,
    dx = 1.0e-3_f64,
    eps_r = 2.2_f64,
    slab_i0 = 50_usize,
    slab_i1 = 55_usize,
    freq_hz = 10.0e9_f64,
    n_steps = 600_usize,
    settle  = 200_usize,
))]
pub fn run_fresnel_tfsf(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    eps_r: f64,
    slab_i0: usize,
    slab_i1: usize,
    freq_hz: f64,
    n_steps: usize,
    settle: usize,
) -> PyFresnelTfsfResult {
    use ndarray::Array3;
    use yee_fdtd::{CpmlParams, PlaneWaveDirection, PlaneWaveSource, WalkingSkeletonSolver};

    let npml: usize = 10;

    // TF box: i₀=12, i₁=nx-npml-1, j₀=1, j₁=ny-2, k₀=1, k₁=nz-2
    let tf_i0: usize = 12;
    let tf_i1: usize = nx.saturating_sub(npml + 1);
    let tf_j0: usize = 1;
    let tf_j1: usize = ny.saturating_sub(2);
    let tf_k0: usize = 1;
    let tf_k1: usize = nz.saturating_sub(2);

    // Probes: incident (vacuum, inside TF, before slab);
    //         transmitted (vacuum, inside TF, after slab).
    let probe_inc = (25_usize.min(nx - 1), ny / 2, nz / 2);
    let probe_trans = (62_usize.min(nx - 1), ny / 2, nz / 2);

    // Helper: run one simulation (vacuum or slab), return E_z trace at probe.
    let run_sim = |with_slab: bool, probe: (usize, usize, usize)| -> Vec<f64> {
        let grid = if with_slab {
            let mut eps_cells = Array3::<f64>::from_elem((nx + 1, ny + 1, nz + 1), 1.0);
            for i in slab_i0..slab_i1 {
                for j in 0..=ny {
                    for k in 0..=nz {
                        eps_cells[(i, j, k)] = eps_r;
                    }
                }
            }
            YeeGrid::vacuum(nx, ny, nz, dx).with_eps_r_cells(eps_cells)
        } else {
            YeeGrid::vacuum(nx, ny, nz, dx)
        };
        let dt = grid.dt;
        let params = CpmlParams::for_grid(&grid, npml);
        let mut solver = WalkingSkeletonSolver::with_cpml(grid, params);
        let mut pw = PlaneWaveSource::new(
            tf_i0,
            tf_i1,
            tf_j0,
            tf_j1,
            tf_k0,
            tf_k1,
            PlaneWaveDirection::PlusX,
            freq_hz,
            50,
            dx,
            dt,
            4,
        );
        let mut trace = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            solver.step_with_plane_wave(&mut pw);
            trace.push(solver.grid().ez[probe]);
        }
        trace
    };

    // Vacuum run: measures incident amplitude at probe_inc.
    let trace_vac = run_sim(false, probe_inc);
    // Slab run: measures transmitted amplitude at probe_trans.
    let trace_slab = run_sim(true, probe_trans);

    // Peak amplitude after settling.
    let settle_clamp = settle.min(n_steps);
    let a_inc = trace_vac[settle_clamp..]
        .iter()
        .cloned()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    let a_trans = trace_slab[settle_clamp..]
        .iter()
        .cloned()
        .fold(0.0_f64, |a, v| a.max(v.abs()));
    let t_measured = if a_inc > 0.0 { a_trans / a_inc } else { 0.0 };

    // Analytic Born-Wolf §1.6.2 transfer-matrix via the canonical yee-validation
    // reference implementation (yee-validation is already a dep of yee-py).
    let t_analytic = yee_validation::fdtd204_t_analytic(eps_r, slab_i1 - slab_i0, dx, freq_hz);

    let rel_err = if t_analytic > 0.0 {
        (t_measured - t_analytic).abs() / t_analytic
    } else {
        0.0
    };

    PyFresnelTfsfResult {
        t_measured,
        t_analytic,
        rel_err,
        passed: rel_err < 0.05,
    }
}
