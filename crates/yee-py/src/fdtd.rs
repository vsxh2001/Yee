//! FDTD driver Python bindings.
//!
//! Exposes [`yee_fdtd::FdtdDriver`], [`yee_fdtd::FdtdDriverConfig`], and
//! [`yee_fdtd::RadiationPattern`] to Python so notebook / scripting callers
//! can run a short-dipole radiation-pattern simulation end-to-end and
//! receive numpy arrays for the θ-cut.
//!
//! The Python wrapper holds only the configuration plus the grid
//! parameters (`nx`, `ny`, `nz`, `dx`); the underlying Rust
//! [`yee_fdtd::FdtdDriver`] is constructed fresh inside every
//! [`PyFdtdDriver::run`] call because the Rust `run` method consumes
//! `self`. As a result a `PyFdtdDriver` can be `.run()` multiple times,
//! each call returning an independent [`PyRadiationPattern`].

use pyo3::prelude::*;
use yee_fdtd::FdtdDriverConfig as RustCfg;

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
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub(crate) theta_deg: Vec<f64>,
    #[allow(dead_code)]
    pub(crate) e_theta_phi0: Vec<f64>,
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
    #[allow(dead_code)]
    pub(crate) nx: usize,
    #[allow(dead_code)]
    pub(crate) ny: usize,
    #[allow(dead_code)]
    pub(crate) nz: usize,
    #[allow(dead_code)]
    pub(crate) dx: f64,
    #[allow(dead_code)]
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
}
