//! Python bindings for the `yee-compute` GPU/CPU engine (ADR-0175..0177).
//!
//! Exposed as the `yee.compute` submodule. The shape is a config-accumulating
//! [`PyFdtdSim`] builder plus an immutable [`PyFdtdResult`]: `run()`
//! constructs a fresh engine each call (the `PyFdtdDriver` idiom), so one
//! sim object can be `.run()` repeatedly and on different backends
//! (`"cpu"`, `"gpu"`, `"auto"`).
//!
//! Material maps arrive as numpy `float64` arrays of shape
//! `(nx+1, ny+1, nz+1)` and PEC masks as `bool` arrays with the component's
//! staggered shape — the same conventions as `yee_fdtd::YeeGrid`. Probe
//! series and final fields come back as numpy arrays.

use numpy::{IntoPyArray, PyArray1, PyArray3, PyReadonlyArray3};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use yee_compute::{
    Boundary, ComputeError, CpmlConfig, Drive, EComponent, FdtdSpec, Fields, Materials, Probe,
    ResistivePort, SoftSource, Waveform,
};

fn parse_component(name: &str) -> PyResult<EComponent> {
    match name {
        "ex" => Ok(EComponent::Ex),
        "ey" => Ok(EComponent::Ey),
        "ez" => Ok(EComponent::Ez),
        other => Err(PyValueError::new_err(format!(
            "component must be 'ex' | 'ey' | 'ez', got {other:?}"
        ))),
    }
}

fn compute_err(e: ComputeError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// Configuration builder + runner for a `yee-compute` FDTD simulation.
#[pyclass(name = "FdtdSim", module = "yee.compute")]
pub struct PyFdtdSim {
    spec: FdtdSpec,
    materials: Materials,
    boundary: Boundary,
    drive: Drive,
}

#[pymethods]
impl PyFdtdSim {
    /// Build a uniform-vacuum sim of `nx × ny × nz` cubic cells of size
    /// `dx` metres, `dt` at 0.9× the Courant limit, no boundary phase.
    #[new]
    fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self {
            spec: FdtdSpec::vacuum(nx, ny, nz, dx),
            materials: Materials::default(),
            boundary: Boundary::None,
            drive: Drive::default(),
        }
    }

    /// Time step in seconds.
    #[getter]
    fn dt(&self) -> f64 {
        self.spec.dt
    }

    /// Grid dims `(nx, ny, nz)`.
    #[getter]
    fn dims(&self) -> (usize, usize, usize) {
        (self.spec.nx, self.spec.ny, self.spec.nz)
    }

    /// Attach a per-cell relative-permittivity map, shape `(nx+1, ny+1, nz+1)`.
    fn set_eps_r_cells(&mut self, eps_r: PyReadonlyArray3<'_, f64>) -> PyResult<()> {
        self.materials.eps_r_cells = Some(check_cell_map(&self.spec, eps_r, "eps_r")?);
        Ok(())
    }

    /// Attach a per-cell relative-permeability map, shape `(nx+1, ny+1, nz+1)`.
    fn set_mu_r_cells(&mut self, mu_r: PyReadonlyArray3<'_, f64>) -> PyResult<()> {
        self.materials.mu_r_cells = Some(check_cell_map(&self.spec, mu_r, "mu_r")?);
        Ok(())
    }

    /// Attach a per-cell conductivity map (S/m), shape `(nx+1, ny+1, nz+1)`.
    fn set_sigma_cells(&mut self, sigma: PyReadonlyArray3<'_, f64>) -> PyResult<()> {
        self.materials.sigma_cells = Some(check_cell_map(&self.spec, sigma, "sigma")?);
        Ok(())
    }

    /// Attach an interior PEC mask for one E component (`bool` array with
    /// that component's staggered shape).
    fn set_pec_mask(&mut self, component: &str, mask: PyReadonlyArray3<'_, bool>) -> PyResult<()> {
        let comp = parse_component(component)?;
        let dims = match comp {
            EComponent::Ex => self.spec.ex_dims(),
            EComponent::Ey => self.spec.ey_dims(),
            EComponent::Ez => self.spec.ez_dims(),
        };
        let shape = mask.as_array().shape().to_vec();
        if shape != vec![dims.0, dims.1, dims.2] {
            return Err(PyValueError::new_err(format!(
                "pec mask for {component:?} must have shape {dims:?}, got {shape:?}"
            )));
        }
        let flat: Vec<bool> = mask.as_array().iter().copied().collect();
        match comp {
            EComponent::Ex => self.materials.pec_mask_ex = Some(flat),
            EComponent::Ey => self.materials.pec_mask_ey = Some(flat),
            EComponent::Ez => self.materials.pec_mask_ez = Some(flat),
        }
        Ok(())
    }

    /// Select the outer boundary: `"none"`, `"pec"`, or `"cpml"` (with
    /// `npml` layers and a per-axis enable mask).
    #[pyo3(signature = (kind, npml = 10, axes = (true, true, true)))]
    fn set_boundary(&mut self, kind: &str, npml: usize, axes: (bool, bool, bool)) -> PyResult<()> {
        self.boundary = match kind {
            "none" => Boundary::None,
            "pec" => Boundary::PecBox,
            "cpml" => Boundary::Cpml(
                CpmlConfig::for_spec(&self.spec, npml).with_axes([axes.0, axes.1, axes.2]),
            ),
            other => {
                return Err(PyValueError::new_err(format!(
                    "boundary must be 'none' | 'pec' | 'cpml', got {other:?}"
                )));
            }
        };
        Ok(())
    }

    /// Add an additive Gaussian point source `exp(−((t−t0)/sigma)²)` on
    /// `component` at cell `(i, j, k)`.
    fn add_gaussian_source(
        &mut self,
        component: &str,
        cell: (usize, usize, usize),
        t0: f64,
        sigma: f64,
    ) -> PyResult<()> {
        self.drive.soft_sources.push(SoftSource {
            component: parse_component(component)?,
            cell,
            waveform: Waveform::Gaussian { t0, sigma },
        });
        Ok(())
    }

    /// Add a lumped resistive drive port on the `E_z` cell `(i, j, k)` with
    /// a modulated-Gaussian EMF (`v0`, carrier `f0` Hz, FWHM bandwidth `bw`
    /// Hz, centre `t0_steps` steps).
    fn add_resistive_port(
        &mut self,
        cell: (usize, usize, usize),
        resistance: f64,
        v0: f64,
        f0: f64,
        bw: f64,
        t0_steps: usize,
    ) {
        self.drive.ports.push(ResistivePort {
            cell,
            resistance,
            waveform: Waveform::GaussianPulse {
                v0,
                f0,
                bw,
                t0_steps,
            },
        });
    }

    /// Record `component` at cell `(i, j, k)` once per step.
    fn add_probe(&mut self, component: &str, cell: (usize, usize, usize)) -> PyResult<()> {
        self.drive.probes.push(Probe {
            component: parse_component(component)?,
            cell,
        });
        Ok(())
    }

    /// Run `n_steps` on `backend` (`"cpu"`, `"gpu"`, or `"auto"` — GPU with
    /// silent CPU fallback). Returns an [`PyFdtdResult`]. Each call builds a
    /// fresh engine from the accumulated configuration.
    #[pyo3(signature = (n_steps, backend = "auto"))]
    fn run(&self, py: Python<'_>, n_steps: usize, backend: &str) -> PyResult<PyFdtdResult> {
        let spec = self.spec;
        let fields = Fields::zero(&spec);
        let materials = self.materials.clone();
        let boundary = self.boundary.clone();
        let drive = self.drive.clone();
        // The solve holds no Python objects — release the GIL for its
        // duration (multi-second runs are the norm).
        let (probes, fields, backend_used) = py
            .detach(move || -> Result<_, ComputeError> {
                match backend {
                    "cpu" => {
                        let mut engine = yee_compute::CpuFdtd::with_drive(
                            spec, fields, materials, boundary, drive,
                        );
                        engine.step_n(n_steps);
                        Ok((
                            engine.probe_series().to_vec(),
                            engine.fields().clone(),
                            "cpu",
                        ))
                    }
                    "gpu" => {
                        let mut engine = yee_compute::GpuFdtd::with_drive(
                            spec, fields, materials, boundary, drive, n_steps,
                        )?;
                        engine.step_n(n_steps)?;
                        Ok((engine.read_probes()?, engine.read_fields()?, "gpu"))
                    }
                    "auto" => {
                        match yee_compute::GpuFdtd::with_drive(
                            spec,
                            fields.clone(),
                            materials.clone(),
                            boundary.clone(),
                            drive.clone(),
                            n_steps,
                        ) {
                            Ok(mut engine) => {
                                engine.step_n(n_steps)?;
                                Ok((engine.read_probes()?, engine.read_fields()?, "gpu"))
                            }
                            Err(ComputeError::NoAdapter) => {
                                let mut engine = yee_compute::CpuFdtd::with_drive(
                                    spec, fields, materials, boundary, drive,
                                );
                                engine.step_n(n_steps);
                                Ok((
                                    engine.probe_series().to_vec(),
                                    engine.fields().clone(),
                                    "cpu",
                                ))
                            }
                            Err(other) => Err(other),
                        }
                    }
                    other => Err(ComputeError::Device(format!(
                        "backend must be 'cpu' | 'gpu' | 'auto', got {other:?}"
                    ))),
                }
            })
            .map_err(compute_err)?;
        Ok(PyFdtdResult {
            spec,
            probes,
            fields,
            backend: backend_used.to_string(),
        })
    }
}

fn check_cell_map(
    spec: &FdtdSpec,
    map: PyReadonlyArray3<'_, f64>,
    name: &str,
) -> PyResult<Vec<f64>> {
    let want = [spec.nx + 1, spec.ny + 1, spec.nz + 1];
    let shape = map.as_array().shape().to_vec();
    if shape != want.to_vec() {
        return Err(PyValueError::new_err(format!(
            "{name} map must have shape {want:?}, got {shape:?}"
        )));
    }
    Ok(map.as_array().iter().copied().collect())
}

/// Result of a [`PyFdtdSim::run`]: probe series + the six final field
/// components as numpy arrays, and which backend actually ran.
#[pyclass(name = "FdtdResult", module = "yee.compute")]
pub struct PyFdtdResult {
    spec: FdtdSpec,
    probes: Vec<Vec<f64>>,
    fields: Fields,
    backend: String,
}

#[pymethods]
impl PyFdtdResult {
    /// Backend that executed the run: `"cpu"` or `"gpu"`.
    #[getter]
    fn backend(&self) -> &str {
        &self.backend
    }

    /// Number of recorded probes.
    #[getter]
    fn n_probes(&self) -> usize {
        self.probes.len()
    }

    /// Probe series `index` as a 1-D `float64` array (one sample per step).
    fn probe<'py>(&self, py: Python<'py>, index: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        self.probes
            .get(index)
            .map(|s| s.clone().into_pyarray(py))
            .ok_or_else(|| PyValueError::new_err(format!("probe index {index} out of range")))
    }

    /// A final field component (`"ex"` … `"hz"`) as a 3-D `float64` array
    /// in its staggered shape.
    fn field<'py>(&self, py: Python<'py>, component: &str) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let s = &self.spec;
        let (dims, data) = match component {
            "ex" => (s.ex_dims(), &self.fields.ex),
            "ey" => (s.ey_dims(), &self.fields.ey),
            "ez" => (s.ez_dims(), &self.fields.ez),
            "hx" => (s.hx_dims(), &self.fields.hx),
            "hy" => (s.hy_dims(), &self.fields.hy),
            "hz" => (s.hz_dims(), &self.fields.hz),
            other => {
                return Err(PyValueError::new_err(format!(
                    "component must be one of ex/ey/ez/hx/hy/hz, got {other:?}"
                )));
            }
        };
        let arr = ndarray::Array3::from_shape_vec(dims, data.clone())
            .expect("field length matches its staggered shape by construction");
        Ok(arr.into_pyarray(py))
    }
}

/// Register the `yee.compute` submodule.
pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFdtdSim>()?;
    m.add_class::<PyFdtdResult>()?;
    Ok(())
}
