//! `yee.touchstone` submodule — `read(path) -> dict` and `write(path, file)`.
//!
//! The dict schema is:
//!
//! ```text
//! {
//!   "z0":        float,
//!   "freq_unit": "Hz" | "kHz" | "MHz" | "GHz",
//!   "format":    "RI" | "MA" | "DB",
//!   "n_ports":   int,
//!   "freq_hz":   np.ndarray[F, f64],
//!   "data":      np.ndarray[F, N, N, c128],
//!   "comments":  list[str],
//! }
//! ```
//!
//! Note: `yee_io::touchstone::Format` variants are `RealImag`, `MagAngle`,
//! `DecibelAngle` (NOT `MagnitudeAngle` as initially noted in the spec).

use num_complex::Complex64;
use numpy::{IntoPyArray, PyReadonlyArray1, PyReadonlyArray3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::PathBuf;
use yee_io::touchstone::{self, File, Format, FreqUnit};

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    Ok(())
}

fn freq_unit_to_str(u: FreqUnit) -> &'static str {
    match u {
        FreqUnit::Hz => "Hz",
        FreqUnit::KHz => "kHz",
        FreqUnit::MHz => "MHz",
        FreqUnit::GHz => "GHz",
    }
}

fn format_to_str(f: Format) -> &'static str {
    match f {
        Format::RealImag => "RI",
        Format::MagAngle => "MA",
        Format::DecibelAngle => "DB",
    }
}

fn parse_freq_unit(s: &str) -> PyResult<FreqUnit> {
    match s.to_ascii_lowercase().as_str() {
        "hz" => Ok(FreqUnit::Hz),
        "khz" => Ok(FreqUnit::KHz),
        "mhz" => Ok(FreqUnit::MHz),
        "ghz" => Ok(FreqUnit::GHz),
        _ => Err(PyValueError::new_err(format!(
            "unknown freq_unit `{s}`; expected one of Hz, kHz, MHz, GHz"
        ))),
    }
}

fn parse_format(s: &str) -> PyResult<Format> {
    match s.to_ascii_uppercase().as_str() {
        "RI" => Ok(Format::RealImag),
        "MA" => Ok(Format::MagAngle),
        "DB" => Ok(Format::DecibelAngle),
        _ => Err(PyValueError::new_err(format!(
            "unknown format `{s}`; expected one of RI, MA, DB"
        ))),
    }
}

#[pyfunction]
fn read<'py>(py: Python<'py>, path: PathBuf) -> PyResult<Bound<'py, PyDict>> {
    let file = touchstone::read(&path).map_err(crate::errors::io_to_py)?;
    let d = PyDict::new(py);
    d.set_item("z0", file.z0)?;
    d.set_item("freq_unit", freq_unit_to_str(file.freq_unit))?;
    d.set_item("format", format_to_str(file.format))?;
    d.set_item("n_ports", file.n_ports)?;
    d.set_item("freq_hz", file.freq_hz.clone().into_pyarray(py))?;

    let f = file.data.len();
    let n = file.n_ports;
    let mut buf: Vec<Complex64> = Vec::with_capacity(f * n * n);
    for row in &file.data {
        buf.extend_from_slice(row);
    }
    let arr = numpy::ndarray::Array3::from_shape_vec((f, n, n), buf)
        .expect("buffer length matches shape by construction")
        .into_pyarray(py);
    d.set_item("data", arr)?;
    d.set_item("comments", file.comments.clone())?;
    Ok(d)
}

fn get_required<'py>(dict: &Bound<'py, PyDict>, key: &str) -> PyResult<Bound<'py, PyAny>> {
    dict.get_item(key)?
        .ok_or_else(|| PyValueError::new_err(format!("touchstone dict missing key `{key}`")))
}

#[pyfunction]
fn write(path: PathBuf, file: &Bound<'_, PyDict>) -> PyResult<()> {
    let z0: f64 = get_required(file, "z0")?.extract()?;
    let freq_unit_str: String = get_required(file, "freq_unit")?.extract()?;
    let freq_unit = parse_freq_unit(&freq_unit_str)?;
    let format_str: String = get_required(file, "format")?.extract()?;
    let format = parse_format(&format_str)?;
    let n_ports: usize = get_required(file, "n_ports")?.extract()?;

    let freq_hz_any = get_required(file, "freq_hz")?;
    let freq_hz_arr: PyReadonlyArray1<'_, f64> = freq_hz_any.extract()?;
    let freq_hz: Vec<f64> = freq_hz_arr.as_array().iter().copied().collect();

    let data_any = get_required(file, "data")?;
    let data_arr: PyReadonlyArray3<'_, Complex64> = data_any.extract()?;
    let data_view = data_arr.as_array();
    let shape = data_view.shape();
    let f = freq_hz.len();
    if shape != [f, n_ports, n_ports] {
        return Err(PyValueError::new_err(format!(
            "data must have shape [{f}, {n_ports}, {n_ports}]; got {shape:?}"
        )));
    }
    let mut data: Vec<Vec<Complex64>> = Vec::with_capacity(f);
    for k in 0..f {
        let mut mat = Vec::with_capacity(n_ports * n_ports);
        for i in 0..n_ports {
            for j in 0..n_ports {
                mat.push(data_view[[k, i, j]]);
            }
        }
        data.push(mat);
    }

    let comments_any = get_required(file, "comments")?;
    let comments_list = comments_any
        .cast::<PyList>()
        .map_err(|e| PyValueError::new_err(format!("`comments` must be a list[str]: {e}")))?;
    let mut comments: Vec<String> = Vec::with_capacity(comments_list.len());
    for item in comments_list.iter() {
        comments.push(item.extract()?);
    }

    let file_struct = File {
        n_ports,
        z0,
        freq_unit,
        format,
        freq_hz,
        data,
        comments,
    };
    touchstone::write(&path, &file_struct).map_err(crate::errors::io_to_py)
}
