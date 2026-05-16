//! Validation-aggregator Python bindings.
//!
//! Exposes [`yee_validation::Report::run_all`] to Python as
//! [`run_validation`], together with read-only wrappers around the
//! underlying [`yee_validation::Report`] / [`yee_validation::CaseResult`]
//! shapes ([`PyValidationReport`] / [`PyValidationCase`]).
//!
//! The wrappers are intentionally narrow: they expose the fields a
//! notebook user typically wants (case id, pass/fail status, notes,
//! wall time, plot paths) plus the two top-level helpers
//! ([`PyValidationReport::has_failures`], [`PyValidationReport::to_json`]).
//! Everything else is reachable via the JSON dump.

use pyo3::prelude::*;
use yee_validation::{CaseResult, Report, Status};

/// Python view of a single [`yee_validation::CaseResult`] row.
///
/// All getters are cheap clones of the underlying Rust struct. The
/// `status` getter returns one of the strings `"Passed"`, `"Failed"`,
/// `"Skipped"` to match the Rust [`yee_validation::Status`] variants
/// without leaking a numeric discriminant.
#[pyclass(module = "yee", name = "ValidationCase")]
pub struct PyValidationCase {
    inner: CaseResult,
}

#[pymethods]
impl PyValidationCase {
    /// Stable case identifier, e.g. `"mom-001"`, `"cpml-001"`.
    #[getter]
    fn id(&self) -> String {
        self.inner.id.clone()
    }

    /// Pass / fail / skip status as `"Passed"`, `"Failed"`, `"Skipped"`.
    #[getter]
    fn status(&self) -> String {
        match self.inner.status {
            Status::Passed => "Passed".to_string(),
            Status::Failed => "Failed".to_string(),
            Status::Skipped => "Skipped".to_string(),
        }
    }

    /// Diagnostic note (measured value on pass, failing assertion on
    /// fail, deferral reason on skip).
    #[getter]
    fn notes(&self) -> String {
        self.inner.notes.clone()
    }

    /// Wall time spent inside this case body, in seconds.
    #[getter]
    fn wall_time_seconds(&self) -> f64 {
        self.inner.wall_time_seconds
    }

    /// Filesystem paths to any plot artifacts emitted by this case.
    #[getter]
    fn plot_paths(&self) -> Vec<String> {
        self.inner
            .plot_paths
            .iter()
            .map(|p| p.display().to_string())
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "ValidationCase(id={}, status={}, wall_time_seconds={:.3})",
            self.inner.id,
            match self.inner.status {
                Status::Passed => "Passed",
                Status::Failed => "Failed",
                Status::Skipped => "Skipped",
            },
            self.inner.wall_time_seconds,
        )
    }
}

/// Python view of a [`yee_validation::Report`].
///
/// `cases` yields a list of [`PyValidationCase`] wrappers in
/// registration order. `has_failures()` matches the Rust helper of the
/// same name and treats `Skipped` cases as non-failing. `to_json()`
/// returns the same pretty-printed JSON the CLI / CI artifact emits.
#[pyclass(module = "yee", name = "ValidationReport")]
pub struct PyValidationReport {
    inner: Report,
}

#[pymethods]
impl PyValidationReport {
    /// All [`PyValidationCase`] rows in registration order.
    #[getter]
    fn cases(&self) -> Vec<PyValidationCase> {
        self.inner
            .cases
            .iter()
            .cloned()
            .map(|c| PyValidationCase { inner: c })
            .collect()
    }

    /// `True` iff at least one case has `status == "Failed"`. Skipped
    /// cases do not count as failures.
    fn has_failures(&self) -> bool {
        self.inner.has_failures()
    }

    /// Pretty-printed JSON dump of the underlying [`yee_validation::Report`].
    fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!("ValidationReport(n_cases={})", self.inner.cases.len())
    }
}

/// Run all known validation cases and return a [`PyValidationReport`].
///
/// **Warning: this is slow.** It invokes the real `mom-001`
/// half-wave-dipole gate (~7-8 min wall time in release builds) plus
/// the `mom-002` microstrip free-space placeholder (~seconds) among
/// others. Intended for offline / notebook use, not for inner-loop
/// iteration.
///
/// The GIL is held for the duration of the call — the wrapped Rust
/// aggregator is single-threaded and does not need to call back into
/// Python, so releasing the GIL would only complicate the wrapper
/// without buying parallelism.
#[pyfunction]
pub fn run_validation(_py: Python<'_>) -> PyValidationReport {
    PyValidationReport {
        inner: Report::run_all(),
    }
}
