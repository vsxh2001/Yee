//! Validation-case aggregator for the Yee workspace.
//!
//! Each known validation case (`mom-001`, `cpml-001`, `ntff-001`, ...) is
//! registered here. This commit lands the report **schema** only:
//! [`Report`] and [`CaseResult`] / [`CaseStatus`]. The runner, the
//! rendering helpers, and the placeholder cases land in subsequent
//! commits.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::Serialize;

/// Outcome of a single validation case.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum CaseStatus {
    /// Case ran to completion and met its tolerances.
    Passed,
    /// Case ran but failed its tolerance check, or its setup errored.
    Failed,
    /// Case was not executed (Phase deferral, missing toolchain, etc.).
    Skipped,
}

/// Result row for a single validation case.
#[derive(Debug, Clone, Serialize)]
pub struct CaseResult {
    /// Stable identifier (e.g. `mom-001-fast`, `cpml-001`).
    pub id: String,
    /// Human-readable one-line description (geometry, reference, tolerance).
    pub description: String,
    /// Pass/fail/skip status.
    pub status: CaseStatus,
    /// Diagnostic message: for `Passed`, the measured quantity; for
    /// `Failed`, the assertion that fired; for `Skipped`, the reason.
    pub message: String,
    /// Wall time spent inside the case body, in seconds.
    pub wall_time_seconds: f64,
}

/// Aggregated report over all registered validation cases.
#[derive(Debug, Clone, Serialize)]
pub struct Report {
    /// ISO-8601-ish UTC timestamp at the moment the runner completed.
    pub generated_at: String,
    /// Git SHA at build time, picked up from the `GIT_SHA` env var if set.
    pub git_sha: Option<String>,
    /// One [`CaseResult`] per registered case, in registration order.
    pub cases: Vec<CaseResult>,
}

/// Error type used by case bodies when their setup or solve step fails.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// I/O failure (mesh fixture missing, etc.).
    #[error("io: {0}")]
    Io(String),
    /// Solver setup or solve failure.
    #[error("solver: {0}")]
    Solver(String),
}
