//! Validation-case aggregator for the Yee workspace.
//!
//! Each known validation case (`mom-001`, `cpml-001`, `ntff-001`, ...) is
//! registered here. The [`run_all`] entry point executes them and produces
//! a structured [`Report`].
//!
//! # Phase 1.validation.0 walking-skeleton scope
//!
//! For the walking-skeleton landing, the aggregator is **mostly skips**.
//! The underlying validation suites currently live in private
//! `#[cfg(test)]` modules and integration tests of `yee-mom` / `yee-fdtd`
//! (`crates/yee-mom/tests/dipole.rs`, `crates/yee-fdtd/tests/cpml_reflection.rs`,
//! `crates/yee-fdtd/tests/ntff_dipole.rs`, `crates/yee-fdtd/tests/dispersive.rs`)
//! and are not callable from a sibling crate.
//!
//! The aggregator's value in Phase 1.validation.0 is:
//!
//! 1. Providing the report **schema** that future phases will populate.
//! 2. Demonstrating the [`run_all`] entry point.
//! 3. Documenting the "private test code can't be reached from a sibling
//!    crate" friction so Phase 1.validation.1 can either (a) move the
//!    test fixtures into public APIs of `yee-mom` / `yee-fdtd`, or
//!    (b) shell out to `cargo test --release --message-format=json` and
//!    parse the structured output.
//!
//! Cases are not faked: every placeholder reports [`CaseStatus::Skipped`]
//! with a non-empty `message` explaining the Phase 1.validation.1 deferral.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use serde::Serialize;
use std::time::Instant;

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
    /// ISO-8601-ish UTC timestamp at the moment [`run_all`] completed.
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

/// Run all known validation cases and return a structured report.
///
/// Phase 1.validation.0 walking-skeleton scope: every case currently
/// reports [`CaseStatus::Skipped`] with an explanatory `message`. The
/// full `mom-001` 24x176 gate (~7-8 min) and the `dipole_full_sweep`
/// continue to live in their respective `cargo test --release` paths
/// and are exercised by workspace CI.
pub fn run_all() -> Report {
    let cases = vec![
        run_mom_001_fast(),
        run_cpml_001(),
        run_ntff_001(),
        run_dispersive_001(),
    ];

    Report {
        generated_at: chrono_iso_now(),
        git_sha: option_env!("GIT_SHA").map(String::from),
        cases,
    }
}

fn chrono_iso_now() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Naive ISO-ish stamp; full `chrono` dep is overkill for this crate.
    format!("epoch+{secs}s")
}

fn run_mom_001_fast() -> CaseResult {
    let t0 = Instant::now();
    // The yee-mom test fixture `thin_cylinder` lives in
    // `crates/yee-mom/tests/fixtures/` and is **not** exported from the
    // public crate API. Two unblock paths for Phase 1.validation.1:
    //   (a) inline a 30-LOC thin-cylinder mesher in this crate, or
    //   (b) promote `yee_mom::fixtures::thin_cylinder` to a public
    //       (possibly feature-gated) module.
    //
    // For the walking-skeleton landing we surface the gap explicitly
    // rather than fake a pass.
    let result: Result<(f64, f64), Error> = Err(Error::Solver(
        "yee_mom::fixtures::thin_cylinder is private; inline mesher or \
         public-API exposure is a Phase 1.validation.1 task"
            .into(),
    ));

    let elapsed = t0.elapsed().as_secs_f64();
    let (status, message) = match result {
        Ok((re, im)) => (
            CaseStatus::Passed,
            format!("Z_in = {re:.3} + j{im:.3} Ohm"),
        ),
        Err(e) => (CaseStatus::Skipped, format!("{e}")),
    };
    CaseResult {
        id: "mom-001-fast".into(),
        description:
            "Half-wave dipole, NEC-4 reference (Z ~= 87+j41), 24x24 mesh (fast variant)"
                .into(),
        status,
        message,
        wall_time_seconds: elapsed,
    }
}

fn run_cpml_001() -> CaseResult {
    CaseResult {
        id: "cpml-001".into(),
        description: "CPML attenuates >= 30 dB vs PEC (FDTD)".into(),
        status: CaseStatus::Skipped,
        message:
            "Phase 1.validation.0: cpml_reflection is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.1"
                .into(),
        wall_time_seconds: 0.0,
    }
}

fn run_ntff_001() -> CaseResult {
    CaseResult {
        id: "ntff-001".into(),
        description: "NTFF broadside/endfire null >= 20 dB".into(),
        status: CaseStatus::Skipped,
        message:
            "Phase 1.validation.0: ntff_dipole is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.1"
                .into(),
        wall_time_seconds: 0.0,
    }
}

fn run_dispersive_001() -> CaseResult {
    CaseResult {
        id: "dispersive-001".into(),
        description: "Drude slab Fresnel reflection within 20%".into(),
        status: CaseStatus::Skipped,
        message:
            "Phase 1.validation.0: drude_slab is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.1"
                .into(),
        wall_time_seconds: 0.0,
    }
}
