//! Validation-case aggregator for the Yee workspace.
//!
//! Each known validation case (`mom-001`, `cpml-001`, `ntff-001`, ...) is
//! registered here. The [`Report::run_all`] entry point executes them and
//! produces a structured [`Report`] that can be serialized to JSON (via
//! [`Report::to_json`]) or rendered as Markdown for CI artifacts.
//!
//! # Phase 1.validation.1 scope
//!
//! The aggregator now actually invokes the `mom-001` solver path
//! (inlined cylinder fixture, [`yee_mom::PlanarMoM`] sweep at the
//! λ/2 = 1 m resonance frequency, NEC-4 reference `Z ≈ 87 + j41 Ω` with
//! 5%/10% tolerance on Re/Im). `mom-002` / `mom-003` remain
//! [`CaseStatus::Skipped`] until the Phase 1.1.1 Sommerfeld /
//! multi-image DCIM extraction lands — see CLAUDE.md §10.
//!
//! The FDTD cases (`cpml-001`, `ntff-001`, `dispersive-001`) continue
//! to report [`CaseStatus::Skipped`] until their test fixtures are
//! promoted out of `#[cfg(test)]` modules; that work is Phase
//! 1.validation.2 territory and explicitly out of lane for this
//! landing.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use num_complex::Complex64;
use serde::Serialize;
use std::path::PathBuf;
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

/// Convenience alias matching the spec-level name `Status`. Aliased
/// to [`CaseStatus`] so downstream consumers (notebooks, CLI bins,
/// integration tests) can write `yee_validation::Status::Passed`
/// without coupling to the historical "Case" prefix.
pub type Status = CaseStatus;

/// Result row for a single validation case.
#[derive(Debug, Clone, Serialize)]
pub struct CaseResult {
    /// Stable identifier (e.g. `mom-001`, `cpml-001`).
    pub id: String,
    /// Human-readable one-line description (geometry, reference, tolerance).
    pub description: String,
    /// Pass/fail/skip status.
    pub status: CaseStatus,
    /// Diagnostic note: for `Passed`, the measured quantity; for
    /// `Failed`, the assertion that fired; for `Skipped`, the reason.
    pub notes: String,
    /// Wall time spent inside the case body, in seconds.
    pub wall_time_seconds: f64,
    /// Filesystem paths to any plot artifacts (PNGs / SVGs) the case
    /// emitted. Empty for `Skipped` cases and for `Failed` cases whose
    /// solver errored before reaching the plot step. Paths are
    /// CWD-relative when written under `validation/results/`, matching
    /// the CI artifact-upload glob.
    #[serde(default)]
    pub plot_paths: Vec<PathBuf>,
}

/// Aggregated report over all registered validation cases.
#[derive(Debug, Clone, Serialize)]
pub struct Report {
    /// ISO-8601-ish UTC timestamp at the moment [`Report::run_all`] completed.
    pub generated_at: String,
    /// Git SHA at build time, picked up from the `GIT_SHA` env var if set.
    pub git_sha: Option<String>,
    /// One [`CaseResult`] per registered case, in registration order.
    pub cases: Vec<CaseResult>,
}

impl Report {
    /// Run all known validation cases and return a structured report.
    ///
    /// `mom-001` runs the real 24x176 thin-cylinder dipole solve
    /// (~7-8 min wall time in `--release`). The remaining cases are
    /// Phase deferrals; see the crate-level documentation for the
    /// full status of each.
    pub fn run_all() -> Report {
        let cases = vec![
            run_mom_001(),
            run_mom_002(),
            run_mom_003(),
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

    /// Return `true` iff any case in the report has status
    /// [`CaseStatus::Failed`]. `Skipped` does **not** count as a
    /// failure — skipped cases carry an explanatory `notes` string
    /// and are expected to remain in that state until the relevant
    /// upstream work lands.
    pub fn has_failures(&self) -> bool {
        self.cases.iter().any(|c| c.status == CaseStatus::Failed)
    }

    /// Render the report as a GitHub-flavoured Markdown table.
    ///
    /// Intended for the CI artifact step: drop the output into a job
    /// summary or a Markdown file for reviewers.
    pub fn to_markdown(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "# Yee Validation Report\n\nGenerated: {}\n\n",
            self.generated_at
        ));
        if let Some(sha) = &self.git_sha {
            out.push_str(&format!("Git SHA: `{sha}`\n\n"));
        }
        out.push_str("| ID | Status | Wall time | Notes |\n");
        out.push_str("|----|--------|-----------|-------|\n");
        for c in &self.cases {
            let icon = match c.status {
                CaseStatus::Passed => "PASS",
                CaseStatus::Failed => "FAIL",
                CaseStatus::Skipped => "SKIP",
            };
            out.push_str(&format!(
                "| `{}` | {} {:?} | {:.2}s | {} |\n",
                c.id, icon, c.status, c.wall_time_seconds, c.notes
            ));
        }
        out
    }

    /// Pretty-print the report as JSON. Intended for machine consumers
    /// (CI annotations, dashboards, regression bots).
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// Module-level free function alias retained for backwards
/// compatibility with the Phase 1.validation.0 skeleton.
///
/// Prefer [`Report::run_all`] in new code.
pub fn run_all() -> Report {
    Report::run_all()
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

/// Error type used by case bodies when their setup or solve step fails.
///
/// Each `run_<case>` body returns a `Result<_, Error>` internally; the
/// outer [`CaseResult`] folds that into a [`CaseStatus`] + notes.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// I/O failure (mesh fixture missing, etc.).
    #[error("io: {0}")]
    Io(String),
    /// Solver setup or solve failure.
    #[error("solver: {0}")]
    Solver(String),
}

// ---------------------------------------------------------------------
// Case bodies
// ---------------------------------------------------------------------

/// `mom-001` reference frequency `f0 = c/(2L)` for `L = 1 m`. This is
/// the same `λ/2 = 1 m` resonance point used in
/// `crates/yee-mom/tests/dipole.rs`; spelled out here to avoid a
/// runtime dependency on `yee_core::units::C0` for the sole purpose of
/// the divide.
const MOM_001_F0_HZ: f64 = 149_896_229.0;

const MOM_001_Z0_REF: f64 = 50.0;
// NEC-4 finite-radius reference Z at a = 5 mm, L = 1 m, half-wave,
// delta-gap on the lateral cylinder surface.
const MOM_001_REFERENCE_RE: f64 = 87.0;
const MOM_001_REFERENCE_IM: f64 = 41.0;
const MOM_001_TOL_RE_REL: f64 = 0.05;
const MOM_001_TOL_IM_REL: f64 = 0.10;
const MOM_001_N_AXIAL: usize = 24;
const MOM_001_N_AROUND: usize = 176;

// ---- Coarse-mesh plot-only sweep ----
//
// The pass/fail Z check above uses the 24x176 fine mesh at a single
// frequency (~8 min wall-time). Repeating that 21 times for a frequency
// sweep would blow past any sane CI budget (PlanarMoM::run rebuilds the
// impedance matrix from scratch per frequency — see
// crates/yee-mom/src/solve.rs::s_parameters_sweep). Instead, the plots
// are generated from a coarse 24x44 mesh sweep across [100, 200] MHz,
// which keeps the resonance dip visible near f0 = 149.896 MHz while
// running in well under one fine-mesh tick. The notes string records
// that this is a coarse-mesh preview, not the same numerics as the
// pass/fail check.
const MOM_001_PLOT_N_AROUND: usize = 44;
const MOM_001_PLOT_F_MIN_HZ: f64 = 100.0e6;
const MOM_001_PLOT_F_MAX_HZ: f64 = 200.0e6;
const MOM_001_PLOT_N_POINTS: usize = 21;

/// Z_in = Z0 · (1 + S11) / (1 − S11). Identical convention to the
/// `z_in_from_s11` helper in `crates/yee-mom/tests/dipole.rs`.
fn z_in_from_s11(s11: Complex64, z0: f64) -> Complex64 {
    Complex64::new(z0, 0.0) * (Complex64::new(1.0, 0.0) + s11) / (Complex64::new(1.0, 0.0) - s11)
}

/// Inlined thin-cylinder mesher. Mirrors
/// `crates/yee-mom/tests/fixtures/cylinder.rs::thin_cylinder` so the
/// aggregator does not depend on test fixtures of `yee-mom` (which
/// are not part of the crate's public API). Any divergence between
/// the two should be treated as a bug in this copy.
fn mom_001_cylinder_mesh(
    length_m: f64,
    radius_m: f64,
    n_axial: usize,
    n_around: usize,
) -> yee_mesh::TriMesh {
    use nalgebra::Vector3;

    assert!(
        n_axial >= 2 && n_axial.is_multiple_of(2),
        "n_axial must be even and >= 2"
    );
    assert!(n_around >= 3, "n_around must be >= 3");

    let mut vertices: Vec<Vector3<f64>> = Vec::with_capacity((n_axial + 1) * n_around);
    let dz = length_m / (n_axial as f64);
    let z0 = -length_m / 2.0;
    let dtheta = std::f64::consts::TAU / (n_around as f64);

    for i in 0..=n_axial {
        let z = z0 + (i as f64) * dz;
        for j in 0..n_around {
            let theta = (j as f64) * dtheta;
            vertices.push(Vector3::new(
                radius_m * theta.cos(),
                radius_m * theta.sin(),
                z,
            ));
        }
    }

    let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(2 * n_axial * n_around);
    let mut tags: Vec<u32> = Vec::with_capacity(2 * n_axial * n_around);
    let central_ring = n_axial / 2;

    for i in 0..n_axial {
        for j in 0..n_around {
            let j_next = (j + 1) % n_around;
            let a = (i * n_around + j) as u32;
            let b = (i * n_around + j_next) as u32;
            let c = ((i + 1) * n_around + j_next) as u32;
            let d = ((i + 1) * n_around + j) as u32;
            triangles.push([a, b, c]);
            triangles.push([a, c, d]);
            let tag = if i == central_ring - 1 {
                1
            } else if i == central_ring {
                2
            } else {
                0
            };
            tags.push(tag);
            tags.push(tag);
        }
    }

    yee_mesh::TriMesh::new(vertices, triangles, tags).expect("cylinder mesh invariants")
}

/// mom-001: half-wave dipole impedance gate.
///
/// Builds the 24x176 lateral-surface cylinder mesh, runs the planar
/// MoM sweep at `f0 = c/(2L) ≈ 149.896 MHz`, and compares the
/// extracted `Z_in = Z0 (1 + S11) / (1 − S11)` to the NEC-4
/// finite-radius reference `Z ≈ 87 + j41 Ω`. Pass requires
/// `|Re-87|/87 < 5%` and `|Im-41|/41 < 10%`.
fn run_mom_001() -> CaseResult {
    use yee_core::{FreqRange, Solver};
    use yee_mom::PlanarMoM;

    let t0 = Instant::now();
    let result: Result<Complex64, Error> = (|| -> Result<Complex64, Error> {
        let mesh = mom_001_cylinder_mesh(1.0, 0.005, MOM_001_N_AXIAL, MOM_001_N_AROUND);
        // Single-point sweep at the λ/2 resonance. FreqRange requires
        // `stop > start` even for a one-point evaluation; this matches
        // the convention used in crates/yee-mom/tests/dipole.rs.
        let freq = FreqRange::new(MOM_001_F0_HZ, MOM_001_F0_HZ + 1.0, 1)
            .map_err(|e| Error::Solver(format!("FreqRange::new: {e}")))?;
        let solver = PlanarMoM::default();
        let s = solver
            .run(&mesh, freq)
            .map_err(|e| Error::Solver(format!("PlanarMoM::run: {e}")))?;
        let s11 = s.data[0][0];
        Ok(z_in_from_s11(s11, MOM_001_Z0_REF))
    })();

    let elapsed = t0.elapsed().as_secs_f64();
    let (status, notes) = match result {
        Ok(z_in) => {
            let err_re = (z_in.re - MOM_001_REFERENCE_RE).abs() / MOM_001_REFERENCE_RE;
            let err_im = (z_in.im - MOM_001_REFERENCE_IM).abs() / MOM_001_REFERENCE_IM;
            let passed = err_re <= MOM_001_TOL_RE_REL && err_im <= MOM_001_TOL_IM_REL;
            let status = if passed {
                CaseStatus::Passed
            } else {
                CaseStatus::Failed
            };
            let notes = format!(
                "Z_in = {:.3} + j{:.3} Ohm (NEC-4 ref 87 + j41); \
                 |dRe|/Re = {:.4} (tol {:.2}), |dIm|/Im = {:.4} (tol {:.2})",
                z_in.re, z_in.im, err_re, MOM_001_TOL_RE_REL, err_im, MOM_001_TOL_IM_REL
            );
            (status, notes)
        }
        Err(e) => (CaseStatus::Failed, format!("{e}")),
    };
    // Plot generation runs regardless of pass/fail so reviewers can
    // eyeball the failure. Solver errors during the coarse sweep, or
    // plotter I/O errors, are surfaced via `plot_notes` and appended
    // to the case notes but do NOT downgrade a Passed case to Failed —
    // the pass/fail decision is the fine-mesh Z above.
    let (plot_paths, plot_notes) = match generate_mom_001_plots() {
        Ok(paths) => (
            paths,
            format!(
                " | plots: coarse {n_ax}x{n_ar} mesh, {n} freqs in [{f0:.1}, {f1:.1}] MHz",
                n_ax = MOM_001_N_AXIAL,
                n_ar = MOM_001_PLOT_N_AROUND,
                n = MOM_001_PLOT_N_POINTS,
                f0 = MOM_001_PLOT_F_MIN_HZ * 1e-6,
                f1 = MOM_001_PLOT_F_MAX_HZ * 1e-6,
            ),
        ),
        Err(e) => (Vec::new(), format!(" | plot generation failed: {e}")),
    };

    CaseResult {
        id: "mom-001".into(),
        description: "Half-wave dipole, NEC-4 reference (Z ~= 87+j41 Ohm), 24x176 cylinder mesh"
            .into(),
        status,
        notes: format!("{notes}{plot_notes}"),
        wall_time_seconds: elapsed,
        plot_paths,
    }
}

/// Resolve the output directory for plot artifacts.
///
/// Resolution order:
/// 1. `$YEE_VALIDATION_OUT_DIR` if set (CI override).
/// 2. `<workspace_root>/validation/results/` — workspace root is
///    discovered by walking up from `CARGO_MANIFEST_DIR` until a
///    directory containing `Cargo.lock` is found. This is the
///    standard convention for locating the workspace root and lets
///    `cargo test` (CWD = crate dir) and `cargo run --bin
///    yee-validate` (CWD = wherever the user invoked it) both write
///    to the same place.
/// 3. `./validation/results/` relative to the current working
///    directory if neither of the above resolves.
fn validation_results_dir() -> PathBuf {
    if let Ok(p) = std::env::var("YEE_VALIDATION_OUT_DIR") {
        return PathBuf::from(p);
    }
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut cur: Option<&std::path::Path> = Some(manifest_dir.as_path());
    while let Some(d) = cur {
        if d.join("Cargo.lock").is_file() {
            return d.join("validation").join("results");
        }
        cur = d.parent();
    }
    PathBuf::from("validation").join("results")
}

/// Generate the S11 dB + Smith chart PNGs for mom-001 under
/// `validation/results/` (CWD-relative). Uses a coarse mesh
/// (`MOM_001_N_AXIAL` x `MOM_001_PLOT_N_AROUND`) and a
/// `MOM_001_PLOT_N_POINTS`-point sweep across
/// `[MOM_001_PLOT_F_MIN_HZ, MOM_001_PLOT_F_MAX_HZ]`. The sweep
/// resolution is intentionally lower than the fine-mesh pass/fail
/// check above — see the constant block for the budget rationale.
///
/// Returns the list of paths written on success, or an [`Error`] if
/// the solver or the plotter failed. The caller folds either into the
/// `CaseResult` notes; plot failures do not flip a Passed status to
/// Failed.
fn generate_mom_001_plots() -> Result<Vec<PathBuf>, Error> {
    use yee_core::{FreqRange, Solver};
    use yee_mom::PlanarMoM;
    use yee_plotters::{PlotConfig, PlotFormat, plot_s11_db, plot_smith_chart};

    let mesh = mom_001_cylinder_mesh(1.0, 0.005, MOM_001_N_AXIAL, MOM_001_PLOT_N_AROUND);
    let freq = FreqRange::new(
        MOM_001_PLOT_F_MIN_HZ,
        MOM_001_PLOT_F_MAX_HZ,
        MOM_001_PLOT_N_POINTS,
    )
    .map_err(|e| Error::Solver(format!("FreqRange::new (plot sweep): {e}")))?;
    let solver = PlanarMoM::default();
    let s = solver
        .run(&mesh, freq)
        .map_err(|e| Error::Solver(format!("PlanarMoM::run (plot sweep): {e}")))?;

    let freq_hz = s.freq_hz.clone();
    let s11: Vec<Complex64> = s.data.iter().map(|row| row[0]).collect();

    let dir = validation_results_dir();
    std::fs::create_dir_all(&dir).map_err(|e| Error::Io(format!("create_dir_all: {e}")))?;

    let s11_db_path = dir.join("mom-001-s11-db.png");
    let smith_path = dir.join("mom-001-smith.png");

    plot_s11_db(
        &freq_hz,
        &s11,
        &s11_db_path,
        &PlotConfig {
            width_px: 800,
            height_px: 600,
            title: "mom-001 |S11| dB (coarse-mesh preview)".to_string(),
            format: PlotFormat::Png,
        },
    )
    .map_err(|e| Error::Io(format!("plot_s11_db: {e}")))?;

    plot_smith_chart(
        &s11,
        &smith_path,
        &PlotConfig {
            width_px: 600,
            height_px: 600,
            title: "mom-001 S11 Smith chart (coarse-mesh preview)".to_string(),
            format: PlotFormat::Png,
        },
    )
    .map_err(|e| Error::Io(format!("plot_smith_chart: {e}")))?;

    Ok(vec![s11_db_path, smith_path])
}

/// mom-002: microstrip Z0 — rides on the `MultilayerGreens`
/// one-image DCIM placeholder per CLAUDE.md §10. Reported as
/// [`CaseStatus::Skipped`] until Phase 1.1.1 ships real Sommerfeld /
/// multi-image extraction.
fn run_mom_002() -> CaseResult {
    CaseResult {
        id: "mom-002".into(),
        description: "Microstrip characteristic impedance Z0 (loose tolerance until Phase 1.1.1)"
            .into(),
        status: CaseStatus::Skipped,
        notes: "Phase 1.1.0 MultilayerGreens placeholder: one-image DCIM only. \
             Awaiting Phase 1.1.1 Sommerfeld-integral / multi-image DCIM extraction \
             before a meaningful tolerance can be asserted."
            .into(),
        wall_time_seconds: 0.0,
        plot_paths: Vec::new(),
    }
}

/// mom-003: 2.4 GHz patch resonance — same `MultilayerGreens`
/// placeholder dependency as mom-002, same deferral.
fn run_mom_003() -> CaseResult {
    CaseResult {
        id: "mom-003".into(),
        description: "2.4 GHz patch antenna resonance (loose tolerance until Phase 1.1.1)".into(),
        status: CaseStatus::Skipped,
        notes: "Phase 1.1.0 MultilayerGreens placeholder: one-image DCIM only. \
             Awaiting Phase 1.1.1 Sommerfeld-integral / multi-image DCIM extraction \
             before a meaningful tolerance can be asserted."
            .into(),
        wall_time_seconds: 0.0,
        plot_paths: Vec::new(),
    }
}

fn run_cpml_001() -> CaseResult {
    CaseResult {
        id: "cpml-001".into(),
        description: "CPML attenuates >= 30 dB vs PEC (FDTD)".into(),
        status: CaseStatus::Skipped,
        notes: "Phase 1.validation.1: cpml_reflection is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.2"
            .into(),
        wall_time_seconds: 0.0,
        plot_paths: Vec::new(),
    }
}

fn run_ntff_001() -> CaseResult {
    CaseResult {
        id: "ntff-001".into(),
        description: "NTFF broadside/endfire null >= 20 dB".into(),
        status: CaseStatus::Skipped,
        notes: "Phase 1.validation.1: ntff_dipole is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.2"
            .into(),
        wall_time_seconds: 0.0,
        plot_paths: Vec::new(),
    }
}

fn run_dispersive_001() -> CaseResult {
    CaseResult {
        id: "dispersive-001".into(),
        description: "Drude slab Fresnel reflection within 20%".into(),
        status: CaseStatus::Skipped,
        notes: "Phase 1.validation.1: drude_slab is a yee-fdtd integration test; \
             aggregator integration deferred to Phase 1.validation.2"
            .into(),
        wall_time_seconds: 0.0,
        plot_paths: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Standalone wall-time probe for the coarse plot sweep, used
    /// once during Phase 1.validation.2 development to confirm the
    /// `MOM_001_PLOT_N_AROUND` / `MOM_001_PLOT_N_POINTS` budget is
    /// reasonable. Ignored because it still does ~21 LU solves on a
    /// 24x44 cylinder mesh.
    #[test]
    #[ignore = "slow: exercises the coarse plot sweep (~minutes)"]
    fn mom_001_plot_sweep_runs() {
        let paths = generate_mom_001_plots().expect("plot sweep");
        assert_eq!(paths.len(), 2, "expected S11 dB + Smith chart");
        for p in &paths {
            assert!(p.exists(), "plot missing: {}", p.display());
        }
    }

    /// Cheap unit test: does NOT call [`Report::run_all`] because
    /// `mom-001` takes 7-8 minutes in `--release`. The full pipeline
    /// is exercised by `tests/integration.rs` under `--include-ignored`.
    #[test]
    fn report_skip_only_subset_renders() {
        let report = Report {
            generated_at: chrono_iso_now(),
            git_sha: None,
            cases: vec![
                run_mom_002(),
                run_mom_003(),
                run_cpml_001(),
                run_ntff_001(),
                run_dispersive_001(),
            ],
        };
        let md = report.to_markdown();
        assert!(md.starts_with("# Yee Validation Report"));
        assert!(md.contains("mom-002"));
        let j = report.to_json().expect("json");
        assert!(j.contains("\"cases\""));
        assert!(!report.has_failures());
    }

    #[test]
    fn skipped_cases_carry_explanatory_notes() {
        for case in [
            run_mom_002(),
            run_mom_003(),
            run_cpml_001(),
            run_ntff_001(),
            run_dispersive_001(),
        ] {
            assert_eq!(case.status, CaseStatus::Skipped);
            assert!(
                !case.notes.is_empty(),
                "skipped case {} has empty notes",
                case.id
            );
        }
    }

    #[test]
    fn status_alias_resolves_to_case_status() {
        let s: Status = Status::Passed;
        assert_eq!(s, CaseStatus::Passed);
    }
}
