//! Full-wave verify flow for the studio (R.5b, ADR-0199): run the
//! designed filter through the engine — reference through-line solve,
//! then the DUT solve — with per-step progress surfaced to the caller,
//! and return the measured directional |S21| next to the design's
//! coupling-matrix curve. The measurement fixture is the shared
//! [`yee_engine::board`] builder, so the studio measures exactly the way
//! the gates do.

use serde::{Deserialize, Serialize};
use yee_engine::board::{TwoPortBoardOptions, reference_through_line, two_port_board_job};
use yee_engine::{JobEvent, sparams};
use yee_filter::{coupling_matrix_s_params, dimension_hairpin_with_fold, synthesize};
use yee_layout::{HairpinSectionParams, Layout, hairpin_bpf_sections};

use crate::design::{FilterDesignRequest, filter_spec_for, substrate_for};

/// Feed length used for verify layouts: long enough for the 3-probe
/// triples with clearance (the R.4 gate value).
const VERIFY_FEED_LEN_M: f64 = 12.0e-3;

/// A verify request: the design spec plus solve fidelity knobs.
#[derive(Debug, Clone, Deserialize)]
pub struct FilterVerifyRequest {
    /// The design being verified.
    pub design: FilterDesignRequest,
    /// Uniform cell size, metres (default 0.2 mm — the R.4 gate value).
    #[serde(default = "default_dx")]
    pub dx_m: f64,
    /// Time steps per solve (default 13000 — the R.4 gate value).
    #[serde(default = "default_steps")]
    pub n_steps: usize,
}

fn default_dx() -> f64 {
    0.2e-3
}
fn default_steps() -> usize {
    13000
}

/// Progress for one phase of the verify flow.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct VerifyProgress {
    /// `"reference"` or `"dut"`.
    pub phase: &'static str,
    /// Steps completed in this phase.
    pub step: usize,
    /// Total steps in this phase.
    pub total: usize,
}

/// The verify result: measured vs designed |S21| over the design band.
#[derive(Debug, Clone, Serialize)]
pub struct FilterVerifyResponse {
    /// Frequency raster, Hz.
    pub freqs_hz: Vec<f64>,
    /// Measured directional |S21|, dB (S.12 3-probe observable).
    pub measured_s21_db: Vec<f64>,
    /// The design's coupling-matrix |S21|, dB, for overlay.
    pub design_s21_db: Vec<f64>,
    /// Backend the solves ran on.
    pub backend: String,
}

/// Run one engine job, forwarding progress; returns the probe series.
fn run_job(
    spec: yee_engine::JobSpec,
    phase: &'static str,
    on_progress: &mut dyn FnMut(VerifyProgress),
) -> Result<(Vec<Vec<f64>>, String), String> {
    let handle = yee_engine::submit(spec);
    for event in handle.events() {
        match event {
            JobEvent::Progress { step, total } => on_progress(VerifyProgress {
                phase,
                step,
                total,
            }),
            JobEvent::Done { result } => return Ok((result.probes, result.backend)),
            JobEvent::Error { message } => return Err(message),
        }
    }
    Err("engine stream ended without a result".into())
}

/// Full-wave verify of an arbitrary two-port board layout (the pipe the
/// filter wrapper rides; headlessly testable at reduced fidelity).
pub fn verify_layout_impl(
    layout: &Layout,
    opts: &TwoPortBoardOptions,
    freqs_hz: &[f64],
    on_progress: &mut dyn FnMut(VerifyProgress),
) -> Result<(Vec<f64>, String), String> {
    let reference = reference_through_line(layout);
    let ref_job = two_port_board_job(&reference, opts)?;
    let dut_job = two_port_board_job(layout, opts)?;
    if (ref_job.spec.nx, ref_job.spec.ny, ref_job.spec.nz)
        != (dut_job.spec.nx, dut_job.spec.ny, dut_job.spec.nz)
    {
        return Err("reference and DUT grids diverged".into());
    }
    let (ref_p, _) = run_job(ref_job.spec, "reference", on_progress)?;
    let (dut_p, backend) = run_job(dut_job.spec, "dut", on_progress)?;
    if !ref_p[3].iter().any(|v| *v != 0.0) || !dut_p[3].iter().any(|v| *v != 0.0) {
        return Err("output probes are silent — measurement fixture broken".into());
    }
    let s21 = sparams::directional_transmission_db(
        [&dut_p[3], &dut_p[4], &dut_p[5]],
        [&ref_p[3], &ref_p[4], &ref_p[5]],
        dut_job.dt_s,
        dut_job.spacing_m,
        freqs_hz,
    );
    Ok((s21, backend))
}

/// The filter verify flow: rebuild the designed hairpin with verify-length
/// feeds, measure it, and return measured vs designed curves.
pub fn verify_filter_impl(
    req: &FilterVerifyRequest,
    on_progress: &mut dyn FnMut(VerifyProgress),
) -> Result<FilterVerifyResponse, String> {
    let project = synthesize(&filter_spec_for(&req.design));
    let substrate = substrate_for(&req.design);
    let dims = dimension_hairpin_with_fold(&project, &substrate, req.design.fold_widths)
        .map_err(|e| e.to_string())?;
    let layout = hairpin_bpf_sections(&HairpinSectionParams {
        substrate,
        arm_length_m: dims.arm_length_m,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        gaps_m: dims.gaps_m.clone(),
        tap_offset_m: dims.tap_offset_m,
        feed_width_m: dims.line_width_m,
        feed_length_m: VERIFY_FEED_LEN_M,
    });

    let f0 = req.design.f0_hz;
    let bw_hz = 4.0 * req.design.fbw * f0; // drive covers the design band
    let mut opts = TwoPortBoardOptions::for_band(f0, bw_hz);
    opts.dx_m = req.dx_m;
    opts.n_steps = req.n_steps;

    // Same raster as the design response: f0·(1 ± 2·FBW), 41 points.
    let n_pts = 41;
    let f_lo = f0 * (1.0 - 2.0 * req.design.fbw);
    let f_hi = f0 * (1.0 + 2.0 * req.design.fbw);
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|k| f_lo + (f_hi - f_lo) * k as f64 / (n_pts - 1) as f64)
        .collect();

    let (measured_s21_db, backend) = verify_layout_impl(&layout, &opts, &freqs_hz, on_progress)?;
    let design_s21_db: Vec<f64> =
        coupling_matrix_s_params(&project.coupling, &freqs_hz, f0, req.design.fbw)
            .iter()
            .map(|(_, s21)| 20.0 * s21.norm().max(1e-12).log10())
            .collect();
    Ok(FilterVerifyResponse {
        freqs_hz,
        measured_s21_db,
        design_s21_db,
        backend,
    })
}
