//! Filter design flow for the studio (R.5, ADR-0198): spec → synthesis →
//! dimensions → layout → design response → export artifacts, all closed
//! form and instant. The full-wave verify/refine loop stays on the
//! existing `run_job` engine path; this module is the *design* side the
//! studio drives interactively, and every artifact it returns is built by
//! the same validated library calls the CLI uses (`yee_io::touchstone`,
//! `yee_export::layout_to_gerber*`, `coupling_matrix_s_params`).

use serde::{Deserialize, Serialize};
use yee_export::{GerberOptions, OutlineOptions, layout_to_gerber, layout_to_gerber_outline};
use yee_filter::{
    Approximation, FilterSpec, Response, SpecMask, coupling_matrix_s_params,
    dimension_hairpin_with_fold, synthesize,
};
use yee_layout::{HairpinSectionParams, Substrate, hairpin_bpf_sections};

/// A filter-design request from the spec form.
#[derive(Debug, Clone, Deserialize)]
pub struct FilterDesignRequest {
    /// Centre frequency, Hz.
    pub f0_hz: f64,
    /// Fractional bandwidth `(f2 − f1)/f0`.
    pub fbw: f64,
    /// Filter order (number of resonators).
    pub order: usize,
    /// Chebyshev passband ripple in dB; `None` → Butterworth.
    #[serde(default)]
    pub ripple_db: Option<f64>,
    /// Substrate relative permittivity.
    pub eps_r: f64,
    /// Substrate height, metres.
    pub height_m: f64,
    /// Hairpin fold pitch in line widths (default 2.0).
    #[serde(default = "default_fold_widths")]
    pub fold_widths: f64,
}

fn default_fold_widths() -> f64 {
    2.0
}

/// The design response: the coupling-matrix S-parameters over the design
/// band plus ready-to-save export artifacts.
#[derive(Debug, Clone, Serialize)]
pub struct FilterDesignResponse {
    /// Frequency raster, Hz.
    pub freqs_hz: Vec<f64>,
    /// |S11| in dB at each raster point (coupling-matrix design response).
    pub s11_db: Vec<f64>,
    /// |S21| in dB at each raster point.
    pub s21_db: Vec<f64>,
    /// The design response as a Touchstone `.s2p` file (Real/Imag, GHz).
    pub s2p: String,
    /// The synthesized layout's copper layer as Gerber RS-274X.
    pub gerber_copper: String,
    /// The board outline as Gerber RS-274X (Edge.Cuts).
    pub gerber_outline: String,
    /// Synthesized dimensions, for display.
    pub line_width_m: f64,
    /// Resonator arm length (fold-corrected), metres.
    pub arm_length_m: f64,
    /// Per-section inter-resonator gaps, metres.
    pub gaps_m: Vec<f64>,
    /// qe→tap feed offset, metres.
    pub tap_offset_m: f64,
}

/// The [`yee_filter::FilterSpec`] a design request describes (shared with
/// the verify flow, R.5b).
pub fn filter_spec_for(req: &FilterDesignRequest) -> FilterSpec {
    FilterSpec {
        response: Response::Bandpass,
        approximation: match req.ripple_db {
            Some(r) => Approximation::Chebyshev { ripple_db: r },
            None => Approximation::Butterworth,
        },
        f0_hz: req.f0_hz,
        fbw: req.fbw,
        order: Some(req.order),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: req.ripple_db.unwrap_or(3.0),
            return_loss_db: 10.0,
            stopband: vec![],
        },
    }
}

/// The substrate a design request describes (shared with the verify flow).
pub fn substrate_for(req: &FilterDesignRequest) -> Substrate {
    Substrate {
        eps_r: req.eps_r,
        height_m: req.height_m,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// Pure design flow (no Tauri types), unit/e2e-testable headlessly.
pub fn design_filter_impl(req: &FilterDesignRequest) -> Result<FilterDesignResponse, String> {
    if !(req.f0_hz > 0.0 && req.fbw > 0.0 && req.fbw < 1.0) {
        return Err("f0 must be positive and 0 < FBW < 1".into());
    }
    if req.order < 2 {
        return Err("order must be at least 2 (inter-resonator coupling)".into());
    }
    let spec = filter_spec_for(req);
    let substrate = substrate_for(req);

    let project = synthesize(&spec);
    let dims = dimension_hairpin_with_fold(&project, &substrate, req.fold_widths)
        .map_err(|e| e.to_string())?;
    let layout = hairpin_bpf_sections(&HairpinSectionParams {
        substrate,
        arm_length_m: dims.arm_length_m,
        line_width_m: dims.line_width_m,
        fold_spacing_m: dims.fold_spacing_m,
        gaps_m: dims.gaps_m.clone(),
        tap_offset_m: dims.tap_offset_m,
        feed_width_m: dims.feed_width_m,
        feed_length_m: dims.arm_length_m,
    });

    // Design response over f0·(1 ± 2·FBW), 81 points.
    let n_pts = 81;
    let f_lo = req.f0_hz * (1.0 - 2.0 * req.fbw);
    let f_hi = req.f0_hz * (1.0 + 2.0 * req.fbw);
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|k| f_lo + (f_hi - f_lo) * k as f64 / (n_pts - 1) as f64)
        .collect();
    let s = coupling_matrix_s_params(&project.coupling, &freqs_hz, req.f0_hz, req.fbw);
    let s11_db: Vec<f64> = s
        .iter()
        .map(|(s11, _)| 20.0 * s11.norm().max(1e-12).log10())
        .collect();
    let s21_db: Vec<f64> = s
        .iter()
        .map(|(_, s21)| 20.0 * s21.norm().max(1e-12).log10())
        .collect();

    // Touchstone: reciprocal symmetric two-port (qe_in = qe_out for the
    // symmetric prototypes synthesis emits); the coupling-matrix model is
    // passive by construction.
    let data: Vec<Vec<num_complex::Complex64>> = s
        .iter()
        .map(|(s11, s21)| vec![*s11, *s21, *s21, *s11])
        .collect();
    let file = yee_io::touchstone::File {
        n_ports: 2,
        z0: 50.0,
        freq_unit: yee_io::touchstone::FreqUnit::GHz,
        format: yee_io::touchstone::Format::RealImag,
        freq_hz: freqs_hz.clone(),
        data,
        comments: vec![format!(
            "yee-studio filter design: N={} f0={:.3} GHz FBW={} ({})",
            req.order,
            req.f0_hz / 1e9,
            req.fbw,
            match req.ripple_db {
                Some(r) => format!("Chebyshev {r} dB"),
                None => "Butterworth".into(),
            }
        )],
    };
    let s2p = yee_io::touchstone::to_string(&file).map_err(|e| e.to_string())?;

    let gerber_copper = layout_to_gerber(&layout, &GerberOptions::default());
    let gerber_outline = layout_to_gerber_outline(&layout, &OutlineOptions::default());

    Ok(FilterDesignResponse {
        freqs_hz,
        s11_db,
        s21_db,
        s2p,
        gerber_copper,
        gerber_outline,
        line_width_m: dims.line_width_m,
        arm_length_m: dims.arm_length_m,
        gaps_m: dims.gaps_m,
        tap_offset_m: dims.tap_offset_m,
    })
}
