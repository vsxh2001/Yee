//! Board import flow for the studio (FS.3.1b, ADR-0209): Gerber copper
//! (+ optional Edge.Cuts outline) → [`Layout`] → preview + round-trip
//! artifacts. Import is the studio's "bring your own board" door; the
//! parse is the gate-certified `yee_export::import` (gerber-rt-001/002),
//! and the response carries a **re-export of what was understood** so the
//! UI (and the e2e gate) can prove losslessness byte-for-byte before the
//! user trusts a verify run on it.
//!
//! Gerber carries no stackup or ports — the request supplies both, same
//! as every commercial import dialog.

use serde::{Deserialize, Serialize};
use yee_export::{GerberOptions, gerber_to_layout, gerber_to_outline, layout_to_gerber};
use yee_layout::{Point2, PortRef, Substrate};

/// A port the user places on the imported board.
#[derive(Debug, Clone, Deserialize)]
pub struct ImportPort {
    /// Port x, metres (layout frame).
    pub x_m: f64,
    /// Port y, metres.
    pub y_m: f64,
    /// Feed width at the port, metres.
    pub width_m: f64,
    /// Reference impedance, ohms.
    pub z0_ohm: f64,
}

/// A board-import request from the import dialog.
#[derive(Debug, Clone, Deserialize)]
pub struct ImportRequest {
    /// The copper-layer Gerber text.
    pub copper_gerber: String,
    /// Optional Edge.Cuts outline Gerber text.
    #[serde(default)]
    pub outline_gerber: Option<String>,
    /// Substrate relative permittivity.
    pub eps_r: f64,
    /// Substrate height, metres.
    pub height_m: f64,
    /// Dielectric loss tangent.
    #[serde(default)]
    pub loss_tangent: f64,
    /// User-placed ports (Gerber carries none).
    pub ports: Vec<ImportPort>,
}

/// The import response: what was understood, provably.
#[derive(Debug, Clone, Serialize)]
pub struct ImportResponse {
    /// Number of copper polygons imported.
    pub trace_count: usize,
    /// Layout bbox width, metres.
    pub bbox_w_m: f64,
    /// Layout bbox height, metres.
    pub bbox_h_m: f64,
    /// SVG preview of the imported layout.
    pub svg: String,
    /// Re-export of the imported copper — byte-identical to the input
    /// when the input came from this suite (the gerber-rt-001 property),
    /// and the exact record of what was understood otherwise.
    pub gerber_copper_echo: String,
    /// Outline corner coordinates (metres), if an outline was supplied.
    pub outline_m: Option<Vec<(f64, f64)>>,
    /// The imported layout as JSON, ready to feed the verify flow.
    pub layout_json: String,
}

/// Core import flow (pure; the Tauri command wraps this).
pub fn import_gerber_impl(req: &ImportRequest) -> Result<ImportResponse, String> {
    if req.ports.is_empty() {
        return Err("place at least one port on the imported board".into());
    }
    if !(req.eps_r >= 1.0 && req.height_m > 0.0) {
        return Err("substrate needs eps_r >= 1 and height > 0".into());
    }
    let substrate = Substrate {
        eps_r: req.eps_r,
        height_m: req.height_m,
        loss_tangent: req.loss_tangent,
        metal_thickness_m: 35e-6,
    };
    let ports: Vec<PortRef> = req
        .ports
        .iter()
        .map(|p| PortRef {
            at: Point2::new(p.x_m, p.y_m),
            width_m: p.width_m,
            ref_impedance_ohm: p.z0_ohm,
        })
        .collect();
    let layout =
        gerber_to_layout(&req.copper_gerber, substrate, ports).map_err(|e| e.to_string())?;
    let outline_m = match &req.outline_gerber {
        None => None,
        Some(g) => Some(
            gerber_to_outline(g)
                .map_err(|e| format!("outline: {e}"))?
                .into_iter()
                .map(|p| (p.x, p.y))
                .collect(),
        ),
    };
    let echo = layout_to_gerber(&layout, &GerberOptions::default());
    let layout_json = serde_json::to_string(&layout).map_err(|e| e.to_string())?;
    Ok(ImportResponse {
        trace_count: layout.traces.len(),
        bbox_w_m: layout.bbox.width(),
        bbox_h_m: layout.bbox.height(),
        svg: layout.to_svg(),
        gerber_copper_echo: echo,
        outline_m,
        layout_json,
    })
}
