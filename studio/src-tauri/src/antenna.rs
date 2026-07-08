//! Antenna design + verify flow for the studio (R.5c, ADR-0203): the
//! filter panel's pattern applied to the A-track — instant closed-form
//! patch design (Balanis dims + inset feed + Gerber artifacts), then a
//! full-wave single-run |S11| verify over the job protocol (the A.1/A.3
//! measurement: directional 3-probe reflection under the A.2 open-top
//! per-face CPML boundary).

use serde::{Deserialize, Serialize};
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_export::{GerberOptions, OutlineOptions, layout_to_gerber, layout_to_gerber_outline};
use yee_layout::{Layout, Substrate, inset_fed_patch_with_depth, patch_antenna_dims};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

use crate::verify::VerifyProgress;

/// An antenna-design request from the spec form.
#[derive(Debug, Clone, Deserialize)]
pub struct AntennaDesignRequest {
    /// Design (resonance) frequency, Hz.
    pub f0_hz: f64,
    /// Substrate relative permittivity.
    pub eps_r: f64,
    /// Substrate height, metres.
    pub height_m: f64,
    /// Inset depth as a fraction of the patch length. The closed-form
    /// seed (`x₀ = (L/π)·acos(√(Z0/R_edge))`) lands near 0.40·L but the
    /// G₁-only slot model overestimates R_edge; the A.3 design loop
    /// measured −25.7 dB at **0.25·L** — the default here.
    #[serde(default = "default_inset_frac")]
    pub inset_frac: f64,
}

fn default_inset_frac() -> f64 {
    0.25
}

/// The design response: dims + ready-to-save Gerber artifacts (instant).
#[derive(Debug, Clone, Serialize)]
pub struct AntennaDesignResponse {
    /// Radiating-edge width, metres (Balanis §14.2).
    pub width_m: f64,
    /// Resonant length, metres.
    pub length_m: f64,
    /// Effective permittivity at the patch width.
    pub eps_eff: f64,
    /// Inset depth used, metres.
    pub inset_m: f64,
    /// Copper Gerber (RS-274X).
    pub gerber_copper: String,
    /// Outline Gerber (Edge.Cuts).
    pub gerber_outline: String,
}

/// The verify response: the measured single-run directional |S11|.
#[derive(Debug, Clone, Serialize)]
pub struct AntennaVerifyResponse {
    /// Frequency raster, Hz.
    pub freqs_hz: Vec<f64>,
    /// Measured |S11|, dB.
    pub s11_db: Vec<f64>,
    /// Deepest-dip frequency, Hz.
    pub f_dip_hz: f64,
    /// Dip depth, dB.
    pub dip_db: f64,
    /// Backend the solve ran on.
    pub backend: String,
}

fn substrate_for(req: &AntennaDesignRequest) -> Substrate {
    Substrate {
        eps_r: req.eps_r,
        height_m: req.height_m,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

fn layout_for(req: &AntennaDesignRequest) -> Result<Layout, String> {
    if !(req.f0_hz > 0.0 && req.eps_r > 1.0 && req.height_m > 0.0) {
        return Err("f0 and h must be positive, eps_r > 1".into());
    }
    if !(0.02..=0.49).contains(&req.inset_frac) {
        return Err("inset_frac must be in [0.02, 0.49] (fraction of the patch length)".into());
    }
    let dims = patch_antenna_dims(req.f0_hz, req.eps_r, req.height_m);
    Ok(inset_fed_patch_with_depth(
        req.f0_hz,
        &substrate_for(req),
        50.0,
        req.inset_frac * dims.length_m,
    ))
}

/// Pure design flow: dims + Gerbers, instant.
pub fn design_antenna_impl(req: &AntennaDesignRequest) -> Result<AntennaDesignResponse, String> {
    let dims = patch_antenna_dims(req.f0_hz, req.eps_r, req.height_m);
    let layout = layout_for(req)?;
    Ok(AntennaDesignResponse {
        width_m: dims.width_m,
        length_m: dims.length_m,
        eps_eff: dims.eps_eff,
        inset_m: req.inset_frac * dims.length_m,
        gerber_copper: layout_to_gerber(&layout, &GerberOptions::default()),
        gerber_outline: layout_to_gerber_outline(&layout, &OutlineOptions::default()),
    })
}

/// A verify request: the design plus solve fidelity knobs.
#[derive(Debug, Clone, Deserialize)]
pub struct AntennaVerifyRequest {
    /// The design being verified.
    pub design: AntennaDesignRequest,
    /// Uniform cell size, metres (default 0.3 mm — the A-gate value).
    #[serde(default = "default_dx")]
    pub dx_m: f64,
    /// Time steps (default 9000 — the A-gate value).
    #[serde(default = "default_steps")]
    pub n_steps: usize,
}

fn default_dx() -> f64 {
    0.3e-3
}
fn default_steps() -> usize {
    9000
}

/// Full-wave verify: one engine solve, the A.1 single-run directional
/// |S11| under the A.2 open-top boundary.
pub fn verify_antenna_impl(
    req: &AntennaVerifyRequest,
    on_progress: &mut dyn FnMut(VerifyProgress),
) -> Result<AntennaVerifyResponse, String> {
    const SPACING_CELLS: usize = 17;
    let f0 = req.design.f0_hz;
    let layout = layout_for(&req.design)?;
    let model = voxelize_microstrip(
        &layout,
        &VoxelOptions {
            dx_m: req.dx_m,
            xy_margin_cells: 34,
            air_above_cells: 34,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - 34.0 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_a = i_for(layout.ports[0].at.x + 12.0e-3);

    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - 34.0 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < w_feed / 2.0 };
    let j_lo = (0..ny)
        .find(|&j| in_band(j))
        .ok_or("feed band rasterized to zero cells")?;
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ex: model
            .grid
            .pec_mask_ex
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        pec_mask_ey: model
            .grid
            .pec_mask_ey
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };

    let bw_hz = 0.8 * f0;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw_hz)) / dt)
            .ceil() as usize;

    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, j_strip, k_probe),
    };
    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: req.dx_m,
        n_steps: req.n_steps,
        // Open top + PEC ground, CPML side walls (A.2, ADR-0192).
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, true],
            faces: Some([[true, true], [true, true], [false, true]]),
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![AperturePortSpec {
            i: model.port_cells[0].0,
            j_lo,
            j_hi,
            k_lo: 0,
            k_top,
            resistance_ohm: 50.0,
            v0: 1.0,
            f0_hz: f0,
            bw_hz,
            t0_steps,
            record: false,
        }],
        probes: vec![
            mk_probe(i_a),
            mk_probe(i_a + SPACING_CELLS),
            mk_probe(i_a + 2 * SPACING_CELLS),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        backend: BackendChoice::Cpu,
    };

    let handle = yee_engine::submit(spec);
    let mut probes = None;
    let mut backend = String::new();
    for event in handle.events() {
        match event {
            JobEvent::Progress { step, total } => on_progress(VerifyProgress {
                phase: "antenna",
                step,
                total,
            }),
            JobEvent::Done { result } => {
                backend = result.backend;
                probes = Some(result.probes);
                break;
            }
            JobEvent::Error { message } => return Err(message),
        }
    }
    let p = probes.ok_or("engine stream ended without a result")?;
    if !p[0].iter().any(|v| *v != 0.0) {
        return Err("feed probe silent — fixture broken".into());
    }

    // f0·(1 ± 0.35), 57 points — the A-gate raster shape.
    let n_pts = 57;
    let f_lo = 0.65 * f0;
    let f_hi = 1.35 * f0;
    let freqs_hz: Vec<f64> = (0..n_pts)
        .map(|k| f_lo + (f_hi - f_lo) * k as f64 / (n_pts - 1) as f64)
        .collect();
    let s11_db = sparams::directional_reflection_db(
        [&p[0], &p[1], &p[2]],
        dt,
        SPACING_CELLS as f64 * dx,
        &freqs_hz,
    );
    let (n_dip, dip_db) = s11_db
        .iter()
        .enumerate()
        .map(|(n, db)| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .ok_or("empty spectrum")?;
    Ok(AntennaVerifyResponse {
        f_dip_hz: freqs_hz[n_dip],
        dip_db,
        freqs_hz,
        s11_db,
        backend,
    })
}
