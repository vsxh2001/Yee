//! Gate `engine-closs-001` (R.0b, ADR-0202): **conductor loss on the
//! board path** — a lossless-dielectric FR-4 line whose masked trace
//! carries the resistive-sheet boundary; the measured attenuation is
//! held against **Pozar's strip conductor loss**:
//!
//! ```text
//! α_c = R_s / (Z₀ · W)   [Np/m],   R_s = √(π f μ₀ / σ)
//! ```
//!
//! `MaterialsSpec::sheet_r_ohm` rides the protocol into the compute
//! kernel's sheet relation `E_tan = R_s·K` on the masked trace edges
//! (compute-017). The R.0 loss-shaped fixture is reused verbatim: two
//! S.12 directional probe triples ~3 λ_g apart; `α = ln(|fwd_A|/|fwd_B|)/d`
//! at f₀ — reflections and the backward wave drop out.
//!
//! **Engineered σ = 5.8e4 S/m** (R_s ≈ 0.583 Ω → α_c ≈ 3.9 Np/m, ~4 dB
//! over the span): real copper's α_c here (~0.1 dB) sits below the
//! fixture's measured ±0.24 dB ripple, so the gate pins the sheet
//! mechanics and the R_s scaling at measurable SNR; real-copper
//! validation needs a high-Q resonator scenario (follow-on).
//!
//! **What the first honest run measured, and what the gate therefore
//! asserts** (two release solves):
//!
//! 1. **Linearity in R_s**: α(R_s)/α(R_s/2) = 2 within ±10 % — the sheet
//!    mechanics (dissipation ∝ R_s in the small-R_s regime), immune to
//!    every closed-form ambiguity below.
//! 2. **Absolute ratio band vs the Pozar total**: measured
//!    α_meas/α_c = 0.415. The decomposition: Pozar's `R_s/(Z₀W)` is the
//!    **strip + ground** first-order total with the ground treated
//!    strip-width (each conductor ~half); this sheet losses the **strip
//!    only** (the ground is a boundary face, not a mask — documented
//!    scope limit), and the single zero-thickness sheet dissipates on the
//!    **net** current `|H_below − H_above|²` where a thick two-faced
//!    conductor dissipates `|H_below|² + |H_above|²` — a further modest
//!    undercount. Gate band [0.30, 0.60], centred on the strip-only
//!    share.
//!
//! `#[ignore]`'d (one multi-minute release FDTD run on the 6 λ_g grid):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_conductor_loss -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, surface_resistance_ohm, voxelize_microstrip};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
/// Engineered conductor conductivity (S/m) — 1000x below copper so the
/// loss is measurable over this span (see the module docs).
const SIGMA_C: f64 = 5.8e4;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 17;

/// One release solve at sheet resistance `r_s`; returns measured α (Np/m).
fn measure_alpha(r_s: f64) -> f64 {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 6.0 * lam_g;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            eps_r: EPS_R,
            height_m: H_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
        ],
        bbox,
    };

    let model = voxelize_microstrip(
        &layout,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    // Two directional triples: A centred near 1.5 λ_g, B near 4.5 λ_g.
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_a0 = i_for(1.5 * lam_g);
    let i_b0 = i_for(4.5 * lam_g);
    let d_m = (i_b0 - i_a0) as f64 * dx;

    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() - W_M / 2.0 < -1e-12 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        sheet_r_ohm: Some(r_s),
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

    let bw = 0.8 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let mk_probe = |i: usize| ProbeSpec {
        component: "ez".into(),
        cell: (i, j_strip, k_probe),
    };
    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![
            AperturePortSpec {
                i: model.port_cells[0].0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: DRIVE_V0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
        ],
        probes: vec![
            mk_probe(i_a0),
            mk_probe(i_a0 + SPACING_CELLS),
            mk_probe(i_a0 + 2 * SPACING_CELLS),
            mk_probe(i_b0),
            mk_probe(i_b0 + SPACING_CELLS),
            mk_probe(i_b0 + 2 * SPACING_CELLS),
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        backend: BackendChoice::Cpu,
    };

    let handle = yee_engine::submit(spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    let p = &result.probes;
    let spacing_m = SPACING_CELLS as f64 * DX_M;

    // Forward-wave amplitude at each plane, f₀ only.
    let fwd_mag = |a: &[f64], b: &[f64], c: &[f64]| -> f64 {
        let va = sparams::single_bin_dft(a, dt, F0_HZ);
        let vb = sparams::single_bin_dft(b, dt, F0_HZ);
        let vc = sparams::single_bin_dft(c, dt, F0_HZ);
        let split = sparams::fit_standing_wave(va, vb, vc, spacing_m);
        (split.fwd.0.powi(2) + split.fwd.1.powi(2)).sqrt()
    };
    let mag_a = fwd_mag(&p[0], &p[1], &p[2]);
    let mag_b = fwd_mag(&p[3], &p[4], &p[5]);
    assert!(mag_a > 0.0 && mag_b > 0.0, "forward wave lost");
    (mag_a / mag_b).ln() / d_m
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs; engine-closs-001 gate (R.0b) — run with --release --ignored"]
fn conductor_attenuation_scales_and_sits_in_the_strip_only_band() {
    let r_s = surface_resistance_ohm(F0_HZ, SIGMA_C);
    let alpha_full = measure_alpha(r_s);
    let alpha_half = measure_alpha(r_s / 2.0);

    // Pozar strip+ground first-order total.
    let alpha_pozar = r_s / (Z0_OHM * W_M);
    let ratio = alpha_full / alpha_pozar;
    let linearity = alpha_full / alpha_half;

    eprintln!(
        "engine-closs-001: α(R_s) = {alpha_full:.3} Np/m, α(R_s/2) = {alpha_half:.3} Np/m \
         (linearity {linearity:.3}, want 2) | vs Pozar total α_c = {alpha_pozar:.3} Np/m → \
         ratio {ratio:.3} (strip-only band [0.30, 0.60]; R_s = {r_s:.3} Ω/sq)"
    );

    // 1. Sheet mechanics: dissipation linear in R_s (measured on the first
    //    honest run pair; ±10 %).
    assert!(
        (1.8..=2.2).contains(&linearity),
        "engine-closs-001 FAILED: α(R_s)/α(R_s/2) = {linearity:.3} not ≈ 2 — the sheet \
         dissipation is not linear in R_s"
    );
    // 2. Absolute band: strip-only share of the strip+ground Pozar total,
    //    single-sheet net-current model (first run: 0.415 — see module docs).
    assert!(
        (0.30..=0.60).contains(&ratio),
        "engine-closs-001 FAILED: α_meas/α_pozar = {ratio:.3} outside the strip-only \
         band [0.30, 0.60] (α_meas = {alpha_full:.3}, α_pozar = {alpha_pozar:.3} Np/m)"
    );
}
