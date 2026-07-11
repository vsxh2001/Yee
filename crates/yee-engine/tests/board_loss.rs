//! Gate `engine-loss-001` (R.0, ADR-0194): **dielectric loss on the board
//! path** — a lossy FR-4 line's measured attenuation matches the
//! **Pozar closed form** (Microwave Engineering §3.199):
//!
//! ```text
//! α_d = k₀ ε_r (ε_eff − 1) tan δ / (2 √ε_eff (ε_r − 1))   [Np/m]
//! ```
//!
//! The substrate `tan δ` maps to per-cell σ at the design frequency
//! (`yee_voxel::substrate_sigma_cells`, σ = 2π f ε₀ ε_r tan δ) riding the
//! S.5 `sigma_cells` protocol field into the E.1 lossy CA/CB update. The
//! measurement is loss-shaped for robustness: two S.12 **directional
//! probe triples** ~3 λ_g apart on the line; the forward-wave amplitude
//! ratio gives `α = ln(|fwd_A|/|fwd_B|) / d` at f₀ — reflections and the
//! backward wave drop out of the observable entirely.
//!
//! `tan δ = 0.05` (lossy but physical) puts ~3.7 dB over the span — far
//! above measurement ripple. Gate ±5 % (measured **0.1 %**).
//!
//! `#[ignore]`'d (one multi-minute release FDTD run on the 6 λ_g grid):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_loss -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, substrate_sigma_cells, voxelize_microstrip};

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
const TAN_D: f64 = 0.05;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 17;

#[test]
#[ignore = "slow: one multi-minute release FDTD run; engine-loss-001 gate (R.0) — run with --release --ignored"]
fn dielectric_attenuation_matches_the_pozar_closed_form() {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 6.0 * lam_g;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            eps_r: EPS_R,
            height_m: H_M,
            loss_tangent: TAN_D,
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

    let sigma = substrate_sigma_cells(&model, TAN_D, F0_HZ);
    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        sigma_cells: Some(sigma),
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
                k_top,
                resistance_ohm: Z0_OHM,
                v0: DRIVE_V0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
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
        spacings: None,
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
    let alpha_meas = (mag_a / mag_b).ln() / d_m;

    // Pozar §3.199 dielectric attenuation.
    let k0 = 2.0 * PI * F0_HZ / C0_M_S;
    let alpha_ref = k0 * EPS_R * (e_eff - 1.0) * TAN_D / (2.0 * e_eff.sqrt() * (EPS_R - 1.0));
    let rel_err = (alpha_meas - alpha_ref).abs() / alpha_ref;

    eprintln!(
        "engine-loss-001: |fwd_A| = {mag_a:.4e}, |fwd_B| = {mag_b:.4e} over d = {:.1} mm \
         | α_meas = {alpha_meas:.3} Np/m vs Pozar α_d = {alpha_ref:.3} Np/m → err {:.1} % \
         ({:.2} dB over the span)",
        d_m * 1e3,
        rel_err * 100.0,
        8.686 * alpha_meas * d_m,
    );

    // Measured 0.1 % (ADR-0194) — the σ map is exact at f_ref and the
    // forward-wave observable is reflection-free; gate at ±5 %.
    assert!(
        rel_err <= 0.05,
        "engine-loss-001 FAILED: α_meas = {alpha_meas:.3} Np/m vs Pozar {alpha_ref:.3} Np/m \
         (err {:.1} % > 5 %; measured 0.1 %)",
        rel_err * 100.0
    );
}
