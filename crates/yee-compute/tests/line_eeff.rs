//! Gate `compute-008`: the engine reproduces `fdtd-line-eeff-001` — the
//! filter pipeline's first full-wave gate (F1.1b.1, ADR-0108) — running the
//! FDTD on `yee-compute`'s CPU backend instead of `yee-fdtd`.
//!
//! Identical scenario to `crates/yee-voxel/tests/fdtd_line_eeff_001.rs`:
//! a dimensioned FR-4 microstrip line (W = 3 mm, h = 1.6 mm, ε_r = 4.4,
//! L ≈ 6 λ_g at 5 GHz) voxelized by `yee_voxel::voxelize_microstrip`, driven
//! by a 50 Ω resistive port with a modulated-Gaussian EMF, hard-PEC box,
//! time-gated two-probe phase-velocity measurement →
//! `ε_eff = (c/v_p)²` vs the **published Hammerstad–Jensen / Pozar closed
//! form** (`yee_layout::eps_eff`), ≤ 15 % relative (the walking-skeleton
//! band of the original gate).
//!
//! `#[ignore]`'d (multi-minute release run on a ~2.9 M-cell grid):
//!
//! ```bash
//! cargo test -p yee-compute --release --test line_eeff -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_compute::{
    Boundary, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, Probe, ResistivePort,
    Waveform,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const PORT_R_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const FREQ_SPAN: f64 = 0.8;

#[test]
#[ignore = "slow: multi-minute release FDTD; compute-008 line-eeff gate (E.2) — run with --release --ignored"]
fn line_eeff_on_engine_matches_hammerstad_jensen() {
    let eps_eff_ref = eps_eff(W_M, H_M, EPS_R);

    // ---- geometry: identical to fdtd-line-eeff-001 ----
    let lam_g = C0_M_S / (F0_HZ * eps_eff_ref.sqrt());
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
        ports: vec![PortRef {
            at: Point2::new(0.5e-3, W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: PORT_R_OHM,
        }],
        bbox,
    };

    // ---- voxelize with the same options, then hand the material grid to
    // the engine ----
    let model = voxelize_microstrip(
        &layout,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
    );
    let (nx, ny, nz) = model.dims;
    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let mut spec = FdtdSpec::vacuum(nx, ny, nz, DX_M);
    spec.dt = model.grid.dt;
    let dt = spec.dt;
    let dx = model.dx_m;

    let materials = Materials {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        mu_r_cells: None,
        sigma_cells: None,
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
        pec_mask_ez: None,
        sheet_r_ohm: None,
    };
    assert!(
        materials.eps_r_cells.is_some() && materials.pec_mask_ex.is_some(),
        "voxelizer no longer attaches eps/PEC — scenario broken"
    );

    // ---- probes: same mapping as run_line_eeff ----
    let k_probe = k_top.saturating_sub(1).max(1);
    let x_a = 2.5 * lam_g;
    let x_b = x_a + lam_g / 3.0;
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| -> usize {
        (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize
    };
    let (i_a, i_b) = (i_for(x_a), i_for(x_b));
    assert!(i_b > i_a, "probe planes collapsed");
    let delta_x = (i_b - i_a) as f64 * dx;

    // ---- time gate: stop before the far-wall reflection reaches probe B ----
    let v_p_ref = C0_M_S / eps_eff_ref.sqrt();
    let x_drive = 0.5e-3;
    let t_refl_b = ((l_m - x_drive) + (l_m - x_b)) / v_p_ref;
    let gate_steps = (0.9 * t_refl_b / dt) as usize;
    let n_steps = gate_steps + 200;

    // ---- drive: 50 Ω resistive port, modulated Gaussian (same recipe) ----
    let bw = FREQ_SPAN * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;
    let drive = Drive {
        soft_sources: vec![],
        ports: vec![ResistivePort {
            cell: model.port_cells[0],
            resistance: PORT_R_OHM,
            waveform: Waveform::GaussianPulse {
                v0: DRIVE_V0,
                f0: F0_HZ,
                bw,
                t0_steps,
            },
        }],
        aperture_ports: vec![],
        probes: vec![
            Probe {
                component: EComponent::Ez,
                cell: (i_a, j_strip, k_probe),
            },
            Probe {
                component: EComponent::Ez,
                cell: (i_b, j_strip, k_probe),
            },
        ],
        h_probes: vec![],
    };

    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        materials,
        Boundary::PecBox,
        drive,
    );
    engine.step_n(n_steps);

    // ---- time-gated single-bin DFT at f0, phase advance A → B ----
    let omega = 2.0 * PI * F0_HZ;
    let series = engine.probe_series();
    let mut acc = [0.0_f64; 4];
    let gate = gate_steps.min(n_steps);
    for (n, (a, b)) in series[0][..gate].iter().zip(&series[1][..gate]).enumerate() {
        let phase = omega * n as f64 * dt;
        let (c, s) = (phase.cos(), phase.sin());
        acc[0] += a * c;
        acc[1] -= a * s;
        acc[2] += b * c;
        acc[3] -= b * s;
    }
    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * PI;
    }
    while delta_phi > 2.0 * PI {
        delta_phi -= 2.0 * PI;
    }
    let v_p = omega * delta_x / delta_phi;
    let eps_eff_fdtd = (C0_M_S / v_p).powi(2);
    let rel_err = (eps_eff_fdtd - eps_eff_ref).abs() / eps_eff_ref;

    eprintln!(
        "compute-008 line-eeff on engine: grid {nx}x{ny}x{nz}, {n_steps} steps \
         (gate {gate_steps}) | Δx = {:.3} mm, Δφ = {:.4} rad, v_p = {:.4e} \
         | ε_eff = {:.4} vs HJ {:.4} → err {:.3} %",
        delta_x * 1e3,
        delta_phi,
        v_p,
        eps_eff_fdtd,
        eps_eff_ref,
        rel_err * 100.0
    );
    assert!(
        rel_err <= 0.15,
        "compute-008 FAILED: engine ε_eff = {eps_eff_fdtd:.4}, HJ = {eps_eff_ref:.4}, \
         err = {:.3} % (> 15 %)",
        rel_err * 100.0
    );
}
