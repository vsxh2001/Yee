//! Gate `engine-eff-001` (FS.2c, ADR-0207): **radiation efficiency** and
//! the full-sphere export — the last FS.2 deliverables. Two solves of the
//! A.1 patch on the FS.2b.1 finite-board fixture, full `sphere_grid(12,
//! 16)` NTFF raster + FS.2a port records:
//!
//! 1. **lossless** (tan δ = 0): η must sit near 1 — the sphere integral
//!    of gain over 4π recovers the accepted power (band pinned from
//!    measurement; the NTFF scale is certified to 3–5 % so η inherits
//!    ~6–10 % plus quadrature + absorber leakage);
//! 2. **lossy** (tan δ = 0.02, the real FR-4 value, via the R.0
//!    `substrate_sigma_cells` map): η must drop by a clear margin — the
//!    direction gate that makes efficiency an actual loss meter.
//!
//! Plus the byte-checkable full-sphere CSV artifact
//! (`farfield::pattern_csv`) — header + 192 rows, stable formatting.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_efficiency -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::sparams::single_bin_dft;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, NtffSpec,
    farfield,
};
use yee_layout::{Substrate, inset_fed_patch};
use yee_voxel::{VoxelOptions, substrate_sigma_cells, voxelize_finite_board};

const F0_HZ: f64 = 2.45e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const AIR_BELOW_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 2.0e9;
const N_STEPS: usize = 9000;
const N_THETA: usize = 12;
const N_PHI: usize = 16;

/// One patch solve → (η, csv) at `F0_HZ`.
fn measure_efficiency(tan_d: f64, label: &str) -> (f64, String) {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: tan_d,
        metal_thickness_m: 35e-6,
    };
    let layout = inset_fed_patch(F0_HZ, &sub, Z0_OHM);

    let dx = auto_dx(&layout, F0_HZ + BW_HZ / 2.0);
    let model = voxelize_finite_board(
        &layout,
        &VoxelOptions {
            dx_m: dx,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
        15.0e-3,
        AIR_BELOW_CELLS,
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let (_i_drive, _j_strip, k_top) = model.port_cells[0];

    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let tap_y = layout.ports[0].at.y;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx - tap_y).abs() < w_feed / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);

    let sigma_cells = if tan_d > 0.0 {
        Some(substrate_sigma_cells(&model, tan_d, F0_HZ))
    } else {
        None
    };
    let materials = MaterialsSpec {
        eps_r_cells: model
            .grid
            .eps_r_cells
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        sigma_cells,
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

    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * BW_HZ)) / dt).ceil() as usize;
    let directions = farfield::sphere_grid(N_THETA, N_PHI);

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: dx,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, true],
            faces: Some([[true, true], [true, true], [true, true]]),
        },
        sources: vec![],
        ports: vec![],
        aperture_ports: vec![AperturePortSpec {
            i: model.port_cells[0].0,
            j_lo,
            j_hi,
            k_lo: model.k_gnd,
            k_top,
            resistance_ohm: Z0_OHM,
            v0: 1.0,
            f0_hz: F0_HZ,
            bw_hz: BW_HZ,
            t0_steps,
            record: true,
        }],
        probes: vec![],
        slice: None,
        ntff: Some(NtffSpec {
            f_hz: F0_HZ,
            margin_cells: 15,
            k_min: None,
            directions: directions.clone(),
        }),
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
    let e_far = result.far_field.expect("no far field");
    let rec = &result.port_records.expect("no port records")[0];

    let v_series: Vec<f64> = rec.iter().map(|&(v_src, _, _)| v_src).collect();
    let i_series: Vec<f64> = rec.iter().map(|&(_, _, i)| i).collect();
    let p_acc = farfield::accepted_power_density(
        single_bin_dft(&v_series, dt, F0_HZ),
        single_bin_dft(&i_series, dt, F0_HZ),
        Z0_OHM,
        dt,
    );
    let eta = farfield::radiation_efficiency(&e_far, N_THETA, N_PHI, p_acc);
    let csv = farfield::pattern_csv(&directions, &e_far, p_acc);
    eprintln!("engine-eff-001 [{label}]: η = {eta:.4} (p_acc = {p_acc:.4e})");
    (eta, csv)
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs + full-sphere NTFF; engine-eff-001 gate (FS.2c) — run with --release --ignored"]
fn lossless_efficiency_is_unity_and_loss_drops_it() {
    let (eta_ll, csv) = measure_efficiency(0.0, "lossless");
    let (eta_lossy, _) = measure_efficiency(0.02, "tan δ = 0.02");

    // Full-sphere CSV artifact: header + one row per direction,
    // byte-stable for identical inputs.
    assert_eq!(csv.lines().count(), 1 + N_THETA * N_PHI, "CSV shape");
    assert!(csv.starts_with("theta_deg,phi_deg,e_far,gain_dbi\n"));

    // Lossless: η near 1. Band loose until measured, then pinned (the
    // NTFF scale alone contributes ~6–10 % in power).
    assert!(
        (0.7..=1.25).contains(&eta_ll),
        "engine-eff-001 FAILED: lossless η = {eta_ll:.4} outside [0.7, 1.25]"
    );
    // Loss is a loss: real FR-4 substrate loss eats a clear share.
    assert!(
        eta_lossy < eta_ll - 0.1,
        "engine-eff-001 FAILED: tan δ = 0.02 did not drop η \
         (lossless {eta_ll:.4} vs lossy {eta_lossy:.4})"
    );
}
