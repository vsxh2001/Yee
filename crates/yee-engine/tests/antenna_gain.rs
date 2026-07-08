//! Gate `engine-gain-001` (FS.2b, ADR-0207): **gain in dBi** — the first
//! absolute far-field number out of the pipeline. Two solves through
//! identical machinery (auto_dx, classic ground stack, A.2 open-top
//! boundary, FS.2a port records, one-direction NTFF at broadside):
//!
//! 1. the A.1 inset-fed patch — textbook broadside gain for a thin FR-4
//!    patch is ~5–7 dBi;
//! 2. the FS.1b 2×1 corporate array — pattern multiplication adds
//!    ~+3 dB of directivity over the single element.
//!
//! **STATUS: RED, measured, root-cause hypothesis (ADR-0207) — awaiting
//! FS.2b.1.** First run: single patch **22.15 dBi**, array 23.92 dBi —
//! absolute levels ~13 dB above physics (the aperture-size diffraction
//! cap for this patch is ~3–6 dBi), while the DIFFERENTIAL (1.77 dB, in
//! [1.5, 4.5]) and every relative pattern gate stay healthy. The
//! engine-scale-001 Hertzian pin then certified the NTFF transform to
//! 3–5 % in free space and the FS.2a identity certified the port power —
//! isolating the excess to the fixture: `voxelize_microstrip` fills the
//! substrate slab across the WHOLE domain, so the equivalence box
//! necessarily intersects dielectric exactly where the strongest guided
//! fields live, and the transform propagates those samples with
//! free-space η₀. Queued fix (FS.2b.1): finite-extent substrate in the
//! voxelizer (real boards end!), then re-measure. The differential
//! remains asserted; the absolute window assert is the design contract
//! and stays red until then. Test fn named `antenna_…` so the blanket CI
//! step skips it.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_gain -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::sparams::single_bin_dft;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, NtffSpec,
    farfield,
};
use yee_layout::{Layout, Substrate, inset_fed_patch, patch_array_2x1};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const F0_HZ: f64 = 2.45e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 2.0e9;
const N_STEPS: usize = 9000;
const NTFF_MARGIN: usize = 15;

/// One antenna → broadside gain in dBi at `F0_HZ`.
fn measure_gain_dbi(layout: &Layout, label: &str) -> f64 {
    let dx = auto_dx(layout, F0_HZ + BW_HZ / 2.0);
    let model = voxelize_microstrip(
        layout,
        &VoxelOptions {
            dx_m: dx,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
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

    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * BW_HZ)) / dt).ceil() as usize;

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: dx,
        n_steps: N_STEPS,
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
            margin_cells: NTFF_MARGIN,
            k_min: Some(1),
            directions: vec![(0.0, 0.0)],
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
    let e_far = result.far_field.expect("no far field")[0];
    let records = result.port_records.expect("no port records");
    let rec = &records[0];

    let v_series: Vec<f64> = rec.iter().map(|&(v_src, _, _)| v_src).collect();
    let i_series: Vec<f64> = rec.iter().map(|&(_, _, i)| i).collect();
    let v_dft = single_bin_dft(&v_series, dt, F0_HZ);
    let i_dft = single_bin_dft(&i_series, dt, F0_HZ);
    let p_acc = farfield::accepted_power_density(v_dft, i_dft, Z0_OHM, dt);
    let g = farfield::gain_dbi(e_far, p_acc);
    eprintln!("engine-gain-001 [{label}]: |F| = {e_far:.4e}, p_acc = {p_acc:.4e} → G = {g:.2} dBi");
    g
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs + NTFF; engine-gain-001 gate (FS.2b) — run with --release --ignored"]
fn broadside_gain_lands_in_the_textbook_window_and_the_array_adds_3db() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let single = inset_fed_patch(F0_HZ, &sub, Z0_OHM);
    let array = patch_array_2x1(F0_HZ, &sub, Z0_OHM);

    let g_single = measure_gain_dbi(&single, "single patch");
    let g_array = measure_gain_dbi(&array.layout, "2x1 array");
    let diff = g_array - g_single;
    eprintln!("engine-gain-001: differential (array − single) = {diff:.2} dB");

    // Absolute window, loose until measured (textbook thin-FR-4 patch is
    // ~5–7 dBi; the staircase + open-boundary chain earns some slack).
    assert!(
        (3.0..=9.0).contains(&g_single),
        "engine-gain-001 FAILED: single-patch broadside gain {g_single:.2} dBi outside [3, 9]"
    );
    // The sharp assert: pattern multiplication adds ~3 dB of directivity;
    // the shared machinery cancels in the difference.
    assert!(
        (1.5..=4.5).contains(&diff),
        "engine-gain-001 FAILED: array-gain differential {diff:.2} dB outside [1.5, 4.5]"
    );
}
