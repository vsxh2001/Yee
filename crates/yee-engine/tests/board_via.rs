//! Gate `engine-via-001` (R.1, ADR-0194): **vias on the board path** — a
//! ground-to-trace via (`yee_voxel::with_via_at_cell`, a vertical PEC
//! column of `E_z` edges riding the existing `pec_mask_ez` protocol
//! field) demonstrably changes the circuit the way transmission-line
//! theory says it must.
//!
//! Differential scenario on ONE grid (the S.6 λ/4 open-stub board):
//!
//! - **Control** (open stub): the stub is λ/4 at 5 GHz → input short →
//!   deep |S21| notch at 5 GHz (the certified S.6 physics).
//! - **DUT** (same stub + via at its far end): a **shorted** λ/4 stub is
//!   an open circuit at its input → the 5 GHz notch must **vanish**.
//!
//! Both runs share the reference (bare line) run for the transmission
//! ratio — three solves total. Asserts: control notch ≤ −8 dB at
//! ~5 GHz; via variant ≥ −3 dB there (the notch is gone).
//!
//! `#[ignore]`'d (three multi-minute release FDTD runs):
//!
//! ```bash
//! cargo test -p yee-engine --release --test board_via -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, voxelize_microstrip, with_via_at_cell};

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
const N_STEPS: usize = 9000;

fn open_end_delta_l(w_m: f64, h_m: f64, e_eff: f64) -> f64 {
    let u = w_m / h_m;
    0.412 * h_m * ((e_eff + 0.3) * (u + 0.264)) / ((e_eff - 0.258) * (u + 0.8))
}

/// Build one run: the stub board (or the bare-line reference), optionally
/// with a via at the stub's far end. Returns `(spec, dt, freqs helper
/// inputs)` — same grid for every variant (shared DUT bbox).
fn job(with_stub: bool, with_via: bool) -> (JobSpec, f64) {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);

    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    let bbox = BBox::from_polygons(&[line.clone(), stub.clone()]);
    let traces = if with_stub {
        vec![line, stub]
    } else {
        vec![line]
    };
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

    let mut model = voxelize_microstrip(
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

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let j_for = |yp: f64| ((yp - y0) / dx).round().clamp(0.0, ny as f64 - 1.0) as usize;

    if with_via {
        // Via at the stub's far-end centre: shorts the open end to ground.
        let i_via = i_for(l_m / 2.0);
        let j_via = j_for(W_M + stub_len - 0.5e-3);
        with_via_at_cell(&mut model, i_via, j_via, k_top);
    }

    let i_m = i_for(l_m - 3.0e-3);

    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < W_M / 2.0 };
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
        pec_mask_ez: model
            .grid
            .pec_mask_ez
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };

    let bw = 0.8 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

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
        probes: vec![ProbeSpec {
            component: "ez".into(),
            cell: (i_m, j_strip, k_probe),
        }],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
        backend: BackendChoice::Cpu,
    };
    (spec, dt)
}

fn run(spec: JobSpec) -> Vec<f64> {
    let handle = yee_engine::submit(spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    result.probes.into_iter().next().expect("no probe")
}

#[test]
#[ignore = "slow: three multi-minute release FDTD runs; engine-via-001 gate (R.1) — run with --release --ignored"]
fn via_shorted_stub_removes_the_quarter_wave_notch() {
    let (ref_spec, dt) = job(false, false);
    let (open_spec, _) = job(true, false);
    let (via_spec, _) = job(true, true);

    let reference = run(ref_spec);
    let open_stub = run(open_spec);
    let via_stub = run(via_spec);

    let freqs: Vec<f64> = (0..=40).map(|n| 4.5e9 + n as f64 * 25.0e6).collect();
    let s21_open = sparams::transmission_db(&open_stub, &reference, dt, &freqs);
    let s21_via = sparams::transmission_db(&via_stub, &reference, dt, &freqs);

    let min_of = |v: &[f64]| {
        v.iter()
            .enumerate()
            .map(|(n, db)| (n, *db))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .expect("empty")
    };
    let (n_open, open_db) = min_of(&s21_open);
    let via_at_notch = s21_via[n_open];

    eprintln!(
        "engine-via-001: open stub notch {:.1} dB at {:.3} GHz | with via at the stub end: \
         {:.1} dB at that frequency (min over band {:.1} dB)",
        open_db,
        freqs[n_open] / 1e9,
        via_at_notch,
        min_of(&s21_via).1,
    );

    assert!(
        open_db <= -8.0,
        "engine-via-001 FAILED: control (open-stub) notch only {open_db:.1} dB — \
         the S.6 physics regressed"
    );
    assert!(
        via_at_notch >= -3.0,
        "engine-via-001 FAILED: with the via the notch frequency still reads \
         {via_at_notch:.1} dB — the via is not shorting the stub"
    );
}
