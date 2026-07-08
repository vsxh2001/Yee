//! Gate `engine-antenna-001` (A.0, ADR-0190): the antenna track's walking
//! skeleton — a rectangular patch **synthesized by closed forms** (Balanis
//! §14.2 transmission-line model, `yee_layout::patch_antenna_dims`) and
//! **verified by the engine**: the measured |S11| dip sits at the designed
//! resonance.
//!
//! Scenario: 2.45 GHz edge-fed patch on FR-4 (ε_r 4.4, h 1.6 mm), the
//! S.9/S.10-certified measurement stack (CPML-xy walls + PEC ground/lid,
//! aperture drive port), S.7 two-run |S11| (reference = the bare feed
//! line on the same bbox/grid). Edge feeding is deliberately unmatched
//! (patch edge resistance is hundreds of ohms), so the dip is shallow but
//! **localized at the cavity resonance** — the assert is on position
//! (±10 % of the designed f₀; the closed-form model itself is good to a
//! few percent) plus a loose ≥ 2 dB prominence below the band median.
//!
//! `#[ignore]`'d (two multi-minute release FDTD runs):
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_patch_s11 -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{Layout, Substrate, edge_fed_patch};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const F0_HZ: f64 = 2.45e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const BW_HZ: f64 = 2.0e9; // drive envelope ≈ 1.45–3.45 GHz at −3 dB
const N_STEPS: usize = 9000;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// Reference: the bare feed line on the DUT's bbox → identical grid.
fn reference_layout(dut: &Layout) -> Layout {
    Layout {
        substrate: dut.substrate,
        traces: vec![dut.traces[0].clone()], // feed line only
        ports: dut.ports.clone(),
        bbox: dut.bbox,
    }
}

/// Voxelize and express one run as a JobSpec; returns `(spec, dt)`.
fn job_for(layout: &Layout) -> (JobSpec, f64) {
    let model = voxelize_microstrip(
        layout,
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
    let k_probe = k_top.saturating_sub(1).max(1);

    // P1 (incident/reflected separation plane) 12 mm down the feed.
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_p1 = i_for(layout.ports[0].at.x + 12.0e-3);

    // Aperture j band: the feed width rasterized with the voxelizer's
    // cell-centre convention (feed centred on y = 0).
    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < w_feed / 2.0 };
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
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        // One-port antenna: aperture drive at the feed end (S.10).
        aperture_ports: vec![AperturePortSpec {
            i: model.port_cells[0].0,
            j_lo,
            j_hi,
            k_lo: 0,
            k_top,
            resistance_ohm: Z0_OHM,
            v0: DRIVE_V0,
            f0_hz: F0_HZ,
            bw_hz: BW_HZ,
            t0_steps,
        }],
        probes: vec![ProbeSpec {
            component: "ez".into(),
            cell: (i_p1, j_strip, k_probe),
        }],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
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
    assert_eq!(result.steps_done, N_STEPS);
    result.probes.into_iter().next().expect("no probe series")
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs; engine-antenna-001 gate (A.0) — run with --release --ignored"]
fn patch_s11_dip_sits_at_the_designed_resonance() {
    let dut = edge_fed_patch(F0_HZ, &fr4(), Z0_OHM);
    let reference = reference_layout(&dut);

    let (dut_spec, dt) = job_for(&dut);
    let (ref_spec, dt2) = job_for(&reference);
    assert_eq!(dt, dt2, "runs must share dt");
    let ref_p1 = run(ref_spec);
    let dut_p1 = run(dut_spec);

    // 1.6–3.4 GHz, 25 MHz raster.
    let freqs: Vec<f64> = (0..=72).map(|n| 1.6e9 + n as f64 * 25.0e6).collect();
    let s11_db = sparams::reflection_db(&dut_p1, &ref_p1, dt, &freqs);

    let (n_dip, dip_db) = s11_db
        .iter()
        .enumerate()
        .map(|(n, db)| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("empty spectrum");
    let f_dip = freqs[n_dip];
    let mut sorted = s11_db.clone();
    sorted.sort_by(f64::total_cmp);
    let median_db = sorted[sorted.len() / 2];

    eprintln!("engine-antenna-001: 2.45 GHz edge-fed patch |S11|(f):");
    for n in (0..freqs.len()).step_by(6) {
        eprintln!("  {:>4.2} GHz: {:>7.2} dB", freqs[n] / 1e9, s11_db[n]);
    }
    let rel_err = (f_dip - F0_HZ).abs() / F0_HZ;
    eprintln!(
        "  dip {:.3} GHz ({:.1} dB) vs designed {:.2} GHz → err {:.1} % | band median {:.1} dB",
        f_dip / 1e9,
        dip_db,
        F0_HZ / 1e9,
        rel_err * 100.0,
        median_db,
    );

    assert!(
        rel_err <= 0.10,
        "engine-antenna-001 FAILED: |S11| dip at {:.3} GHz, designed {:.2} GHz (err {:.1} % > 10 %)",
        f_dip / 1e9,
        F0_HZ / 1e9,
        rel_err * 100.0
    );
    assert!(
        dip_db <= median_db - 2.0,
        "engine-antenna-001 FAILED: dip {:.1} dB not ≥ 2 dB below the band median {:.1} dB",
        dip_db,
        median_db
    );
}
