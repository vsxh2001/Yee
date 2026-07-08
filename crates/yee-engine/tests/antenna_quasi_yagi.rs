//! Gate `engine-antenna-005` (FS.1a.1, ADR-0205): the **quasi-Yagi** — the
//! first non-broadside, truncated-ground antenna through the pipeline.
//! `yee_layout::quasi_yagi` synthesizes the Kaneda/Deal topology from
//! scaling rules; the FS.1a.0 ground truncation makes the ground edge the
//! reflector; grid seeded by the FS.0a rulebook (no hand-set dx); measured
//! with the single-run directional |S11| (the A.1 machinery).
//!
//! **FS.1a.1b applied**: the lifted stack (`voxelize_microstrip_open`,
//! air + CPML below the mid-domain ground sheet) + `AperturePortSpec::k_lo`
//! replace the floor-ground fixture whose PEC bottom face acted as an
//! infinite image plane — the measured root cause of the first two runs'
//! |S11| ≈ 0 dB (dipole never driven; recorded in ADR-0205):
//! Two instrumented runs (133/140 s) measured |S11| ≈ 0 dB across the
//! band: the dipole is not driven. Root cause is the z-stack, not the
//! layout: `voxelize_microstrip` puts the ground at `k = 0` **on the
//! domain floor**, and the antenna boundary keeps that bottom face PEC —
//! so past the truncation the *boundary* still provides an infinite image
//! plane 1.6 mm under the dipole, annihilating its radiation resistance.
//! Fixing it needs open space below the substrate: a lifted voxel stack
//! (air + CPML below a mid-domain ground sheet) and a `k_lo` on
//! `AperturePortSpec` (which currently drives `k = 0 .. k_top`,
//! ground-at-floor hard-wired) across the protocol and both backends —
//! queued as FS.1a.1b. The 7.15–7.5 GHz dip both runs saw is a feed-
//! structure resonance (it moved only 5 % when the dipole grew 29 %).
//!
//! **GREEN (measured 2026-07-08)**: dip **5.950 GHz / −20.9 dB** vs the
//! designed 5.80 GHz → **2.6 %** error, broadband |S11| baseline −3…−4 dB
//! (real radiation), matched band ≈ 5.8–6.4 GHz below −10 dB. The ε = 1.61
//! dipole calibration verified blind. The test fn is named `antenna_…` so
//! the blanket CI engine-gates step skips it; it runs in the antenna CI
//! job.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_quasi_yagi -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{Substrate, quasi_yagi};
use yee_voxel::{VoxelOptions, truncate_ground_at_cell, voxelize_microstrip_open};

const F0_HZ: f64 = 5.8e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const AIR_BELOW_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 4.0e9;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 10;

#[test]
#[ignore = "slow: one multi-minute release FDTD run; engine-antenna-005 gate (FS.1a.1) — run with --release --ignored"]
fn antenna_quasi_yagi_resonates_at_the_designed_frequency() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let qy = quasi_yagi(F0_HZ, &sub, Z0_OHM);
    let layout = &qy.layout;
    let dims = qy.dims;

    // FS.0a rulebook seeds the grid — no hand-set dx anywhere.
    let dx = auto_dx(layout, F0_HZ + BW_HZ / 2.0);
    eprintln!(
        "engine-antenna-005: auto_dx = {:.3} mm; dipole {:.1} mm, director {:.1} mm, \
         reflector gap {:.1} mm, x_gnd {:.1} mm",
        dx * 1e3,
        dims.dipole_len_m * 1e3,
        dims.director_len_m * 1e3,
        dims.reflector_gap_m * 1e3,
        dims.x_gnd_m * 1e3,
    );

    // Lifted stack (FS.1a.1b): air + absorber BELOW the ground sheet, so
    // past the truncation the antenna sees open space underneath instead
    // of the domain floor's image plane.
    let mut model = voxelize_microstrip_open(
        layout,
        &VoxelOptions {
            dx_m: dx,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: AIR_ABOVE_CELLS,
        },
        AIR_BELOW_CELLS,
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;

    // FS.1a.0: the ground plane ends at the reflector edge.
    truncate_ground_at_cell(&mut model, i_for(dims.x_gnd_m));

    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let k_probe = k_top.saturating_sub(1).max(1);

    // Directional probe triple on the feed, inside the feed length
    // (6 mm past the port; the T-junction starts at 0.75 λg ≈ 21 mm).
    let i_a = i_for(layout.ports[0].at.x + 6.0e-3);
    let i_b = i_a + SPACING_CELLS;
    let i_c = i_b + SPACING_CELLS;

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
        dx_m: dx,
        n_steps: N_STEPS,
        // All six faces absorb: the lifted stack frees the bottom face,
        // so the only PEC in the domain is the masked (truncated) ground
        // sheet and the traces — true open-space radiation.
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
            record: false,
        }],
        probes: vec![
            ProbeSpec {
                component: "ez".into(),
                cell: (i_a, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_b, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_c, j_strip, k_probe),
            },
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
    assert_eq!(result.steps_done, N_STEPS);
    let p = &result.probes;

    // 4.0–7.6 GHz sweep, 50 MHz raster.
    let freqs: Vec<f64> = (0..=72).map(|n| 4.0e9 + n as f64 * 50.0e6).collect();
    let spacing_m = SPACING_CELLS as f64 * dx;
    let s11_db = sparams::directional_reflection_db([&p[0], &p[1], &p[2]], dt, spacing_m, &freqs);

    let (n_dip, dip_db) = s11_db
        .iter()
        .enumerate()
        .map(|(n, db)| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("empty spectrum");
    let f_dip = freqs[n_dip];

    eprintln!("engine-antenna-005: quasi-Yagi 5.8 GHz, directional |S11|(f):");
    for n in (0..freqs.len()).step_by(4) {
        eprintln!("  {:>4.2} GHz: {:>7.2} dB", freqs[n] / 1e9, s11_db[n]);
    }
    let rel_err = (f_dip - F0_HZ).abs() / F0_HZ;
    eprintln!(
        "  dip {:.3} GHz ({:.1} dB) vs designed {:.2} GHz → err {:.1} %",
        f_dip / 1e9,
        dip_db,
        F0_HZ / 1e9,
        rel_err * 100.0,
    );

    assert!(
        rel_err <= 0.10,
        "engine-antenna-005 FAILED: dip at {:.3} GHz, designed {:.2} GHz (err {:.1} % > 10 %)",
        f_dip / 1e9,
        F0_HZ / 1e9,
        rel_err * 100.0
    );
    // Measured-then-pinned (house pattern): the first green run read
    // −20.9 dB; −10 dB (≈ 2× linear margin) also IS the practical
    // "usable match" bar for an antenna.
    assert!(
        dip_db <= -10.0,
        "engine-antenna-005 FAILED: dip only {dip_db:.1} dB (measured −20.9 at ship time)"
    );
}
