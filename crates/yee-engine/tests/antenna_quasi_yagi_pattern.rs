//! Gate `engine-antenna-006` (FS.1a.2, ADR-0205): the quasi-Yagi's
//! **end-fire pattern** — the antenna's purpose made machine-checkable,
//! and the wrong-phase-balun detector (a common-mode-fed dipole still
//! matches but radiates split/broadside instead of end-fire).
//!
//! Same lifted-stack fixture as `engine-antenna-005`; the NTFF cut runs
//! at the A.5-measured resonance (5.95 GHz). With the lifted stack the
//! NTFF box's **bottom face sits in free air below the ground sheet** —
//! better than the A.2 patch box (which hugged the ground); the −x face
//! still crosses the ground sheet + feed (unavoidable: the feed enters
//! there), which is the standard grounded-antenna approximation, hence
//! qualitative + measured-then-pinned asserts.
//!
//! Azimuth cut (θ = 90°, the substrate plane): the beam must point
//! toward **+x** (φ = 0, the director side): forward beats backward
//! (φ = 180°, over the ground/reflector) and both broadside directions
//! (φ = ±90°); F/B pinned ≥ 6 dB.
//!
//! **GREEN first run (2026-07-08)**: F/B **12.3 dB**, main lobe −1.9 dB
//! at ±30°, −6.4 dB at 60°, pattern minimum exactly at φ = 180° — the
//! balun's 180° split works and the ground edge reflects.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_quasi_yagi_pattern -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, NtffSpec,
};
use yee_layout::{Substrate, quasi_yagi};
use yee_voxel::{VoxelOptions, truncate_ground_at_cell, voxelize_microstrip_open};

const F0_HZ: f64 = 5.8e9;
const F_MEAS_HZ: f64 = 5.95e9; // the engine-antenna-005 measured dip
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const AIR_BELOW_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const BW_HZ: f64 = 4.0e9;
const N_STEPS: usize = 9000;
const NTFF_MARGIN: usize = 15; // npml (10) + 5, clear of the absorber

#[test]
#[ignore = "slow: one multi-minute release FDTD run + per-step NTFF; engine-antenna-006 gate (FS.1a.2) — run with --release --ignored"]
fn antenna_quasi_yagi_beam_fires_end_on() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let qy = quasi_yagi(F0_HZ, &sub, Z0_OHM);
    let layout = &qy.layout;
    let dims = qy.dims;

    let dx = auto_dx(layout, F0_HZ + BW_HZ / 2.0);
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
    truncate_ground_at_cell(&mut model, i_for(dims.x_gnd_m));

    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let _ = j_strip;

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

    // Azimuth cut θ = 90° (the substrate plane), φ = 0 → +x (director).
    let directions: Vec<(f64, f64)> = (0..12).map(|n| (PI / 2.0, n as f64 * PI / 6.0)).collect();

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
            record: false,
        }],
        probes: vec![],
        slice: None,
        ntff: Some(NtffSpec {
            f_hz: F_MEAS_HZ,
            margin_cells: NTFF_MARGIN,
            k_min: None,
            directions: directions.clone(),
        }),
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
    assert_eq!(result.steps_done, N_STEPS);
    let ff = result.far_field.expect("no far field returned");
    assert_eq!(ff.len(), directions.len());

    let fwd = ff[0]; // φ = 0, +x
    let bwd = ff[6]; // φ = 180°, −x (over the reflector)
    eprintln!(
        "engine-antenna-006: quasi-Yagi azimuth cut at {:.2} GHz:",
        F_MEAS_HZ / 1e9
    );
    for ((_, phi), mag) in directions.iter().zip(&ff) {
        eprintln!(
            "  φ = {:>4.0}°: |E| = {:.4e} ({:+.1} dB vs forward)",
            phi.to_degrees(),
            mag,
            20.0 * (mag / fwd).log10(),
        );
    }
    eprintln!("  front-to-back: {:.1} dB", 20.0 * (fwd / bwd).log10());

    assert!(
        fwd.is_finite() && fwd > 0.0,
        "engine-antenna-006 FAILED: no forward radiation captured"
    );
    // End-fire: forward beats backward and both broadside directions.
    // Measured-then-pinned: the first green run read F/B = 12.3 dB (real
    // quasi-Yagi territory); 6 dB is the 2× floor.
    assert!(
        20.0 * (fwd / bwd).log10() >= 6.0,
        "engine-antenna-006 FAILED: front-to-back only {:.1} dB (measured 12.3 at ship time)",
        20.0 * (fwd / bwd).log10()
    );
    for &n in &[3usize, 9] {
        assert!(
            fwd > ff[n],
            "engine-antenna-006 FAILED: forward {fwd:.3e} not above broadside φ = {:.0}° ({:.3e})",
            directions[n].1.to_degrees(),
            ff[n],
        );
    }
}
