//! Gate `engine-antenna-008` (FS.1b.2, ADR-0206): the 2×1 array's
//! **pattern multiplication** — the reason arrays exist, made
//! machine-checkable. Two in-phase elements at 0.5 λ₀ along y multiply
//! the patch element pattern by AF(θ) = cos(π/2·sin θ) in the **array
//! plane** (the y-z cut): −3 dB at θ = 30°, −13.6 dB at θ = 60°, null at
//! the horizon — a much faster rolloff than the single patch's gentle
//! H-plane taper. The E-plane (x-z) cut stays patch-like (broadside
//! beats θ ≥ 60°, the A.2 assert).
//!
//! Same fixture as `engine-antenna-007` (classic ground stack, A.2
//! open-top boundary, ground-hugging NTFF box `k_min = 1`), NTFF at the
//! measured-on-design 2.45 GHz.
//!
//! **GREEN first run (2026-07-08)** — measured vs AF theory in the array
//! plane: θ = 30°: −3.4 dB (AF −3.0); θ = 60°: −14.2/−13.9 dB (AF −13.6);
//! θ = 75°: −21.7/−21.2 dB; symmetric to 0.3 dB between the two sides.
//! E-plane: −0.6…−6.7 dB gentle patch taper. Pattern multiplication
//! confirmed quantitatively (~0.6 dB of AF), not just qualitatively.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_patch_array_pattern -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, NtffSpec,
};
use yee_layout::{Substrate, patch_array_2x1};
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

#[test]
#[ignore = "slow: one multi-minute release FDTD run + per-step NTFF; engine-antenna-008 gate (FS.1b.2) — run with --release --ignored"]
fn antenna_patch_array_narrows_the_beam_in_the_array_plane() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let pa = patch_array_2x1(F0_HZ, &sub, Z0_OHM);
    let layout = &pa.layout;

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

    // Index map: 0 = broadside; 1–5 = array-plane (y-z, φ = 90°) θ =
    // 15°..75°; 6–10 = mirrored (φ = 270°); 11–14 = E-plane (x-z) θ =
    // 30°/60° at φ = 0°/180°.
    let deg = |d: f64| d * PI / 180.0;
    let mut directions: Vec<(f64, f64)> = vec![(0.0, 0.0)];
    for n in 1..=5 {
        directions.push((deg(15.0 * n as f64), deg(90.0)));
    }
    for n in 1..=5 {
        directions.push((deg(15.0 * n as f64), deg(270.0)));
    }
    directions.extend([
        (deg(30.0), 0.0),
        (deg(60.0), 0.0),
        (deg(30.0), deg(180.0)),
        (deg(60.0), deg(180.0)),
    ]);

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
            record: false,
        }],
        thin_wires: vec![],
        probes: vec![],
        slice: None,
        ntff: Some(NtffSpec {
            f_hz: F0_HZ,
            margin_cells: NTFF_MARGIN,
            k_min: Some(1),
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

    let broadside = ff[0];
    eprintln!(
        "engine-antenna-008: 2×1 array cuts at {:.2} GHz:",
        F0_HZ / 1e9
    );
    for ((theta, phi), mag) in directions.iter().zip(&ff) {
        eprintln!(
            "  θ = {:>3.0}°, φ = {:>3.0}°: |E| = {:.4e} ({:+.1} dB vs broadside)",
            theta.to_degrees(),
            phi.to_degrees(),
            mag,
            20.0 * (mag / broadside).log10(),
        );
    }

    assert!(
        broadside.is_finite() && broadside > 0.0,
        "engine-antenna-008 FAILED: no broadside radiation captured"
    );
    // Array-plane narrowing: AF alone is −13.6 dB at θ = 60°; the element
    // pattern only helps. Assert θ = 60° and 75° at least 10 dB down on
    // both sides of the cut (indices 4, 5, 9, 10).
    for &n in &[4usize, 5, 9, 10] {
        let db = 20.0 * (ff[n] / broadside).log10();
        assert!(
            db <= -10.0,
            "engine-antenna-008 FAILED: array-plane θ = {:.0}°, φ = {:.0}° only {db:.1} dB \
             below broadside (need ≤ −10; AF alone predicts −13.6 at 60°)",
            directions[n].0.to_degrees(),
            directions[n].1.to_degrees(),
        );
    }
    // E-plane stays patch-like: broadside beats θ ≥ 60° (A.2 idiom).
    for &n in &[12usize, 14] {
        assert!(
            broadside > ff[n],
            "engine-antenna-008 FAILED: broadside not above E-plane θ = 60° (φ = {:.0}°)",
            directions[n].1.to_degrees(),
        );
    }
}
