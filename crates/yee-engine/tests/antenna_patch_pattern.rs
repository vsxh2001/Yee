//! Gate `engine-antenna-003` (A.2, ADR-0192): the far-field **radiation
//! pattern over the job protocol** — the inset-fed patch radiates a
//! broadside beam, measured by the validated `yee_fdtd::NtffState`
//! transform running inside the engine job (`JobSpec.ntff` →
//! `JobResult.far_field`), under the A.2 per-face open boundary
//! (absorbing side walls + **open top**, PEC ground).
//!
//! Physics asserted (walking-skeleton qualitative facts of a patch
//! pattern): the upper-hemisphere E-plane cut peaks at/near broadside
//! (θ = 0, +z above the patch) and falls toward the horizon — broadside
//! must beat every θ ≥ 60° direction on both sides of the cut. Documented
//! approximation: the NTFF box's bottom face (`k_min = 1`) hugs the
//! ground plane and crosses the substrate (the equivalence surface is not
//! fully in homogeneous air) — standard practice for grounded antennas,
//! hence qualitative asserts.
//!
//! `#[ignore]`'d (one multi-minute release FDTD run with per-step NTFF
//! sampling):
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_patch_pattern -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, NtffSpec,
};
use yee_layout::{Substrate, inset_fed_patch};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const F0_HZ: f64 = 2.425e9; // the A.1-measured resonance of this exact layout
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const BW_HZ: f64 = 1.0e9;
const N_STEPS: usize = 9000;
const NTFF_MARGIN: usize = 15; // npml (10) + 5, clear of the absorber

#[test]
#[ignore = "slow: one multi-minute release FDTD run + per-step NTFF; engine-antenna-003 gate (A.2) — run with --release --ignored"]
fn patch_radiates_a_broadside_beam() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let layout = inset_fed_patch(2.45e9, &sub, Z0_OHM);

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

    // E-plane cut (φ = 0 / φ = π): θ = 0° (broadside, +z) out to 80°.
    let mut directions: Vec<(f64, f64)> = vec![(0.0, 0.0)];
    for deg in [20.0_f64, 40.0, 60.0, 80.0] {
        directions.push((deg.to_radians(), 0.0));
        directions.push((deg.to_radians(), PI));
    }

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
            v0: DRIVE_V0,
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

    eprintln!(
        "engine-antenna-003: patch E-plane cut at {:.3} GHz:",
        F0_HZ / 1e9
    );
    for ((theta, phi), mag) in directions.iter().zip(&ff) {
        eprintln!(
            "  θ = {:>4.0}°, φ = {:>4.0}°: |E| = {:.4e} ({:+.1} dB vs broadside)",
            theta.to_degrees(),
            phi.to_degrees(),
            mag,
            20.0 * (mag / ff[0]).log10(),
        );
    }

    let broadside = ff[0];
    assert!(
        broadside.is_finite() && broadside > 0.0,
        "engine-antenna-003 FAILED: no broadside radiation captured"
    );
    // The beam points up: broadside beats every θ ≥ 60° direction on both
    // sides of the cut.
    for ((theta, phi), mag) in directions.iter().zip(&ff) {
        if theta.to_degrees() >= 59.0 {
            assert!(
                broadside > *mag,
                "engine-antenna-003 FAILED: |E(broadside)| = {broadside:.3e} not above \
                 |E(θ={:.0}°, φ={:.0}°)| = {mag:.3e} — beam not pointing up",
                theta.to_degrees(),
                phi.to_degrees(),
            );
        }
    }
}
