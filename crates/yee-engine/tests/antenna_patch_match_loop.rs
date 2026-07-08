//! Gate `engine-antenna-004` (A.3, ADR-0193): the **antenna design loop**
//! — the engine tunes the inset depth against the measured return loss,
//! closing the closed-form model gap that A.1 exposed (the G₁-only slot
//! model badly misestimates the edge resistance on thick / high-ε_r
//! FR-4, so the formula inset realizes only a −1.2 dB match).
//!
//! Loop: coarse scan of the inset depth `x ∈ {0.10, 0.20, 0.30, 0.40}·L`
//! (0.40·L ≈ the closed-form seed), then one parabolic/neighbour
//! refinement around the best point — each evaluation is one engine job
//! measuring the single-run directional |S11| dip (A.1 observable) under
//! the open-top boundary (A.2). Asserts: the tuned inset reaches a
//! well-matched ≤ −15 dB (measured **−25.7 dB at 0.25·L**) and improves
//! on the closed-form seed by ≥ 10 dB (measured 24.8 dB).
//!
//! `#[ignore]`'d (five release FDTD solves, ~25 min):
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_patch_match_loop -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{Substrate, inset_fed_patch_with_depth, patch_antenna_dims};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const F0_HZ: f64 = 2.45e9;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const Z0_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const BW_HZ: f64 = 2.0e9;
const N_STEPS: usize = 9000;
const SPACING_CELLS: usize = 17;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// One engine job at inset depth `x_inset_m`; returns `(f_dip, dip_db)`
/// from the single-run directional |S11| (the A.1 observable, under the
/// A.2 open-top boundary).
fn measure_dip(x_inset_m: f64) -> (f64, f64) {
    let layout = inset_fed_patch_with_depth(F0_HZ, &fr4(), Z0_OHM, x_inset_m);

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
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_a = i_for(layout.ports[0].at.x + 12.0e-3);
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
    let p = &result.probes;
    let freqs: Vec<f64> = (0..=72).map(|n| 1.6e9 + n as f64 * 25.0e6).collect();
    let s11_db = sparams::directional_reflection_db(
        [&p[0], &p[1], &p[2]],
        dt,
        SPACING_CELLS as f64 * DX_M,
        &freqs,
    );
    let (n_dip, dip_db) = s11_db
        .iter()
        .enumerate()
        .map(|(n, db)| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("empty spectrum");
    (freqs[n_dip], dip_db)
}

#[test]
#[ignore = "slow: five release FDTD solves (~25 min); engine-antenna-004 gate (A.3) — run with --release --ignored"]
fn design_loop_tunes_the_inset_to_a_real_match() {
    let dims = patch_antenna_dims(F0_HZ, EPS_R, H_M);
    let l = dims.length_m;

    // Coarse scan; 0.40·L ≈ the closed-form seed depth (0.396·L).
    let fracs = [0.10_f64, 0.20, 0.30, 0.40];
    let mut points: Vec<(f64, f64, f64)> = Vec::new(); // (frac, f_dip, dip_db)
    for &fr in &fracs {
        let (f_dip, dip_db) = measure_dip(fr * l);
        eprintln!(
            "  inset {:.2}·L: dip {:.3} GHz, {:.1} dB",
            fr,
            f_dip / 1e9,
            dip_db
        );
        points.push((fr, f_dip, dip_db));
    }
    let seed_db = points.last().expect("scan empty").2; // 0.40·L ≈ closed form

    // One refinement: midpoint between the best point and its better
    // neighbour.
    let best_idx = (0..points.len())
        .min_by(|&a, &b| points[a].2.total_cmp(&points[b].2))
        .expect("scan empty");
    let neighbour = match best_idx {
        0 => 1,
        i if i == points.len() - 1 => points.len() - 2,
        i => {
            if points[i - 1].2 < points[i + 1].2 {
                i - 1
            } else {
                i + 1
            }
        }
    };
    let fr_mid = 0.5 * (points[best_idx].0 + points[neighbour].0);
    let (f_dip, dip_db) = measure_dip(fr_mid * l);
    eprintln!(
        "  inset {:.2}·L (refined): dip {:.3} GHz, {:.1} dB",
        fr_mid,
        f_dip / 1e9,
        dip_db
    );
    points.push((fr_mid, f_dip, dip_db));

    let best = points
        .iter()
        .min_by(|a, b| a.2.total_cmp(&b.2))
        .expect("no points");
    eprintln!(
        "engine-antenna-004: best inset {:.2}·L → dip {:.3} GHz at {:.1} dB \
         (closed-form seed depth 0.40·L: {:.1} dB)",
        best.0,
        best.1 / 1e9,
        best.2,
        seed_db,
    );

    // Measured (ADR-0193): 0.10·L −6.5 dB, 0.20·L −13.2 dB, 0.30·L
    // −9.3 dB, 0.40·L (closed-form seed) −0.9 dB, refined 0.25·L
    // **−25.7 dB** — the loop turns an unusable closed-form match into a
    // well-matched antenna. Thresholds set with ~10 dB headroom.
    assert!(
        best.2 <= -15.0,
        "engine-antenna-004 FAILED: best return loss {:.1} dB (need ≤ −15 dB; measured −25.7)",
        best.2
    );
    assert!(
        best.2 <= seed_db - 10.0,
        "engine-antenna-004 FAILED: loop did not beat the closed-form seed by ≥ 10 dB \
         (seed {seed_db:.1} dB → best {:.1} dB; measured improvement 24.8 dB)",
        best.2
    );
    let f_err = (best.1 - F0_HZ).abs() / F0_HZ;
    assert!(
        f_err <= 0.10,
        "engine-antenna-004 FAILED: matched resonance at {:.3} GHz drifted {:.1} % from design",
        best.1 / 1e9,
        f_err * 100.0
    );
}
