//! Gate `engine-antenna-007` (FS.1b.1, ADR-0206): the **2×1 corporate-fed
//! patch array** — two A.0/A.3-certified inset patches at 0.5 λ₀ H-plane
//! spacing behind a symmetric corporate tree (50 Ω spine → junction →
//! λg/4 70.7 Ω transformers → 50 Ω branches → the A.3-measured 0.25·L
//! insets). The tree is the DUT: phase balance is exact by mirror
//! symmetry, so what this gate measures is the corporate match (the
//! transformer pair presenting 50 Ω at the junction) plus the mutual-
//! coupling detune at 0.5 λ₀. Grid seeded by the FS.0a rulebook; classic
//! floor-ground stack (patches radiate up), A.2 open-top boundary,
//! single-run directional |S11| (the A.1 machinery).
//!
//! **GREEN first run (2026-07-08)**: dip **2.450 GHz / −21.1 dB — 0.0 %**
//! from design (the 25 MHz raster hit f₀ exactly); the 0.5 λ₀ mutual-
//! coupling detune is below the raster. Asserts pinned: err ≤ 5 %,
//! depth ≤ −10 dB.
//!
//! ```bash
//! cargo test -p yee-engine --release --test antenna_patch_array -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::auto_dx;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
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
const SPACING_CELLS: usize = 10;

#[test]
#[ignore = "slow: one multi-minute release FDTD run; engine-antenna-007 gate (FS.1b.1) — run with --release --ignored"]
fn antenna_patch_array_matches_at_the_designed_resonance() {
    let sub = Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let pa = patch_array_2x1(F0_HZ, &sub, Z0_OHM);
    let layout = &pa.layout;

    let dx = auto_dx(layout, F0_HZ + BW_HZ / 2.0);
    eprintln!(
        "engine-antenna-007: auto_dx = {:.3} mm; patch {:.1}×{:.1} mm, spacing {:.1} mm, \
         xfmr {:.2} mm wide × {:.1} mm",
        dx * 1e3,
        pa.dims.patch.length_m * 1e3,
        pa.dims.patch.width_m * 1e3,
        pa.dims.spacing_m * 1e3,
        pa.dims.xfmr_width_m * 1e3,
        pa.dims.xfmr_len_m * 1e3,
    );

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
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;

    let (_i_drive, j_strip, k_top) = model.port_cells[0];
    let k_probe = k_top.saturating_sub(1).max(1);

    // Probe triple on the spine (6 mm past the port; the junction is
    // λg/2 ≈ 34 mm downstream).
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
        // A.2 open-top boundary: side walls + top absorb, ground PEC.
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

    // 1.6–3.4 GHz sweep, 25 MHz raster.
    let freqs: Vec<f64> = (0..=72).map(|n| 1.6e9 + n as f64 * 25.0e6).collect();
    let spacing_m = SPACING_CELLS as f64 * dx;
    let s11_db = sparams::directional_reflection_db([&p[0], &p[1], &p[2]], dt, spacing_m, &freqs);

    let (n_dip, dip_db) = s11_db
        .iter()
        .enumerate()
        .map(|(n, db)| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("empty spectrum");
    let f_dip = freqs[n_dip];

    eprintln!("engine-antenna-007: 2×1 array directional |S11|(f):");
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
        rel_err <= 0.05,
        "engine-antenna-007 FAILED: dip at {:.3} GHz, designed {:.2} GHz (err {:.1} % > 5 %; \
         measured 0.0 % at ship time)",
        f_dip / 1e9,
        F0_HZ / 1e9,
        rel_err * 100.0
    );
    // Measured-then-pinned: first green run read −21.1 dB; −10 dB is the
    // 2× floor and the practical usable-match bar.
    assert!(
        dip_db <= -10.0,
        "engine-antenna-007 FAILED: dip only {dip_db:.1} dB (measured −21.1 at ship time)"
    );
}
