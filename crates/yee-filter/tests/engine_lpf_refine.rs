//! Gate `engine-refine-001` (F1.2.1.0 / S.11-S.12, ADR-0188/0189): the
//! first **closed design loop** on the engine — the EM measurement
//! corrects the synthesis and the corrected design re-verifies on
//! target. Convergence required S.12's **directional** |S21| (3-probe
//! standing-wave fit, `sparams::directional_transmission_db`): on the
//! plain single-probe ratio the port-to-port reflected wave rippled the
//! skirt and the secant oscillated (full history in ADR-0188/0189).
//! Measured with the directional observable: seed −27.0 % →
//! +25.5 % → +10.5 % → **+1.0 %** in four map points.
//!
//! Loop (one scalar knob — the synthesis frequency — driven by a secant
//! iteration on the measured map; the loop structure, not an optimizer):
//! synthesize the N = 5 Butterworth stepped-impedance LPF at
//! `f_c = 2 GHz` → verify on the engine (S.8 machinery: CPML-xy walls,
//! S.10 aperture ports) → the measured response deviates because
//! closed-form dimensions are seeds (staircased high-Z widths,
//! un-de-embedded junctions) → correct the synthesis frequency and
//! repeat (up to 4 map points). The observable is the **whole-curve
//! fitted Butterworth cutoff** (see `verify_cutoff`). Asserts:
//!
//! 1. the seed error is real (≥ 5 % — there is something to fix);
//! 2. the final cutoff error is at most half the seed error;
//! 3. the final error is within ±10 % of the 2 GHz design target.
//!
//! `#[ignore]`'d (up to eight release FDTD solves, typically six):
//!
//! ```bash
//! cargo test -p yee-filter --release --test engine_lpf_refine -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_filter::{Approximation, dimension_stepped_impedance_layout, ideal_response_lowpass};
use yee_layout::{Layout, Polygon, Substrate};
use yee_synth::prototype;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

// The S.8-certified N = 5 scenario, unchanged: an N = 3 mini-board was
// tried first to save solve time and rejected — its short 3-section body
// leaks a parasitic over-substrate air path to the probe, putting
// +3..+5 dB transmission bumps in the stopband that corrupt the cutoff
// observable (measured, recorded in ADR-0188). The N = 5 board's skirt
// is clean and monotone.
const ORDER: usize = 5;
const F_TARGET_HZ: f64 = 2.0e9;
const Z0_OHM: f64 = 50.0;
const Z_HIGH_OHM: f64 = 120.0;
const Z_LOW_OHM: f64 = 20.0;
const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const DRIVE_V0: f64 = 1.0;
const F_DRIVE_HZ: f64 = 2.4e9;
const BW_HZ: f64 = 3.4e9;
const N_STEPS: usize = 9000;
// Probe-triple spacing for the directional fit: 17 cells = 5.1 mm keeps
// beta*d in ~(0.15, 0.85) rad over the 0.8-4.2 GHz band (away from the
// 0 / pi degeneracies).
const SPACING_CELLS: usize = 17;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// Synthesize the LPF at `f_c_synth_hz` (the refinement knob).
fn dut_layout(f_c_synth_hz: f64) -> Layout {
    let proto = prototype(Approximation::Butterworth, ORDER);
    dimension_stepped_impedance_layout(&proto, f_c_synth_hz, Z0_OHM, Z_HIGH_OHM, Z_LOW_OHM, &fr4())
        .expect("stepped-impedance synthesis failed")
}

fn reference_layout(dut: &Layout) -> Layout {
    let p0 = dut.ports[0].at;
    let p1 = dut.ports[1].at;
    let w = dut.ports[0].width_m;
    Layout {
        substrate: dut.substrate,
        traces: vec![Polygon::rect(p0.x, -w / 2.0, p1.x - p0.x, w)],
        ports: dut.ports.clone(),
        bbox: dut.bbox,
    }
}

/// Voxelize and express one run as a JobSpec (S.9 CPML-xy walls + S.10
/// aperture ports); returns `(spec, dt)`.
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
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    // Directional probe triple (S.12): three planes on the output feed,
    // SPACING_CELLS apart, outermost 3 mm before the load; the
    // standing-wave fit needs them equally spaced along +x.
    let i_p2c = i_for(layout.ports[1].at.x - 3.0e-3);
    let i_p2b = i_p2c - SPACING_CELLS;
    let i_p2a = i_p2b - SPACING_CELLS;

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
        aperture_ports: vec![
            AperturePortSpec {
                i: model.port_cells[0].0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: DRIVE_V0,
                f0_hz: F_DRIVE_HZ,
                bw_hz: BW_HZ,
                t0_steps,
                record: false,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F_DRIVE_HZ,
                bw_hz: BW_HZ,
                t0_steps,
                record: false,
            },
        ],
        probes: vec![
            ProbeSpec {
                component: "ez".into(),
                cell: (i_p2a, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_p2b, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_p2c, j_strip, k_probe),
            },
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        backend: BackendChoice::Cpu,
    };
    (spec, dt)
}

fn run(spec: JobSpec) -> Vec<Vec<f64>> {
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
    assert_eq!(result.probes.len(), 3, "expected the directional triple");
    result.probes
}

/// One full verify: synthesize at `f_c_synth_hz`, run DUT + reference,
/// return the measured **fitted cutoff**: the cutoff of the ideal
/// Butterworth |S21| curve that best matches the whole measured spectrum
/// in least squares (dB domain, samples where the ideal is above −25 dB
/// so the deep stopband cannot dominate).
///
/// A threshold-crossing detector was tried first and rejected (recorded
/// in ADR-0188): the measured skirts carry local bumps (a spurious
/// stopband response near 2.8 GHz), so "first/last −3 dB crossing" is
/// metric-dependent — two defensible detectors read the SAME board 1.7
/// and 2.9 GHz. The whole-curve fit uses all 69 points and one bump
/// cannot dominate it.
fn verify_cutoff(f_c_synth_hz: f64) -> f64 {
    let dut = dut_layout(f_c_synth_hz);
    let reference = reference_layout(&dut);
    let (dut_spec, dt) = job_for(&dut);
    let (ref_spec, dt2) = job_for(&reference);
    assert_eq!(dt, dt2, "runs must share dt");
    let ref_p = run(ref_spec);
    let dut_p = run(dut_spec);

    let freqs: Vec<f64> = (0..=68).map(|n| 0.8e9 + n as f64 * 50.0e6).collect();
    // Directional |S21| (S.12): forward-wave ratio from the three-probe
    // standing-wave fit — immune to the port-to-port reflected wave that
    // oscillated the ADR-0188 secant.
    let spacing_m = SPACING_CELLS as f64 * DX_M;
    let s21_db = sparams::directional_transmission_db(
        [&dut_p[0], &dut_p[1], &dut_p[2]],
        [&ref_p[0], &ref_p[1], &ref_p[2]],
        dt,
        spacing_m,
        &freqs,
    );
    eprintln!(
        "  spectrum for synthesis f_c = {:.3} GHz:",
        f_c_synth_hz / 1e9
    );
    for n in (0..freqs.len()).step_by(8) {
        eprintln!("    {:>4.2} GHz: {:>7.2} dB", freqs[n] / 1e9, s21_db[n]);
    }

    // Grid-search the fitted cutoff over 0.8–3.6 GHz in 10 MHz steps.
    let mut best = (f64::INFINITY, 0.0);
    let mut c = 0.8e9;
    while c <= 3.6e9 {
        let ideal = ideal_response_lowpass(Approximation::Butterworth, ORDER, c, &freqs);
        let mut sse = 0.0;
        let mut n_used = 0usize;
        for (m, i) in s21_db.iter().zip(&ideal) {
            let ideal_db = 20.0 * i.norm().log10();
            if ideal_db > -25.0 {
                let d = m - ideal_db;
                sse += d * d;
                n_used += 1;
            }
        }
        let score = sse / n_used as f64;
        if score < best.0 {
            best = (score, c);
        }
        c += 10.0e6;
    }
    let (score, f_fit) = best;
    eprintln!(
        "    → fitted Butterworth cutoff: {:.3} GHz (rms residual {:.2} dB)",
        f_fit / 1e9,
        score.sqrt()
    );
    f_fit
}

#[test]
#[ignore = "slow: up to eight release FDTD solves; engine-refine-001 gate (F1.2.1.0/S.11) — run with --release --ignored"]
fn em_in_the_loop_secant_moves_the_cutoff_to_target() {
    // Seed: closed-form dimensions at the design frequency.
    let mut x1 = F_TARGET_HZ; // synthesis frequency (the knob)
    let mut y1 = verify_cutoff(x1); // measured fitted cutoff
    let err_seed = (y1 - F_TARGET_HZ).abs() / F_TARGET_HZ;

    // First correction: multiplicative guess (lengths scale as 1/f_c).
    // The synth→measured map is NOT a constant ratio (measured 0.73 at
    // 2.0 GHz vs 0.92 at 2.74 GHz on this scenario), so subsequent steps
    // are proper secant updates on the two most recent map points.
    let mut x2 = x1 * F_TARGET_HZ / y1;
    let mut y2 = verify_cutoff(x2);
    let mut history = vec![(x1, y1), (x2, y2)];

    while (y2 - F_TARGET_HZ).abs() / F_TARGET_HZ > 0.08 && history.len() < 4 {
        let x3 = x2 + (F_TARGET_HZ - y2) * (x2 - x1) / (y2 - y1);
        let y3 = verify_cutoff(x3);
        history.push((x3, y3));
        (x1, y1) = (x2, y2);
        (x2, y2) = (x3, y3);
    }
    let err_final = (y2 - F_TARGET_HZ).abs() / F_TARGET_HZ;

    eprintln!(
        "engine-refine-001: N={ORDER} Butterworth LPF, target f_c = {:.1} GHz",
        F_TARGET_HZ / 1e9
    );
    for (n, (x, y)) in history.iter().enumerate() {
        eprintln!(
            "  iter {n}: synth {:.3} GHz → measured {:.3} GHz (err {:+.1} %)",
            x / 1e9,
            y / 1e9,
            (y - F_TARGET_HZ) / F_TARGET_HZ * 100.0
        );
    }

    assert!(
        err_seed >= 0.05,
        "engine-refine-001: seed error only {:.1} % — the premise (closed-form seeds \
         are off) no longer holds, re-scope the gate",
        err_seed * 100.0
    );
    assert!(
        err_final <= 0.5 * err_seed,
        "engine-refine-001 FAILED: refinement did not halve the error \
         (seed {:.1} % → final {:.1} %)",
        err_seed * 100.0,
        err_final * 100.0
    );
    assert!(
        err_final <= 0.10,
        "engine-refine-001 FAILED: refined cutoff still {:.1} % off target (> 10 %)",
        err_final * 100.0
    );
}
