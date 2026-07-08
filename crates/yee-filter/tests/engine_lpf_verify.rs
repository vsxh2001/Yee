//! Gate `engine-filter-verify-001` (F1.3.0 / S.8, ADR-0185): the first
//! filter **synthesized by the pipeline** and **verified by the engine**
//! against its own design intent — the closed loop the project exists for.
//!
//! Design side (all shipped closed forms): `yee_synth::prototype`
//! (Butterworth N = 5) → `dimension_stepped_impedance_layout`
//! (f_c = 2 GHz, Z₀ = 50 Ω, Z_high = 120 Ω, Z_low = 20 Ω, FR-4) → a
//! two-port `Layout`. Verify side (S.5–S.7 machinery): voxelize, run the
//! DUT and a same-grid Z₀ through-line reference as two `JobSpec`s over
//! the job protocol, extract |S21|(f) (transmission ratio) and |S11|(f)
//! (incident/reflected separation), and hold the measurement against the
//! `ideal_response_lowpass` design targets:
//!
//! 1. measured −3 dB cutoff within ±20 % of the designed 2 GHz;
//! 2. passband(1 GHz) − stopband(4 GHz) rejection ≥ 20 dB (ideal 30.1 dB);
//! 3. passband mean within ±3 dB of 0 dB and return loss ≤ −6 dB — the
//!    absolute bounds the S.9 CPML-xy walls + S.10 aperture ports earn
//!    (measured −0.49 dB / −9.2 dB).
//!
//! Known measurement limits (documented, accepted at walking-skeleton
//! tolerance): dx = 0.3 mm staircases the ~0.55 mm high-Z sections to
//! ~2 cells (impedance error → the measured 15 % cutoff shift — a
//! design-side error for F1.2.1's EM-in-the-loop refinement to close);
//! feeds/junctions are not de-embedded.
//!
//! `#[ignore]`'d (two multi-minute release FDTD runs):
//!
//! ```bash
//! cargo test -p yee-filter --release --test engine_lpf_verify -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_filter::{Approximation, dimension_stepped_impedance_layout, ideal_response_lowpass};
use yee_layout::{Layout, Polygon, Substrate};
use yee_synth::prototype;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_voxel::{MicrostripModel, VoxelOptions, voxelize_microstrip};

const ORDER: usize = 5;
const F_C_HZ: f64 = 2.0e9;
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
const BW_HZ: f64 = 3.4e9; // −3 dB drive envelope ≈ 0.7–4.1 GHz
const N_STEPS: usize = 9000;

fn fr4() -> Substrate {
    Substrate {
        eps_r: EPS_R,
        height_m: H_M,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// The synthesized DUT layout (design side of the loop).
fn dut_layout() -> Layout {
    let proto = prototype(Approximation::Butterworth, ORDER);
    dimension_stepped_impedance_layout(&proto, F_C_HZ, Z0_OHM, Z_HIGH_OHM, Z_LOW_OHM, &fr4())
        .expect("stepped-impedance synthesis failed")
}

/// The reference: a straight Z₀ through line spanning the same
/// port-to-port extent, on the DUT's bbox → the identical voxel grid.
fn reference_layout(dut: &Layout) -> Layout {
    let p0 = dut.ports[0].at;
    let p1 = dut.ports[1].at;
    let w = dut.ports[0].width_m;
    Layout {
        substrate: dut.substrate,
        traces: vec![Polygon::rect(p0.x, -w / 2.0, p1.x - p0.x, w)],
        ports: dut.ports.clone(),
        bbox: dut.bbox, // shared bbox = shared grid dims and origin
    }
}

/// Voxelize and express one run as a JobSpec; returns (spec, dt).
fn job_for(layout: &Layout) -> (JobSpec, f64) {
    let model: MicrostripModel = voxelize_microstrip(
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

    // P2 (transmission) 3 mm before the load; P1 (reflection reference
    // plane) 12 mm after the drive, on the input feed.
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_p2 = i_for(layout.ports[1].at.x - 3.0e-3);
    let i_p1 = i_for(layout.ports[0].at.x + 12.0e-3);

    // Aperture j band: the feed-trace width, rasterized with the same
    // cell-centre convention as the voxelizer (feed centred on y = 0).
    let w_feed = layout.ports[0].width_m;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx).abs() < w_feed / 2.0 };
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("feed band empty");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);
    assert!(j_hi > j_lo, "aperture band empty");

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
        // Side-wall CPML, PEC ground/lid (S.9, ADR-0186). Three boundary
        // options were measured on this scenario:
        // - PEC box: correct LPF shape but every DUT reflection becomes a
        //   cavity mode — single-probe ripple up to +17.8 dB.
        // - All-face CPML (npml = 10): |S21| collapsed below −3 dB across
        //   the whole band. Root cause (ADR-0186): the substrate is ~5
        //   cells tall, so the ENTIRE line sat inside the 10-layer z-min
        //   absorber and propagated ~76 mm through it.
        // - CPML on x/y only (this setting): side-wall cavity modes are
        //   absorbed while the ground plane and lid stay PEC — the
        //   correct board-level open boundary for a thin stack.
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, false],
            faces: None,
        },
        sources: vec![],
        ports: vec![],
        // Aperture ports (S.10, ADR-0187): the validated
        // LumpedRlcPort::aperture scheme — one aggregate 50 Ω branch over
        // the modal (y, z) port face (trace width × substrate height) —
        // the ADR-0186 residual-ripple fix. The j band is the feed-trace
        // width rasterized with the voxelizer's cell-centre convention.
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
            },
            // Passive matched load (v0 = 0 → pure resistor branch).
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
            },
        ],
        probes: vec![
            ProbeSpec {
                component: "ez".into(),
                cell: (i_p2, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_p1, j_strip, k_probe),
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

/// Run one job; returns `(p2_series, p1_series)`.
fn run(spec: JobSpec) -> (Vec<f64>, Vec<f64>) {
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
    let mut probes = result.probes.into_iter();
    let p2 = probes.next().expect("no transmission probe");
    let p1 = probes.next().expect("no P1 probe");
    (p2, p1)
}

/// Mean of the dB samples whose frequency lies within ±window of centre.
fn band_mean_db(freqs: &[f64], db: &[f64], centre_hz: f64, window_hz: f64) -> f64 {
    let vals: Vec<f64> = freqs
        .iter()
        .zip(db)
        .filter(|(f, _)| (**f - centre_hz).abs() <= window_hz)
        .map(|(_, v)| *v)
        .collect();
    assert!(!vals.is_empty(), "no samples near {centre_hz}");
    vals.iter().sum::<f64>() / vals.len() as f64
}

#[test]
#[ignore = "slow: two multi-minute release FDTD runs; engine-filter-verify-001 gate (F1.3.0/S.8) — run with --release --ignored"]
fn synthesized_lpf_verifies_against_its_design_on_the_engine() {
    let dut = dut_layout();
    let reference = reference_layout(&dut);

    let (dut_spec, dt) = job_for(&dut);
    let (ref_spec, dt2) = job_for(&reference);
    assert_eq!(dt, dt2, "runs must share dt");
    assert_eq!(
        (dut_spec.nx, dut_spec.ny, dut_spec.nz),
        (ref_spec.nx, ref_spec.ny, ref_spec.nz),
        "runs must share the grid"
    );

    let (ref_p2, ref_p1) = run(ref_spec);
    let (dut_p2, dut_p1) = run(dut_spec);
    assert!(ref_p2.iter().any(|v| *v != 0.0), "reference probe silent");
    assert!(dut_p2.iter().any(|v| *v != 0.0), "dut probe silent");

    // 0.8–4.2 GHz, 50 MHz raster.
    let freqs: Vec<f64> = (0..=68).map(|n| 0.8e9 + n as f64 * 50.0e6).collect();
    let s21_db = sparams::transmission_db(&dut_p2, &ref_p2, dt, &freqs);
    let s11_db = sparams::reflection_db(&dut_p1, &ref_p1, dt, &freqs);

    // Measured −3 dB cutoff: the lowest frequency from which |S21| stays
    // below −3.01 dB for at least 5 consecutive samples (250 MHz), so a
    // single ripple dip does not fake the crossing.
    let below: Vec<bool> = s21_db.iter().map(|db| *db < -3.01).collect();
    let f_3db = freqs
        .iter()
        .enumerate()
        .find(|(n, _)| below[*n..].iter().take(5).filter(|b| **b).count() == 5.min(below.len() - n))
        .map(|(_, f)| *f)
        .expect("|S21| never sustains a −3 dB crossing — no cutoff in band");

    let passband_db = band_mean_db(&freqs, &s21_db, 1.0e9, 0.1e9);
    let stopband_db = band_mean_db(&freqs, &s21_db, 4.0e9, 0.1e9);
    let s11_passband_db = band_mean_db(&freqs, &s11_db, 1.0e9, 0.1e9);

    // Record: measured vs the closed-form design response.
    let table_freqs = [0.8e9, 1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9, 4.0e9];
    let ideal = ideal_response_lowpass(Approximation::Butterworth, ORDER, F_C_HZ, &table_freqs);
    eprintln!(
        "engine-filter-verify-001: N={ORDER} Butterworth stepped-impedance LPF, designed f_c = {:.1} GHz",
        F_C_HZ / 1e9
    );
    for (i, f) in table_freqs.iter().enumerate() {
        let measured = band_mean_db(&freqs, &s21_db, *f, 0.026e9);
        eprintln!(
            "  {:>4.2} GHz: measured {:>7.2} dB | ideal {:>7.2} dB",
            f / 1e9,
            measured,
            20.0 * ideal[i].norm().log10()
        );
    }
    let rel_err = (f_3db - F_C_HZ).abs() / F_C_HZ;
    eprintln!(
        "  cutoff: measured −3 dB at {:.3} GHz vs designed {:.1} GHz → err {:.1} % \
         | passband {:.2} dB @1 GHz (S11 {:.1} dB) | stopband {:.2} dB @4 GHz",
        f_3db / 1e9,
        F_C_HZ / 1e9,
        rel_err * 100.0,
        passband_db,
        s11_passband_db,
        stopband_db,
    );

    assert!(
        rel_err <= 0.20,
        "engine-filter-verify-001 FAILED: −3 dB cutoff at {:.3} GHz, designed {:.1} GHz \
         (err {:.1} % > 20 %)",
        f_3db / 1e9,
        F_C_HZ / 1e9,
        rel_err * 100.0
    );
    // Relative rejection is ripple-robust: whatever the standing waves do
    // to absolute levels, the designed filter must separate its passband
    // from its stopband by a wide margin (ideal N=5 Butterworth: 30.1 dB
    // between 1 GHz and 4 GHz).
    let rejection_db = passband_db - stopband_db;
    assert!(
        rejection_db >= 20.0,
        "engine-filter-verify-001 FAILED: passband−stopband rejection only \
         {rejection_db:.1} dB (need ≥ 20 dB)"
    );
    // Absolute passband level: with CPML-xy walls (S.9) + aperture ports
    // (S.10) the measurement is clean enough to bound absolutely — the
    // lossless through path must sit within ±3 dB of 0 dB (measured
    // −0.49 dB; the pre-S.10 single-cell ports read +3.4 dB and the
    // PEC-box +17.8 dB band-edge ripple is gone entirely).
    assert!(
        passband_db.abs() <= 3.0,
        "engine-filter-verify-001 FAILED: passband |S21| mean {passband_db:.2} dB @1 GHz — \
         beyond the ±3 dB aperture-port budget (S.10 regression?)"
    );
    // Passband return loss must be physical and reasonably matched
    // (measured −9.2 dB with aperture ports; the single-cell ports read a
    // non-physical +7 dB).
    assert!(
        s11_passband_db <= -6.0,
        "engine-filter-verify-001 FAILED: passband |S11| = {s11_passband_db:.1} dB @1 GHz \
         (need ≤ −6 dB return loss)"
    );
}
