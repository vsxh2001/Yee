//! Gate `engine-sparams-001` (S.6, ADR-0183): |S21|(f) out of engine jobs —
//! the two-run transmission method on the job protocol recovers the
//! **textbook λ/4 open-stub bandstop notch** (Pozar, *Microwave
//! Engineering*: a shunt open-circuited stub shorts the line where its
//! electrical length is a quarter wave).
//!
//! Scenario: the S.5-certified FR-4 microstrip stack (W = 3 mm,
//! h = 1.6 mm, ε_r = 4.4, dx = 0.3 mm) as a 3 λ_g feed line with a 50 Ω
//! drive port at one end and a **passive** 50 Ω port (`v0 = 0` — a lumped
//! resistor load) at the other. Two jobs over `submit()`:
//! reference = bare line, DUT = line + open stub sized
//! `L_s = λ_g/4 − ΔL` (Hammerstad open-end correction), so closed forms
//! alone predict the notch at f₀ = 5 GHz. `sparams::transmission_db`
//! divides the runs; assert the notch lands within ±15 % of the
//! prediction (the filter pipeline's walking-skeleton band), is ≥ 8 dB
//! deep, and the band-edge ripple stays bounded. S.7 (ADR-0184) adds a
//! port-1 reference-plane probe to both runs: `sparams::reflection_db`
//! separates incident/reflected, and the gate additionally asserts
//! |S11| peaks near 0 dB at the notch with a physical |S11|²+|S21|².
//!
//! `#[ignore]`'d (two multi-minute release FDTD runs on a ~1.7 M-cell grid):
//!
//! ```bash
//! cargo test -p yee-engine --release --test sparams_stub_notch -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, PortSpec, ProbeSpec, sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff};
use yee_voxel::{VoxelOptions, voxelize_microstrip};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.3e-3;
const MARGIN_CELLS: usize = 34;
const AIR_ABOVE_CELLS: usize = 34;
const PORT_R_OHM: f64 = 50.0;
const DRIVE_V0: f64 = 1.0;
const FREQ_SPAN: f64 = 0.8;
const N_STEPS: usize = 9000;

/// Hammerstad microstrip open-end length correction ΔL (the fringing field
/// makes an open stub electrically longer than its metal):
/// `ΔL = 0.412·h·(ε_eff + 0.3)(W/h + 0.264) / ((ε_eff − 0.258)(W/h + 0.8))`.
fn open_end_delta_l(w_m: f64, h_m: f64, e_eff: f64) -> f64 {
    let u = w_m / h_m;
    0.412 * h_m * ((e_eff + 0.3) * (u + 0.264)) / ((e_eff - 0.258) * (u + 0.8))
}

/// Build the JobSpec for one run: the shared feed line, plus the stub when
/// `with_stub`. Both runs use the DUT's bbox so the grids are identical.
fn stub_job(with_stub: bool) -> (JobSpec, f64) {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);

    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    // The bbox (and therefore the voxel grid) comes from the FULL DUT
    // geometry in both runs, so reference and DUT share dims, origin, and
    // every port/probe cell — only the stub metal differs.
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
                ref_impedance_ohm: PORT_R_OHM,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: PORT_R_OHM,
            },
        ],
        bbox,
    };

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
    let load_cell = model.port_cells[1];
    let k_probe = k_top.saturating_sub(1).max(1);

    // Transmission probe under the trace, 3 mm before the load port; S11
    // probe (S.7) at the port-1 reference plane, between the launch and
    // the stub.
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let i_m = i_for(l_m - 3.0e-3);
    let i_p1 = i_for(12.0e-3);

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

    let bw = FREQ_SPAN * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Pec,
        sources: vec![],
        ports: vec![
            PortSpec {
                cell: model.port_cells[0],
                resistance_ohm: PORT_R_OHM,
                v0: DRIVE_V0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
            },
            // Passive matched load: zero EMF leaves the pure-resistor arm
            // of the lumped-port update — the far-end termination.
            PortSpec {
                cell: load_cell,
                resistance_ohm: PORT_R_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
            },
        ],
        aperture_ports: vec![],
        probes: vec![
            ProbeSpec {
                component: "ez".into(),
                cell: (i_m, j_strip, k_probe),
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

/// Run one job; returns `(transmission_probe, p1_probe)` series.
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

#[test]
#[ignore = "slow: two multi-minute release FDTD runs; engine-sparams-001 gate (S.6) — run with --release --ignored"]
fn open_stub_notch_matches_transmission_line_theory() {
    let (ref_spec, dt) = stub_job(false);
    let (dut_spec, dt2) = stub_job(true);
    assert_eq!(dt, dt2, "runs must share dt");
    assert_eq!(
        (ref_spec.nx, ref_spec.ny, ref_spec.nz),
        (dut_spec.nx, dut_spec.ny, dut_spec.nz),
        "runs must share the grid"
    );

    let (reference, ref_p1) = run(ref_spec);
    let (dut, dut_p1) = run(dut_spec);
    assert!(
        reference.iter().any(|v| *v != 0.0),
        "reference probe silent"
    );
    assert!(dut.iter().any(|v| *v != 0.0), "dut probe silent");

    // |S21|(f) and |S11|(f) over the drive band, 50 MHz raster.
    let freqs: Vec<f64> = (0..=80).map(|n| 3.0e9 + n as f64 * 50.0e6).collect();
    let s21_db = sparams::transmission_db(&dut, &reference, dt, &freqs);
    let s11_db = sparams::reflection_db(&dut_p1, &ref_p1, dt, &freqs);

    // Deepest point of the notch, searched inside 3.5–6.5 GHz.
    let (n_min, db_min) = freqs
        .iter()
        .zip(&s21_db)
        .enumerate()
        .filter(|(_, (f, _))| (3.5e9..=6.5e9).contains(*f))
        .map(|(n, (_, db))| (n, *db))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("no in-band samples");
    let f_notch = freqs[n_min];

    // Passband sanity: the band edges must stay within ripple range in BOTH
    // directions — a genuine narrow dip, not broadband loss and not a
    // pathological ratio. Out-of-band the stub still partially reflects,
    // and those reflections exist only in the DUT run, so the single-probe
    // ratio carries standing-wave ripple of either sign (measured
    // +8.7 dB @3 GHz / +5.2 dB @7 GHz on the shipped scenario); the bound
    // is |edge| ≤ 12 dB.
    let edge_db = s21_db[0].abs().max(s21_db[s21_db.len() - 1].abs());

    // S.7: at resonance the stub reflects nearly everything — |S11| must
    // peak near 0 dB at the notch, and the lossless-DUT energy balance
    // |S11|² + |S21|² must stay physical (within the known band-edge
    // ripple budget).
    let s11_notch_db = s11_db[n_min];
    let energy_notch = 10.0_f64.powf(s11_notch_db / 10.0) + 10.0_f64.powf(s21_db[n_min] / 10.0);

    let rel_err = (f_notch - F0_HZ).abs() / F0_HZ;
    eprintln!(
        "engine-sparams-001 open-stub notch over the job protocol: \
         notch {:.3} GHz ({:.1} dB) vs TL-theory {:.1} GHz → err {:.2} % \
         | band edges: {:.1} dB @3 GHz, {:.1} dB @7 GHz \
         | S11 at notch: {:.2} dB, |S11|²+|S21|² = {:.3}",
        f_notch / 1e9,
        db_min,
        F0_HZ / 1e9,
        rel_err * 100.0,
        s21_db[0],
        s21_db[s21_db.len() - 1],
        s11_notch_db,
        energy_notch,
    );

    assert!(
        rel_err <= 0.15,
        "engine-sparams-001 FAILED: notch at {:.3} GHz, predicted {:.1} GHz (err {:.2} % > 15 %)",
        f_notch / 1e9,
        F0_HZ / 1e9,
        rel_err * 100.0
    );
    assert!(
        db_min <= -8.0,
        "engine-sparams-001 FAILED: notch only {db_min:.1} dB deep (need ≤ −8 dB)"
    );
    assert!(
        edge_db <= 12.0,
        "engine-sparams-001 FAILED: band-edge |ripple| {edge_db:.1} dB > 12 dB — not a clean notch measurement"
    );
    assert!(
        s11_notch_db >= -4.0,
        "engine-sparams-001 FAILED: |S11| at the notch only {s11_notch_db:.2} dB — \
         a λ/4 open stub must reflect nearly everything at resonance"
    );
    assert!(
        (0.5..=1.3).contains(&energy_notch),
        "engine-sparams-001 FAILED: |S11|²+|S21|² = {energy_notch:.3} at the notch — \
         outside the physical band for a lossless DUT (+ ripple budget)"
    );
}
