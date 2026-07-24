//! Gate `engine-stackup-via-001` (FS.4.1, ADR-0221): a **through-via in
//! a multilayer stackup** demonstrably changes the circuit the way
//! transmission-line theory says it must — the R.1 `engine-via-001`
//! differential fixture moved onto the FS.4.0 symmetric stripline.
//!
//! Differential scenario on ONE grid (three runs):
//!
//! - **Control** (open stub): a mid-line stub that is λ/4 at ~5 GHz →
//!   input short → deep |S21| notch there.
//! - **DUT** (same stub + a `with_through_via_at_cell` barrel at its far
//!   end, ground → trace → lid): a **shorted** λ/4 stub is an input
//!   *open* → the notch must vanish at that frequency (the shorted
//!   stub's own notch is the λ/2 condition near 2 f₀, outside the band).
//! - **Reference** (bare line), shared by both transmission ratios
//!   (`sparams::transmission_db` — notch-location/depth-shaped asserts,
//!   the ADR-0204-sanctioned use of the single ratio).
//!
//! Grid/window hygiene, mapped from `stripline_eeff.rs` (ADR-0215):
//!
//! - **Resolution**: the stripline mode is confined (transverse decay
//!   scale b/π); lidded confined modes need ≥ ~16 cells across the
//!   ground-to-ground gap b. Here b = 3.2 mm at dx = 0.2 mm = 16 cells.
//! - **Box modes**: the ADR-0215 PEC box turned into a lateral waveguide
//!   whose TE₁₀ cutoff had to clear the band. This fixture instead uses
//!   **CPML side walls** (axes x + y; z stays PEC = ground and lid), so
//!   there is no lateral cavity to resonate: parallel-plate waves
//!   launched sideways by the tee/via are absorbed. The stub end sits
//!   14 working cells (2.8 mm) inside the CPML — lateral evanescent
//!   fringing exp(−πx/b) is down to ~6 % there.
//! - **Windowing**: no time gate — absorbing terminations (resistive
//!   aperture ports at both ends + CPML beyond) replace the PEC-box
//!   pulse-tail rule; N_STEPS covers launch + transit + stub ring-down
//!   (12000 steps ≈ 4.6 ns, ~23 f₀ periods after the pulse peak).
//!
//! `#[ignore]`'d (three multi-minute release FDTD runs):
//!
//! ```bash
//! cargo test -p yee-engine --release --test stackup_via -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    sparams,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, Substrate};
use yee_voxel::{VoxelOptions, voxelize_stackup, with_through_via_at_cell};

const EPS_R: f64 = 4.4;
/// Ground-to-ground spacing b (trace at b/2) — 16 cells at DX (ADR-0215).
const B_M: f64 = 3.2e-3;
const W_M: f64 = 1.5e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.2e-3;
/// 10 CPML cells + 14 working cells between the metal and the absorber.
const MARGIN_CELLS: usize = 24;
const NPML: usize = 10;
const Z0_OHM: f64 = 50.0;
const N_STEPS: usize = 12_000;

/// Stripline open-end length correction Δl ≈ b·ln2/π — pre-compensates
/// the stub so the open-stub notch lands near F0 (the measured notch
/// pins the asserts either way; this only centres it in the scan band).
fn stripline_open_end_delta_l(b_m: f64) -> f64 {
    b_m * std::f64::consts::LN_2 / PI
}

/// Build one run on the shared grid: bare line (reference), open stub
/// (control), or stub + through-via (DUT).
fn job(with_stub: bool, with_via: bool) -> (JobSpec, f64) {
    // Exact TEM: ε_eff = ε_r (homogeneous fill between the plates).
    let lam_g = C0_M_S / (F0_HZ * EPS_R.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - stripline_open_end_delta_l(B_M);

    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    // Shared bbox → identical grid for all three variants.
    let bbox = BBox::from_polygons(&[line.clone(), stub.clone()]);
    let traces = if with_stub {
        vec![line, stub]
    } else {
        vec![line]
    };
    let layout = Layout {
        substrate: Substrate {
            // Unused by the stackup path; kept for the Layout contract.
            eps_r: EPS_R,
            height_m: B_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
        ],
        bbox,
    };

    let stack = Stackup::symmetric_stripline(EPS_R, B_M);
    let mut model = voxelize_stackup(
        &layout,
        &stack,
        0,
        &VoxelOptions {
            dx_m: DX_M,
            xy_margin_cells: MARGIN_CELLS,
            air_above_cells: 0, // lidded: ignored
        },
    );
    let (nx, ny, nz) = model.dims;
    let dt = model.grid.dt;
    let dx = model.dx_m;
    let (_i_drive, j_strip, k_trace) = model.port_cells[0];
    let load_cell = model.port_cells[1];
    let k_probe = k_trace - 1;

    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let y0 = layout.bbox.min.y - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| ((xp - x0) / dx).round().clamp(0.0, nx as f64 - 1.0) as usize;
    let j_for = |yp: f64| ((yp - y0) / dx).round().clamp(0.0, ny as f64 - 1.0) as usize;

    if with_via {
        // Through-via at the stub's far-end centre: ground → trace → lid
        // barrel; the stub end is shorted to the bottom ground (and lid).
        let i_via = i_for(l_m / 2.0);
        let j_via = j_for(W_M + stub_len - 0.5e-3);
        with_through_via_at_cell(&mut model, i_via, j_via);
    }

    // Probe past the stub, λ_g/2-ish before the load port.
    let i_m = i_for(l_m - 3.0e-3);

    let in_band = |j: usize| -> bool { (y0 + (j as f64 + 0.5) * dx - W_M / 2.0).abs() < W_M / 2.0 };
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
        pec_mask_ez: model
            .grid
            .pec_mask_ez
            .as_ref()
            .map(|a| a.as_slice().unwrap().to_vec()),
        ..MaterialsSpec::default()
    };

    let bw = 0.8 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: NPML,
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
                k_top: k_trace,
                resistance_ohm: Z0_OHM,
                v0: 1.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
            AperturePortSpec {
                i: load_cell.0,
                j_lo,
                j_hi,
                k_lo: 0,
                k_top: k_trace,
                resistance_ohm: Z0_OHM,
                v0: 0.0,
                f0_hz: F0_HZ,
                bw_hz: bw,
                t0_steps,
                record: false,
            },
        ],
        thin_wires: vec![],
        probes: vec![ProbeSpec {
            component: "ez".into(),
            cell: (i_m, j_strip, k_probe),
        }],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
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
    result.probes.into_iter().next().expect("no probe")
}

#[test]
#[ignore = "slow: three multi-minute release FDTD runs; engine-stackup-via-001 gate (FS.4.1) — run with --release --ignored"]
fn through_via_shorted_stripline_stub_removes_the_open_stub_notch() {
    let (ref_spec, dt) = job(false, false);
    let (open_spec, _) = job(true, false);
    let (via_spec, _) = job(true, true);
    eprintln!(
        "engine-stackup-via-001: grid {}x{}x{} (b = 16 cells), {} steps x 3 runs",
        ref_spec.nx, ref_spec.ny, ref_spec.nz, N_STEPS
    );

    let reference = run(ref_spec);
    let open_stub = run(open_spec);
    let via_stub = run(via_spec);

    // Scan 4.0–5.6 GHz: wide enough to catch the open-end-shifted notch,
    // and the shorted stub's λ/2 notch (~2 f₀) stays far outside.
    let freqs: Vec<f64> = (0..=64).map(|n| 4.0e9 + n as f64 * 25.0e6).collect();
    let s21_open = sparams::transmission_db(&open_stub, &reference, dt, &freqs);
    let s21_via = sparams::transmission_db(&via_stub, &reference, dt, &freqs);

    let min_of = |v: &[f64]| {
        v.iter()
            .enumerate()
            .map(|(n, db)| (n, *db))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .expect("empty")
    };
    let (n_open, open_db) = min_of(&s21_open);
    let via_at_notch = s21_via[n_open];
    let (n_via_min, via_min_db) = min_of(&s21_via);

    eprintln!(
        "engine-stackup-via-001: open-stub notch {:.2} dB at {:.3} GHz | through-via variant \
         at that frequency: {:.2} dB (its own band min {:.2} dB at {:.3} GHz)",
        open_db,
        freqs[n_open] / 1e9,
        via_at_notch,
        via_min_db,
        freqs[n_via_min] / 1e9,
    );

    // Measured 2026-07-13 (release, CPU, grid 477x88x16, 3 x 12000
    // steps, ~5.5 min total): open-stub notch −39.81 dB at 5.075 GHz
    // (design f₀ = 5 GHz + Δl pre-compensation: 1.5 % off — the open-end
    // model is good); through-via variant at 5.075 GHz: +0.62 dB, and
    // its whole-band min is −1.18 dB at 4.575 GHz — no notch anywhere in
    // the scan (the small positive dB is the ADR-0204 single-ratio
    // launch artifact, fine for notch-shaped asserts). Pinned with
    // ~2.5× dB headroom: control ≤ −15 dB, via ≥ −3 dB everywhere.
    assert!(
        open_db <= -15.0,
        "engine-stackup-via-001 FAILED: control (open-stub) notch only {open_db:.1} dB — \
         the stripline stub physics regressed"
    );
    assert!(
        via_at_notch >= -3.0,
        "engine-stackup-via-001 FAILED: with the through-via the notch frequency still reads \
         {via_at_notch:.1} dB — the via is not shorting the stub through the stack"
    );
    assert!(
        via_min_db >= -3.0,
        "engine-stackup-via-001 FAILED: the through-via variant has a {via_min_db:.1} dB dip \
         at {:.3} GHz — the shorted stub should not notch anywhere in this band",
        freqs[n_via_min] / 1e9
    );
}
