//! Gate `engine-stripline-eeff-001` (FS.4.0, ADR-0215): a symmetric
//! stripline built by the multilayer stackup voxelizer propagates at the
//! **exact TEM velocity** — ε_eff = ε_r, no approximation involved
//! (homogeneous dielectric between two ground planes), which makes this
//! the buried-line analog of S.5's microstrip ε_eff gate with a
//! *stronger* reference than Hammerstad–Jensen.
//!
//! Method: identical to `verify_line_eeff.rs` (S.5) — resistive-port
//! drive, hard-PEC box, time-gated two-probe single-bin DFT phase
//! advance → v_p → ε_eff.
//!
//! `#[ignore]`'d (multi-minute release run):
//!
//! ```bash
//! cargo test -p yee-engine --release --test stripline_eeff -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::{
    BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, PortSpec, ProbeSpec,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, Substrate};
use yee_voxel::{VoxelOptions, voxelize_stackup};

const EPS_R: f64 = 4.4;
/// Ground-to-ground spacing b (trace at b/2).
const B_M: f64 = 3.2e-3;
// Box-mode hygiene: the lidded, homogeneously-filled box is itself a
// waveguide. Its TE10 cutoff across the box width w = W + 2*margin is
// f_c = c/(2 w sqrt(eps_r)); margin = 10 cells (4 mm) puts f_c at
// ~7.5 GHz, above the drive band (5 +/- 2 GHz half-width via bw = 0.4 f0),
// so every box mode is evanescent and the two-probe phase fit sees the
// TEM line mode only. (First attempt: margin 26 -> w = 22.3 mm ->
// f_c = 3.2 GHz PROPAGATING; measured eps_eff 5.007 vs 4.4 = 13.8 % high
// from mode mixing. Lateral field decay exp(-pi x/b) at 4 mm: ~2 %.)
// Window hygiene: the narrow band's LONGER pulse (t0 ~ 950 steps at
// bw = 2 GHz) must clear probe B before the far-wall-reflection gate;
// L = 8 lambda_g gives gate ~ 3400 steps vs the pulse tail at ~2650.
// (Second attempt at L = 6 lambda_g: gate 2375 clipped the tail ->
// 14.5 % high. Both hygiene bounds must hold at once.)
// Resolution hygiene (the FS.4.0 measured lesson, ADR-0215): the
// stripline mode is CONFINED — transverse decay scale b/pi ~ 1 mm — and
// under-resolving it breaks the discrete TEM Laplacian cancellation,
// inflating beta. Measured: b = 8 cells (dx = 0.4 mm) -> eps_eff 5.03,
// +14.3 % (probe-separation doubling scaled the phase excess 2x: a real
// beta error, not an artifact); b = 16 cells (dx = 0.2 mm) -> 0.065 %.
// Confined lidded modes need b >= ~16 cells.
const W_M: f64 = 1.5e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
const DX_M: f64 = 0.2e-3;
const MARGIN_CELLS: usize = 20;
const PORT_R_OHM: f64 = 50.0;

#[test]
#[ignore = "slow: multi-minute release FDTD; engine-stripline-eeff-001 gate (FS.4.0) — run with --release --ignored"]
fn stripline_eeff_matches_the_exact_tem_value() {
    // TEM in homogeneous dielectric: exact, not an approximation.
    let eps_eff_ref = EPS_R;

    let lam_g = C0_M_S / (F0_HZ * eps_eff_ref.sqrt());
    let l_m = 8.0 * lam_g;
    let traces = vec![Polygon::rect(0.0, 0.0, l_m, W_M)];
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: Substrate {
            // Unused by the stackup path; kept for the Layout contract.
            eps_r: EPS_R,
            height_m: B_M,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![PortRef {
            at: Point2::new(0.5e-3, W_M / 2.0),
            width_m: W_M,
            ref_impedance_ohm: PORT_R_OHM,
        }],
        bbox,
    };

    let stack = Stackup::symmetric_stripline(EPS_R, B_M);
    let model = voxelize_stackup(
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
    let (_i_drive, j_strip, k_trace) = model.port_cells[0];
    let dt = model.grid.dt;
    let dx = model.dx_m;
    eprintln!(
        "engine-stripline-eeff-001: grid {nx}x{ny}x{nz}, trace at k = {k_trace} \
         (b = {} cells), L = {:.1} mm",
        nz,
        l_m * 1e3
    );

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

    // Probes: two planes λ_g/3 apart, past the launch transient.
    let k_probe = k_trace.saturating_sub(1).max(1);
    let x_a = 2.5 * lam_g;
    let x_b = x_a + lam_g / 3.0;
    let x0 = layout.bbox.min.x - MARGIN_CELLS as f64 * dx;
    let i_for = |xp: f64| (((xp - x0) / dx).round() as isize).clamp(0, nx as isize - 1) as usize;
    let (i_a, i_b) = (i_for(x_a), i_for(x_b));
    assert!(i_b > i_a, "probe planes collapsed");
    let delta_x = (i_b - i_a) as f64 * dx;

    // Time gate: stop before the far-wall reflection reaches probe B.
    let v_p_ref = C0_M_S / eps_eff_ref.sqrt();
    let x_drive = 0.5e-3;
    let t_refl_b = ((l_m - x_drive) + (l_m - x_b)) / v_p_ref;
    let gate_steps = (0.9 * t_refl_b / dt) as usize;
    let n_steps = gate_steps + 200;

    let bw = 0.4 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: DX_M,
        n_steps,
        boundary: BoundarySpec::Pec,
        sources: vec![],
        ports: vec![PortSpec {
            cell: model.port_cells[0],
            resistance_ohm: PORT_R_OHM,
            v0: 1.0,
            f0_hz: F0_HZ,
            bw_hz: bw,
            t0_steps,
        }],
        aperture_ports: vec![],
        probes: vec![
            ProbeSpec {
                component: "ez".into(),
                cell: (i_a, j_strip, k_probe),
            },
            ProbeSpec {
                component: "ez".into(),
                cell: (i_b, j_strip, k_probe),
            },
        ],
        slice: None,
        ntff: None,
        materials: Some(materials),
        dt_s: Some(dt),
        spacings: None,
        backend: BackendChoice::Cpu,
    };

    let handle = yee_engine::submit(spec);
    let mut result = None;
    for event in handle.events() {
        match event {
            JobEvent::Done { result: r } => result = Some(r),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => {}
        }
    }
    let result = result.expect("no Done event");

    // Time-gated single-bin DFT at f0, phase advance A → B.
    let omega = 2.0 * PI * F0_HZ;
    let series = &result.probes;
    let mut acc = [0.0_f64; 4];
    let gate = gate_steps.min(n_steps);
    for (n, (a, b)) in series[0][..gate].iter().zip(&series[1][..gate]).enumerate() {
        let phase = omega * n as f64 * dt;
        let (c, s) = (phase.cos(), phase.sin());
        acc[0] += a * c;
        acc[1] -= a * s;
        acc[2] += b * c;
        acc[3] -= b * s;
    }
    let phi_a = acc[1].atan2(acc[0]);
    let phi_b = acc[3].atan2(acc[2]);
    let mut delta_phi = phi_a - phi_b;
    while delta_phi <= 0.0 {
        delta_phi += 2.0 * PI;
    }
    while delta_phi > 2.0 * PI {
        delta_phi -= 2.0 * PI;
    }
    let v_p = omega * delta_x / delta_phi;
    let eps_eff_fdtd = (C0_M_S / v_p).powi(2);
    let rel_err = (eps_eff_fdtd - eps_eff_ref).abs() / eps_eff_ref;

    eprintln!(
        "  eps_eff = {eps_eff_fdtd:.4} vs exact TEM {eps_eff_ref:.4} -> err {:.3} % \
         (dphi = {delta_phi:.4} rad over {:.3} mm, {n_steps} steps, gate {gate_steps})",
        rel_err * 100.0,
        delta_x * 1e3
    );
    // Measured 2026-07-12: 0.065 % at b = 16 cells (grid 1184x48x16,
    // 7027 steps, ~40 s release). Pinned at 2 % — 30x headroom while
    // staying an order under the coarse-resolution failure mode.
    assert!(
        rel_err <= 0.02,
        "engine-stripline-eeff-001 FAILED: eps_eff = {eps_eff_fdtd:.4} vs exact {eps_eff_ref} \
         (err {:.3} % > 2 %)",
        rel_err * 100.0
    );
}
