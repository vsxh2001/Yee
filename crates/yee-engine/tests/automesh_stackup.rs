//! Gate `engine-automesh-stackup-001` (FS.4.2c, ADR-0227): the
//! `engine-stripline-eeff-001` symmetric-stripline fixture (ADR-0215) with
//! **no hand-set dx anywhere** — [`auto_dx_stackup`] alone seeds the grid
//! from the [`Stackup`] + drive band, and the measured ε_eff still lands
//! inside the same ≤ 2 % bar as the hand-tuned gate. The point: the
//! ADR-0215 "b ≥ 16 cells across a lidded confined mode" lesson is now a
//! rule the rulebook enforces automatically, not a fact a fixture author
//! has to already know.
//!
//! Method: identical to `stripline_eeff.rs` (resistive-port drive,
//! hard-PEC box, time-gated two-probe single-bin DFT phase advance → v_p
//! → ε_eff); the only difference is the grid comes from
//! [`auto_dx_stackup`] instead of a hand-picked `DX_M`.
//!
//! `#[ignore]`'d (multi-minute release run):
//!
//! ```bash
//! cargo test -p yee-engine --release --test automesh_stackup -- --ignored --nocapture
//! ```

use std::f64::consts::PI;

use yee_engine::automesh::{auto_dx_stackup, min_feature_m};
use yee_engine::{
    BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, PortSpec, ProbeSpec,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Stackup, Substrate};
use yee_voxel::{VoxelOptions, voxelize_stackup};

const EPS_R: f64 = 4.4;
/// Ground-to-ground spacing b (trace at b/2) — identical fixture numbers
/// to `engine-stripline-eeff-001` (ADR-0215).
const B_M: f64 = 3.2e-3;
const W_M: f64 = 1.5e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;
/// Physical CPML margin, held constant in **metres** (ADR-0204 hygiene)
/// and converted to cells at whatever dx the rulebook returns. Matches
/// the 4 mm the hand-tuned gate picked as 20 cells at its hand-set
/// dx = 0.2 mm.
const MARGIN_M: f64 = 4.0e-3;
const PORT_R_OHM: f64 = 50.0;

#[test]
#[ignore = "slow: multi-minute release FDTD; engine-automesh-stackup-001 gate (FS.4.2c) — run with --release --ignored"]
fn automesh_stackup_matches_the_exact_tem_value() {
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
    let bw = 0.4 * F0_HZ;
    // Drive-band ceiling for the wavelength rule: centre frequency plus
    // the Gaussian pulse's half-bandwidth (the same "f_max from the
    // drive" idiom other automesh gates use, e.g. board_automesh.rs).
    let f_max_hz = F0_HZ + bw / 2.0;

    // No hand-set dx anywhere: the rulebook alone seeds the grid.
    let dx_seed = auto_dx_stackup(&layout, &stack, f_max_hz);

    // Which rule binds — computed here, not assumed, so a regression that
    // silently loosens a rule shows up as an assertion failure naming the
    // wrong rule, not a quietly-passing coincidence.
    let eps_r_max = stack.layers.iter().fold(0.0_f64, |m, l| m.max(l.eps_r));
    let lambda_min = C0_M_S / (f_max_hz * eps_r_max.sqrt());
    let by_wavelength = lambda_min / 20.0;
    let by_layer = stack
        .layers
        .iter()
        .fold(f64::INFINITY, |m, l| m.min(l.height_m / 3.0));
    let by_feature = min_feature_m(&layout) / 2.0;
    let by_lid = stack.total_height_m() / 16.0;
    let terms = [
        ("wavelength lambda/20", by_wavelength),
        ("per-layer h/3", by_layer),
        ("feature/2", by_feature),
        ("lid b/16 (ADR-0215)", by_lid),
    ];
    let &(binding_name, binding_val) = terms
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("terms is non-empty");
    eprintln!(
        "engine-automesh-stackup-001: dx = {:.4} mm, binding rule = {binding_name} \
         ({:.4} mm) [wavelength {:.4}, h/3 {:.4}, feature/2 {:.4}, lid b/16 {:.4}]",
        dx_seed * 1e3,
        binding_val * 1e3,
        by_wavelength * 1e3,
        by_layer * 1e3,
        by_feature * 1e3,
        by_lid * 1e3
    );
    assert_eq!(
        binding_name,
        "lid b/16 (ADR-0215)",
        "engine-automesh-stackup-001 FAILED: expected the ADR-0215 lid rule (b/16) to \
         bind for the symmetric-stripline fixture (its confined transverse mode is the \
         reason that rule exists) but '{binding_name}' bound instead at {:.4} mm -- a \
         different binding rule here means the rulebook is no longer exercising the \
         b >= 16-cells-across-a-confined-mode lesson this gate exists to check",
        binding_val * 1e3
    );
    assert!(
        (dx_seed - by_lid).abs() < 1e-12,
        "auto_dx_stackup returned {dx_seed} but the recomputed lid term is {by_lid}"
    );

    let margin_cells = (MARGIN_M / dx_seed).round() as usize;
    let model = voxelize_stackup(
        &layout,
        &stack,
        0,
        &VoxelOptions {
            dx_m: dx_seed,
            xy_margin_cells: margin_cells,
            air_above_cells: 0, // lidded: ignored
        },
    );
    let (nx, ny, nz) = model.dims;
    let (_i_drive, j_strip, k_trace) = model.port_cells[0];
    let dt = model.grid.dt;
    let dx = model.dx_m;
    eprintln!(
        "  grid {nx}x{ny}x{nz}, trace at k = {k_trace} (b = {} cells), L = {:.1} mm",
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
    let x0 = layout.bbox.min.x - margin_cells as f64 * dx;
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

    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (PI * bw)) / dt).ceil() as usize;

    let spec = JobSpec {
        nx,
        ny,
        nz,
        dx_m: dx,
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
    // Measured 2026-07-24: 0.065 % at the rulebook's dx = 0.200 mm (grid
    // 1184x48x16, 7027 steps, ~23 s release) — bit-for-bit the hand-tuned
    // engine-stripline-eeff-001 gate's own numbers, this time with no
    // fixture author needing to already know the b >= 16 cells lesson.
    // Pinned at the same 2 % bar as that gate (30x headroom while staying
    // an order under the b = 8 cells failure mode ADR-0215 measured,
    // +14.3 %).
    assert!(
        rel_err <= 0.02,
        "engine-automesh-stackup-001 FAILED: eps_eff = {eps_eff_fdtd:.4} vs exact {eps_eff_ref} \
         (err {:.3} % > 2 %)",
        rel_err * 100.0
    );
}
