//! Gate `engine-graded-001` (FS.0b.1, ADR-0210): the FS.0a S.6 λ/4
//! open-stub notch board solved ONCE on an [`auto_spacings`] graded grid
//! must reproduce the uniform-converged answer — notch within ±2 % of the
//! ADR-0204-measured 4.850 GHz at ≤ −20 dB — at a measured cell-count
//! reduction vs the uniform dx = 0.267 mm convergence pass. This is the
//! FS.0b payoff assertion: refine only the staircase-limited feature
//! regions, not everywhere.
//!
//! The JobSpec is built directly on the graded voxelizer (the uniform
//! `board.rs` fixture is not rewired — that integration is FS.0b.2). Both
//! runs (DUT + through-line reference) share the DUT-derived grid — the
//! ADR-0204 same-physical-problem lesson — and the measurement is the
//! launch-normalized double ratio (`sparams::forward_transfer`), never
//! the single ratio (ADR-0204).
//!
//! **Measured 2026-07-11 (ADR-0210):** grid 282×110×41 = 1,271,820 cells
//! (coarse 0.533 mm, fine 0.267 mm, guard 1.6 mm, k_top 6) vs the uniform
//! pass-2 506×176×75 = 6,679,200 — **ratio 0.190**; notch **4.900 GHz @
//! −37.2 dB** (err 1.03 % of the uniform-converged 4.850 GHz); both
//! solves 306 s release on 4 cores.
//!
//! `#[ignore]`'d (2 release FDTD solves, ~1.3 M cells × ~10 k steps):
//!
//! ```bash
//! cargo test -p yee-engine --release --test engine_graded_notch -- --ignored --nocapture
//! ```

use std::time::Instant;

use yee_engine::automesh::{GradedMeshOptions, auto_dx, auto_spacings};
use yee_engine::board::{TwoPortBoardOptions, reference_through_line, two_port_board_job};
use yee_engine::sparams;
use yee_engine::{
    AperturePortSpec, BackendChoice, BoundarySpec, JobEvent, JobSpec, MaterialsSpec, ProbeSpec,
    submit,
};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff, open_end_delta_l};
use yee_voxel::{GradedMicrostripModel, GradedVoxelGrid, voxelize_microstrip_graded};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const F_MAX_HZ: f64 = 6.0e9;
const C0_M_S: f64 = 299_792_458.0;
const Z0_OHM: f64 = 50.0;
/// The uniform-converged notch (ADR-0204, dx = 0.267 mm pass: 4.850 GHz
/// at −34.2 dB) — the answer the graded grid must reproduce.
const F_NOTCH_UNIFORM_HZ: f64 = 4.85e9;

/// The S.6 scenario: a through line with a λ/4 open stub at mid-length
/// (byte-for-byte the `board_automesh.rs` fixture).
fn stub_layout() -> Layout {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);
    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    let traces = vec![line, stub];
    let bbox = BBox::from_polygons(&traces);
    Layout {
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
                ref_impedance_ohm: Z0_OHM,
            },
            PortRef {
                at: Point2::new(l_m - 0.5e-3, W_M / 2.0),
                width_m: W_M,
                ref_impedance_ohm: Z0_OHM,
            },
        ],
        bbox,
    }
}

/// First index `i` such that a full `2·sp`-cell probe span starting at
/// node `i` lies on bit-equal-coarse cells with `nodes[i] ≥ x_min`.
fn coarse_run_at(
    nodes: &[f64],
    widths: &[f64],
    coarse: f64,
    sp: usize,
    x_min: f64,
    from_right_of: f64,
) -> usize {
    (0..widths.len() - 2 * sp)
        .find(|&i| {
            nodes[i] >= x_min
                && nodes[i + 2 * sp] <= from_right_of
                && widths[i..i + 2 * sp].iter().all(|d| *d == coarse)
        })
        .expect("no uniform-coarse probe stretch found")
}

/// Run one graded job to completion; returns (probes, dt).
fn run(spec: JobSpec) -> (Vec<Vec<f64>>, f64) {
    let handle = submit(spec);
    for event in handle.events() {
        match event {
            JobEvent::Done { result } => return (result.probes, result.dt_s),
            JobEvent::Error { message } => panic!("graded job failed: {message}"),
            _ => {}
        }
    }
    panic!("engine stream ended without a result");
}

#[test]
#[ignore = "slow: 2 release FDTD solves (~1.4M cells x ~10k steps); engine-graded-001 gate (FS.0b.1) — run with --release --ignored"]
fn graded_rules_reproduce_the_converged_notch_at_a_fraction_of_the_cells() {
    let layout = stub_layout();
    let reference = reference_through_line(&layout);

    // --- The graded rulebook, no hand-set spacing anywhere. -------------
    let opts = GradedMeshOptions::for_board(&layout, F_MAX_HZ);
    let spac = auto_spacings(&layout, F_MAX_HZ, &opts).expect("auto_spacings failed");
    let (nx, ny, nz) = (spac.dx.len(), spac.dy.len(), spac.dz.len());
    eprintln!(
        "engine-graded-001: coarse = {:.4} mm, fine = {:.4} mm, guard = {:.2} mm, \
         grid {nx}x{ny}x{nz} = {} cells (k_top = {})",
        spac.coarse_m * 1e3,
        spac.fine_m * 1e3,
        opts.guard_m * 1e3,
        spac.cell_count(),
        spac.k_top,
    );

    // --- Both layouts voxelized on the SAME (DUT-derived) grid. ---------
    let grid = GradedVoxelGrid {
        dx_m: spac.dx.clone(),
        dy_m: spac.dy.clone(),
        dz_m: spac.dz.clone(),
        x0_m: spac.x0_m,
        y0_m: spac.y0_m,
        k_gnd: spac.k_gnd,
        k_top: spac.k_top,
    };
    let dut_model = voxelize_microstrip_graded(&layout, &grid);
    let ref_model = voxelize_microstrip_graded(&reference, &grid);
    assert_eq!(dut_model.dims, ref_model.dims);

    // --- Fixture geometry (all sizes in METRES, ADR-0204 hygiene). ------
    // Aperture j-band: the feed width centred on the port height, from
    // the true graded cell centres.
    let tap_y = layout.ports[0].at.y;
    let yc = |j: usize| (dut_model.y_nodes_m[j] + dut_model.y_nodes_m[j + 1]) / 2.0;
    let in_band = |j: usize| (yc(j) - tap_y).abs() < W_M / 2.0;
    let j_lo = (0..ny).find(|&j| in_band(j)).expect("empty feed band");
    let j_hi = (j_lo..ny).find(|&j| !in_band(j)).unwrap_or(ny);
    let j_strip = (j_lo + j_hi) / 2;
    let k_top = spac.k_top;
    let k_probe = k_top - 1;

    // Probe triples on uniform-coarse stretches (fit_standing_wave needs
    // equal spacing), 12 coarse cells apart = the FS.0a 6.4 mm.
    let sp = 12usize;
    let spacing_m = sp as f64 * spac.coarse_m;
    let clearance = 2.4e-3;
    let i_a0 = coarse_run_at(
        &dut_model.x_nodes_m,
        &spac.dx,
        spac.coarse_m,
        sp,
        layout.ports[0].at.x + clearance,
        layout.bbox.max.x,
    );
    let i_b0 = {
        // Last coarse run whose far probe stays clear of the output port.
        let limit = layout.ports[1].at.x - clearance;
        (0..nx - 2 * sp)
            .rev()
            .find(|&i| {
                dut_model.x_nodes_m[i + 2 * sp] <= limit
                    && spac.dx[i..i + 2 * sp].iter().all(|d| *d == spac.coarse_m)
            })
            .expect("no coarse stretch for triple B")
    };
    assert!(i_b0 > i_a0 + 2 * sp, "probe triples overlap");

    // Time base: the FS.0a physical window at the fine spacing (the
    // graded Courant dt equals the uniform-at-fine dt: every axis
    // minimum is the fine spacing).
    let min_d = |a: &[f64]| a.iter().copied().fold(f64::INFINITY, f64::min);
    let (mx, my, mz) = (min_d(&spac.dx), min_d(&spac.dy), min_d(&spac.dz));
    let dt = 0.9 / (C0_M_S * (1.0 / (mx * mx) + 1.0 / (my * my) + 1.0 / (mz * mz)).sqrt());
    let n_steps = (9000.0 * 0.3e-3 / mx).round() as usize;
    let bw_hz = 0.8 * F0_HZ;
    let t0_steps =
        ((3.5 * (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * bw_hz)) / dt)
            .ceil() as usize;

    let job_for = |model: &GradedMicrostripModel| -> JobSpec {
        let materials = MaterialsSpec {
            eps_r_cells: Some(model.eps_r_cells.as_slice().unwrap().to_vec()),
            pec_mask_ex: Some(model.pec_mask_ex.as_slice().unwrap().to_vec()),
            pec_mask_ey: Some(model.pec_mask_ey.as_slice().unwrap().to_vec()),
            ..MaterialsSpec::default()
        };
        let mk_probe = |i: usize| ProbeSpec {
            component: "ez".into(),
            cell: (i, j_strip, k_probe),
        };
        let mk_port = |i: usize, v0: f64| AperturePortSpec {
            i,
            j_lo,
            j_hi,
            k_lo: 0,
            k_top,
            resistance_ohm: Z0_OHM,
            v0,
            f0_hz: F0_HZ,
            bw_hz,
            t0_steps,
            record: false,
        };
        JobSpec {
            nx,
            ny,
            nz,
            // The nominal spacing feeding the CPML sigma_max recipe: the
            // absorbing layers are exactly coarse (ADR-0208 scope rule).
            dx_m: spac.coarse_m,
            n_steps,
            boundary: BoundarySpec::Cpml {
                npml: opts.npml,
                axes: [true, true, false],
                faces: None,
            },
            sources: vec![],
            ports: vec![],
            aperture_ports: vec![
                mk_port(model.port_cells[0].0, 1.0),
                mk_port(model.port_cells[1].0, 0.0),
            ],
            probes: vec![
                mk_probe(i_a0),
                mk_probe(i_a0 + sp),
                mk_probe(i_a0 + 2 * sp),
                mk_probe(i_b0),
                mk_probe(i_b0 + sp),
                mk_probe(i_b0 + 2 * sp),
            ],
            slice: None,
            ntff: None,
            materials: Some(materials),
            dt_s: Some(dt),
            spacings: Some(spac.to_spacings()),
            backend: BackendChoice::Cpu,
        }
    };

    eprintln!(
        "  fixture: ports i = {} / {}, j band [{j_lo}, {j_hi}), triples at \
         i = {i_a0}/{i_b0} (+{sp}/+{}), spacing {:.3} mm, {n_steps} steps, \
         dt = {:.4e} s",
        dut_model.port_cells[0].0,
        dut_model.port_cells[1].0,
        2 * sp,
        spacing_m * 1e3,
        dt,
    );

    // --- Solve (reference first, then DUT), timed. -----------------------
    let t_start = Instant::now();
    let (ref_probes, dt_ref) = run(job_for(&ref_model));
    let t_ref = t_start.elapsed();
    let (dut_probes, dt_dut) = run(job_for(&dut_model));
    let t_total = t_start.elapsed();
    assert_eq!(dt_ref, dt_dut, "runs diverged in dt");
    eprintln!(
        "  solves: reference {:.1} s, total {:.1} s",
        t_ref.as_secs_f64(),
        t_total.as_secs_f64()
    );

    // --- Launch-normalized double-ratio |S21| (ADR-0204). ----------------
    let freqs: Vec<f64> = (0..=50).map(|n| 3.5e9 + n as f64 * 50.0e6).collect();
    let transfer = |p: &[Vec<f64>]| {
        sparams::forward_transfer(
            [&p[0], &p[1], &p[2]],
            [&p[3], &p[4], &p[5]],
            dt_ref,
            spacing_m,
            &freqs,
        )
    };
    let t_dut = transfer(&dut_probes);
    let t_ref_w = transfer(&ref_probes);
    let s21_db: Vec<f64> = t_dut
        .iter()
        .zip(&t_ref_w)
        .map(|(d, r)| 20.0 * (d.0.hypot(d.1) / r.0.hypot(r.1)).log10())
        .collect();
    for (f, db) in freqs.iter().zip(&s21_db) {
        eprintln!("    {:.2} GHz: {db:8.2} dB", f / 1e9);
    }

    let (i_min, &db_min) = s21_db
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .unwrap();
    let f_notch = freqs[i_min];
    let err = (f_notch - F_NOTCH_UNIFORM_HZ).abs() / F_NOTCH_UNIFORM_HZ;

    // --- Cell budget vs the uniform converged pass (built, not solved):
    // the exact options the FS.0a loop used at its third pass
    // (dx0 -> dx0/sqrt2 -> dx0/2, cells rescaled to constant metres).
    let dx0 = auto_dx(&layout, F_MAX_HZ);
    let mut dx_u = dx0;
    dx_u /= std::f64::consts::SQRT_2;
    dx_u /= std::f64::consts::SQRT_2;
    let mut uopts = TwoPortBoardOptions::for_band(F0_HZ, bw_hz);
    uopts.dx_m = dx_u;
    uopts.n_steps = (9000.0 * 0.3e-3 / dx_u).round() as usize;
    uopts.margin_cells = (34.0 * dx0 / dx_u).round() as usize;
    uopts.air_above_cells = (34.0 * dx0 / dx_u).round() as usize;
    uopts.npml = (10.0 * dx0 / dx_u).round() as usize;
    uopts.spacing_cells = (12.0 * dx0 / dx_u).round() as usize;
    let ujob = two_port_board_job(&layout, &uopts).expect("uniform job build failed");
    let cells_uniform = ujob.spec.nx * ujob.spec.ny * ujob.spec.nz;
    let cells_graded = spac.cell_count();
    let ratio = cells_graded as f64 / cells_uniform as f64;

    eprintln!(
        "  notch {:.3} GHz @ {db_min:.1} dB vs uniform-converged {:.2} GHz -> err {:.2} %",
        f_notch / 1e9,
        F_NOTCH_UNIFORM_HZ / 1e9,
        err * 100.0
    );
    eprintln!(
        "  cells: graded {cells_graded} ({nx}x{ny}x{nz}) vs uniform pass-2 \
         {cells_uniform} ({}x{}x{}) -> ratio {ratio:.3}",
        ujob.spec.nx, ujob.spec.ny, ujob.spec.nz
    );

    assert!(
        err <= 0.02,
        "engine-graded-001 FAILED: notch {:.3} GHz vs uniform-converged {:.2} GHz \
         (err {:.2} % > 2 %)",
        f_notch / 1e9,
        F_NOTCH_UNIFORM_HZ / 1e9,
        err * 100.0
    );
    assert!(
        db_min <= -20.0,
        "engine-graded-001 FAILED: notch only {db_min:.1} dB deep (need <= -20)"
    );
    // Measured 2026-07-11 (ADR-0210): ratio 0.190 (1,271,820 graded vs
    // 6,679,200 uniform cells; notch 4.900 GHz @ −37.2 dB, err 1.03 %;
    // both solves 306 s release). Pinned at 0.25 with margin — the graded
    // grid must stay at a quarter of the uniform pass's cells.
    assert!(
        ratio < 0.25,
        "engine-graded-001 FAILED: cells_graded/cells_uniform = {ratio:.3} >= 0.25 \
         (measured 0.190) — the graded mesh lost its payoff"
    );
}
