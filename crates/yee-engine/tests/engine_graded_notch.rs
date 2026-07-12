//! Gate `engine-graded-001` (FS.0b.1, ADR-0210): the FS.0a S.6 λ/4
//! open-stub notch board solved ONCE on an [`auto_spacings`] graded grid
//! must reproduce the uniform-converged answer — notch within ±2 % of the
//! ADR-0204-measured 4.850 GHz at ≤ −20 dB — at a measured cell-count
//! reduction vs the uniform dx = 0.267 mm convergence pass. This is the
//! FS.0b payoff assertion: refine only the staircase-limited feature
//! regions, not everywhere.
//!
//! As of FS.0b.2a the fixture is the library builder
//! [`two_port_board_jobs_graded`] (both runs share the DUT-derived grid;
//! the measurement is the launch-normalized double ratio — the ADR-0204
//! lessons live in the fixture API now). This gate therefore certifies
//! the shared builder end-to-end.
//!
//! **Measured 2026-07-11 (ADR-0210):** grid 282×110×41 = 1,271,820 cells
//! (coarse 0.533 mm, fine 0.267 mm, guard 1.6 mm, k_top 6) vs the uniform
//! pass-2 506×176×75 = 6,679,200 — **ratio 0.190**; notch **4.900 GHz @
//! −37.2 dB** (err 1.03 % of the uniform-converged 4.850 GHz); both
//! solves 306 s release on 4 cores (merged-tree re-run: identical
//! numbers, 221 s).
//!
//! `#[ignore]`'d (2 release FDTD solves, ~1.3 M cells × ~10 k steps):
//!
//! ```bash
//! cargo test -p yee-engine --release --test engine_graded_notch -- --ignored --nocapture
//! ```

use std::time::Instant;

use yee_engine::automesh::auto_dx;
use yee_engine::board::{
    GradedBoardOptions, TwoPortBoardOptions, two_port_board_job, two_port_board_jobs_graded,
};
use yee_engine::sparams;
use yee_engine::{JobEvent, JobSpec, submit};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff, open_end_delta_l};

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
    let bw_hz = 0.8 * F0_HZ;

    // --- The FS.0b.2a fixture: rulebook grid + shared-grid job pair, no
    // hand-set spacing anywhere.
    let opts = GradedBoardOptions::for_board(&layout, F_MAX_HZ, F0_HZ, bw_hz);
    let (dut_job, ref_job) =
        two_port_board_jobs_graded(&layout, F_MAX_HZ, &opts).expect("graded fixture failed");
    let (nx, ny, nz) = (dut_job.spec.nx, dut_job.spec.ny, dut_job.spec.nz);
    eprintln!(
        "engine-graded-001: coarse = {:.4} mm (spec dx_m), grid {nx}x{ny}x{nz} = {} cells, \
         {} steps, dt = {:.4e} s, probe spacing {:.3} mm",
        dut_job.spec.dx_m * 1e3,
        dut_job.cells,
        dut_job.spec.n_steps,
        dut_job.dt_s,
        dut_job.spacing_m * 1e3,
    );

    // --- Solve (reference first, then DUT), timed. -----------------------
    let t_start = Instant::now();
    let (ref_probes, dt_ref) = run(ref_job.spec.clone());
    let t_ref = t_start.elapsed();
    let (dut_probes, dt_dut) = run(dut_job.spec.clone());
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
            dut_job.spacing_m,
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
    let cells_graded = dut_job.cells;
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
