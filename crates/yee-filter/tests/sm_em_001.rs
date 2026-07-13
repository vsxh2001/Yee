//! Gate `sm-em-001` (FS.5b.1, ADR-0218): **aggressive space mapping with
//! the engine as the fine model** — the FS.5b.0 ASM machinery (ADR-0213,
//! proven on a closed-form warp) closes the loop against real EM solves.
//!
//! Scenario: the S.6 open-stub notch board, one knob (stub length), one
//! response (measured notch frequency). The coarse model is the TL
//! formula the S.6 synthesis uses — known biased vs the engine by ~1 %
//! (the ADR-0216 trajectory) plus staircase effects; ASM must land the
//! *measured* notch on an off-design 5.3 GHz target in a handful of fine
//! evaluations, beating the coarse-only seed.
//!
//! `#[ignore]`'d (~2 release FDTD solves per fine eval):
//!
//! ```bash
//! cargo test -p yee-filter --release --test sm_em_001 -- --ignored --nocapture
//! ```

use std::cell::RefCell;

use yee_engine::board::{GradedBoardOptions, two_port_board_jobs_graded};
use yee_engine::{JobEvent, sparams, submit};
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate, eps_eff, open_end_delta_l};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const C0_M_S: f64 = 299_792_458.0;
const Z0_OHM: f64 = 50.0;
/// Mesh design frequency: keeps the measured band (≤ 6.2 GHz) inside the
/// ADR-0216-safe 0.96·f_max.
const F_MAX_HZ: f64 = 6.5e9;
/// Drive centre; the notch target sits above the S.6 design point.
const F0_HZ: f64 = 5.0e9;
const F_TARGET_HZ: f64 = 5.3e9;

/// The S.6 board with an explicit stub length (the ASM design knob).
fn stub_layout(stub_len_m: f64) -> Layout {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len_m);
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

/// The coarse model: TL notch frequency (GHz) for a stub length (mm).
/// Units chosen so both design and response are O(1)–O(10).
fn coarse_notch_ghz(l_mm: f64) -> f64 {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let dl = open_end_delta_l(W_M, H_M, e_eff);
    C0_M_S / (4.0 * (l_mm * 1e-3 + dl) * e_eff.sqrt()) / 1e9
}

/// The fine model: measured notch frequency (GHz) of the graded two-port
/// double-ratio |S21|, refined sub-bin by a parabola through the three
/// dB bins around the minimum.
fn measured_notch_ghz(stub_len_mm: f64) -> f64 {
    let dut = stub_layout(stub_len_mm * 1e-3);
    let mut opts = GradedBoardOptions::for_board(&dut, F_MAX_HZ, F0_HZ, 0.8 * F0_HZ);
    // Edge-snapped meshing (ADR-0218): without it the stub length
    // quantizes to fine-cell multiples and the notch response is a
    // staircase in the design variable — the first run measured three
    // lengths spanning 34 µm reading the identical 5.3530 GHz, and
    // Broyden oscillated inside the step. Snapping makes the rasterized
    // stub end track the requested length continuously.
    opts.mesh.snap_edges = true;
    let (dut_job, ref_job) =
        two_port_board_jobs_graded(&dut, F_MAX_HZ, &opts).expect("graded job build failed");
    let (dt, spacing) = (dut_job.dt_s, dut_job.spacing_m);
    let run = |spec: yee_engine::JobSpec| -> Vec<Vec<f64>> {
        let handle = submit(spec);
        for event in handle.events() {
            match event {
                JobEvent::Done { result } => return result.probes,
                JobEvent::Error { message } => panic!("job failed: {message}"),
                _ => {}
            }
        }
        panic!("engine stream ended without a result");
    };
    // 25 MHz bins over 4.6–6.2 GHz (≤ 0.96·f_max per ADR-0216).
    let freqs: Vec<f64> = (0..=64).map(|n| 4.6e9 + n as f64 * 25.0e6).collect();
    let ref_p = run(ref_job.spec);
    let dut_p = run(dut_job.spec);
    let transfer = |p: &[Vec<f64>]| {
        sparams::forward_transfer(
            [&p[0], &p[1], &p[2]],
            [&p[3], &p[4], &p[5]],
            dt,
            spacing,
            &freqs,
        )
    };
    let t_dut = transfer(&dut_p);
    let t_ref = transfer(&ref_p);
    let s21_db: Vec<f64> = t_dut
        .iter()
        .zip(&t_ref)
        .map(|(d, r)| 20.0 * (d.0.hypot(d.1) / r.0.hypot(r.1)).log10())
        .collect();
    let i_min = s21_db
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .unwrap()
        .0;
    assert!(
        i_min > 0 && i_min < freqs.len() - 1,
        "notch at the band edge (bin {i_min}) — widen the sweep"
    );
    // Parabolic sub-bin refinement through (i−1, i, i+1) in dB.
    let (a, b, c) = (s21_db[i_min - 1], s21_db[i_min], s21_db[i_min + 1]);
    let denom = a - 2.0 * b + c;
    let shift = if denom.abs() > 1e-12 {
        0.5 * (a - c) / denom
    } else {
        0.0
    };
    (freqs[i_min] + shift.clamp(-0.5, 0.5) * 25.0e6) / 1e9
}

#[test]
#[ignore = "slow: ~2 release FDTD solves per fine eval; sm-em-001 gate (FS.5b.1) — run with --release --ignored"]
fn graded_space_mapping_lands_the_measured_notch_on_target() {
    let f_t_ghz = F_TARGET_HZ / 1e9;
    // Coarse-optimal design: the exact TL inverse of the target.
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let dl = open_end_delta_l(W_M, H_M, e_eff);
    let z_star_mm = (C0_M_S / (4.0 * F_TARGET_HZ * e_eff.sqrt()) - dl) * 1e3;
    eprintln!(
        "sm-em-001: target {f_t_ghz:.3} GHz, coarse-optimal stub {z_star_mm:.4} mm \
         (coarse check: {:.4} GHz)",
        coarse_notch_ghz(z_star_mm)
    );

    // The fine closure logs every (length, frequency) pair.
    let log: RefCell<Vec<(f64, f64)>> = RefCell::new(Vec::new());
    let fine = |x: &[f64]| -> Vec<f64> {
        let f = measured_notch_ghz(x[0]);
        eprintln!(
            "  fine eval {}: stub {:.4} mm → {:.4} GHz",
            log.borrow().len(),
            x[0],
            f
        );
        log.borrow_mut().push((x[0], f));
        vec![f]
    };
    let coarse = |z: &[f64]| -> Vec<f64> { vec![coarse_notch_ghz(z[0])] };

    let cfg = yee_surrogate::spacemap::SpaceMapConfig {
        max_fine_evals: 5,
        // 0.005 scaled by the nominal length ≈ 0.5 % in frequency via
        // df/dl ≈ −f/l.
        tol: 0.005,
        scale: vec![z_star_mm],
        extract: yee_surrogate::spacemap::ExtractConfig::default(),
    };
    let result =
        yee_surrogate::spacemap::space_map(&fine, &coarse, &[z_star_mm], &[z_star_mm], &cfg);

    let log = log.into_inner();
    let seed_err = (log[0].1 - f_t_ghz).abs() / f_t_ghz;
    let (l_fin, f_fin) = *log.last().unwrap();
    let final_err = (f_fin - f_t_ghz).abs() / f_t_ghz;
    eprintln!(
        "  ASM: {} fine evals, converged = {}, misalignment {:.5}; \
         seed err {:.3} % → final err {:.3} % (stub {:.4} mm, notch {:.4} GHz)",
        result.n_fine_evals,
        result.converged,
        result.misalignment,
        seed_err * 100.0,
        final_err * 100.0,
        l_fin,
        f_fin
    );

    assert!(
        result.converged,
        "sm-em-001 FAILED: ASM not converged in {} fine evals (misalignment {:.5})",
        result.n_fine_evals, result.misalignment
    );
    assert!(
        result.n_fine_evals <= 4,
        "sm-em-001 FAILED: needed {} fine evals (> 4)",
        result.n_fine_evals
    );
    assert!(
        final_err <= 0.0075,
        "sm-em-001 FAILED: final measured notch {:.4} GHz vs target {f_t_ghz:.3} \
         (err {:.3} % > 0.75 %)",
        f_fin,
        final_err * 100.0
    );
    assert!(
        final_err < seed_err,
        "sm-em-001 FAILED: mapping did not beat the coarse seed \
         ({:.3} % vs {:.3} %)",
        final_err * 100.0,
        seed_err * 100.0
    );
}
