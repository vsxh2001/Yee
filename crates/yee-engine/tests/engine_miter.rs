//! Gate `engine-miter-001` (FS.3.2a, ADR-0217): the repo's first
//! **non-axis-aligned polygon edges under full-wave test**. Two
//! double-jog through lines (four 90° bends each) — square corners vs
//! 45°-mitered outer corners (Douville & James) — measured through the
//! certified graded two-port fixture. The miter physics is the reference:
//! chopping the outer corner removes the bend's excess capacitance, so
//! the mitered line must transmit at least as well as the square one —
//! at every bin. (The first run measured the advantage U-shaped across
//! the band — +0.021 linear at the edges, ~0 near 4.45 GHz — because the
//! four bends' reflections interfere with path-length-dependent phase;
//! a naive \"advantage grows with f\" excess-C assert is mis-specified
//! for a multi-bend line.)
//!
//! `#[ignore]`'d (4 release FDTD solves):
//!
//! ```bash
//! cargo test -p yee-engine --release --test engine_miter -- --ignored --nocapture
//! ```

use yee_engine::board::{GradedBoardOptions, two_port_board_jobs_graded};
use yee_engine::{JobEvent, sparams, submit};
use yee_layout::{Layout, MiterStyle, Substrate, double_jog};

const F0_HZ: f64 = 5.0e9;
const F_MAX_HZ: f64 = 6.0e9;
const W_M: f64 = 3.0e-3;

fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

fn jog(style: MiterStyle) -> Layout {
    double_jog(&substrate(), W_M, 24.0e-3, 9.0e-3, 9.0e-3, style)
}

/// Launch-normalized double-ratio |S21| (linear) over `freqs` for one
/// DUT: both jobs on the DUT-derived graded grid (the ADR-0204 lesson).
fn s21_lin(dut: &Layout, freqs: &[f64]) -> Vec<f64> {
    let opts = GradedBoardOptions::for_board(dut, F_MAX_HZ, F0_HZ, 0.8 * F0_HZ);
    let (dut_job, ref_job) =
        two_port_board_jobs_graded(dut, F_MAX_HZ, &opts).expect("graded job build failed");
    let (dt, spacing) = (dut_job.dt_s, dut_job.spacing_m);
    eprintln!(
        "  grid {} cells, dt {:.3e} s, {} steps",
        dut_job.cells, dt, dut_job.spec.n_steps
    );
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
    let ref_p = run(ref_job.spec);
    let dut_p = run(dut_job.spec);
    let transfer = |p: &[Vec<f64>]| {
        sparams::forward_transfer(
            [&p[0], &p[1], &p[2]],
            [&p[3], &p[4], &p[5]],
            dt,
            spacing,
            freqs,
        )
    };
    let t_dut = transfer(&dut_p);
    let t_ref = transfer(&ref_p);
    t_dut
        .iter()
        .zip(&t_ref)
        .map(|(d, r)| d.0.hypot(d.1) / r.0.hypot(r.1))
        .collect()
}

#[test]
#[ignore = "slow: 4 release FDTD solves; engine-miter-001 gate (FS.3.2a) — run with --release --ignored"]
fn graded_miter_bends_transmit_better_than_square() {
    // The ADR-0216 criterion band: stop at 0.96·f_max.
    let freqs: Vec<f64> = (0..=45).map(|n| 3.5e9 + n as f64 * 50.0e6).collect();

    eprintln!("engine-miter-001: square corners");
    let sq = s21_lin(&jog(MiterStyle::Square), &freqs);
    eprintln!("engine-miter-001: mitered corners (f = 0.7)");
    let mi = s21_lin(&jog(MiterStyle::Mitered { f: 0.7 }), &freqs);

    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let db = |x: f64| 20.0 * x.log10();
    for (i, f) in freqs.iter().enumerate() {
        eprintln!(
            "    {:.2} GHz: square {:7.2} dB, mitered {:7.2} dB (Δ {:+.3} lin)",
            f / 1e9,
            db(sq[i]),
            db(mi[i]),
            mi[i] - sq[i]
        );
    }
    let (mean_sq, mean_mi) = (mean(&sq), mean(&mi));
    let min_gap = mi
        .iter()
        .zip(&sq)
        .map(|(m, s)| m - s)
        .fold(f64::INFINITY, f64::min);
    eprintln!(
        "  mean |S21|: square {:.4} ({:.2} dB), mitered {:.4} ({:.2} dB); \
         min per-bin miter advantage {:+.4} lin",
        mean_sq,
        db(mean_sq),
        mean_mi,
        db(mean_mi),
        min_gap
    );

    // Measured 2026-07-12 (first run, 621 s release): square mean 0.9665
    // (−0.30 dB), mitered 0.9738 (−0.23 dB) — gap 0.0073; min per-bin
    // advantage +0.000 (4.45 GHz, the interference null); mitered worst
    // in-band −0.27 dB.
    // (a) The miter physics, band-mean in linear magnitude, at half the
    // measured margin.
    assert!(
        mean_mi >= mean_sq + 0.004,
        "engine-miter-001 FAILED: mitered mean |S21| {mean_mi:.4} vs square {mean_sq:.4} \
         (need +0.004 linear margin)"
    );
    // (b) The mitered line stays a clean through-line across the band.
    let worst_mi = mi.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(
        db(worst_mi) >= -1.0,
        "engine-miter-001 FAILED: mitered worst in-band |S21| {:.2} dB < -1 dB",
        db(worst_mi)
    );
    // (c) Mitering never hurts, bin by bin (ripple slack: the measured
    // minimum sits at the four-bend interference null, ~0.000).
    assert!(
        min_gap >= -0.005,
        "engine-miter-001 FAILED: mitered worse than square somewhere in-band \
         (min per-bin gap {min_gap:+.4})"
    );
}
