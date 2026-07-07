//! Gate `studio-verify-e2e-001` (R.5b, ADR-0199): the studio's full-wave
//! verify **pipe**, driven headlessly at reduced fidelity — the physics of
//! the measurement fixture is gated in `yee-engine`/`yee-filter`; this
//! gate pins the studio-side flow:
//!
//! 1. progress streams for both phases, reference strictly before DUT,
//!    monotone within each phase;
//! 2. the identity case measures exactly: verifying a straight line whose
//!    "reference" is the same line runs two bit-identical solves, so the
//!    directional |S21| must read 0 dB to numerical identity;
//! 3. a broken layout (feeds too short for the probe triples) surfaces
//!    the shared `yee_engine::board` builder's error, not a panic.

use yee_engine::board::TwoPortBoardOptions;
use yee_layout::{BBox, Layout, Point2, Polygon, PortRef, Substrate};
use yee_studio_app::verify::verify_layout_impl;

fn line_layout(len_m: f64) -> Layout {
    let w = 1.5e-3;
    let traces = vec![Polygon::rect(0.0, -w / 2.0, len_m, w)];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: Substrate {
            eps_r: 4.4,
            height_m: 0.8e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        },
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.0, 0.0),
                width_m: w,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(len_m, 0.0),
                width_m: w,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

/// Reduced-fidelity options: coarse grid, short run — a pipe exercise,
/// not a physics measurement.
fn fast_opts() -> TwoPortBoardOptions {
    let mut opts = TwoPortBoardOptions::for_band(5.0e9, 4.0e9);
    opts.dx_m = 0.4e-3;
    opts.n_steps = 700;
    opts
}

#[test]
fn verify_pipe_streams_phases_and_measures_the_identity_case_exactly() {
    let layout = line_layout(30.0e-3);
    let freqs: Vec<f64> = (0..5).map(|n| 4.0e9 + n as f64 * 0.5e9).collect();
    let mut events: Vec<(String, usize, usize)> = vec![];
    let (s21_db, backend) = verify_layout_impl(&layout, &fast_opts(), &freqs, &mut |p| {
        events.push((p.phase.to_string(), p.step, p.total));
    })
    .expect("verify pipe failed");

    // 1. Both phases streamed, reference strictly before dut, monotone.
    let first_dut = events
        .iter()
        .position(|(ph, _, _)| ph == "dut")
        .expect("no dut progress");
    assert!(
        events[..first_dut].iter().all(|(ph, _, _)| ph == "reference"),
        "phases interleaved"
    );
    assert!(first_dut > 0, "no reference progress");
    for w in events.windows(2) {
        if w[0].0 == w[1].0 {
            assert!(w[1].1 >= w[0].1, "progress went backwards");
        }
    }

    // 2. Identity case: the "DUT" is the reference line itself, so the two
    //    solves are bit-identical and the directional ratio is exactly 1.
    assert_eq!(s21_db.len(), freqs.len());
    for (f, db) in freqs.iter().zip(&s21_db) {
        assert!(
            db.abs() < 1e-9,
            "identity |S21| at {:.1} GHz reads {db} dB (want exactly 0)",
            f / 1e9
        );
    }
    assert!(!backend.is_empty());
}

#[test]
fn short_feeds_surface_the_builder_error() {
    let layout = line_layout(5.0e-3);
    let err = verify_layout_impl(&layout, &fast_opts(), &[5.0e9], &mut |_| {})
        .expect_err("short feeds must be rejected");
    assert!(err.contains("too short"), "{err}");
}
