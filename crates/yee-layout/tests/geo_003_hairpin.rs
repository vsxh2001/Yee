//! `geo-003` — hairpin generator geometry smoke gate.
//!
//! The `edge_coupled_bpf` generator is gated by `geo-001`; this file gives
//! `hairpin_bpf` the same minimum coverage (per the crate's validation
//! convention — every generator carries a gate). Builds a 3-resonator hairpin
//! and asserts: trace count `n*3 + 2` (two arms + a bend per resonator, plus
//! the two feed lines), exactly two ports, every polygon non-degenerate
//! (≥ 4 vertices, positive shoelace signed area), and a serde JSON round-trip.
//! Pure geometry — deterministic ground truth.

use yee_layout::{HairpinParams, Layout, Substrate, hairpin_bpf};

fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

fn params(n: usize) -> HairpinParams {
    HairpinParams {
        substrate: fr4(),
        n,
        arm_length_m: 10.0e-3,
        line_width_m: 1.0e-3,
        fold_spacing_m: 2.0e-3,
        coupling_gap_m: 0.3e-3,
        tap_offset_m: 2.0e-3,
        feed_width_m: 1.0e-3,
        feed_length_m: 5.0e-3,
    }
}

#[test]
fn hairpin_trace_and_port_counts() {
    let n = 3;
    let layout = hairpin_bpf(&params(n));
    assert_eq!(
        layout.traces.len(),
        n * 3 + 2,
        "hairpin should emit two arms + a bend per resonator plus two feeds"
    );
    assert_eq!(layout.ports.len(), 2, "a 2-port band-pass filter");
}

#[test]
fn hairpin_polygons_non_degenerate() {
    let layout = hairpin_bpf(&params(3));
    for (i, poly) in layout.traces.iter().enumerate() {
        assert!(
            poly.verts.len() >= 4,
            "trace {i} must have >= 4 vertices, got {}",
            poly.verts.len()
        );
        assert!(
            poly.signed_area() > 0.0,
            "trace {i} must have positive (CCW) signed area, got {}",
            poly.signed_area()
        );
    }
}

#[test]
fn hairpin_serde_round_trip() {
    let layout = hairpin_bpf(&params(4));
    let json = serde_json::to_string(&layout).expect("serialize");
    let back: Layout = serde_json::from_str(&json).expect("deserialize");

    // Round-trip invariant: structural shape preserved and every coordinate
    // recovered to within 1 nm. We use an approximate (not bit-exact) compare:
    // the hairpin generator's coordinates carry harmless floating-point
    // accumulation (e.g. 0.0086 + 0.001) whose decimal serialization is not
    // bit-stable across a parse — 1e-9 m is ~1e-6 of the resonator dimensions,
    // far tighter than any manufacturing tolerance.
    const TOL_M: f64 = 1e-9;
    assert_eq!(back.traces.len(), layout.traces.len(), "trace count");
    assert_eq!(back.ports.len(), layout.ports.len(), "port count");
    for (t, (pa, pb)) in layout.traces.iter().zip(&back.traces).enumerate() {
        assert_eq!(pa.verts.len(), pb.verts.len(), "trace {t} vertex count");
        for (va, vb) in pa.verts.iter().zip(&pb.verts) {
            assert!(
                (va.x - vb.x).abs() < TOL_M && (va.y - vb.y).abs() < TOL_M,
                "trace {t} vertex moved across round-trip: ({}, {}) vs ({}, {})",
                va.x,
                va.y,
                vb.x,
                vb.y
            );
        }
    }
}
