//! `geo-001` — edge-coupled generator geometry gate.
//!
//! Builds a 2-section edge-coupled layout with known dimensions and asserts:
//! the trace and port counts are as expected; the bounding box width/height is
//! within 1% of the hand-computed extents; every polygon is non-degenerate
//! (≥ 4 vertices, positive shoelace signed area); the `Layout` round-trips
//! through serde JSON unchanged; and `to_svg()` returns a non-empty SVG
//! document. Pure geometry — these numbers are deterministic ground truth.

use yee_layout::{
    BBox, EdgeCoupledParams, EdgeCoupledSection, Layout, Substrate, edge_coupled_bpf,
};

/// FR-4 substrate fixture.
fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    }
}

// Hand-chosen, hand-computable section dims (metres).
const L: f64 = 20.0e-3; // strip length
const W: f64 = 2.0e-3; // strip width
const G: f64 = 0.3e-3; // coupling gap
const FEED_W: f64 = 1.5e-3; // feed width (< W so it does not extend min.y)
const FEED_L: f64 = 5.0e-3; // feed length

fn two_section_layout() -> Layout {
    let params = EdgeCoupledParams {
        substrate: fr4(),
        sections: vec![
            EdgeCoupledSection {
                length_m: L,
                width_m: W,
                gap_m: G,
            },
            EdgeCoupledSection {
                length_m: L,
                width_m: W,
                gap_m: G,
            },
        ],
        feed_width_m: FEED_W,
        feed_length_m: FEED_L,
    };
    edge_coupled_bpf(&params)
}

#[test]
fn trace_and_port_counts() {
    let layout = two_section_layout();
    // 2 coupled strips + 2 feed lines.
    assert_eq!(
        layout.traces.len(),
        4,
        "expected 2 strips + 2 feeds = 4 traces, got {}",
        layout.traces.len()
    );
    assert_eq!(
        layout.ports.len(),
        2,
        "edge-coupled BPF has exactly 2 ports"
    );
    for port in &layout.ports {
        assert_eq!(port.ref_impedance_ohm, 50.0, "ports default to 50 Ω");
    }
}

#[test]
fn bbox_within_1pct_of_hand_computed() {
    let layout = two_section_layout();

    // x: input feed starts at -FEED_L; strip 1 is staggered to start at L/2 and
    // ends at L/2 + L = 1.5 L; the output feed extends FEED_L beyond that.
    let expect_min_x = -FEED_L;
    let expect_max_x = 1.5 * L + FEED_L;
    let expect_w = expect_max_x - expect_min_x; // 1.5 L + 2 FEED_L

    // y: strip 0 spans [0, W], strip 1 spans [W+G, 2W+G]; feeds (width FEED_W <
    // W) stay inside that band, so the extent is [0, 2W+G].
    let expect_h = 2.0 * W + G;

    let w = layout.bbox.width();
    let h = layout.bbox.height();
    assert!(
        (w - expect_w).abs() / expect_w < 0.01,
        "bbox width {w:.6e} vs expected {expect_w:.6e}"
    );
    assert!(
        (h - expect_h).abs() / expect_h < 0.01,
        "bbox height {h:.6e} vs expected {expect_h:.6e}"
    );
    assert!(
        (layout.bbox.min.x - expect_min_x).abs() < 1e-9,
        "bbox min.x {:.6e} vs expected {expect_min_x:.6e}",
        layout.bbox.min.x
    );
}

#[test]
fn polygons_non_degenerate() {
    let layout = two_section_layout();
    for (i, poly) in layout.traces.iter().enumerate() {
        assert!(
            poly.verts.len() >= 4,
            "trace {i} has {} verts (< 4)",
            poly.verts.len()
        );
        let area = poly.signed_area();
        assert!(
            area > 0.0,
            "trace {i} signed area {area:.6e} is not positive (degenerate or CW)"
        );
    }
}

#[test]
fn serde_json_round_trip() {
    let layout = two_section_layout();
    let json = serde_json::to_string(&layout).expect("serialize Layout to JSON");
    let back: Layout = serde_json::from_str(&json).expect("deserialize Layout from JSON");
    assert_eq!(layout, back, "Layout did not survive a JSON round-trip");
}

#[test]
fn to_svg_is_well_formed() {
    let layout = two_section_layout();
    let svg = layout.to_svg();
    assert!(svg.contains("<svg"), "SVG missing opening <svg tag");
    assert!(svg.contains("</svg>"), "SVG missing closing </svg> tag");
    // Every trace should appear as a polygon element.
    assert_eq!(
        svg.matches("<polygon").count(),
        layout.traces.len(),
        "SVG polygon count should match trace count"
    );
}

#[test]
fn bbox_from_polygons_matches_layout() {
    // Sanity: the stored bbox equals one recomputed from the traces.
    let layout = two_section_layout();
    let recomputed = BBox::from_polygons(&layout.traces);
    assert_eq!(layout.bbox, recomputed);
}
