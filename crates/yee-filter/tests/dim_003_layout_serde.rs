//! dim-003 (Filter Phase F1.2.0): layout assembly + serde round-trip.
//!
//! `dimension_edge_coupled_layout` returns a non-degenerate `yee_layout::Layout`
//! (≥ 1 polygon, the expected two ports), and `EdgeCoupledDimensions` serde
//! round-trips byte-identically through JSON.

use yee_filter::{
    Approximation, EdgeCoupledDimensions, FilterSpec, Response, SpecMask, dimension_edge_coupled,
    dimension_edge_coupled_layout, synthesize,
};
use yee_layout::Substrate;

fn fixture() -> (FilterSpec, Substrate) {
    let spec = FilterSpec {
        response: Response::Bandpass,
        approximation: Approximation::Chebyshev { ripple_db: 0.5 },
        f0_hz: 2.0e9,
        fbw: 0.10,
        order: Some(5),
        z0_ohm: 50.0,
        mask: SpecMask {
            passband_ripple_db: 0.5,
            return_loss_db: 10.0,
            stopband: vec![(2.4e9, 30.0)],
        },
    };
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    (spec, substrate)
}

#[test]
fn dim_003_layout_nondegenerate() {
    let (spec, substrate) = fixture();
    let proj = synthesize(&spec);
    let layout = dimension_edge_coupled_layout(&proj, &substrate).expect("fixture layout");

    // N=5 resonators + 2 feed lines = 7 trace polygons; assert at least one,
    // and that the count matches the edge-coupled generator's contract.
    assert!(
        !layout.traces.is_empty(),
        "layout must have at least one trace polygon"
    );
    assert_eq!(
        layout.traces.len(),
        5 + 2,
        "N=5 → 5 coupled strips + 2 feed lines"
    );

    // Edge-coupled BPF has exactly two ports (input + output feed).
    assert_eq!(layout.ports.len(), 2, "edge-coupled BPF has 2 ports");

    // Non-degenerate bounding box.
    assert!(
        layout.bbox.width() > 0.0 && layout.bbox.height() > 0.0,
        "layout bounding box must have positive extent"
    );
}

#[test]
fn dim_003_dimensions_serde_roundtrip() {
    let (spec, substrate) = fixture();
    let proj = synthesize(&spec);
    let dims = dimension_edge_coupled(&proj, &substrate).expect("fixture dimensions");

    let json = serde_json::to_string(&dims).expect("serialize EdgeCoupledDimensions");
    let back: EdgeCoupledDimensions =
        serde_json::from_str(&json).expect("deserialize EdgeCoupledDimensions");
    let json2 = serde_json::to_string(&back).expect("re-serialize EdgeCoupledDimensions");

    assert_eq!(
        dims, back,
        "round-tripped dimensions must equal the original"
    );
    assert_eq!(json, json2, "serde round-trip must be byte-identical");
}
