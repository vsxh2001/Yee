//! `kicad-001` — KiCad `.kicad_pcb` structural-validity gate.
//!
//! For a known small `Layout` with two rectangular traces, assert the emitted
//! board: (a) starts with `(kicad_pcb`; (b) has balanced parentheses; (c)
//! contains a `(layers` block naming both `F.Cu` and `Edge.Cuts`; (d) has
//! exactly one `(gr_poly` carrying `(layer "F.Cu")` per trace polygon; (e) has
//! exactly one `(gr_poly` carrying `(layer "Edge.Cuts")`.
//!
//! Like `gerber-001`/`gerber-003`, this is the I/O-structural analogue of a
//! physics gate — `.kicad_pcb` is an interchange format, so its gate is
//! structural validity (the CI box has no KiCad to open it with).

use yee_export::{KicadPcbOptions, layout_to_kicad_pcb};
use yee_layout::{BBox, Layout, Polygon, Substrate};

/// A small, deterministic two-rectangle layout with a known bbox.
fn sample_layout() -> Layout {
    let traces = vec![
        Polygon::rect(2.0e-3, 1.0e-3, 10.0e-3, 0.5e-3),
        Polygon::rect(2.0e-3, 4.0e-3, 10.0e-3, 0.5e-3),
    ];
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    Layout {
        substrate,
        bbox: BBox::from_polygons(&traces),
        traces,
        ports: vec![],
    }
}

#[test]
fn kicad_001_structure() {
    let layout = sample_layout();
    let n_traces = layout.traces.len();
    let pcb = layout_to_kicad_pcb(&layout, &KicadPcbOptions::default());

    // (a) Top-level form is a (kicad_pcb …).
    assert!(
        pcb.starts_with("(kicad_pcb"),
        "board must start with (kicad_pcb, got:\n{}",
        &pcb[..pcb.len().min(80)]
    );

    // (b) Parentheses are balanced (count of '(' equals count of ')').
    let opens = pcb.matches('(').count();
    let closes = pcb.matches(')').count();
    assert_eq!(
        opens, closes,
        "parentheses must balance: {opens} '(' vs {closes} ')'"
    );

    // (c) A (layers block naming both F.Cu and Edge.Cuts.
    assert!(
        pcb.contains("(layers"),
        "board must contain a (layers table"
    );
    assert!(pcb.contains("\"F.Cu\""), "layer table must name F.Cu");
    assert!(
        pcb.contains("\"Edge.Cuts\""),
        "layer table must name Edge.Cuts"
    );

    // (d) Exactly one (gr_poly carrying (layer "F.Cu") per trace polygon.
    let n_fcu = pcb.matches("(layer \"F.Cu\")").count();
    assert_eq!(
        n_fcu, n_traces,
        "expected one F.Cu gr_poly per trace ({n_traces}), found {n_fcu}"
    );

    // (e) Exactly one (gr_poly carrying (layer "Edge.Cuts").
    let n_edge = pcb.matches("(layer \"Edge.Cuts\")").count();
    assert_eq!(
        n_edge, 1,
        "expected exactly one Edge.Cuts gr_poly, found {n_edge}"
    );

    // Cross-check: total gr_poly count is traces + 1 (the outline).
    let n_grpoly = pcb.matches("(gr_poly").count();
    assert_eq!(
        n_grpoly,
        n_traces + 1,
        "expected {} gr_poly ({n_traces} traces + 1 outline), found {n_grpoly}",
        n_traces + 1
    );
}
