//! `kicad-002` — KiCad `.kicad_pcb` coordinate-geometry gate.
//!
//! Parse the `(xy X Y)` pairs out of the copper (`F.Cu`) `gr_poly` lines and
//! confirm they equal the trace vertices in millimetres (metres × 1e3) within
//! `1e-6`; then parse the `Edge.Cuts` `gr_poly` and confirm its rectangle
//! equals `bbox ± outline_margin_mm`. This confirms the emitted board carries
//! the layout geometry faithfully — the part that actually matters when the
//! file is opened in KiCad.

use yee_export::{KicadPcbOptions, layout_to_kicad_pcb};
use yee_layout::{BBox, Layout, Polygon, Substrate};

/// Pull every `(xy X Y)` pair out of one line, as `(x_mm, y_mm)` floats.
fn parse_xy(line: &str) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    let mut rest = line;
    while let Some(start) = rest.find("(xy ") {
        let after = &rest[start + 4..];
        let end = after.find(')').expect("(xy … missing closing paren");
        let body = &after[..end];
        let mut it = body.split_whitespace();
        let x: f64 = it.next().expect("missing x").parse().expect("bad x");
        let y: f64 = it.next().expect("missing y").parse().expect("bad y");
        out.push((x, y));
        rest = &after[end + 1..];
    }
    out
}

#[test]
fn kicad_002_geometry() {
    // Deliberately non-round metre coordinates so the float conversion is
    // exercised (not just integers).
    let traces = vec![
        Polygon::rect(3.059e-3, 1.234e-3, 10.0e-3, 0.5e-3),
        Polygon::rect(3.059e-3, 4.0e-3, 8.0e-3, 0.5e-3),
    ];
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate,
        bbox,
        traces: traces.clone(),
        ports: vec![],
    };

    // Non-default margin (0.75 mm) to prove outline_margin_mm is honoured.
    let opts = KicadPcbOptions {
        copper_layer: "F.Cu".into(),
        outline_margin_mm: 0.75,
        generator: "yee-export".into(),
    };
    let pcb = layout_to_kicad_pcb(&layout, &opts);

    const TOL_MM: f64 = 1.0e-6;

    // --- Copper polygons: each F.Cu gr_poly equals its trace's vertices -----
    let fcu_lines: Vec<&str> = pcb
        .lines()
        .filter(|l| l.contains("(gr_poly") && l.contains("(layer \"F.Cu\")"))
        .collect();
    assert_eq!(
        fcu_lines.len(),
        traces.len(),
        "expected one F.Cu gr_poly per trace"
    );

    for (line, poly) in fcu_lines.iter().zip(traces.iter()) {
        let parsed = parse_xy(line);
        assert_eq!(
            parsed.len(),
            poly.verts.len(),
            "F.Cu gr_poly must carry every trace vertex (no implicit-close repeat)"
        );
        for (got, v) in parsed.iter().zip(poly.verts.iter()) {
            let exp = (v.x * 1.0e3, v.y * 1.0e3);
            assert!(
                (got.0 - exp.0).abs() <= TOL_MM && (got.1 - exp.1).abs() <= TOL_MM,
                "vertex mismatch: emitted {got:?} mm vs expected {exp:?} mm"
            );
        }
    }

    // --- Board outline: Edge.Cuts rectangle equals bbox ± margin ------------
    let edge_line = pcb
        .lines()
        .find(|l| l.contains("(gr_poly") && l.contains("(layer \"Edge.Cuts\")"))
        .expect("missing Edge.Cuts gr_poly");
    let outline = parse_xy(edge_line);
    assert_eq!(
        outline.len(),
        4,
        "Edge.Cuts gr_poly must have exactly 4 corners (closed implicitly)"
    );

    let m = opts.outline_margin_mm * 1.0e-3; // metres.
    let expected_mm = [
        ((bbox.min.x - m) * 1.0e3, (bbox.min.y - m) * 1.0e3),
        ((bbox.max.x + m) * 1.0e3, (bbox.min.y - m) * 1.0e3),
        ((bbox.max.x + m) * 1.0e3, (bbox.max.y + m) * 1.0e3),
        ((bbox.min.x - m) * 1.0e3, (bbox.max.y + m) * 1.0e3),
    ];
    for (got, exp) in outline.iter().zip(expected_mm.iter()) {
        assert!(
            (got.0 - exp.0).abs() <= TOL_MM && (got.1 - exp.1).abs() <= TOL_MM,
            "outline corner mismatch: emitted {got:?} mm vs expected {exp:?} mm"
        );
    }
}
