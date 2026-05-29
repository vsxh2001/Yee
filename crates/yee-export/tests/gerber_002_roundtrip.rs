//! `gerber-002` — coordinate round-trip gate.
//!
//! Parse the `X<int>Y<int>` coordinate words back out of one polygon's
//! `G36*…G37*` region and assert they reproduce that polygon's vertices after
//! inverting the metres → mm → 4.6-fixed-point conversion, to within the
//! `1e-6 mm` format quantisation. This validates the coordinate emission — the
//! part that actually matters for fabrication.

use yee_export::{GerberOptions, layout_to_gerber};
use yee_layout::{BBox, Layout, Point2, Polygon, Substrate};

/// Inverse of the emitter's conversion: a `X<int>Y<int>` 4.6-fixed-point pair
/// (millimetres) back to metres.
fn fixed46_to_metres(fixed: i64) -> f64 {
    // int = round(mm * 1e6)  ⇒  mm = int / 1e6  ⇒  m = mm / 1e3.
    (fixed as f64) / 1.0e6 / 1.0e3
}

/// Pull every `X<int>Y<int>` coordinate word out of the *first* `G36*…G37*`
/// region of a Gerber string, as `(x_metres, y_metres)` pairs (drops the
/// trailing `Dnn*` operation code).
fn parse_first_region_coords(gerber: &str) -> Vec<(f64, f64)> {
    let start = gerber.find("G36*").expect("no region in output");
    let rest = &gerber[start..];
    let end = rest.find("G37*").expect("region not terminated");
    let region = &rest[..end];

    let mut coords = Vec::new();
    for line in region.lines() {
        let line = line.trim();
        if !line.starts_with('X') {
            continue;
        }
        // Form: X<int>Y<int>D0n*  — split on the markers.
        let after_x = &line[1..];
        let y_pos = after_x.find('Y').expect("coordinate word missing Y");
        let x_str = &after_x[..y_pos];
        let after_y = &after_x[y_pos + 1..];
        // The Y integer runs up to the next non-digit (the 'D' op code).
        let d_pos = after_y
            .find(|c: char| !c.is_ascii_digit() && c != '-')
            .unwrap_or(after_y.len());
        let y_str = &after_y[..d_pos];

        let xi: i64 = x_str.parse().expect("bad X integer");
        let yi: i64 = y_str.parse().expect("bad Y integer");
        coords.push((fixed46_to_metres(xi), fixed46_to_metres(yi)));
    }
    coords
}

#[test]
fn gerber_002_coordinate_roundtrip() {
    // A single rectangle with deliberately non-round metre coordinates so the
    // round-trip exercises the fixed-point rounding, not just integers.
    // 3.0590 mm = 3.059e-3 m → X3059000.
    let rect = Polygon::rect(3.059e-3, 1.234e-3, 10.0e-3, 0.500e-3);
    let expected = rect.verts.clone();

    let traces = vec![rect];
    let substrate = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.02,
        metal_thickness_m: 35e-6,
    };
    let layout = Layout {
        substrate,
        bbox: BBox::from_polygons(&traces),
        traces,
        ports: vec![],
    };

    let gerber = layout_to_gerber(&layout, &GerberOptions::default());
    let parsed = parse_first_region_coords(&gerber);

    // The emitter draws each vertex once, then closes back to the first
    // vertex (the rect is not explicitly closed), so we expect N+1 coords.
    assert_eq!(
        parsed.len(),
        expected.len() + 1,
        "expected {} drawn vertices (N + closing), got {}",
        expected.len() + 1,
        parsed.len()
    );

    // The first N parsed coordinates must reproduce the polygon vertices.
    // Tolerance = 1e-6 mm = 1e-9 m (the 4.6 fixed-point quantisation).
    const TOL_M: f64 = 1.0e-9;
    for (i, (exp, &(gx, gy))) in expected.iter().zip(&parsed).enumerate() {
        assert!(
            (exp.x - gx).abs() <= TOL_M && (exp.y - gy).abs() <= TOL_M,
            "vertex {i} mismatch: expected ({}, {}), got ({}, {})",
            exp.x,
            exp.y,
            gx,
            gy
        );
    }

    // The closing coordinate must equal the first vertex.
    let first = Point2::new(parsed[0].0, parsed[0].1);
    let last = parsed.last().unwrap();
    assert!(
        (first.x - last.0).abs() <= TOL_M && (first.y - last.1).abs() <= TOL_M,
        "region must close back to the first vertex"
    );
}
