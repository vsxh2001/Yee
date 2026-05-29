//! `gerber-004` — board-outline coordinate-geometry gate.
//!
//! Parse the `X<int>Y<int>` coordinate words out of the outline Gerber and
//! assert the four distinct corners equal the layout `bbox.min/max` expanded by
//! the requested margin, after inverting the metres → mm → 4.6-fixed-point
//! conversion, to within the `1e-6 mm` format quantisation. This confirms the
//! emitted profile encloses the layout with exactly the requested margin — the
//! part that actually matters for fabrication.

use yee_export::{OutlineOptions, layout_to_gerber_outline};
use yee_layout::{BBox, Layout, Polygon, Substrate};

/// Inverse of the emitter's conversion: a `X<int>Y<int>` 4.6-fixed-point pair
/// (millimetres) back to metres.
fn fixed46_to_metres(fixed: i64) -> f64 {
    // int = round(mm * 1e6)  ⇒  mm = int / 1e6  ⇒  m = mm / 1e3.
    (fixed as f64) / 1.0e6 / 1.0e3
}

/// Pull every `X<int>Y<int>` coordinate word out of a Gerber string, as
/// `(x_metres, y_metres)` pairs (drops the trailing `Dnn*` operation code).
fn parse_coords(gerber: &str) -> Vec<(f64, f64)> {
    let mut coords = Vec::new();
    for line in gerber.lines() {
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

/// Assert `(x, y)` matches one of the `expected` corners within `TOL_M`.
fn assert_is_corner(x: f64, y: f64, expected: &[(f64, f64)]) {
    const TOL_M: f64 = 1.0e-9; // 1e-6 mm format quantisation.
    let hit = expected
        .iter()
        .any(|&(ex, ey)| (ex - x).abs() <= TOL_M && (ey - y).abs() <= TOL_M);
    assert!(
        hit,
        "emitted coordinate ({x}, {y}) is not one of the expected corners {expected:?}"
    );
}

#[test]
fn gerber_004_outline_geometry() {
    // A single rectangle with deliberately non-round metre coordinates so the
    // bbox ± margin corners exercise the fixed-point rounding.
    let traces = vec![Polygon::rect(3.059e-3, 1.234e-3, 10.0e-3, 5.0e-3)];
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
        traces,
        ports: vec![],
    };

    // A non-default margin (0.75 mm) to prove margin_mm is honoured.
    let opts = OutlineOptions {
        layer_name: "Edge.Cuts".into(),
        margin_mm: 0.75,
    };
    let m = opts.margin_mm * 1.0e-3; // metres.

    // Expected corners = bbox.min/max ± margin.
    let expected = [
        (bbox.min.x - m, bbox.min.y - m),
        (bbox.max.x + m, bbox.min.y - m),
        (bbox.max.x + m, bbox.max.y + m),
        (bbox.min.x - m, bbox.max.y + m),
    ];

    let gerber = layout_to_gerber_outline(&layout, &opts);
    let parsed = parse_coords(&gerber);

    // The emitter draws corner 0 (move), corners 1/2/3, then closes back to
    // corner 0 → five coordinate words.
    assert_eq!(
        parsed.len(),
        5,
        "expected 5 coordinate words (4 corners + explicit close), got {}",
        parsed.len()
    );

    // Every emitted coordinate must be one of the four expected corners.
    for (gx, gy) in &parsed {
        assert_is_corner(*gx, *gy, &expected);
    }

    // The four distinct corners must all be present (the contour visits each).
    const TOL_M: f64 = 1.0e-9;
    for &(ex, ey) in &expected {
        let present = parsed
            .iter()
            .any(|&(gx, gy)| (ex - gx).abs() <= TOL_M && (ey - gy).abs() <= TOL_M);
        assert!(
            present,
            "expected corner ({ex}, {ey}) not found among emitted coordinates {parsed:?}"
        );
    }

    // The contour must close: the move (first) and the final draw coincide.
    let first = parsed[0];
    let last = *parsed.last().unwrap();
    assert!(
        (first.0 - last.0).abs() <= TOL_M && (first.1 - last.1).abs() <= TOL_M,
        "outline must close back to the first corner: {first:?} vs {last:?}"
    );
}
