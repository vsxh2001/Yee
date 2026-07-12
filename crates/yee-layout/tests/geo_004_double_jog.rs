//! Structural tests for the FS.3.2a [`double_jog`] generator: the square
//! variant is all axis-aligned rects; the mitered variant carries exactly
//! four 45° edges (one per bend, on the outer corner), the segment
//! overlaps stay clear of every cut, and both variants share ports/bbox.

use yee_layout::{MiterStyle, Point2, Substrate, double_jog};

const W: f64 = 3.0e-3;
const RUN: f64 = 24.0e-3;
const GAP: f64 = 9.0e-3;
const DY: f64 = 9.0e-3;

fn substrate() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

/// Count edges at exactly 45° (|dx| == |dy|, both nonzero).
fn diagonal_edges(l: &yee_layout::Layout) -> usize {
    let mut n = 0;
    for p in &l.traces {
        let v = &p.verts;
        for i in 0..v.len() {
            let a = v[i];
            let b = v[(i + 1) % v.len()];
            let (dx, dy) = ((b.x - a.x).abs(), (b.y - a.y).abs());
            if dx > 1e-12 && dy > 1e-12 && (dx - dy).abs() < 1e-12 {
                n += 1;
            }
        }
    }
    n
}

#[test]
fn square_variant_is_all_axis_aligned() {
    let l = double_jog(&substrate(), W, RUN, GAP, DY, MiterStyle::Square);
    assert_eq!(l.traces.len(), 9, "5 straights + 4 corners");
    assert_eq!(diagonal_edges(&l), 0);
    // Every corner polygon is a 4-vert rect.
    for p in &l.traces[5..] {
        assert_eq!(p.verts.len(), 4);
    }
}

#[test]
fn mitered_variant_has_exactly_four_45_degree_cuts() {
    let l = double_jog(
        &substrate(),
        W,
        RUN,
        GAP,
        DY,
        MiterStyle::Mitered { f: 0.7 },
    );
    assert_eq!(l.traces.len(), 9);
    assert_eq!(diagonal_edges(&l), 4, "one 45° cut per bend");
    for p in &l.traces[5..] {
        assert_eq!(p.verts.len(), 5, "mitered corner is a pentagon");
    }
}

#[test]
fn cut_positions_sit_on_the_outer_corners() {
    let f = 0.7;
    let c = f * W;
    let l = double_jog(&substrate(), W, RUN, GAP, DY, MiterStyle::Mitered { f });
    let (xa, xb, dy) = (RUN, RUN + W + GAP, DY);
    // The outer corner point of each bend must NOT be a vertex of its
    // pentagon; the two cut endpoints must be, each `c` along an edge.
    let outers = [
        (xa + W, 0.0, xa + W - c, 0.0, xa + W, c),
        (xa, dy + W, xa + c, dy + W, xa, dy + W - c),
        (xb + W, dy + W, xb + W, dy + W - c, xb + W - c, dy + W),
        (xb, 0.0, xb, c, xb + c, 0.0),
    ];
    let has_vert = |p: &yee_layout::Polygon, x: f64, y: f64| {
        p.verts
            .iter()
            .any(|v| (v.x - x).abs() < 1e-12 && (v.y - y).abs() < 1e-12)
    };
    for (i, &(ox, oy, ax, ay, bx, by)) in outers.iter().enumerate() {
        let p = &l.traces[5 + i];
        assert!(!has_vert(p, ox, oy), "corner {i}: outer corner not cut");
        assert!(has_vert(p, ax, ay), "corner {i}: cut endpoint A missing");
        assert!(has_vert(p, bx, by), "corner {i}: cut endpoint B missing");
    }
}

#[test]
fn ports_face_x_at_equal_y_and_variants_share_bbox() {
    let sq = double_jog(&substrate(), W, RUN, GAP, DY, MiterStyle::Square);
    let mi = double_jog(
        &substrate(),
        W,
        RUN,
        GAP,
        DY,
        MiterStyle::Mitered { f: 0.7 },
    );
    assert_eq!(sq.ports.len(), 2);
    assert!((sq.ports[0].at.y - sq.ports[1].at.y).abs() < 1e-15);
    assert!((sq.ports[0].at.y - W / 2.0).abs() < 1e-15);
    // Mitering removes outer-corner metal only; the extremes are set by
    // the straights, so the bbox (and thus the voxel grid) is identical.
    assert_eq!(sq.bbox, mi.bbox);
    assert_eq!(sq.ports, mi.ports);
    // Port 2 sits just inside the far end.
    let end = 2.0 * RUN + 2.0 * W + GAP;
    assert!((sq.ports[1].at.x - (end - 0.5e-3)).abs() < 1e-15);
}

#[test]
fn segment_overlaps_stay_clear_of_the_cuts() {
    // The overlap rule: straights reach 0.2·w into each corner; the cut
    // spans f·w from the outer edges. For f = 0.7 the clearance is
    // (1 − 0.7 − 0.2)·w = 0.1·w > 0. Verify geometrically: no straight
    // rect contains any point strictly inside any cut triangle.
    let f = 0.7;
    let c = f * W;
    let l = double_jog(&substrate(), W, RUN, GAP, DY, MiterStyle::Mitered { f });
    let (xa, xb, dy) = (RUN, RUN + W + GAP, DY);
    // A representative interior point of each cut triangle (its centroid).
    let centroids = [
        Point2::new(xa + W - c / 3.0, c / 3.0),
        Point2::new(xa + c / 3.0, dy + W - c / 3.0),
        Point2::new(xb + W - c / 3.0, dy + W - c / 3.0),
        Point2::new(xb + c / 3.0, c / 3.0),
    ];
    for (i, pt) in centroids.iter().enumerate() {
        for (j, poly) in l.traces[..5].iter().enumerate() {
            let xs: Vec<f64> = poly.verts.iter().map(|v| v.x).collect();
            let ys: Vec<f64> = poly.verts.iter().map(|v| v.y).collect();
            let inside = pt.x > xs.iter().cloned().fold(f64::INFINITY, f64::min)
                && pt.x < xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                && pt.y > ys.iter().cloned().fold(f64::INFINITY, f64::min)
                && pt.y < ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            assert!(!inside, "straight {j} intrudes into cut {i}");
        }
    }
}
