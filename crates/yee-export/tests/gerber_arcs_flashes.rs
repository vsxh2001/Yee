//! Gate `gerber-rt-003` (FS.3.2b, ADR-0220): circular-arc region
//! segments (G02/G03 under G75) and flashed C/R apertures (D03) import
//! to **analytically predictable** polygons:
//!
//! - arc endpoints and rect-flash corners are exact (≤ 0.5 nm, the 4.6
//!   fixed-point quantum idiom of `gerber-rt-001`);
//! - tessellation counts are pinned by the documented chord-tolerance
//!   formula (`ARC_CHORD_TOL_M` = 1 µm): quarter arc r = 1 mm → **18**
//!   segments, full-circle region r = 1 mm → **71** vertices, circle
//!   flash ⌀1 mm → **50** vertices;
//! - every synthesized vertex lies on the circle and the measured
//!   sagitta of every chord stays ≤ the tolerance;
//! - everything still outside the subset is rejected **by name**
//!   (single-quadrant G74, arcs without G75, flashes in regions,
//!   unknown/unsupported apertures, broken arc geometry — plus the
//!   FS.3.0 rejections re-asserted: inches, polarity, stroked draws).
//!
//! Instant, non-ignored.

use yee_export::{ARC_CHORD_TOL_M, GerberImportError, gerber_to_polygons};
use yee_layout::Point2;

/// The standard writer-dialect header the importer accepts.
const HEADER: &str = "%FSLAX46Y46*%\n%MOMM*%\nG04 F.Cu*\n%ADD10C,0.010*%\nD10*\n";

fn parse(body: &str) -> Vec<yee_layout::Polygon> {
    let file = format!("{HEADER}{body}M02*\n");
    gerber_to_polygons(&file).unwrap_or_else(|e| panic!("{e}\n---\n{file}"))
}

/// Exact-coordinate assert at the 4.6 quantum (same idiom as
/// `gerber-rt-001`).
fn assert_exact(p: &Point2, x: f64, y: f64, what: &str) {
    assert!(
        (p.x - x).abs() < 0.5e-9 && (p.y - y).abs() < 0.5e-9,
        "{what}: ({}, {}) vs expected ({x}, {y})",
        p.x,
        p.y
    );
}

/// Check consecutive arc vertices: on-circle to ≤ 1 nm and chord
/// sagitta ≤ `ARC_CHORD_TOL_M` (the pinned tessellation contract).
fn assert_on_circle_and_sagitta(verts: &[Point2], cx: f64, cy: f64, r: f64, what: &str) {
    for v in verts {
        let d = (v.x - cx).hypot(v.y - cy);
        assert!((d - r).abs() < 1.0e-9, "{what}: off-circle by {}", d - r);
    }
    for w in verts.windows(2) {
        let (mx, my) = (0.5 * (w[0].x + w[1].x), 0.5 * (w[0].y + w[1].y));
        let sag = r - (mx - cx).hypot(my - cy);
        assert!(
            sag <= ARC_CHORD_TOL_M + 1.0e-12,
            "{what}: sagitta {sag} > {ARC_CHORD_TOL_M}"
        );
    }
}

#[test]
fn quarter_arc_ccw_pinned_tessellation() {
    // r = 1 mm CCW quarter about the origin: (1, 0) → (0, 1) mm, then
    // two straight edges through the centre close the wedge.
    // NB: interpolation mode is modal — the straight closing edges need
    // an explicit G01 after the arc.
    let polys = parse(
        "G75*\nG36*\nX1000000Y0D02*\nG03X0Y1000000I-1000000J0D01*\nG01*\nX0Y0D01*\nX1000000Y0D01*\nG37*\n",
    );
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    // 1 start + 18 pinned arc segments (17 interior + exact end) + the
    // centre vertex; the explicit closing vertex is dropped.
    assert_eq!(v.len(), 20, "quarter arc: pinned n = 18 segments");
    assert_exact(&v[0], 1.0e-3, 0.0, "quarter arc start");
    assert_exact(&v[18], 0.0, 1.0e-3, "quarter arc end (exact)");
    assert_exact(&v[19], 0.0, 0.0, "wedge centre");
    assert_on_circle_and_sagitta(&v[0..=18], 0.0, 0.0, 1.0e-3, "quarter arc");
    // CCW: polar angle strictly increases along the arc.
    for w in v[0..=18].windows(2) {
        assert!(
            w[1].y.atan2(w[1].x) > w[0].y.atan2(w[0].x),
            "quarter arc: not CCW"
        );
    }
}

#[test]
fn quarter_arc_cw_mirrors_ccw() {
    // Same wedge traced the other way: (0, 1) → (1, 0) mm clockwise.
    let polys = parse(
        "G75*\nG36*\nX0Y1000000D02*\nG02X1000000Y0I0J-1000000D01*\nG01*\nX0Y0D01*\nX0Y1000000D01*\nG37*\n",
    );
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    assert_eq!(v.len(), 20, "CW quarter arc: pinned n = 18 segments");
    assert_exact(&v[0], 0.0, 1.0e-3, "CW arc start");
    assert_exact(&v[18], 1.0e-3, 0.0, "CW arc end (exact)");
    assert_on_circle_and_sagitta(&v[0..=18], 0.0, 0.0, 1.0e-3, "CW quarter arc");
    // CW: polar angle strictly decreases along the arc.
    for w in v[0..=18].windows(2) {
        assert!(
            w[1].y.atan2(w[1].x) < w[0].y.atan2(w[0].x),
            "CW quarter arc: not CW"
        );
    }
}

#[test]
fn full_circle_region_start_equals_end() {
    // start == end in multi-quadrant mode is a full 360° circle
    // (Ucamco G75 rule): r = 1 mm about (2, 0) mm.
    let polys = parse("G75*\nG36*\nX3000000Y0D02*\nG03X3000000Y0I-1000000J0D01*\nG37*\n");
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    // Pinned n = 71 segments; the coincident closing vertex is dropped.
    assert_eq!(v.len(), 71, "full circle: pinned n = 71");
    assert_exact(&v[0], 3.0e-3, 0.0, "full circle start");
    let mut closed = v.clone();
    closed.push(v[0]);
    assert_on_circle_and_sagitta(&closed, 2.0e-3, 0.0, 1.0e-3, "full circle");
}

#[test]
fn rect_flash_is_exact() {
    // R,2X1 flashed at (5, 5) mm: exactly its 4 half-extent corners,
    // CCW from lower-left.
    let polys = parse("%ADD11R,2X1*%\nD11*\nX5000000Y5000000D03*\n");
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    assert_eq!(v.len(), 4);
    assert_exact(&v[0], 4.0e-3, 4.5e-3, "rect corner 0");
    assert_exact(&v[1], 6.0e-3, 4.5e-3, "rect corner 1");
    assert_exact(&v[2], 6.0e-3, 5.5e-3, "rect corner 2");
    assert_exact(&v[3], 4.0e-3, 5.5e-3, "rect corner 3");
}

#[test]
fn circle_flash_pinned_tessellation() {
    // C,1 (r = 0.5 mm) flashed at the origin: pinned n = 50 vertices,
    // vertex 0 on the +x axis, CCW, on-circle + sagitta bound (with
    // the wrap-around chord included).
    let polys = parse("%ADD12C,1*%\nD12*\nX0Y0D03*\n");
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    assert_eq!(v.len(), 50, "circle flash: pinned n = 50");
    assert_exact(&v[0], 0.5e-3, 0.0, "circle flash vertex 0");
    let mut closed = v.clone();
    closed.push(v[0]);
    assert_on_circle_and_sagitta(&closed, 0.0, 0.0, 0.5e-3, "circle flash");
    assert!(v[1].y > 0.0, "circle flash: not CCW");
}

#[test]
fn regions_and_flashes_interleave_in_file_order() {
    let polys = parse(
        "%ADD11R,1X1*%\n\
         G36*\nX0Y0D02*\nX1000000Y0D01*\nX1000000Y1000000D01*\nX0Y0D01*\nG37*\n\
         D11*\nX5000000Y5000000D03*\n\
         G36*\nX0Y2000000D02*\nX1000000Y2000000D01*\nX1000000Y3000000D01*\nX0Y2000000D01*\nG37*\n",
    );
    assert_eq!(polys.len(), 3, "region, flash, region — in file order");
    assert_eq!(polys[0].verts.len(), 3);
    assert_eq!(polys[1].verts.len(), 4, "the middle polygon is the flash");
    assert_exact(&polys[1].verts[0], 4.5e-3, 4.5e-3, "flash corner 0");
    assert_eq!(polys[2].verts.len(), 3);
}

#[test]
fn unsupported_constructs_stay_named() {
    let g = |body: &str| gerber_to_polygons(&format!("{HEADER}{body}M02*\n"));

    // Single-quadrant arc mode: legacy, ambiguous — rejected by name.
    assert!(matches!(
        g("G74*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
    // An arc draw before G75 would have single-quadrant semantics.
    assert!(matches!(
        g("G36*\nX1000000Y0D02*\nG03X0Y1000000I-1000000J0D01*\nG37*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
    // Flashes are forbidden inside regions by the Gerber spec.
    assert_eq!(
        g("%ADD11R,1X1*%\nD11*\nG36*\nX0Y0D02*\nX1000000Y0D03*\nG37*\n"),
        Err(GerberImportError::FlashInRegion)
    );
    // Flash with no aperture ever selected (the shared HEADER selects
    // D10, so this case uses a bare header) / an undefined D-code.
    assert!(matches!(
        gerber_to_polygons("%FSLAX46Y46*%\n%MOMM*%\nX0Y0D03*\nM02*\n"),
        Err(GerberImportError::UnknownAperture(_))
    ));
    assert!(matches!(
        g("D99*\nX0Y0D03*\n"),
        Err(GerberImportError::UnknownAperture(_))
    ));
    // Obround and holed apertures are bookkept but not flashable.
    assert!(matches!(
        g("%ADD13O,1X0.5*%\nD13*\nX0Y0D03*\n"),
        Err(GerberImportError::UnsupportedAperture(_))
    ));
    assert!(matches!(
        g("%ADD14C,1X0.5*%\nD14*\nX0Y0D03*\n"),
        Err(GerberImportError::UnsupportedAperture(_))
    ));
    // Broken arc geometry: zero radius (I0J0) and mismatched radii.
    assert!(matches!(
        g("G75*\nG36*\nX0Y0D02*\nG03X1000000Y0I0J0D01*\nG37*\n"),
        Err(GerberImportError::BadArc(_))
    ));
    assert!(matches!(
        g("G75*\nG36*\nX1000000Y0D02*\nG03X2000000Y0I-1000000J0D01*\nG37*\n"),
        Err(GerberImportError::BadArc(_))
    ));
    // A stroked ARC draw outside a region is still a stroked draw.
    assert!(matches!(
        g("G75*\nX1000000Y0D02*\nG02X0Y1000000I-1000000J0D01*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));

    // FS.3.0 rejections re-asserted so the subset boundary cannot
    // silently widen: inches, polarity, stroked linear draws.
    assert_eq!(
        gerber_to_polygons("%MOIN*%\nM02*\n"),
        Err(GerberImportError::ImperialUnits)
    );
    assert!(matches!(
        gerber_to_polygons("%LPD*%\nM02*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
    assert!(matches!(
        gerber_to_polygons("X1Y1D02*\nX2Y2D01*\nM02*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
}
