//! Gate `dxf-rt-001` (FS.3.3, ADR-0230): the DXF importer reproduces the
//! S.6 stub-notch trace geometry vertex-exactly (same 0.5 nm tolerance
//! idiom as `gerber-rt-001`), plus bulge-arc tessellation pinned to the
//! **identical** tessellation the Gerber importer already gates
//! (`gerber-rt-003`'s quarter-arc wedge, r = 1 mm, n = 18) — expected,
//! since [`yee_export::dxf`] reuses
//! [`crate::import::arc_vertices`][arcv]'s angle-stepping loop, not a
//! reimplementation.
//!
//! Plus the full named-rejection matrix: open polylines, every
//! unsupported entity kind, nonzero elevation, and `$INSUNITS` outside
//! mm/inch (including missing). Instant, non-ignored.
//!
//! [arcv]: yee_export::import

use yee_export::{ARC_CHORD_TOL_M, DxfImportError, DxfOptions, dxf_to_outline};
use yee_layout::{Point2, Polygon, eps_eff, open_end_delta_l};

const EPS_R: f64 = 4.4;
const H_M: f64 = 1.6e-3;
const W_M: f64 = 3.0e-3;
const F0_HZ: f64 = 5.0e9;
const C0_M_S: f64 = 299_792_458.0;

/// The native S.6 stub-notch trace geometry (feed line + open stub),
/// identical to `import_twin.rs`'s `native_stub_layout()` traces —
/// reusing `yee_layout`'s own public Hammerstad helpers rather than
/// re-deriving them.
fn native_stub_traces() -> (Polygon, Polygon) {
    let e_eff = eps_eff(W_M, H_M, EPS_R);
    let lam_g = C0_M_S / (F0_HZ * e_eff.sqrt());
    let l_m = 3.0 * lam_g;
    let stub_len = lam_g / 4.0 - open_end_delta_l(W_M, H_M, e_eff);
    let line = Polygon::rect(0.0, 0.0, l_m, W_M);
    let stub = Polygon::rect(l_m / 2.0 - W_M / 2.0, W_M, W_M, stub_len);
    (line, stub)
}

/// Emit a closed `LWPOLYLINE` DXF entity for an axis-aligned rectangle
/// (metres in, millimetres out — `$INSUNITS 4`), corners in the same CCW
/// order [`Polygon::rect`] uses so vertex-for-vertex comparison needs no
/// reordering.
fn lwpolyline_rect(x0_m: f64, y0_m: f64, w_m: f64, h_m: f64) -> String {
    let corners = [
        (x0_m, y0_m),
        (x0_m + w_m, y0_m),
        (x0_m + w_m, y0_m + h_m),
        (x0_m, y0_m + h_m),
    ];
    let mut s = String::from("0\nLWPOLYLINE\n8\n0\n90\n4\n70\n1\n");
    for (x, y) in corners {
        s.push_str(&format!("10\n{:.15e}\n20\n{:.15e}\n", x * 1.0e3, y * 1.0e3));
    }
    s
}

/// Wrap a raw `ENTITIES`-section body in the minimal `HEADER` (mm units)
/// + `ENTITIES` + `EOF` DXF envelope this importer accepts.
fn dxf_file(entities: &str) -> String {
    format!(
        "0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n4\n0\nENDSEC\n\
         0\nSECTION\n2\nENTITIES\n{entities}0\nENDSEC\n0\nEOF\n"
    )
}

fn assert_exact(p: &Point2, x: f64, y: f64, what: &str) {
    assert!(
        (p.x - x).abs() < 0.5e-9 && (p.y - y).abs() < 0.5e-9,
        "{what}: ({}, {}) vs expected ({x}, {y})",
        p.x,
        p.y
    );
}

/// Same idiom as `gerber_arcs_flashes.rs`: every vertex within
/// [`ARC_CHORD_TOL_M`] of the true circle, sagitta of every chord no
/// larger than the tolerance.
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
fn dxf_rt_001_vertex_exact_vs_native_stub() {
    let (line, stub) = native_stub_traces();
    let Point2 { x: lx0, y: ly0 } = line.verts[0];
    let l_m = line.verts[1].x - lx0;
    let stub_x0 = stub.verts[0].x;
    let stub_y0 = stub.verts[0].y;
    let stub_w = stub.verts[1].x - stub_x0;
    let stub_h = stub.verts[3].y - stub_y0;

    let entities = format!(
        "{}{}",
        lwpolyline_rect(lx0, ly0, l_m, W_M),
        lwpolyline_rect(stub_x0, stub_y0, stub_w, stub_h),
    );
    let dxf = dxf_file(&entities);

    let polys = dxf_to_outline(&dxf, &DxfOptions::default())
        .unwrap_or_else(|e| panic!("dxf-rt-001: import failed: {e}\n---\n{dxf}"));

    assert_eq!(polys.len(), 2, "dxf-rt-001: polygon count");
    for (p_in, p_native) in polys.iter().zip([&line, &stub]) {
        assert_eq!(
            p_in.verts.len(),
            p_native.verts.len(),
            "dxf-rt-001: vertex count"
        );
        for (a, b) in p_in.verts.iter().zip(&p_native.verts) {
            assert!(
                (a.x - b.x).abs() < 0.5e-9 && (a.y - b.y).abs() < 0.5e-9,
                "dxf-rt-001: vertex ({}, {}) vs native ({}, {})",
                a.x,
                a.y,
                b.x,
                b.y
            );
        }
    }
}

/// Same wedge as `gerber_arcs_flashes.rs`'s `quarter_arc_ccw_pinned_tessellation`
/// (r = 1 mm quarter about the origin, closed through the centre) but built
/// from a DXF bulge instead of a Gerber G03 arc: proves the two importers'
/// tessellation is bit-for-bit identical (both call
/// [`yee_export::import::arc_vertices`]).
#[test]
fn bulge_ccw_quarter_matches_gerber_pinned_tessellation() {
    let bulge = (std::f64::consts::FRAC_PI_8).tan();
    // (1,0) -mm-> bulge arc -> (0,1) -straight-> (0,0) -straight-> close.
    let entity = format!(
        "0\nLWPOLYLINE\n8\n0\n90\n3\n70\n1\n\
         10\n1.0\n20\n0.0\n42\n{bulge}\n\
         10\n0.0\n20\n1.0\n\
         10\n0.0\n20\n0.0\n"
    );
    let dxf = dxf_file(&entity);
    let polys =
        dxf_to_outline(&dxf, &DxfOptions::default()).unwrap_or_else(|e| panic!("{e}\n{dxf}"));
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    // Identical structure to the Gerber wedge: 1 start + 17 interior + 1
    // exact end + 1 centre = 20.
    assert_eq!(v.len(), 20, "CCW bulge quarter: pinned n = 18 segments");
    assert_exact(&v[0], 1.0e-3, 0.0, "CCW bulge start");
    assert_exact(&v[18], 0.0, 1.0e-3, "CCW bulge end (exact)");
    assert_exact(&v[19], 0.0, 0.0, "wedge centre");
    assert_on_circle_and_sagitta(&v[0..=18], 0.0, 0.0, 1.0e-3, "CCW bulge quarter");
    for w in v[0..=18].windows(2) {
        assert!(
            w[1].y.atan2(w[1].x) > w[0].y.atan2(w[0].x),
            "CCW bulge quarter: not CCW"
        );
    }
}

/// Same wedge traced the other way (CW): `(0,1) -> (1,0)` with a negative
/// bulge, mirroring `gerber_arcs_flashes.rs`'s `quarter_arc_cw_mirrors_ccw`.
#[test]
fn bulge_cw_quarter_matches_gerber_pinned_tessellation() {
    let bulge = -(std::f64::consts::FRAC_PI_8).tan();
    let entity = format!(
        "0\nLWPOLYLINE\n8\n0\n90\n3\n70\n1\n\
         10\n0.0\n20\n1.0\n42\n{bulge}\n\
         10\n1.0\n20\n0.0\n\
         10\n0.0\n20\n0.0\n"
    );
    let dxf = dxf_file(&entity);
    let polys =
        dxf_to_outline(&dxf, &DxfOptions::default()).unwrap_or_else(|e| panic!("{e}\n{dxf}"));
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    assert_eq!(v.len(), 20, "CW bulge quarter: pinned n = 18 segments");
    assert_exact(&v[0], 0.0, 1.0e-3, "CW bulge start");
    assert_exact(&v[18], 1.0e-3, 0.0, "CW bulge end (exact)");
    assert_exact(&v[19], 0.0, 0.0, "wedge centre");
    assert_on_circle_and_sagitta(&v[0..=18], 0.0, 0.0, 1.0e-3, "CW bulge quarter");
    for w in v[0..=18].windows(2) {
        assert!(
            w[1].y.atan2(w[1].x) < w[0].y.atan2(w[0].x),
            "CW bulge quarter: not CW"
        );
    }
}

#[test]
fn polyline_vertex_chain_parses_closed_rectangle() {
    // R12 fallback: POLYLINE + VERTEX + SEQEND, no bulges.
    let entity = "0\nPOLYLINE\n8\n0\n70\n1\n\
         0\nVERTEX\n8\n0\n10\n0.0\n20\n0.0\n\
         0\nVERTEX\n8\n0\n10\n2.0\n20\n0.0\n\
         0\nVERTEX\n8\n0\n10\n2.0\n20\n1.0\n\
         0\nVERTEX\n8\n0\n10\n0.0\n20\n1.0\n\
         0\nSEQEND\n";
    let dxf = dxf_file(entity);
    let polys =
        dxf_to_outline(&dxf, &DxfOptions::default()).unwrap_or_else(|e| panic!("{e}\n{dxf}"));
    assert_eq!(polys.len(), 1);
    let v = &polys[0].verts;
    assert_eq!(v.len(), 4);
    assert_exact(&v[0], 0.0, 0.0, "POLYLINE corner 0");
    assert_exact(&v[1], 2.0e-3, 0.0, "POLYLINE corner 1");
    assert_exact(&v[2], 2.0e-3, 1.0e-3, "POLYLINE corner 2");
    assert_exact(&v[3], 0.0, 1.0e-3, "POLYLINE corner 3");
}

#[test]
fn layer_filter_skips_non_matching_layers() {
    let entities = format!(
        "{}{}",
        {
            let mut e = lwpolyline_rect(0.0, 0.0, 1.0e-3, 1.0e-3);
            e = e.replacen("8\n0\n", "8\nkeep\n", 1);
            e
        },
        {
            let mut e = lwpolyline_rect(5.0e-3, 5.0e-3, 1.0e-3, 1.0e-3);
            e = e.replacen("8\n0\n", "8\ndrop\n", 1);
            e
        }
    );
    let dxf = dxf_file(&entities);
    let opts = DxfOptions {
        layer: Some("keep".into()),
    };
    let polys = dxf_to_outline(&dxf, &opts).unwrap_or_else(|e| panic!("{e}\n{dxf}"));
    assert_eq!(polys.len(), 1, "layer filter must drop the other layer");
    assert_exact(&polys[0].verts[0], 0.0, 0.0, "kept layer corner 0");
}

#[test]
fn out_of_subset_inputs_are_rejected_explicitly() {
    // $INSUNITS missing entirely.
    let no_units = "0\nSECTION\n2\nENTITIES\n0\nENDSEC\n0\nEOF\n";
    assert!(matches!(
        dxf_to_outline(no_units, &DxfOptions::default()),
        Err(DxfImportError::UnsupportedUnits(v)) if v == "missing"
    ));

    // $INSUNITS explicitly unitless (0) or an unsupported value (2 = feet).
    let bad_units = |code: i32| {
        format!(
            "0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n{code}\n0\nENDSEC\n\
             0\nSECTION\n2\nENTITIES\n0\nENDSEC\n0\nEOF\n"
        )
    };
    assert!(matches!(
        dxf_to_outline(&bad_units(0), &DxfOptions::default()),
        Err(DxfImportError::UnsupportedUnits(v)) if v == "0"
    ));
    assert!(matches!(
        dxf_to_outline(&bad_units(2), &DxfOptions::default()),
        Err(DxfImportError::UnsupportedUnits(v)) if v == "2"
    ));

    // Open polyline (closed flag unset).
    let open = "0\nLWPOLYLINE\n8\n0\n90\n3\n70\n0\n\
         10\n0.0\n20\n0.0\n10\n1.0\n20\n0.0\n10\n1.0\n20\n1.0\n";
    assert_eq!(
        dxf_to_outline(&dxf_file(open), &DxfOptions::default()),
        Err(DxfImportError::OpenPolyline)
    );

    // Nonzero elevation (LWPOLYLINE constant elevation, group 38).
    let elevated = "0\nLWPOLYLINE\n8\n0\n90\n3\n70\n1\n38\n1.0\n\
         10\n0.0\n20\n0.0\n10\n1.0\n20\n0.0\n10\n1.0\n20\n1.0\n";
    assert_eq!(
        dxf_to_outline(&dxf_file(elevated), &DxfOptions::default()),
        Err(DxfImportError::NonzeroElevation)
    );

    // Nonzero Z on a POLYLINE/VERTEX chain.
    let elevated_vertex = "0\nPOLYLINE\n8\n0\n70\n1\n\
         0\nVERTEX\n8\n0\n10\n0.0\n20\n0.0\n30\n0.0\n\
         0\nVERTEX\n8\n0\n10\n1.0\n20\n0.0\n30\n2.0\n\
         0\nVERTEX\n8\n0\n10\n1.0\n20\n1.0\n30\n0.0\n\
         0\nSEQEND\n";
    assert_eq!(
        dxf_to_outline(&dxf_file(elevated_vertex), &DxfOptions::default()),
        Err(DxfImportError::NonzeroElevation)
    );

    // POLYLINE never reaching SEQEND.
    let unclosed = "0\nPOLYLINE\n8\n0\n70\n1\n\
         0\nVERTEX\n8\n0\n10\n0.0\n20\n0.0\n";
    assert_eq!(
        dxf_to_outline(&dxf_file(unclosed), &DxfOptions::default()),
        Err(DxfImportError::UnclosedPolyline)
    );

    // No closed polylines at all.
    let empty = "0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n4\n0\nENDSEC\n\
         0\nSECTION\n2\nENTITIES\n0\nENDSEC\n0\nEOF\n";
    assert_eq!(
        dxf_to_outline(empty, &DxfOptions::default()),
        Err(DxfImportError::NoOutline)
    );

    // The named-entity rejection matrix: CIRCLE, ARC, ELLIPSE, SPLINE,
    // TEXT, INSERT — one typed error each, all funnelling into
    // `UnsupportedEntity`.
    for (name, body) in [
        ("CIRCLE", "0\nCIRCLE\n8\n0\n10\n0.0\n20\n0.0\n40\n1.0\n"),
        (
            "ARC",
            "0\nARC\n8\n0\n10\n0.0\n20\n0.0\n40\n1.0\n50\n0.0\n51\n90.0\n",
        ),
        (
            "ELLIPSE",
            "0\nELLIPSE\n8\n0\n10\n0.0\n20\n0.0\n11\n1.0\n21\n0.0\n40\n0.5\n",
        ),
        ("SPLINE", "0\nSPLINE\n8\n0\n70\n0\n"),
        ("TEXT", "0\nTEXT\n8\n0\n10\n0.0\n20\n0.0\n40\n1.0\n1\nhi\n"),
        ("INSERT", "0\nINSERT\n8\n0\n2\nBLOCK1\n10\n0.0\n20\n0.0\n"),
    ] {
        let dxf = dxf_file(body);
        assert_eq!(
            dxf_to_outline(&dxf, &DxfOptions::default()),
            Err(DxfImportError::UnsupportedEntity(name.into())),
            "{name} must be a named rejection"
        );
    }
}
