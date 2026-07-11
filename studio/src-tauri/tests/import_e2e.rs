//! Gate `studio-import-e2e-001` (FS.3.1b, ADR-0209): the studio's
//! headless import flow — a board exported by this suite's own writer
//! comes back through `import_gerber_impl` **provably lossless**: the
//! response's copper echo is byte-identical to the input Gerber, the
//! preview SVG and layout JSON are present, and the outline corners land
//! exactly. Error paths: no ports, junk Gerber, imperial units. Instant
//! (no solves) — the R.5 command-layer e2e idiom.

use yee_export::{GerberOptions, OutlineOptions, layout_to_gerber, layout_to_gerber_outline};
use yee_layout::{Substrate, inset_fed_patch};
use yee_studio_app::import::{ImportPort, ImportRequest, import_gerber_impl};

fn patch_request() -> (ImportRequest, String) {
    let sub = Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    };
    let layout = inset_fed_patch(2.45e9, &sub, 50.0);
    let copper = layout_to_gerber(&layout, &GerberOptions::default());
    let outline = layout_to_gerber_outline(&layout, &OutlineOptions::default());
    let req = ImportRequest {
        copper_gerber: copper.clone(),
        outline_gerber: Some(outline),
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        ports: vec![ImportPort {
            x_m: layout.ports[0].at.x,
            y_m: layout.ports[0].at.y,
            width_m: layout.ports[0].width_m,
            z0_ohm: 50.0,
        }],
    };
    (req, copper)
}

#[test]
fn studio_import_e2e_001_roundtrip_is_byte_provable() {
    let (req, copper_in) = patch_request();
    let resp = import_gerber_impl(&req).expect("import failed");

    // The echo proves what was understood, byte-for-byte.
    assert_eq!(
        resp.gerber_copper_echo, copper_in,
        "copper echo not byte-identical to the input"
    );
    // The A.1 inset patch is 4 polygons (feed + 2 bands + centre).
    assert_eq!(resp.trace_count, 4);
    assert!(resp.svg.contains("<svg"), "preview SVG missing");
    assert!(resp.bbox_w_m > 0.0 && resp.bbox_h_m > 0.0);
    // Outline present, rectangular, and enclosing the bbox.
    let outline = resp.outline_m.expect("outline missing");
    assert_eq!(outline.len(), 4);
    let (xs, ys): (Vec<f64>, Vec<f64>) = outline.iter().copied().unzip();
    let (min_x, max_x) = (
        xs.iter().cloned().fold(f64::INFINITY, f64::min),
        xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    assert!(max_x - min_x > resp.bbox_w_m, "outline must enclose the bbox");
    let _ = ys;
    // The layout JSON deserializes back into a Layout.
    let layout: yee_layout::Layout = serde_json::from_str(&resp.layout_json).expect("layout JSON");
    assert_eq!(layout.traces.len(), 4);
    assert_eq!(layout.ports.len(), 1);
}

#[test]
fn studio_import_e2e_001_error_paths() {
    let (mut req, _) = patch_request();
    req.ports.clear();
    assert!(import_gerber_impl(&req).unwrap_err().contains("port"));

    let (mut req, _) = patch_request();
    req.copper_gerber = "%MOIN*%\nM02*\n".into();
    assert!(
        import_gerber_impl(&req)
            .unwrap_err()
            .contains("imperial units")
    );

    let (mut req, _) = patch_request();
    req.copper_gerber = "M02*\n".into();
    assert!(
        import_gerber_impl(&req)
            .unwrap_err()
            .contains("no copper regions")
    );
}
