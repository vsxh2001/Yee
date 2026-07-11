//! Gate `gerber-rt-001` (FS.3.0, ADR-0209): the Gerber importer is
//! lossless over the writer dialect, proven two ways on REAL generator
//! layouts (hairpin BPF, inset patch, quasi-Yagi, 2×1 array):
//!
//! 1. `import(export(L))` reproduces every polygon vertex-exactly (the
//!    4.6 fixed-point quantum is 1 nm; generator coordinates are far
//!    coarser, so equality is exact);
//! 2. **byte-stability**: `export(import(export(L))) == export(L)` — the
//!    house artifact philosophy (ADR-0198), and the strongest cheap proof
//!    the parse round-trips the dialect.
//!
//! Plus the explicit-rejection paths: imperial units, out-of-subset
//! commands, draw-before-move, unclosed regions. Instant, non-ignored.

use yee_export::{GerberImportError, GerberOptions, gerber_to_polygons, layout_to_gerber};
use yee_layout::{HairpinSectionParams, Layout, Substrate, hairpin_bpf_sections};
use yee_layout::{inset_fed_patch, patch_array_2x1, quasi_yagi};

fn fr4() -> Substrate {
    Substrate {
        eps_r: 4.4,
        height_m: 1.6e-3,
        loss_tangent: 0.0,
        metal_thickness_m: 35e-6,
    }
}

fn cases() -> Vec<(&'static str, Layout)> {
    let sub = fr4();
    vec![
        ("inset patch", inset_fed_patch(2.45e9, &sub, 50.0)),
        ("quasi-Yagi", quasi_yagi(5.8e9, &sub, 50.0).layout),
        ("2x1 array", patch_array_2x1(2.45e9, &sub, 50.0).layout),
        (
            "hairpin BPF",
            hairpin_bpf_sections(&HairpinSectionParams {
                substrate: sub,
                arm_length_m: 8e-3,
                line_width_m: 1e-3,
                fold_spacing_m: 3e-3,
                gaps_m: vec![0.4e-3, 0.9e-3],
                tap_offset_m: 2e-3,
                feed_width_m: 1e-3,
                feed_length_m: 8e-3,
            }),
        ),
    ]
}

#[test]
fn gerber_rt_001_vertex_exact_and_byte_stable() {
    for (name, layout) in cases() {
        let opts = GerberOptions::default();
        let gerber = layout_to_gerber(&layout, &opts);
        let imported = gerber_to_polygons(&gerber).unwrap_or_else(|e| panic!("{name}: {e}"));

        // 1. Vertex-exact reproduction.
        assert_eq!(imported.len(), layout.traces.len(), "{name}: polygon count");
        for (p_in, p_orig) in imported.iter().zip(&layout.traces) {
            assert_eq!(p_in.verts.len(), p_orig.verts.len(), "{name}: vertex count");
            for (a, b) in p_in.verts.iter().zip(&p_orig.verts) {
                assert!(
                    (a.x - b.x).abs() < 0.5e-9 && (a.y - b.y).abs() < 0.5e-9,
                    "{name}: vertex ({}, {}) vs original ({}, {})",
                    a.x,
                    a.y,
                    b.x,
                    b.y
                );
            }
        }

        // 2. Byte-stability through a full re-export.
        let relayout = Layout {
            substrate: layout.substrate,
            traces: imported,
            ports: layout.ports.clone(),
            bbox: layout.bbox,
        };
        let gerber2 = layout_to_gerber(&relayout, &opts);
        assert_eq!(
            gerber, gerber2,
            "{name}: export∘import∘export not byte-identical"
        );
    }
}

#[test]
fn out_of_subset_inputs_are_rejected_explicitly() {
    assert_eq!(
        gerber_to_polygons("%MOIN*%\nM02*\n"),
        Err(GerberImportError::ImperialUnits)
    );
    assert!(matches!(
        gerber_to_polygons("%LPD*%\nM02*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
    assert!(matches!(
        gerber_to_polygons("G36*\nX1Y1D01*\nG37*\nM02*\n"),
        Err(GerberImportError::DrawBeforeMove)
    ));
    assert_eq!(
        gerber_to_polygons("G36*\nX1Y1D02*\nX2Y2D01*\n"),
        Err(GerberImportError::UnclosedRegion)
    );
    // Arcs are FS.3.1+ — explicit rejection, not mis-parse.
    assert!(matches!(
        gerber_to_polygons("G02*\nM02*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
    // Stroked draws outside a region (the outline layer) are FS.3.1.
    assert!(matches!(
        gerber_to_polygons("X1Y1D02*\nX2Y2D01*\nM02*\n"),
        Err(GerberImportError::UnsupportedCommand(_))
    ));
}
