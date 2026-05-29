//! # yee-export
//!
//! Manufacturing-file emitters for the Yee filter-design studio (Filter Phase
//! F1.4). This first brick is a **single-copper-layer RS-274X Gerber** emitter:
//! it turns a [`yee_layout::Layout`]'s top-metal polygons into filled
//! [Gerber][ucamco] regions (`G36*`/`G37*`).
//!
//! Pure text, no EM, no native dependency — **WASM-safe** so the studio can
//! export client-side (ADR-0089). The walking skeleton is deliberately minimal:
//! one copper layer, one aperture, no drill / board-outline / soldermask /
//! multi-layer (those are F1.4.1+).
//!
//! ## Coordinate model
//!
//! `yee-layout` stores all coordinates in **metres**. RS-274X here uses
//! millimetre units (`%MOMM*%`) with a 4-integer / 6-decimal fixed-point format
//! (`%FSLAX46Y46*%`). The conversion is therefore
//!
//! ```text
//! metres → mm:                  mm  = m * 1e3
//! mm → 4.6 fixed-point integer: int = round(mm * 1e6)
//! ```
//!
//! For example `3.0590 mm → 3059000 → X3059000`. See [`mm_to_fixed46`].
//!
//! ## Example
//!
//! ```
//! use yee_export::{layout_to_gerber, GerberOptions};
//! use yee_layout::{BBox, Layout, Polygon, Substrate};
//!
//! let traces = vec![Polygon::rect(0.0, 0.0, 1.0e-3, 0.5e-3)];
//! let layout = Layout {
//!     substrate: Substrate { eps_r: 4.4, height_m: 1.6e-3, loss_tangent: 0.02, metal_thickness_m: 35e-6 },
//!     bbox: BBox::from_polygons(&traces),
//!     traces,
//!     ports: vec![],
//! };
//! let gerber = layout_to_gerber(&layout, &GerberOptions::default());
//! assert!(gerber.starts_with("%FSLAX46Y46*%\n%MOMM*%\n"));
//! assert!(gerber.trim_end().ends_with("M02*"));
//! ```
//!
//! [ucamco]: https://www.ucamco.com/en/gerber — Ucamco, *The Gerber Layer
//! Format Specification* (RS-274X): region (`G36`/`G37`), format (`%FS…%`) and
//! units (`%MO…%`) statements.

use yee_layout::{Layout, Point2};

/// Options controlling the RS-274X Gerber emission.
///
/// The skeleton emits millimetre units in 4.6 fixed-point; the only tunable is
/// the layer name, which is written as a `G04` comment (Gerber attributes /
/// `%TF.FileFunction` are an F1.4.1+ refinement).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GerberOptions {
    /// Human-readable copper-layer name, emitted as a `G04 <layer_name>*`
    /// comment. Defaults to `"F.Cu"` (KiCad front-copper convention).
    pub layer_name: String,
}

impl Default for GerberOptions {
    fn default() -> Self {
        Self {
            layer_name: "F.Cu".into(),
        }
    }
}

/// Convert a coordinate in **millimetres** to a 4.6 fixed-point Gerber integer.
///
/// The `%FSLAX46Y46*%` format statement declares 4 integer + 6 decimal digits,
/// i.e. a resolution of `1e-6 mm` (1 nm). A value `mm` is therefore emitted as
/// the signed integer `round(mm * 1e6)`: `3.0590 mm → 3059000`. Leading-zero
/// omission (`LA` in the format word) means the integer is printed plainly with
/// no zero padding — Rust's default integer formatting already does this.
///
/// This is the exact inverse of the round-trip the `gerber-002` gate checks:
/// `fixed / 1e6` recovers the millimetre value to within the `1e-6 mm`
/// quantisation.
fn mm_to_fixed46(mm: f64) -> i64 {
    (mm * 1.0e6).round() as i64
}

/// Format one [`Point2`] (in metres) as a Gerber `X<int>Y<int>` coordinate
/// word, doing the metres → mm → 4.6 fixed-point conversion (see
/// [`mm_to_fixed46`]).
fn coord_word(p: &Point2) -> String {
    let ix = mm_to_fixed46(p.x * 1.0e3);
    let iy = mm_to_fixed46(p.y * 1.0e3);
    format!("X{ix}Y{iy}")
}

/// Emit a single-copper-layer **RS-274X Gerber** for a [`Layout`]'s top-metal
/// polygons, each as a filled region (`G36*` … `G37*`).
///
/// Structure (see the module docs for the coordinate model):
///
/// - header `%FSLAX46Y46*%` then `%MOMM*%`, a `G04 <layer_name>*` comment, one
///   aperture `%ADD10C,0.010*%` and `D10*` (a region needs a current aperture
///   selected even though the fill ignores its size);
/// - per polygon: `G36*`, a `D02*` move to the first vertex, a `D01*` draw to
///   each subsequent vertex, a closing `D01*` draw back to the first vertex if
///   the vertex list is not already explicitly closed, then `G37*`;
/// - footer `M02*`.
///
/// Coordinates come straight from the polygon vertices (metres), converted to
/// millimetre 4.6 fixed-point. Polygons with fewer than three vertices are
/// skipped (they have no fillable area).
///
/// The returned `String` is a complete, standalone Gerber file.
pub fn layout_to_gerber(layout: &Layout, opts: &GerberOptions) -> String {
    let mut s = String::new();

    // --- Header -------------------------------------------------------------
    // Format: absolute (A), leading-zero omission (L), X/Y = 4 integer + 6
    // decimal digits.
    s.push_str("%FSLAX46Y46*%\n");
    // Units = millimetres.
    s.push_str("%MOMM*%\n");
    // Layer/function note (informational comment).
    s.push_str(&format!("G04 {}*\n", opts.layer_name));
    // One circular aperture, 0.010 mm — regions ignore the aperture size but a
    // current aperture must be selected.
    s.push_str("%ADD10C,0.010*%\n");
    s.push_str("D10*\n");

    // --- Regions: one per polygon ------------------------------------------
    for poly in &layout.traces {
        let verts = &poly.verts;
        if verts.len() < 3 {
            // Degenerate polygon: no area to fill, skip.
            continue;
        }

        s.push_str("G36*\n");
        // Move (pen up) to the first vertex.
        s.push_str(&format!("{}D02*\n", coord_word(&verts[0])));
        // Draw (pen down) to each subsequent vertex.
        for v in &verts[1..] {
            s.push_str(&format!("{}D01*\n", coord_word(v)));
        }
        // Close the contour back to the first vertex unless it is already
        // explicitly closed (last vertex coincident with the first).
        if verts[verts.len() - 1] != verts[0] {
            s.push_str(&format!("{}D01*\n", coord_word(&verts[0])));
        }
        s.push_str("G37*\n");
    }

    // --- Footer -------------------------------------------------------------
    s.push_str("M02*\n");

    s
}
