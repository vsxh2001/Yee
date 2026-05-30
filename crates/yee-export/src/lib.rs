//! # yee-export
//!
//! Manufacturing-file emitters for the Yee filter-design studio (Filter Phase
//! F1.4). Two RS-274X [Gerber][ucamco] emitters and one KiCad board emitter are
//! provided:
//!
//! - [`layout_to_gerber`] — single-copper-layer: a [`yee_layout::Layout`]'s
//!   top-metal polygons as filled regions (`G36*`/`G37*`) — F1.4.0.
//! - [`layout_to_gerber_outline`] — board-outline (Edge.Cuts): a single
//!   **stroked** rectangular contour around the layout `bbox`, expanded by a
//!   margin — F1.4.1a.
//! - [`layout_to_kicad_pcb`] — a KiCad 7 `.kicad_pcb` S-expression board file:
//!   the top-metal polygons as filled `gr_poly` on the copper layer plus the
//!   board outline on `Edge.Cuts`, so the user can open the layout directly in
//!   the KiCad PCB editor — F1.4.1b.
//!
//! Pure text, no EM, no native dependency — **WASM-safe** so the studio can
//! export client-side (ADR-0089). Drill, soldermask, silkscreen, footprints,
//! pads, vias, zones, and multi-layer stack-ups are F1.4.1c+.
//!
//! ## Coordinate model
//!
//! `yee-layout` stores all coordinates in **metres**. The two coordinate
//! encodings used by the emitters differ:
//!
//! - **Gerber** (`layout_to_gerber` / `layout_to_gerber_outline`) uses
//!   millimetre units (`%MOMM*%`) with a 4-integer / 6-decimal fixed-point
//!   format (`%FSLAX46Y46*%`):
//!
//!   ```text
//!   metres → mm:                  mm  = m * 1e3
//!   mm → 4.6 fixed-point integer: int = round(mm * 1e6)
//!   ```
//!
//!   For example `3.0590 mm → 3059000 → X3059000`. See [`mm_to_fixed46`].
//!
//! - **KiCad** (`layout_to_kicad_pcb`) uses millimetre **floats** directly in
//!   the S-expression (`(xy 3.059000 1.234000)`) — *not* fixed-point integers.
//!   See [`xy_mm`].
//!
//! Neither emitter flips the `y` axis in this walking-skeleton scope; the
//! studio, Gerber, and KiCad all share the layout's metre frame scaled to mm.
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
    coord_word_xy(p.x, p.y)
}

/// Format a raw `(x_m, y_m)` coordinate pair (metres) as a Gerber
/// `X<int>Y<int>` word, doing the metres → mm → 4.6 fixed-point conversion
/// (see [`mm_to_fixed46`]).
///
/// This is the `Point2`-free variant used by the board-outline emitter, whose
/// corners are computed from `bbox ± margin` rather than read from existing
/// vertices. [`coord_word`] delegates to it so both paths share one
/// conversion.
fn coord_word_xy(x_m: f64, y_m: f64) -> String {
    let ix = mm_to_fixed46(x_m * 1.0e3);
    let iy = mm_to_fixed46(y_m * 1.0e3);
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

/// Options controlling the board-outline (Edge.Cuts) Gerber emission.
///
/// The outline is a single rectangular contour around the layout bounding box
/// expanded by [`margin_mm`](Self::margin_mm) on each side. The
/// [`layer_name`](Self::layer_name) is written as a `G04` comment.
#[derive(Debug, Clone, PartialEq)]
pub struct OutlineOptions {
    /// Human-readable outline-layer name, emitted as a `G04 <layer_name>*`
    /// comment. Defaults to `"Edge.Cuts"` (KiCad board-profile convention).
    pub layer_name: String,
    /// Margin added on *each* side of the layout bounding box, in
    /// millimetres. The emitted rectangle therefore spans `bbox` grown by this
    /// amount in `±x` and `±y`. Defaults to `1.0` mm.
    pub margin_mm: f64,
}

impl Default for OutlineOptions {
    fn default() -> Self {
        Self {
            layer_name: "Edge.Cuts".into(),
            margin_mm: 1.0,
        }
    }
}

/// Emit a board-outline (Edge.Cuts) **RS-274X Gerber** for a [`Layout`]: a
/// single closed rectangular contour around the layout `bbox`, expanded by
/// `opts.margin_mm` on each side, **stroked** with a thin aperture (it is a
/// cut path / profile, *not* a region fill — there is no `G36*`/`G37*`).
///
/// Structure (see the module docs for the coordinate model):
///
/// - header `%FSLAX46Y46*%` then `%MOMM*%`, a `G04 <layer_name>*` comment, one
///   thin circular aperture `%ADD10C,0.100*%` (0.1 mm stroke) and `D10*`;
/// - the four rectangle corners derived from `layout.bbox` expanded by
///   `margin_m = margin_mm * 1e-3` metres on each side, in CCW order starting
///   from the lower-left:
///   `(min.x−m, min.y−m)`, `(max.x+m, min.y−m)`, `(max.x+m, max.y+m)`,
///   `(min.x−m, max.y+m)`;
/// - a `D02*` pen-up move to corner 0, then `D01*` pen-down draws to corners
///   1, 2, 3, and a final `D01*` draw back to corner 0 (explicit close);
/// - footer `M02*`.
///
/// Coordinates are converted from metres to millimetre 4.6 fixed-point with
/// the same conversion as [`layout_to_gerber`] (see [`mm_to_fixed46`]).
///
/// The returned `String` is a complete, standalone Gerber file.
pub fn layout_to_gerber_outline(layout: &Layout, opts: &OutlineOptions) -> String {
    let mut s = String::new();

    // --- Header -------------------------------------------------------------
    // Format: absolute (A), leading-zero omission (L), X/Y = 4 integer + 6
    // decimal digits.
    s.push_str("%FSLAX46Y46*%\n");
    // Units = millimetres.
    s.push_str("%MOMM*%\n");
    // Layer/function note (informational comment).
    s.push_str(&format!("G04 {}*\n", opts.layer_name));
    // One thin circular aperture, 0.100 mm — this is a stroked cut path, so
    // the aperture *width* is meaningful (it is the routing tool's trace).
    s.push_str("%ADD10C,0.100*%\n");
    s.push_str("D10*\n");

    // --- Outline contour ----------------------------------------------------
    // Rectangle corners = bbox expanded by `margin_m` on each side, CCW from
    // the lower-left corner.
    let m = opts.margin_mm * 1.0e-3;
    let (min, max) = (layout.bbox.min, layout.bbox.max);
    let corners = [
        (min.x - m, min.y - m),
        (max.x + m, min.y - m),
        (max.x + m, max.y + m),
        (min.x - m, max.y + m),
    ];

    // Move (pen up) to corner 0.
    s.push_str(&format!(
        "{}D02*\n",
        coord_word_xy(corners[0].0, corners[0].1)
    ));
    // Draw (pen down) to corners 1, 2, 3.
    for &(x, y) in &corners[1..] {
        s.push_str(&format!("{}D01*\n", coord_word_xy(x, y)));
    }
    // Close the contour back to corner 0.
    s.push_str(&format!(
        "{}D01*\n",
        coord_word_xy(corners[0].0, corners[0].1)
    ));

    // --- Footer -------------------------------------------------------------
    s.push_str("M02*\n");

    s
}

/// Options controlling the KiCad `.kicad_pcb` board emission.
///
/// The skeleton emits the top-metal trace polygons as filled `gr_poly` on
/// [`copper_layer`](Self::copper_layer) plus a board outline on `Edge.Cuts`
/// expanded by [`outline_margin_mm`](Self::outline_margin_mm). The
/// [`generator`](Self::generator) string is written into the board's
/// `(generator …)` token.
///
/// `PartialEq` (not `Eq`) is derived because [`outline_margin_mm`] is an `f64`,
/// which is not `Eq` — this mirrors [`OutlineOptions`].
#[derive(Debug, Clone, PartialEq)]
pub struct KicadPcbOptions {
    /// Copper layer the trace polygons are placed on, written verbatim into
    /// each `(layer "…")` token. Defaults to `"F.Cu"` (KiCad front copper).
    pub copper_layer: String,
    /// Margin added on *each* side of the layout bounding box for the
    /// `Edge.Cuts` board outline, in millimetres. The outline rectangle spans
    /// `bbox` grown by this amount in `±x` and `±y`. Defaults to `1.0` mm.
    pub outline_margin_mm: f64,
    /// Value written into the board's `(generator "…")` token, identifying the
    /// tool that produced the file. Defaults to `"yee-export"`.
    pub generator: String,
}

impl Default for KicadPcbOptions {
    fn default() -> Self {
        Self {
            copper_layer: "F.Cu".into(),
            outline_margin_mm: 1.0,
            generator: "yee-export".into(),
        }
    }
}

/// Format a raw `(x_m, y_m)` coordinate pair (metres) as a KiCad `(xy X Y)`
/// S-expression element, converting metres → millimetres as a **float**.
///
/// KiCad `.kicad_pcb` coordinates are plain millimetre floating-point numbers,
/// e.g. `(xy 3.059000 1.234000)` — *not* the integer 4.6 fixed-point words used
/// by the Gerber emitters. This is why `xy_mm` is a separate helper from
/// [`mm_to_fixed46`] / [`coord_word_xy`]: the two output formats encode
/// coordinates differently and must not share a converter. Six decimal places
/// give nanometre resolution (`1e-6 mm`), matching the Gerber 4.6 fixed-point
/// precision while staying a human-readable float.
fn xy_mm(x_m: f64, y_m: f64) -> String {
    format!("(xy {:.6} {:.6})", x_m * 1.0e3, y_m * 1.0e3)
}

/// Emit a **KiCad 7 `.kicad_pcb`** S-expression board for a [`Layout`]: the
/// top-metal trace polygons as filled `gr_poly` on `opts.copper_layer`, plus
/// the board outline on `Edge.Cuts`.
///
/// Structure (see the module docs for the coordinate model):
///
/// - the `(kicad_pcb (version 20221018) (generator "<gen>") …)` header, a
///   `(general (thickness <h>))` block carrying the substrate height in mm, a
///   `(paper "A4")` token, a `(layers …)` table declaring `F.Cu` / `B.Cu` /
///   `Edge.Cuts`, and an empty `(setup)` block;
/// - per trace polygon (≥ 3 vertices): one filled
///   `(gr_poly (pts (xy …) …) (layer "<copper_layer>") (width 0) (fill solid))`
///   carrying the polygon's vertices in millimetres (see [`xy_mm`]). KiCad
///   closes a `gr_poly` implicitly, so the first vertex is **not** repeated;
/// - the board outline as a single
///   `(gr_poly (pts (xy …)(xy …)(xy …)(xy …)) (layer "Edge.Cuts") (width 0.1)
///   (fill none))` rectangle = `layout.bbox` expanded by
///   `opts.outline_margin_mm` on each side, in CCW order from the lower-left
///   corner (the same geometry as [`layout_to_gerber_outline`]);
/// - the closing `)` for the top-level `(kicad_pcb …)` form.
///
/// Coordinates are KiCad millimetre floats (metres × 1e3, six decimals) via
/// [`xy_mm`]. Polygons with fewer than three vertices are skipped (they have no
/// fillable area).
///
/// The returned `String` is a complete, standalone `.kicad_pcb` file body.
pub fn layout_to_kicad_pcb(layout: &Layout, opts: &KicadPcbOptions) -> String {
    let mut s = String::new();

    // --- Header -------------------------------------------------------------
    s.push_str("(kicad_pcb\n");
    s.push_str("  (version 20221018)\n");
    s.push_str(&format!("  (generator \"{}\")\n", opts.generator));
    // Board thickness from the substrate height, metres → mm.
    s.push_str(&format!(
        "  (general (thickness {:.6}))\n",
        layout.substrate.height_m * 1.0e3
    ));
    s.push_str("  (paper \"A4\")\n");
    // Layer table: front/back copper signal layers plus the board-profile
    // (Edge.Cuts) user layer.
    s.push_str("  (layers\n");
    s.push_str("    (0 \"F.Cu\" signal)\n");
    s.push_str("    (31 \"B.Cu\" signal)\n");
    s.push_str("    (44 \"Edge.Cuts\" user)\n");
    s.push_str("  )\n");
    s.push_str("  (setup)\n");

    // --- Copper polygons: one filled gr_poly per trace ----------------------
    for poly in &layout.traces {
        let verts = &poly.verts;
        if verts.len() < 3 {
            // Degenerate polygon: no area to fill, skip.
            continue;
        }
        s.push_str("  (gr_poly (pts");
        for v in verts {
            s.push(' ');
            s.push_str(&xy_mm(v.x, v.y));
        }
        s.push_str(&format!(
            ") (layer \"{}\") (width 0) (fill solid))\n",
            opts.copper_layer
        ));
    }

    // --- Board outline on Edge.Cuts -----------------------------------------
    // Rectangle corners = bbox expanded by `margin_m` on each side, CCW from
    // the lower-left corner. KiCad closes the gr_poly implicitly, so corner 0
    // is NOT repeated.
    let m = opts.outline_margin_mm * 1.0e-3;
    let (min, max) = (layout.bbox.min, layout.bbox.max);
    let corners = [
        (min.x - m, min.y - m),
        (max.x + m, min.y - m),
        (max.x + m, max.y + m),
        (min.x - m, max.y + m),
    ];
    s.push_str("  (gr_poly (pts");
    for &(x, y) in &corners {
        s.push(' ');
        s.push_str(&xy_mm(x, y));
    }
    s.push_str(") (layer \"Edge.Cuts\") (width 0.1) (fill none))\n");

    // --- Close the top-level form -------------------------------------------
    s.push_str(")\n");

    s
}
