//! # yee-layout
//!
//! Parametric planar-filter **geometry** for the Yee electromagnetic-simulation
//! studio (Filter Phase F1.0). Pure geometry: no EM, no meshing, no I/O beyond a
//! dependency-free SVG preview.
//!
//! This crate turns explicit physical dimensions (lengths, widths, gaps, tap
//! positions, substrate) into a [`Layout`] — a top-metal-on-substrate footprint
//! (`Vec<Polygon>` of traces plus `PortRef`s and a bounding box) that the later
//! F1.x sub-steps mesh on the FDTD back-end and export to manufacturing
//! formats. The dims→geometry direction only; the coupling-matrix→dimensions
//! mapping is the later F1.2 dimensional-synthesis step (ADR-0086).
//!
//! Two generators ship:
//!
//! - [`edge_coupled_bpf`] — N parallel half-wave coupled half-wavelength
//!   resonators plus end feed lines.
//! - [`hairpin_bpf`] — N U-folded resonators plus a tapped feed.
//!
//! Microstrip line sizing uses the Hammerstad-Jensen closed form
//! ([`microstrip_width`] / [`eps_eff`]). Edge-coupled-line electrical
//! parameters use the Kirschning-Jansen static even/odd model
//! ([`coupled_microstrip`] / [`coupling_coefficient`], in the [`coupled`]
//! module).
//!
//! ## References
//!
//! - Hammerstad & Jensen, "Accurate Models for Microstrip Computer-Aided
//!   Design," *IEEE MTT-S Int. Microwave Symp. Digest*, 1980 (width / ε_eff
//!   synthesis).
//! - Pozar, *Microwave Engineering* 4e, §3.8 (microstrip design equations).
//! - Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*,
//!   chs. 5 (edge-coupled) and 6 (hairpin).
//! - Kirschning & Jansen, "Accurate Wide-Range Design Equations for the
//!   Frequency-Dependent Characteristic of Parallel Coupled Microstrip Lines,"
//!   *IEEE Trans. MTT*, vol. 32, no. 1, pp. 83–90, 1984 (coupled even/odd
//!   model; see the [`coupled`] module).
//!
//! ## Example
//!
//! ```
//! use yee_layout::{microstrip_width, eps_eff};
//!
//! // FR-4 50 Ω microstrip on a 1.6 mm substrate: W ≈ 3.0 mm, ε_eff ≈ 3.3.
//! let h = 1.6e-3;
//! let w = microstrip_width(50.0, 4.4, h);
//! assert!((w - 3.0e-3).abs() / 3.0e-3 < 0.05);
//! let ee = eps_eff(w, h, 4.4);
//! assert!((ee - 3.3).abs() / 3.3 < 0.05);
//! ```

use serde::{Deserialize, Serialize};

pub mod coupled;
pub use coupled::{CoupledMicrostrip, coupled_microstrip, coupling_coefficient};

/// Free-space wave impedance `η₀ = 120π` Ω, used by the Hammerstad-Jensen
/// `B` term (`377π / (2·z0·√εr)` with `377 = 120π`).
const ETA0: f64 = 120.0 * std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Core geometry types
// ---------------------------------------------------------------------------

/// A dielectric substrate with a metal top layer and an implied ground plane.
///
/// All lengths are in metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Substrate {
    /// Relative permittivity `ε_r` of the dielectric.
    pub eps_r: f64,
    /// Substrate height `h` (metal-to-ground spacing), metres.
    pub height_m: f64,
    /// Dielectric loss tangent `tan δ` (dimensionless).
    pub loss_tangent: f64,
    /// Conductor (top-metal) thickness `t`, metres.
    pub metal_thickness_m: f64,
}

/// A point in the substrate top-plane, metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Point2 {
    /// `x` coordinate, metres.
    pub x: f64,
    /// `y` coordinate, metres.
    pub y: f64,
}

impl Point2 {
    /// Construct a point from `(x, y)` in metres.
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

/// A top-metal footprint: a closed polygon given by its vertices in order.
///
/// Generators in this crate emit counter-clockwise (CCW) windings, so
/// [`Polygon::signed_area`] is positive for every trace they produce.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Polygon {
    /// Vertices in order; the polygon is implicitly closed (last → first).
    pub verts: Vec<Point2>,
}

impl Polygon {
    /// Build an axis-aligned rectangle from its lower-left corner and size,
    /// wound counter-clockwise.
    ///
    /// `(x0, y0)` is the minimum corner; `w` and `h` are the extents along
    /// `x` and `y` respectively (metres).
    pub fn rect(x0: f64, y0: f64, w: f64, h: f64) -> Self {
        Self {
            verts: vec![
                Point2::new(x0, y0),
                Point2::new(x0 + w, y0),
                Point2::new(x0 + w, y0 + h),
                Point2::new(x0, y0 + h),
            ],
        }
    }

    /// The signed area of the polygon via the shoelace formula.
    ///
    /// Positive for a counter-clockwise winding, negative for clockwise. A
    /// magnitude near zero indicates a degenerate (collinear / zero-area)
    /// polygon.
    pub fn signed_area(&self) -> f64 {
        let n = self.verts.len();
        if n < 3 {
            return 0.0;
        }
        let mut acc = 0.0;
        for i in 0..n {
            let a = self.verts[i];
            let b = self.verts[(i + 1) % n];
            acc += a.x * b.y - b.x * a.y;
        }
        acc / 2.0
    }
}

/// A lumped port reference: a feed point, its trace width, and the reference
/// impedance the later EM solve de-embeds to.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PortRef {
    /// Port location (centre of the feed-line end), metres.
    pub at: Point2,
    /// Feed-line width at the port, metres.
    pub width_m: f64,
    /// Reference (system) impedance, ohms — typically `50`.
    pub ref_impedance_ohm: f64,
}

/// An axis-aligned bounding box over a layout's metal, metres.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BBox {
    /// Minimum (lower-left) corner.
    pub min: Point2,
    /// Maximum (upper-right) corner.
    pub max: Point2,
}

impl BBox {
    /// The box width (extent along `x`), metres.
    pub fn width(&self) -> f64 {
        self.max.x - self.min.x
    }

    /// The box height (extent along `y`), metres.
    pub fn height(&self) -> f64 {
        self.max.y - self.min.y
    }

    /// Compute the bounding box enclosing every vertex of every polygon.
    ///
    /// # Panics
    ///
    /// Panics if `polys` is empty or contains only empty polygons (no extent
    /// is defined).
    pub fn from_polygons(polys: &[Polygon]) -> Self {
        let mut min = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for poly in polys {
            for v in &poly.verts {
                min.x = min.x.min(v.x);
                min.y = min.y.min(v.y);
                max.x = max.x.max(v.x);
                max.y = max.y.max(v.y);
            }
        }
        assert!(
            min.x.is_finite() && max.x.is_finite(),
            "BBox::from_polygons needs at least one vertex"
        );
        Self { min, max }
    }
}

/// A complete planar-filter geometry: top metal on a substrate (ground plane
/// implied), with port references and a precomputed bounding box.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Layout {
    /// The dielectric substrate the metal sits on.
    pub substrate: Substrate,
    /// Top-metal trace footprints (resonators + feed lines).
    pub traces: Vec<Polygon>,
    /// Port references (feed terminations).
    pub ports: Vec<PortRef>,
    /// Axis-aligned bounding box over all traces.
    pub bbox: BBox,
}

impl Layout {
    /// Render a dependency-free top-view SVG of the layout.
    ///
    /// Metal traces are filled rectangles/polygons; ports are marked with a
    /// small circle. The `viewBox` is the layout's bounding box scaled to
    /// millimetres so the document is human-legible. The returned string is a
    /// complete, standalone `<svg>…</svg>` document.
    pub fn to_svg(&self) -> String {
        // Work in millimetres for a legible coordinate space, with a margin.
        const MM: f64 = 1.0e3;
        let margin_mm = 1.0;
        let min_x = self.bbox.min.x * MM - margin_mm;
        let min_y = self.bbox.min.y * MM - margin_mm;
        let w_mm = self.bbox.width() * MM + 2.0 * margin_mm;
        let h_mm = self.bbox.height() * MM + 2.0 * margin_mm;

        let mut s = String::new();
        s.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{min_x:.3} {min_y:.3} {w_mm:.3} {h_mm:.3}\" width=\"{w_mm:.3}mm\" height=\"{h_mm:.3}mm\">\n"
        ));
        // Background = substrate footprint.
        s.push_str(&format!(
            "  <rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" fill=\"#0a3d2e\" />\n",
            self.bbox.min.x * MM,
            self.bbox.min.y * MM,
            self.bbox.width() * MM,
            self.bbox.height() * MM
        ));
        // Traces = copper polygons.
        for poly in &self.traces {
            let pts: Vec<String> = poly
                .verts
                .iter()
                .map(|p| format!("{:.4},{:.4}", p.x * MM, p.y * MM))
                .collect();
            s.push_str(&format!(
                "  <polygon points=\"{}\" fill=\"#d4942a\" stroke=\"#b87814\" stroke-width=\"0.02\" />\n",
                pts.join(" ")
            ));
        }
        // Ports = markers.
        for port in &self.ports {
            s.push_str(&format!(
                "  <circle cx=\"{:.4}\" cy=\"{:.4}\" r=\"{:.4}\" fill=\"#e23b3b\" />\n",
                port.at.x * MM,
                port.at.y * MM,
                (port.width_m * MM * 0.25).max(0.05)
            ));
        }
        s.push_str("</svg>\n");
        s
    }
}

// ---------------------------------------------------------------------------
// Hammerstad-Jensen microstrip synthesis (Pozar §3.8)
// ---------------------------------------------------------------------------

/// Synthesize the microstrip line width `W` (metres) for a target characteristic
/// impedance.
///
/// Hammerstad-Jensen closed form (spec §"Hammerstad-Jensen synthesis"):
///
/// ```text
/// A = z0/60·√((εr+1)/2) + (εr−1)/(εr+1)·(0.23 + 0.11/εr)
/// B = 377π/(2·z0·√εr)
/// W/h = 8·e^A/(e^{2A}−2)                                          if that ratio < 2
/// W/h = (2/π)[B−1−ln(2B−1) + (εr−1)/(2εr)·(ln(B−1)+0.39−0.61/εr)] otherwise
/// ```
///
/// The wide-line (`>2`) branch is selected only when the thin-line (`<2`)
/// branch yields `W/h ≥ 2`, matching the regime each formula is valid in. For
/// FR-4 50 Ω (`εr = 4.4`, `h = 1.6 mm`) the ratio is ≈ 1.9, so the thin-line
/// branch applies and `W ≈ 3.0 mm`.
///
/// # Arguments
///
/// - `z0_ohm` — target characteristic impedance, ohms.
/// - `eps_r` — substrate relative permittivity.
/// - `h_m` — substrate height, metres.
pub fn microstrip_width(z0_ohm: f64, eps_r: f64, h_m: f64) -> f64 {
    let a = z0_ohm / 60.0 * ((eps_r + 1.0) / 2.0).sqrt()
        + (eps_r - 1.0) / (eps_r + 1.0) * (0.23 + 0.11 / eps_r);
    // ETA0 = 120π = 377; `377π / (2·z0·√εr)`.
    let b = ETA0 * std::f64::consts::PI / (2.0 * z0_ohm * eps_r.sqrt());

    // Thin-line branch first.
    let w_over_h_thin = 8.0 * a.exp() / ((2.0 * a).exp() - 2.0);
    let w_over_h = if w_over_h_thin < 2.0 {
        w_over_h_thin
    } else {
        // Wide-line branch.
        (2.0 / std::f64::consts::PI)
            * (b - 1.0 - (2.0 * b - 1.0).ln()
                + (eps_r - 1.0) / (2.0 * eps_r) * ((b - 1.0).ln() + 0.39 - 0.61 / eps_r))
    };
    w_over_h * h_m
}

/// Effective permittivity `ε_eff` of a microstrip line of width `W` on a
/// substrate of height `h` and permittivity `εr`.
///
/// ```text
/// ε_eff = (εr+1)/2 + (εr−1)/2 · (1 + 12·h/W)^(−1/2)
/// ```
///
/// No conductor-thickness correction is applied at this fidelity (spec
/// §"Math notes").
///
/// # Arguments
///
/// - `w_m` — line width, metres.
/// - `h_m` — substrate height, metres.
/// - `eps_r` — substrate relative permittivity.
pub fn eps_eff(w_m: f64, h_m: f64, eps_r: f64) -> f64 {
    (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * (1.0 + 12.0 * h_m / w_m).powf(-0.5)
}

// ---------------------------------------------------------------------------
// Generators
// ---------------------------------------------------------------------------

/// One section of an edge-coupled band-pass filter: a coupled half-wave strip.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EdgeCoupledSection {
    /// Resonator (coupled-strip) length along `x`, metres.
    pub length_m: f64,
    /// Resonator width along `y`, metres.
    pub width_m: f64,
    /// Edge-coupling gap to the next strip, metres.
    pub gap_m: f64,
}

/// Parameters for [`edge_coupled_bpf`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EdgeCoupledParams {
    /// The dielectric substrate.
    pub substrate: Substrate,
    /// The coupled-strip sections, in order.
    pub sections: Vec<EdgeCoupledSection>,
    /// Feed-line width, metres.
    pub feed_width_m: f64,
    /// Feed-line length, metres.
    pub feed_length_m: f64,
}

/// Generate an **edge-coupled** band-pass filter [`Layout`].
///
/// Lays the `N` coupled half-wave strips along `x`; adjacent strips are offset
/// in `y` by `width + gap` and alternately staggered by a half length to form
/// the edge-coupled overlap of each section (Hong & Lancaster ch. 5). Feed
/// lines of `feed_width_m × feed_length_m` attach at the two outer ends, and a
/// `PortRef` (`50` Ω) sits at each outer feed-line end. The bounding box is
/// computed from all polygons.
///
/// # Panics
///
/// Panics if `p.sections` is empty.
pub fn edge_coupled_bpf(p: &EdgeCoupledParams) -> Layout {
    assert!(
        !p.sections.is_empty(),
        "edge_coupled_bpf needs at least one section"
    );

    let mut traces: Vec<Polygon> = Vec::with_capacity(p.sections.len() + 2);

    // Place the coupled strips. Each strip starts at the running `y`; adjacent
    // strips alternate a half-length stagger in `x` so consecutive sections
    // overlap over half their length (the edge-coupled coupling region).
    let mut y = 0.0_f64;
    let mut strip_records: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(p.sections.len());
    for (i, sec) in p.sections.iter().enumerate() {
        let x0 = if i % 2 == 0 { 0.0 } else { sec.length_m / 2.0 };
        traces.push(Polygon::rect(x0, y, sec.length_m, sec.width_m));
        strip_records.push((x0, y, sec.length_m, sec.width_m));
        // Advance to the next strip's lower edge: this strip's width + gap.
        y += sec.width_m + sec.gap_m;
    }

    // Feed line at the input: attaches to the left end of the first strip,
    // extending leftward (−x).
    let (fx0, fy0, _flen, fwid) = strip_records[0];
    let in_feed_x = fx0 - p.feed_length_m;
    let in_feed_y = fy0 + fwid / 2.0 - p.feed_width_m / 2.0;
    traces.push(Polygon::rect(
        in_feed_x,
        in_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(in_feed_x, fy0 + fwid / 2.0),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };

    // Feed line at the output: attaches to the right end of the last strip,
    // extending rightward (+x).
    let (lx0, ly0, llen, lwid) = *strip_records.last().unwrap();
    let out_feed_x = lx0 + llen;
    let out_feed_y = ly0 + lwid / 2.0 - p.feed_width_m / 2.0;
    traces.push(Polygon::rect(
        out_feed_x,
        out_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(out_feed_x + p.feed_length_m, ly0 + lwid / 2.0),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };

    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: p.substrate,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    }
}

/// Parameters for [`hairpin_bpf`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HairpinParams {
    /// The dielectric substrate.
    pub substrate: Substrate,
    /// Number of hairpin (U-folded) resonators.
    pub n: usize,
    /// Length of each resonator arm along `y`, metres.
    pub arm_length_m: f64,
    /// Line width of each arm / bend, metres.
    pub line_width_m: f64,
    /// Centre-to-centre spacing of the two arms of one hairpin, metres.
    pub fold_spacing_m: f64,
    /// Edge-coupling gap between adjacent resonators, metres.
    pub coupling_gap_m: f64,
    /// Tap offset of the feed up each end-resonator arm, metres.
    pub tap_offset_m: f64,
    /// Feed-line width, metres.
    pub feed_width_m: f64,
    /// Feed-line length, metres.
    pub feed_length_m: f64,
}

/// Generate a **hairpin** band-pass filter [`Layout`].
///
/// Places `n` U-folded resonators (two arms of `arm_length_m` joined by a
/// bend) side by side, spaced by `coupling_gap_m`, and taps a feed line onto
/// the outer arm of each end resonator at height `tap_offset_m` (Hong &
/// Lancaster ch. 6). Each hairpin is built from three rectangles (two arms +
/// the connecting bend) so the metal is a single connected footprint. Ports
/// (`50` Ω) sit at the outer end of each tap feed. The bounding box is computed
/// from all polygons.
///
/// # Panics
///
/// Panics if `p.n == 0`.
pub fn hairpin_bpf(p: &HairpinParams) -> Layout {
    assert!(p.n >= 1, "hairpin_bpf needs at least one resonator");

    let mut traces: Vec<Polygon> = Vec::with_capacity(p.n * 3 + 2);

    // One hairpin occupies, in x, the span of its two arms plus the fold:
    // left arm at x_base, right arm at x_base + fold_spacing. Each resonator's
    // x-pitch is fold_spacing + line_width + coupling_gap.
    let pitch = p.fold_spacing_m + p.line_width_m + p.coupling_gap_m;
    let arm_h = p.arm_length_m;
    let lw = p.line_width_m;

    for i in 0..p.n {
        let x_base = i as f64 * pitch;
        // Left arm.
        traces.push(Polygon::rect(x_base, 0.0, lw, arm_h));
        // Right arm.
        let right_x = x_base + p.fold_spacing_m;
        traces.push(Polygon::rect(right_x, 0.0, lw, arm_h));
        // Connecting bend across the top, joining the two arms.
        let bend_w = (right_x + lw) - x_base;
        traces.push(Polygon::rect(x_base, arm_h - lw, bend_w, lw));
    }

    // Input feed: taps the left arm of the first resonator at tap_offset_m,
    // extending leftward (−x).
    let in_feed_x = -p.feed_length_m;
    let in_feed_y = p.tap_offset_m - p.feed_width_m / 2.0;
    traces.push(Polygon::rect(
        in_feed_x,
        in_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(in_feed_x, p.tap_offset_m),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };

    // Output feed: taps the right arm of the last resonator at tap_offset_m,
    // extending rightward (+x).
    let last_x_base = (p.n - 1) as f64 * pitch;
    let last_right_arm_x = last_x_base + p.fold_spacing_m + lw;
    let out_feed_y = p.tap_offset_m - p.feed_width_m / 2.0;
    traces.push(Polygon::rect(
        last_right_arm_x,
        out_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(last_right_arm_x + p.feed_length_m, p.tap_offset_m),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };

    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: p.substrate,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    }
}

/// Parameters for [`hairpin_bpf_sections`] — the **per-section-gap** hairpin
/// generator (R.4/F1.2.1). Unlike [`HairpinParams`], which bakes one
/// `coupling_gap_m` into a uniform pitch, this carries the N − 1 *distinct*
/// inter-resonator gaps synthesis produces (one per coupling `k_{i,i+1}`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HairpinSectionParams {
    /// The dielectric substrate.
    pub substrate: Substrate,
    /// Length of each resonator arm along `y`, metres.
    pub arm_length_m: f64,
    /// Line width of each arm / bend, metres.
    pub line_width_m: f64,
    /// Centre-to-centre spacing of the two arms of one hairpin, metres.
    pub fold_spacing_m: f64,
    /// Per-section edge-coupling gaps, metres — `gaps_m[i]` separates
    /// resonators `i` and `i + 1`, so the resonator count is
    /// `gaps_m.len() + 1`.
    pub gaps_m: Vec<f64>,
    /// Tap offset of the feed up each end-resonator arm, metres, measured
    /// from the open (bottom) end.
    pub tap_offset_m: f64,
    /// Feed-line width, metres.
    pub feed_width_m: f64,
    /// Feed-line length, metres.
    pub feed_length_m: f64,
}

/// Generate a **per-section-gap hairpin** band-pass filter [`Layout`]
/// (R.4/F1.2.1 — "gap option (a)" from the F1.2.2 dimensioning notes).
///
/// Identical resonator geometry to [`hairpin_bpf`] — `n` U-folded resonators
/// of three rectangles each, tapped feeds on the outer arms of the end
/// resonators — except each adjacent pair sits at its **own** solved gap:
/// resonator `i + 1`'s x-base is resonator `i`'s base advanced by
/// `fold_spacing + line_width + gaps_m[i]`. [`hairpin_bpf`] (and its
/// committed `geo-003` gate) is untouched; a uniform `gaps_m` reproduces its
/// resonator placement exactly.
///
/// # Panics
///
/// Panics if `gaps_m` is empty (a one-resonator "filter" has no coupling
/// section; use [`hairpin_bpf`] for degenerate single-resonator layouts).
pub fn hairpin_bpf_sections(p: &HairpinSectionParams) -> Layout {
    assert!(
        !p.gaps_m.is_empty(),
        "hairpin_bpf_sections needs at least one coupling gap (two resonators)"
    );
    assert!(
        p.fold_spacing_m > p.line_width_m,
        "fold_spacing_m is centre-to-centre: <= line_width_m merges the two arms of the U"
    );
    let n = p.gaps_m.len() + 1;
    let mut traces: Vec<Polygon> = Vec::with_capacity(n * 3 + 2);
    let arm_h = p.arm_length_m;
    let lw = p.line_width_m;

    let mut x_base = 0.0;
    let mut last_x_base = 0.0;
    for i in 0..n {
        traces.push(Polygon::rect(x_base, 0.0, lw, arm_h));
        let right_x = x_base + p.fold_spacing_m;
        traces.push(Polygon::rect(right_x, 0.0, lw, arm_h));
        let bend_w = (right_x + lw) - x_base;
        traces.push(Polygon::rect(x_base, arm_h - lw, bend_w, lw));
        last_x_base = x_base;
        if i < n - 1 {
            x_base += p.fold_spacing_m + lw + p.gaps_m[i];
        }
    }

    // Tapped feeds on the outer arms of the end resonators, exactly as
    // hairpin_bpf places them.
    let in_feed_x = -p.feed_length_m;
    let in_feed_y = p.tap_offset_m - p.feed_width_m / 2.0;
    traces.push(Polygon::rect(
        in_feed_x,
        in_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let in_port = PortRef {
        at: Point2::new(in_feed_x, p.tap_offset_m),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };
    let last_right_arm_x = last_x_base + p.fold_spacing_m + lw;
    traces.push(Polygon::rect(
        last_right_arm_x,
        in_feed_y,
        p.feed_length_m,
        p.feed_width_m,
    ));
    let out_port = PortRef {
        at: Point2::new(last_right_arm_x + p.feed_length_m, p.tap_offset_m),
        width_m: p.feed_width_m,
        ref_impedance_ohm: 50.0,
    };

    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: p.substrate,
        traces,
        ports: vec![in_port, out_port],
        bbox,
    }
}

// ---------------------------------------------------------------------------
// Rectangular patch antenna (A.0, ADR-0190) — Balanis closed forms.
// ---------------------------------------------------------------------------

/// Closed-form rectangular-patch dimensions (Balanis, *Antenna Theory*
/// §14.2, transmission-line model), returned by [`patch_antenna_dims`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PatchDims {
    /// Radiating-edge width `W = c/(2f0) * sqrt(2/(eps_r+1))` (metres).
    pub width_m: f64,
    /// Resonant length `L = c/(2 f0 sqrt(eps_eff)) - 2*dL` (metres).
    pub length_m: f64,
    /// Effective permittivity of the patch section at width `W`.
    pub eps_eff: f64,
    /// Hammerstad open-end length extension `dL` (metres).
    pub delta_l_m: f64,
}

/// Hammerstad microstrip open-end length correction (also used by the
/// engine-sparams stub gate): the fringing field makes an open microstrip
/// edge electrically longer than its metal:
/// `dL = 0.412*h*(eps_eff+0.3)(W/h+0.264) / ((eps_eff-0.258)(W/h+0.8))`.
pub fn open_end_delta_l(w_m: f64, h_m: f64, e_eff: f64) -> f64 {
    let u = w_m / h_m;
    0.412 * h_m * ((e_eff + 0.3) * (u + 0.264)) / ((e_eff - 0.258) * (u + 0.8))
}

/// Design a rectangular patch for resonance at `f0_hz` on a substrate of
/// `eps_r`, height `h_m` (Balanis §14.2):
/// `W = c/(2 f0) * sqrt(2/(eps_r+1))`, `eps_eff` from the Hammerstad
/// closed form ([`eps_eff`]) at that width, and
/// `L = c/(2 f0 sqrt(eps_eff)) - 2*dL` with the open-end extension `dL`.
///
/// # Panics
///
/// Panics if `f0_hz`, `eps_r`, or `h_m` is non-positive.
pub fn patch_antenna_dims(f0_hz: f64, eps_r: f64, h_m: f64) -> PatchDims {
    assert!(
        f0_hz > 0.0 && eps_r > 0.0 && h_m > 0.0,
        "non-physical patch inputs"
    );
    const C: f64 = 299_792_458.0;
    let width_m = C / (2.0 * f0_hz) * (2.0 / (eps_r + 1.0)).sqrt();
    let e_eff = eps_eff(width_m, h_m, eps_r);
    let delta_l_m = open_end_delta_l(width_m, h_m, e_eff);
    let length_m = C / (2.0 * f0_hz * e_eff.sqrt()) - 2.0 * delta_l_m;
    PatchDims {
        width_m,
        length_m,
        eps_eff: e_eff,
        delta_l_m,
    }
}

/// Assemble an **edge-fed** rectangular-patch [`Layout`]: a `feed_z0_ohm`
/// microstrip feed line joining the centre of the patch's radiating edge,
/// one `PortRef` at the outer feed end. The resonant length `L` runs
/// along `x` (the feed direction); the radiating edges are the two
/// `y`-parallel ends. Edge feeding is deliberately unmatched (the edge
/// resistance of a patch is hundreds of ohms) — resonance shows as a
/// localized |S11| dip, and matching (inset feed) is the A.1 follow-on.
///
/// The feed is half a guided wavelength long so measurement probes fit
/// on it upstream of the patch.
pub fn edge_fed_patch(f0_hz: f64, substrate: &Substrate, feed_z0_ohm: f64) -> Layout {
    const C: f64 = 299_792_458.0;
    let dims = patch_antenna_dims(f0_hz, substrate.eps_r, substrate.height_m);
    let feed_w = microstrip_width(feed_z0_ohm, substrate.eps_r, substrate.height_m);
    let feed_e_eff = eps_eff(feed_w, substrate.height_m, substrate.eps_r);
    let feed_len = C / (2.0 * f0_hz * feed_e_eff.sqrt());

    let traces = vec![
        // Feed line along x, centred on y = 0, ending at the patch edge.
        Polygon::rect(-feed_len, -feed_w / 2.0, feed_len, feed_w),
        // The patch: resonant length L along x, width W along y.
        Polygon::rect(0.0, -dims.width_m / 2.0, dims.length_m, dims.width_m),
    ];
    let port = PortRef {
        at: Point2::new(-feed_len, 0.0),
        width_m: feed_w,
        ref_impedance_ohm: feed_z0_ohm,
    };
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: *substrate,
        traces,
        ports: vec![port],
        bbox,
    }
}

/// Quasi-Yagi seed dimensions (FS.1a.1, ADR-0205): the Kaneda/Deal/Qian/
/// Itoh scaling rules evaluated for a given band and stack. All lengths in
/// metres; the layout frame puts the feed start at `x = 0, y = 0` with the
/// beam firing toward `+x`.
#[derive(Debug, Clone, Copy)]
pub struct QuasiYagiDims {
    /// 50 Ω feed width.
    pub feed_width_m: f64,
    /// Balun-branch / element strip width (~70.7 Ω branches).
    pub branch_width_m: f64,
    /// CPS centre-to-centre separation of the two branches.
    pub cps_sep_m: f64,
    /// Driven-dipole tip-to-tip length (≈ 0.46 λ₀/√((ε_r+1)/2)).
    pub dipole_len_m: f64,
    /// Director length (≈ 0.65 × dipole).
    pub director_len_m: f64,
    /// Ground-edge (reflector) → driven-element spacing (≈ 0.25 λ_diel).
    pub reflector_gap_m: f64,
    /// Driven-element → director spacing (≈ 0.20 λ_diel).
    pub director_gap_m: f64,
    /// x of the driven-dipole centreline.
    pub x_dipole_m: f64,
    /// x of the ground-plane truncation edge — pass through
    /// `yee_voxel::truncate_ground_at_cell` (the ground *is* the
    /// reflector; everything at `x < x_gnd_m` is microstrip over ground,
    /// everything beyond is the CPS/antenna region).
    pub x_gnd_m: f64,
}

/// A generated quasi-Yagi: the trace [`Layout`] plus the ground-truncation
/// plane the voxelizer must apply (the layout alone cannot express a
/// partial ground — see `yee_voxel::truncate_ground_at_cell`).
#[derive(Debug, Clone)]
pub struct QuasiYagi {
    /// Traces + the single feed port.
    pub layout: Layout,
    /// The seed dimensions used.
    pub dims: QuasiYagiDims,
}

/// Generate a microstrip-fed **quasi-Yagi** antenna [`Layout`] (FS.1a.1):
/// 50 Ω feed → T-junction balun whose lower branch takes a U-detour
/// adding half a branch guided wavelength (the 180° split) → the two
/// branches continue as a CPS pair past the ground truncation → driven
/// dipole arms (±y) + one director, the truncated ground edge acting as
/// the reflector. All axis-aligned rectangles with overlapped joints
/// (the voxelizer point-samples cell centres; shared-edge-only joints
/// rasterize inconsistently).
///
/// Walking-skeleton fidelity notes (spec 2026-07-08): plain rectangular
/// corners (no miters — the R.6 corner-correction lesson applies if
/// measurement shows balun detune), and the published X-band reference
/// stack (ε_r 10.2, 0.09 mm gaps) is deliberately NOT the target — the
/// scaling rules are evaluated on the caller's stack and the gate holds
/// the measured S11 dip against the design `f0`.
pub fn quasi_yagi(f0_hz: f64, substrate: &Substrate, feed_z0_ohm: f64) -> QuasiYagi {
    const C: f64 = 299_792_458.0;
    assert!(
        f0_hz > 0.0 && substrate.eps_r > 1.0 && substrate.height_m > 0.0,
        "non-physical quasi-Yagi inputs"
    );
    let h = substrate.height_m;
    let eps_r = substrate.eps_r;
    let lambda0 = C / f0_hz;
    // Effective permittivity of the dipole/director elements. The
    // half-space quasi-static value ε_avg = (ε_r+1)/2 assumes the
    // substrate fills a half-space; a resonant dipole on a THIN substrate
    // (h ≪ λ) is loaded far less. FDTD-calibrated (ADR-0205, the R.6
    // single-point-calibration pattern): the (ε_r+1)/2 seed on FR-4
    // 1.6 mm measured its λ/2 resonance 29 % HIGH (7.5 vs 5.8 GHz),
    // giving ε_eff = 2.7·(5.8/7.5)² = 1.61; the linear-in-(ε_r−1) form
    // below reproduces that measurement and degrades gracefully toward
    // ε → 1 for air. Re-verified blind after calibration (the gate).
    let eps_dipole = 1.0 + 0.18 * (eps_r - 1.0);
    let lambda_diel = lambda0 / eps_dipole.sqrt();

    let feed_w = microstrip_width(feed_z0_ohm, eps_r, h);
    let feed_e_eff = eps_eff(feed_w, h, eps_r);
    let lambda_g_feed = C / (f0_hz * feed_e_eff.sqrt());
    // ~70.7 Ω branches: reasonable T-split compromise on mm-scale stacks.
    let branch_w = microstrip_width(feed_z0_ohm * std::f64::consts::SQRT_2, eps_r, h);
    let branch_e_eff = eps_eff(branch_w, h, eps_r);
    let lambda_g_branch = C / (f0_hz * branch_e_eff.sqrt());

    let feed_len = 0.75 * lambda_g_feed;
    let cps_gap = 1.0e-3;
    let cps_sep = cps_gap + branch_w; // centre-to-centre
    let dipole_len = 0.46 * lambda_diel;
    let director_len = 0.65 * dipole_len;
    let reflector_gap = 0.25 * lambda_diel;
    let director_gap = 0.20 * lambda_diel;
    // The balun region (T-bar to ground edge) hosts the U-detour; the
    // detour adds 2·dy = λg_branch/2 of centreline path.
    let balun_len = 0.5 * lambda_g_branch + 3.0 * branch_w;
    let dy = 0.25 * lambda_g_branch;

    let x_t = feed_len;
    let x_gnd = x_t + branch_w + balun_len;
    let x_dip = x_gnd + reflector_gap;
    let x_dir = x_dip + director_gap;
    let w_arm = branch_w;

    let y_top = cps_sep / 2.0; // top branch centreline
    let y_bot = -cps_sep / 2.0;

    // Branch rows and the T-bar.
    let mut traces = vec![
        // Feed along x, centred on y = 0.
        Polygon::rect(0.0, -feed_w / 2.0, feed_len + branch_w, feed_w),
        // T-bar (vertical) spanning both branch centrelines.
        Polygon::rect(x_t, y_bot - branch_w / 2.0, branch_w, cps_sep + branch_w),
        // Top branch (A): straight T-bar → dipole plane.
        Polygon::rect(
            x_t,
            y_top - branch_w / 2.0,
            x_dip + w_arm / 2.0 - x_t,
            branch_w,
        ),
    ];
    // Bottom branch (B): seg1 → down → bottom run → up → seg3, the
    // centreline detour adding exactly 2·dy of path.
    let l1 = 0.25 * balun_len;
    let l2 = 0.35 * balun_len;
    let x_a = x_t + branch_w + l1; // seg1 outer end
    let x_b = x_a + l2; // bottom-run outer end
    let y_low = y_bot - dy;
    traces.extend([
        // seg1 on the branch row.
        Polygon::rect(x_t, y_bot - branch_w / 2.0, branch_w + l1, branch_w),
        // down bar.
        Polygon::rect(
            x_a - branch_w,
            y_low - branch_w / 2.0,
            branch_w,
            dy + branch_w,
        ),
        // bottom run.
        Polygon::rect(
            x_a - branch_w,
            y_low - branch_w / 2.0,
            l2 + branch_w,
            branch_w,
        ),
        // up bar.
        Polygon::rect(
            x_b - branch_w,
            y_low - branch_w / 2.0,
            branch_w,
            dy + branch_w,
        ),
        // seg3: back on the branch row, to the dipole plane.
        Polygon::rect(
            x_b - branch_w,
            y_bot - branch_w / 2.0,
            x_dip + w_arm / 2.0 - (x_b - branch_w),
            branch_w,
        ),
        // Driven arms: from each branch row (overlapping it) out to the
        // dipole tips at ±dipole_len/2. Both heights equal by symmetry:
        // dipole_len/2 − cps_sep/2 + branch_w/2.
        Polygon::rect(
            x_dip - w_arm / 2.0,
            y_top - branch_w / 2.0,
            w_arm,
            dipole_len / 2.0 - (y_top - branch_w / 2.0),
        ),
        Polygon::rect(
            x_dip - w_arm / 2.0,
            -dipole_len / 2.0,
            w_arm,
            (y_bot + branch_w / 2.0) - (-dipole_len / 2.0),
        ),
        // Director.
        Polygon::rect(
            x_dir - w_arm / 2.0,
            -director_len / 2.0,
            w_arm,
            director_len,
        ),
    ]);

    let port = PortRef {
        at: Point2::new(0.0, 0.0),
        width_m: feed_w,
        ref_impedance_ohm: feed_z0_ohm,
    };
    let bbox = BBox::from_polygons(&traces);
    let layout = Layout {
        substrate: *substrate,
        traces,
        ports: vec![port],
        bbox,
    };
    QuasiYagi {
        layout,
        dims: QuasiYagiDims {
            feed_width_m: feed_w,
            branch_width_m: branch_w,
            cps_sep_m: cps_sep,
            dipole_len_m: dipole_len,
            director_len_m: director_len,
            reflector_gap_m: reflector_gap,
            director_gap_m: director_gap,
            x_dipole_m: x_dip,
            x_gnd_m: x_gnd,
        },
    }
}

/// 2×1 patch-array seed dimensions (FS.1b, ADR-0206). Lengths in metres;
/// layout frame: patch left edges at `x = 0`, array symmetric about
/// `y = 0`, port at the −x end of the feed spine.
#[derive(Debug, Clone, Copy)]
pub struct PatchArrayDims {
    /// Single-element dims (Balanis, the A.0-certified closed form).
    pub patch: PatchDims,
    /// Centre-to-centre element spacing (0.5 λ₀).
    pub spacing_m: f64,
    /// 50 Ω line width.
    pub feed_width_m: f64,
    /// λg/4 impedance-transformer width (≈ 70.7 Ω).
    pub xfmr_width_m: f64,
    /// Transformer length (λg at 70.7 Ω / 4).
    pub xfmr_len_m: f64,
    /// Inset depth used on each element (the A.3-measured 0.25·L).
    pub inset_m: f64,
    /// x of the corporate junction / vertical tree centreline.
    pub x_junction_m: f64,
}

/// A generated 2×1 array: the layout plus its seed dims.
#[derive(Debug, Clone)]
pub struct PatchArray {
    /// Traces + the single feed port.
    pub layout: Layout,
    /// The seed dimensions used.
    pub dims: PatchArrayDims,
}

/// Corner style for [`double_jog`] (FS.3.2a, ADR-0217).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MiterStyle {
    /// Plain square corners (the full w×w corner metal).
    Square,
    /// The outer corner of each bend chopped by a 45° edge (Douville &
    /// James): the cut legs run `f·w` along each outer edge. `f` must be
    /// in `(0, 1)`; ~0.7 is near the published optimum for w/h ≈ 1.9.
    Mitered {
        /// Cut-leg fraction of the trace width.
        f: f64,
    },
}

/// Generate a **double-jog through line** (FS.3.2a): port 1 → x-feed →
/// 90° up-bend → y-riser (`jog_dy_m`) → x-run (`gap_x_m` between the
/// risers) → 90° down-bend → x-feed → port 2. Both ports face ±x at the
/// same y, so the standard two-port fixtures (uniform and graded) apply
/// unchanged, and [`MiterStyle::Mitered`] corners carry the repo's first
/// **non-axis-aligned polygon edges** into a full-wave measurement.
///
/// Construction notes:
/// - Straight segments are rects that **overlap** each corner region by
///   `0.2·w` (the voxelizer point-samples cell centres; shared-edge-only
///   joints rasterize inconsistently — the quasi-Yagi lesson). The
///   overlap is capped well inside the un-cut corner metal: the miter
///   cut spans `f·w` from the outer edges, and the overlaps stay `≥
///   (1−f−0.2)·w` clear of it for the default f = 0.7.
/// - Each corner is its **own polygon** (square rect or 5-vertex mitered
///   pentagon), so the automesh rulebook's per-polygon AABBs put fine
///   bands at every bend (a single-outline polygon would present one
///   blob AABB and leave the bends unrefined).
pub fn double_jog(
    substrate: &Substrate,
    w_m: f64,
    run_x_m: f64,
    gap_x_m: f64,
    jog_dy_m: f64,
    style: MiterStyle,
) -> Layout {
    assert!(
        w_m > 0.0 && run_x_m > 0.0 && gap_x_m > 0.0 && jog_dy_m > w_m,
        "non-physical double-jog inputs (need jog_dy > w)"
    );
    if let MiterStyle::Mitered { f } = style {
        assert!(f > 0.0 && f < 1.0, "miter fraction {f} outside (0, 1)");
        assert!(
            f <= 0.8,
            "miter fraction {f} leaves < 0.2·w of corner metal — the segment \
             overlaps would intrude into the cut"
        );
    }
    let w = w_m;
    let ov = 0.2 * w;
    let xa = run_x_m; // riser-A left edge
    let xb = xa + w + gap_x_m; // riser-B left edge
    let dy = jog_dy_m; // mid-run bottom edge
    let end = xb + w + run_x_m;

    // The four corner regions, outer-corner cut per style. Bends: A turns
    // +x→+y (outer corner bottom-right), B +y→+x (top-left), C +x→−y
    // (top-right), D −y→+x (bottom-left).
    let corner = |x0: f64, y0: f64, outer: (f64, f64)| -> Polygon {
        match style {
            MiterStyle::Square => Polygon::rect(x0, y0, w, w),
            MiterStyle::Mitered { f } => {
                let c = f * w;
                let (cx, cy) = outer;
                // Start from the rect's verts and replace the outer
                // corner with the two cut points, each `c` along an
                // adjacent edge toward the corner's rect interior.
                let mut verts = Vec::with_capacity(5);
                for v in Polygon::rect(x0, y0, w, w).verts {
                    if (v.x - cx).abs() < 1e-15 && (v.y - cy).abs() < 1e-15 {
                        // The two edge directions from the outer corner
                        // back into the rect: ±x and ±y toward centre.
                        let sx = if cx > x0 + w / 2.0 { -1.0 } else { 1.0 };
                        let sy = if cy > y0 + w / 2.0 { -1.0 } else { 1.0 };
                        // Order the pair to preserve the rect winding:
                        // the vertex reached first along the incoming
                        // edge keeps that edge's axis.
                        if sy > 0.0 && sx < 0.0 || sy < 0.0 && sx > 0.0 {
                            verts.push(Point2::new(cx + sx * c, cy));
                            verts.push(Point2::new(cx, cy + sy * c));
                        } else {
                            verts.push(Point2::new(cx, cy + sy * c));
                            verts.push(Point2::new(cx + sx * c, cy));
                        }
                    } else {
                        verts.push(v);
                    }
                }
                Polygon { verts }
            }
        }
    };

    let traces = vec![
        // Straights, each overlapping its corner(s) by `ov`.
        Polygon::rect(0.0, 0.0, xa + ov, w),
        Polygon::rect(xa, w - ov, w, dy - w + 2.0 * ov),
        Polygon::rect(xa + w - ov, dy, gap_x_m + 2.0 * ov, w),
        Polygon::rect(xb, w - ov, w, dy - w + 2.0 * ov),
        Polygon::rect(xb + w - ov, 0.0, run_x_m + ov, w),
        // Corners A, B, C, D.
        corner(xa, 0.0, (xa + w, 0.0)),
        corner(xa, dy, (xa, dy + w)),
        corner(xb, dy, (xb + w, dy + w)),
        corner(xb, 0.0, (xb, 0.0)),
    ];
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: *substrate,
        traces,
        ports: vec![
            PortRef {
                at: Point2::new(0.5e-3, w / 2.0),
                width_m: w,
                ref_impedance_ohm: 50.0,
            },
            PortRef {
                at: Point2::new(end - 0.5e-3, w / 2.0),
                width_m: w,
                ref_impedance_ohm: 50.0,
            },
        ],
        bbox,
    }
}

/// Generate a **2×1 corporate-fed patch array** (FS.1b): two inset-fed
/// patches side by side along y (H-plane pair) at 0.5 λ₀ spacing, fed in
/// phase through a symmetric corporate tree — 50 Ω spine → junction → two
/// λg/4 **70.7 Ω transformers** along ±y (each transforms the branch's
/// 50 Ω to 100 Ω; the two in parallel present 50 Ω at the junction) →
/// 50 Ω branches → each element's inset (the A.3-measured 0.25·L depth).
/// Phase balance is exact by mirror symmetry, so the tree's corner path
/// errors (the R.6 κ effect) can only detune the common match — which the
/// S11 gate measures.
pub fn patch_array_2x1(f0_hz: f64, substrate: &Substrate, feed_z0_ohm: f64) -> PatchArray {
    const C: f64 = 299_792_458.0;
    let patch = patch_antenna_dims(f0_hz, substrate.eps_r, substrate.height_m);
    let (w, l) = (patch.width_m, patch.length_m);
    let lambda0 = C / f0_hz;
    let d = 0.5 * lambda0; // centre-to-centre spacing
    let inset = 0.25 * l; // the A.3-measured 50 Ω depth on this stack

    let h = substrate.height_m;
    let w50 = microstrip_width(feed_z0_ohm, substrate.eps_r, h);
    let w_x = microstrip_width(feed_z0_ohm * std::f64::consts::SQRT_2, substrate.eps_r, h);
    let e_x = eps_eff(w_x, h, substrate.eps_r);
    let l_x = C / (f0_hz * e_x.sqrt()) / 4.0; // λg/4 at 70.7 Ω
    let e50 = eps_eff(w50, h, substrate.eps_r);
    let spine_len = C / (2.0 * f0_hz * e50.sqrt()); // λg/2: probe room

    // Vertical tree centreline, left of the patches with room for the
    // horizontal per-element runs.
    let x_v = -(4.0 * w50);
    let x_port = x_v - spine_len;
    let gap = w50; // inset notch gap each side of a feed (A.1 convention)

    let mut traces = vec![
        // Feed spine along x, centred on y = 0, overlapping the tree.
        Polygon::rect(x_port, -w50 / 2.0, spine_len + w_x, w50),
    ];
    for sgn in [1.0_f64, -1.0] {
        let yc = sgn * d / 2.0; // this element's centreline
        // λg/4 transformer from the junction outward.
        let y_lo = if sgn > 0.0 { -w50 / 2.0 } else { -l_x };
        traces.push(Polygon::rect(x_v - w_x / 2.0, y_lo, w_x, l_x + w50 / 2.0));
        // 50 Ω vertical from the transformer end to the element row.
        let seg = d / 2.0 - l_x + w50;
        let y_lo = if sgn > 0.0 {
            l_x - w50 / 2.0
        } else {
            -(l_x - w50 / 2.0) - seg
        };
        traces.push(Polygon::rect(x_v - w50 / 2.0, y_lo, w50, seg));
        // Horizontal 50 Ω feed into the element's inset.
        traces.push(Polygon::rect(
            x_v - w50 / 2.0,
            yc - w50 / 2.0,
            inset - (x_v - w50 / 2.0),
            w50,
        ));
        // The element (the 3 non-feed rects of the A.1 inset construction,
        // translated to yc).
        traces.push(Polygon::rect(
            0.0,
            yc + w50 / 2.0 + gap,
            l,
            w / 2.0 - w50 / 2.0 - gap,
        ));
        traces.push(Polygon::rect(
            0.0,
            yc - w / 2.0,
            l,
            w / 2.0 - w50 / 2.0 - gap,
        ));
        traces.push(Polygon::rect(
            inset,
            yc - (w50 / 2.0 + gap),
            l - inset,
            w50 + 2.0 * gap,
        ));
    }

    let port = PortRef {
        at: Point2::new(x_port, 0.0),
        width_m: w50,
        ref_impedance_ohm: feed_z0_ohm,
    };
    let bbox = BBox::from_polygons(&traces);
    PatchArray {
        layout: Layout {
            substrate: *substrate,
            traces,
            ports: vec![port],
            bbox,
        },
        dims: PatchArrayDims {
            patch,
            spacing_m: d,
            feed_width_m: w50,
            xfmr_width_m: w_x,
            xfmr_len_m: l_x,
            inset_m: inset,
            x_junction_m: x_v,
        },
    }
}

#[cfg(test)]
mod patch_array_tests {
    use super::*;

    fn fr4() -> Substrate {
        Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        }
    }

    fn aabb(p: &Polygon) -> (f64, f64, f64, f64) {
        let (mut x0, mut y0, mut x1, mut y1) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for v in &p.verts {
            x0 = x0.min(v.x);
            y0 = y0.min(v.y);
            x1 = x1.max(v.x);
            y1 = y1.max(v.y);
        }
        (x0, y0, x1, y1)
    }

    #[test]
    fn array_is_mirror_symmetric_and_connected() {
        let pa = patch_array_2x1(2.45e9, &fr4(), 50.0);
        let boxes: Vec<_> = pa.layout.traces.iter().map(aabb).collect();
        // Mirror symmetry about y = 0: every box's mirror is present.
        for b in &boxes {
            let mirrored = (b.0, -b.3, b.2, -b.1);
            assert!(
                boxes.iter().any(|c| (c.0 - mirrored.0).abs() < 1e-12
                    && (c.1 - mirrored.1).abs() < 1e-12
                    && (c.2 - mirrored.2).abs() < 1e-12
                    && (c.3 - mirrored.3).abs() < 1e-12),
                "no mirror partner for {b:?}"
            );
        }
        // Single connected component. Unlike the quasi-Yagi joints (drawn
        // segments that must genuinely overlap), the patch bands reuse the
        // certified inset construction whose rects share EXACT edges (same
        // arithmetic expressions) — contiguous metal under rasterization.
        // So connectivity here = closed-interval touch with positive
        // shared length (corner-point contact does not count).
        let overlaps = |a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)| {
            let x_touch = a.0 <= b.2 + 1e-12 && b.0 <= a.2 + 1e-12;
            let y_touch = a.1 <= b.3 + 1e-12 && b.1 <= a.3 + 1e-12;
            let x_len = a.2.min(b.2) - a.0.max(b.0);
            let y_len = a.3.min(b.3) - a.1.max(b.1);
            x_touch && y_touch && (x_len > 1e-12 || y_len > 1e-12)
        };
        let n = boxes.len();
        let mut comp: Vec<usize> = (0..n).collect();
        fn root(comp: &mut [usize], mut i: usize) -> usize {
            while comp[i] != i {
                comp[i] = comp[comp[i]];
                i = comp[i];
            }
            i
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if overlaps(boxes[i], boxes[j]) {
                    let (ri, rj) = (root(&mut comp, i), root(&mut comp, j));
                    comp[ri] = rj;
                }
            }
        }
        let r0 = root(&mut comp, 0);
        for i in 1..n {
            assert_eq!(root(&mut comp, i), r0, "trace {i} disconnected");
        }
    }

    #[test]
    fn seed_dims_match_the_closed_forms() {
        let pa = patch_array_2x1(2.45e9, &fr4(), 50.0);
        let d = pa.dims;
        let single = patch_antenna_dims(2.45e9, 4.4, 1.6e-3);
        assert!((d.patch.length_m - single.length_m).abs() < 1e-12);
        assert!((d.spacing_m - 0.5 * 299_792_458.0 / 2.45e9).abs() < 1e-12);
        // λg/4 at 70.7 Ω.
        let w_x = microstrip_width(50.0 * std::f64::consts::SQRT_2, 4.4, 1.6e-3);
        let lg = 299_792_458.0 / (2.45e9 * eps_eff(w_x, 1.6e-3, 4.4).sqrt());
        assert!((d.xfmr_len_m - lg / 4.0).abs() < 1e-12);
        assert!((d.inset_m - 0.25 * single.length_m).abs() < 1e-12);
        // Elements clear the tree: junction left of the patches.
        assert!(d.x_junction_m < 0.0);
        // The transformer must end before the element row.
        assert!(d.xfmr_len_m < d.spacing_m / 2.0);
    }
}

#[cfg(test)]
mod quasi_yagi_tests {
    use super::*;

    fn fr4() -> Substrate {
        Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        }
    }

    fn aabb(p: &Polygon) -> (f64, f64, f64, f64) {
        let (mut x0, mut y0, mut x1, mut y1) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for v in &p.verts {
            x0 = x0.min(v.x);
            y0 = y0.min(v.y);
            x1 = x1.max(v.x);
            y1 = y1.max(v.y);
        }
        (x0, y0, x1, y1)
    }

    fn overlaps(a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)) -> bool {
        // Positive-area overlap (shared edges alone rasterize
        // inconsistently, so joints must genuinely overlap).
        a.0 < b.2 - 1e-12 && b.0 < a.2 - 1e-12 && a.1 < b.3 - 1e-12 && b.1 < a.3 - 1e-12
    }

    #[test]
    fn feed_chain_is_connected_and_director_is_parasitic() {
        let qy = quasi_yagi(5.8e9, &fr4(), 50.0);
        let boxes: Vec<_> = qy.layout.traces.iter().map(aabb).collect();
        let n = boxes.len();
        // Union-find over positive-area overlaps.
        let mut comp: Vec<usize> = (0..n).collect();
        fn root(comp: &mut [usize], mut i: usize) -> usize {
            while comp[i] != i {
                comp[i] = comp[comp[i]];
                i = comp[i];
            }
            i
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if overlaps(boxes[i], boxes[j]) {
                    let (ri, rj) = (root(&mut comp, i), root(&mut comp, j));
                    comp[ri] = rj;
                }
            }
        }
        let feed_root = root(&mut comp, 0);
        for i in 1..n - 1 {
            assert_eq!(
                root(&mut comp, i),
                feed_root,
                "trace {i} disconnected from the feed chain"
            );
        }
        // The director (last trace) is deliberately parasitic.
        assert_ne!(
            root(&mut comp, n - 1),
            feed_root,
            "director must NOT touch the driven structure"
        );
    }

    #[test]
    fn geometry_ordering_and_dipole_length() {
        let qy = quasi_yagi(5.8e9, &fr4(), 50.0);
        let d = qy.dims;
        // Feed → T → ground edge → dipole → director along +x.
        assert!(d.x_gnd_m > 0.0 && d.x_dipole_m > d.x_gnd_m);
        assert!((d.x_dipole_m - d.x_gnd_m - d.reflector_gap_m).abs() < 1e-12);
        // Every balun/feed rect stays over the ground (microstrip needs
        // its return plane); only arms + director may cross the edge.
        let boxes: Vec<_> = qy.layout.traces.iter().map(aabb).collect();
        let n = boxes.len();
        for (i, b) in boxes.iter().enumerate().take(n - 3) {
            assert!(
                b.2 <= d.x_gnd_m + d.branch_width_m / 2.0 + d.reflector_gap_m + 1e-12,
                "trace {i} extends past the dipole plane"
            );
        }
        // Balun detour fully over ground.
        for (i, b) in boxes.iter().enumerate().take(8).skip(3) {
            if i < 7 {
                assert!(
                    b.2 <= d.x_gnd_m + 1e-12,
                    "balun trace {i} crosses the ground edge (x1 = {}, x_gnd = {})",
                    b.2,
                    d.x_gnd_m
                );
            }
        }
        // Dipole tip-to-tip = designed length (arms are traces n-3, n-2).
        let up = boxes[n - 3];
        let dn = boxes[n - 2];
        assert!((up.3 - d.dipole_len_m / 2.0).abs() < 1e-12, "upper tip");
        assert!((dn.1 + d.dipole_len_m / 2.0).abs() < 1e-12, "lower tip");
        // Director shorter than the driven element, further out.
        assert!(d.director_len_m < d.dipole_len_m);
        assert!(boxes[n - 1].0 > d.x_dipole_m);
        // Port sits on the feed.
        let p = &qy.layout.ports[0];
        assert!(p.at.x.abs() < 1e-12 && p.at.y.abs() < 1e-12);
    }
}

#[cfg(test)]
mod hairpin_sections_tests {
    use super::*;

    fn substrate() -> Substrate {
        Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        }
    }

    #[test]
    fn per_section_gaps_set_each_resonator_pitch() {
        let p = HairpinSectionParams {
            substrate: substrate(),
            arm_length_m: 8e-3,
            line_width_m: 1e-3,
            fold_spacing_m: 3e-3,
            gaps_m: vec![0.4e-3, 0.9e-3],
            tap_offset_m: 2e-3,
            feed_width_m: 1e-3,
            feed_length_m: 8e-3,
        };
        let layout = hairpin_bpf_sections(&p);
        // 3 resonators × 3 rects + 2 feeds.
        assert_eq!(layout.traces.len(), 11);
        // Resonator i's left arm is trace 3i; x-bases accumulate each gap.
        let min_x = |poly: &Polygon| poly.verts.iter().map(|v| v.x).fold(f64::INFINITY, f64::min);
        let base = |i: usize| min_x(&layout.traces[3 * i]);
        let pitch0 = p.fold_spacing_m + p.line_width_m + p.gaps_m[0];
        let pitch1 = p.fold_spacing_m + p.line_width_m + p.gaps_m[1];
        assert!((base(1) - base(0) - pitch0).abs() < 1e-12);
        assert!((base(2) - base(1) - pitch1).abs() < 1e-12);
        // Ports at tap height, outside the end resonators.
        assert!((layout.ports[0].at.y - p.tap_offset_m).abs() < 1e-12);
        assert!((layout.ports[1].at.y - p.tap_offset_m).abs() < 1e-12);
        assert!(layout.ports[0].at.x < base(0));
        assert!(layout.ports[1].at.x > base(2) + p.fold_spacing_m + p.line_width_m);
    }

    #[test]
    fn uniform_gaps_reproduce_hairpin_bpf_resonator_placement() {
        let gap = 0.6e-3;
        let sections = HairpinSectionParams {
            substrate: substrate(),
            arm_length_m: 8e-3,
            line_width_m: 1e-3,
            fold_spacing_m: 3e-3,
            gaps_m: vec![gap, gap],
            tap_offset_m: 8e-3 / 3.0,
            feed_width_m: 1e-3,
            feed_length_m: 8e-3,
        };
        let uniform = HairpinParams {
            substrate: substrate(),
            n: 3,
            arm_length_m: 8e-3,
            line_width_m: 1e-3,
            fold_spacing_m: 3e-3,
            coupling_gap_m: gap,
            tap_offset_m: 8e-3 / 3.0,
            feed_width_m: 1e-3,
            feed_length_m: 8e-3,
        };
        let a = hairpin_bpf_sections(&sections);
        let b = hairpin_bpf(&uniform);
        assert_eq!(a.traces.len(), b.traces.len());
        for (ra, rb) in a.traces.iter().zip(&b.traces) {
            assert_eq!(ra.verts.len(), rb.verts.len());
            for (va, vb) in ra.verts.iter().zip(&rb.verts) {
                assert!((va.x - vb.x).abs() < 1e-15);
                assert!((va.y - vb.y).abs() < 1e-15);
            }
        }
    }
}

#[cfg(test)]
mod patch_tests {
    use super::*;

    #[test]
    fn patch_dims_match_hand_computed_balanis_values() {
        // 2.45 GHz on FR-4 (eps_r 4.4, h 1.6 mm) — Balanis worked-example
        // arithmetic: W = c/(2 f0) sqrt(2/5.4) = 37.26 mm;
        // eps_eff(W) ~ 4.09; L = c/(2 f0 sqrt(eps_eff)) - 2 dL ~ 28.8 mm.
        let d = patch_antenna_dims(2.45e9, 4.4, 1.6e-3);
        assert!((d.width_m - 37.26e-3).abs() < 0.05e-3, "W = {}", d.width_m);
        assert!((d.eps_eff - 4.09).abs() < 0.03, "eps_eff = {}", d.eps_eff);
        assert!((d.length_m - 28.8e-3).abs() < 0.3e-3, "L = {}", d.length_m);
        assert!(d.delta_l_m > 0.5e-3 && d.delta_l_m < 1.0e-3);
    }

    #[test]
    fn edge_fed_patch_layout_is_connected_and_ported() {
        let sub = Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        };
        let layout = edge_fed_patch(2.45e9, &sub, 50.0);
        assert_eq!(layout.traces.len(), 2);
        assert_eq!(layout.ports.len(), 1);
        // Feed end and patch start meet at x = 0 (connected footprint).
        assert!(layout.bbox.min.x < 0.0 && layout.bbox.max.x > 20e-3);
        assert!(layout.ports[0].at.x == layout.bbox.min.x);
    }
}

/// Balanis single-slot radiating conductance `G1` (siemens) for a patch of
/// width `w_m` at free-space wavelength `lambda0_m` (Antenna Theory
/// eq. 14-8a, piecewise in `W/lambda0`).
pub fn patch_slot_conductance(w_m: f64, lambda0_m: f64) -> f64 {
    let u = w_m / lambda0_m;
    if u < 0.35 {
        u * u / 90.0
    } else if u < 2.0 {
        u / 120.0 - 1.0 / (60.0 * std::f64::consts::PI * std::f64::consts::PI)
    } else {
        u / 120.0
    }
}

/// Assemble an **inset-fed** rectangular-patch [`Layout`] matched to
/// `feed_z0_ohm` (A.1, ADR-0191): edge resistance `R_edge = 1/(2 G1)`
/// from the Balanis slot-conductance model (the mutual term `G12` is
/// neglected — documented walking-skeleton approximation, it lowers
/// `R_edge` by ~20 % and the inset only moves as acos(sqrt(...))), inset
/// depth `x0 = (L/pi) * acos(sqrt(Z0/R_edge))` from the cos^2 cavity
/// current profile.
///
/// Metal is a union of four rectangles: two outer patch bands flanking
/// the notch, the centre band beyond the inset depth, and the feed line
/// running through the notch (gap = one feed width each side). One
/// `PortRef` at the outer feed end; the feed is half a guided wavelength
/// so measurement probes fit upstream.
pub fn inset_fed_patch(f0_hz: f64, substrate: &Substrate, feed_z0_ohm: f64) -> Layout {
    let dims = patch_antenna_dims(f0_hz, substrate.eps_r, substrate.height_m);
    let lambda0 = 299_792_458.0 / f0_hz;
    let g1 = patch_slot_conductance(dims.width_m, lambda0);
    let r_edge = 1.0 / (2.0 * g1);
    assert!(
        feed_z0_ohm < r_edge,
        "inset_fed_patch: Z0 {feed_z0_ohm} >= edge resistance {r_edge} — no inset match exists"
    );
    let x_inset = (dims.length_m / std::f64::consts::PI) * (feed_z0_ohm / r_edge).sqrt().acos();
    inset_fed_patch_with_depth(f0_hz, substrate, feed_z0_ohm, x_inset)
}

/// [`inset_fed_patch`] with an **explicit** inset depth (metres) — the
/// design-loop knob (A.3): the closed-form seed's depth comes from the
/// G1-only slot-conductance model, which overestimates the edge
/// resistance on thick / high-ε_r substrates (measured in ADR-0191), so
/// the engine tunes this depth against the measured return loss.
pub fn inset_fed_patch_with_depth(
    f0_hz: f64,
    substrate: &Substrate,
    feed_z0_ohm: f64,
    x_inset_m: f64,
) -> Layout {
    const C: f64 = 299_792_458.0;
    let dims = patch_antenna_dims(f0_hz, substrate.eps_r, substrate.height_m);
    let (w, l) = (dims.width_m, dims.length_m);
    assert!(
        x_inset_m > 0.0 && x_inset_m < l / 2.0,
        "inset depth {x_inset_m} m outside (0, L/2)"
    );
    let x_inset = x_inset_m;

    let feed_w = microstrip_width(feed_z0_ohm, substrate.eps_r, substrate.height_m);
    let feed_e_eff = eps_eff(feed_w, substrate.height_m, substrate.eps_r);
    let feed_len = C / (2.0 * f0_hz * feed_e_eff.sqrt());
    let gap = feed_w; // notch gap each side of the feed

    let traces = vec![
        // Feed line: from the port, through the notch, joining the patch
        // interior at the inset depth.
        Polygon::rect(-feed_len, -feed_w / 2.0, feed_len + x_inset, feed_w),
        // Outer patch band, +y side of the notch.
        Polygon::rect(0.0, feed_w / 2.0 + gap, l, w / 2.0 - feed_w / 2.0 - gap),
        // Outer patch band, -y side of the notch.
        Polygon::rect(0.0, -w / 2.0, l, w / 2.0 - feed_w / 2.0 - gap),
        // Centre band beyond the inset (between the notch slots' far end
        // and the patch's far edge).
        Polygon::rect(
            x_inset,
            -(feed_w / 2.0 + gap),
            l - x_inset,
            feed_w + 2.0 * gap,
        ),
    ];
    let port = PortRef {
        at: Point2::new(-feed_len, 0.0),
        width_m: feed_w,
        ref_impedance_ohm: feed_z0_ohm,
    };
    let bbox = BBox::from_polygons(&traces);
    Layout {
        substrate: *substrate,
        traces,
        ports: vec![port],
        bbox,
    }
}

#[cfg(test)]
mod inset_patch_tests {
    use super::*;

    #[test]
    fn inset_depth_matches_hand_computed_balanis_arithmetic() {
        // 2.45 GHz FR-4: W = 37.26 mm, lambda0 = 122.36 mm, W/lambda0 =
        // 0.3045 < 0.35 -> G1 = 0.3045^2/90 = 1.030e-3 S ->
        // R_edge = 485.4 ohm; x0 = (L/pi) acos(sqrt(50/485.4))
        //        = (28.8e-3/pi) * acos(0.3210) = 9.17 mm * 1.2440 = 11.41 mm.
        let sub = Substrate {
            eps_r: 4.4,
            height_m: 1.6e-3,
            loss_tangent: 0.0,
            metal_thickness_m: 35e-6,
        };
        let g1 = patch_slot_conductance(37.26e-3, 299_792_458.0 / 2.45e9);
        assert!((g1 - 1.030e-3).abs() < 0.01e-3, "G1 = {g1}");
        let layout = inset_fed_patch(2.45e9, &sub, 50.0);
        // Recover the inset from the geometry: the centre band starts at
        // x_inset (its rect is traces[3]).
        // Polygon::rect stores vertices; use the bbox of that polygon via
        // its first vertex x.
        let x_inset = layout.traces[3].verts[0].x;
        assert!((x_inset - 11.4e-3).abs() < 0.3e-3, "x_inset = {x_inset}");
        assert_eq!(layout.ports.len(), 1);
        assert_eq!(layout.traces.len(), 4);
    }
}

// ===========================================================================
// Multilayer stackups — FS.4.0, ADR-0215
// ===========================================================================

/// One dielectric layer of a [`Stackup`], metres/dimensionless.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StackupLayer {
    /// Relative permittivity `ε_r` of this layer.
    pub eps_r: f64,
    /// Layer thickness, metres.
    pub height_m: f64,
    /// Dielectric loss tangent `tan δ` (dimensionless).
    pub loss_tangent: f64,
}

/// An N-layer dielectric stackup (FS.4.0, ADR-0215): layers bottom-up,
/// ground plane below `layers[0]`, optional PEC lid directly on top of
/// the last layer (stripline / shielded boards). `lid: false` leaves the
/// top open to air — the microstrip case.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stackup {
    /// Dielectric layers, bottom (on the ground plane) to top.
    pub layers: Vec<StackupLayer>,
    /// PEC sheet directly above the last layer.
    pub lid: bool,
}

impl Stackup {
    /// Symmetric stripline: two identical `b/2` layers with the trace at
    /// the mid interface (`trace_layer = 0`), lid on. The canonical
    /// FS.4.0 validation stack — TEM in homogeneous dielectric, so
    /// `ε_eff = ε_r` exactly.
    pub fn symmetric_stripline(eps_r: f64, b_m: f64) -> Self {
        let half = StackupLayer {
            eps_r,
            height_m: b_m / 2.0,
            loss_tangent: 0.0,
        };
        Self {
            layers: vec![half, half],
            lid: true,
        }
    }

    /// Total dielectric thickness, metres.
    pub fn total_height_m(&self) -> f64 {
        self.layers.iter().map(|l| l.height_m).sum()
    }
}

#[cfg(test)]
mod stackup_tests {
    use super::*;

    #[test]
    fn symmetric_stripline_is_two_equal_layers_with_lid() {
        let s = Stackup::symmetric_stripline(4.4, 3.2e-3);
        assert_eq!(s.layers.len(), 2);
        assert_eq!(s.layers[0], s.layers[1]);
        assert!(s.lid);
        assert!((s.total_height_m() - 3.2e-3).abs() < 1e-18);
    }

    #[test]
    fn stackup_round_trips_through_json() {
        let s = Stackup {
            layers: vec![
                StackupLayer {
                    eps_r: 4.4,
                    height_m: 0.2e-3,
                    loss_tangent: 0.02,
                },
                StackupLayer {
                    eps_r: 3.5,
                    height_m: 0.8e-3,
                    loss_tangent: 0.001,
                },
            ],
            lid: false,
        };
        let json = serde_json::to_string(&s).unwrap();
        assert_eq!(serde_json::from_str::<Stackup>(&json).unwrap(), s);
    }
}
