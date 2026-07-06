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
