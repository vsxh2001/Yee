//! RS-274X Gerber **import** (FS.3.0, ADR-0209; extended FS.3.2b,
//! ADR-0220): the region-fill subset our own writer
//! ([`crate::layout_to_gerber`]) emits, plus circular-arc region
//! segments (`G02`/`G03` under `G75` multi-quadrant mode) and flashed
//! `C`/`R` apertures (`D03`), parsed back into [`Polygon`]s — the
//! "bring your own board" door. Everything outside the subset is
//! rejected with an explicit error, never silently mis-parsed (still
//! out: polarity, macro/obround apertures, step-repeat, inches,
//! stroked-draw copper, single-quadrant `G74` arcs).
//!
//! Linear coordinates are exact: the writer's 4.6 fixed-point
//! millimetre words map back to metres as `n·1e-9`, so
//! `export → import → export` is byte-identical (gate `gerber-rt-001`).
//! Arc *endpoints* and rect-flash corners are equally exact; arc and
//! circle-flash interiors are tessellated at the pinned chord tolerance
//! [`ARC_CHORD_TOL_M`] (gate `gerber-rt-003`).

use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

use yee_layout::{Layout, Point2, Polygon, PortRef, Substrate};

/// Chord tolerance for tessellating circular geometry (arcs and circle
/// flashes), in metres: the maximum sagitta — perpendicular distance
/// from any chord to the true arc — is bounded by this value (1 µm, two
/// orders below any λ/20 cell the suite meshes).
///
/// The sagitta of a chord subtending angle `φ` on a circle of radius
/// `r` is `r·(1 − cos(φ/2))`, so the maximum angular step is
/// `φ_max = 2·acos(1 − tol/r)` and an arc sweeping `θ` is split into
/// `n = ceil(θ/φ_max)` uniform segments (min 1; full circles min 4).
pub const ARC_CHORD_TOL_M: f64 = 1.0e-6;

/// Maximum angular step (radians) keeping the chord sagitta of a circle
/// of radius `r_m` below [`ARC_CHORD_TOL_M`]. See the const docs for
/// the formula; degenerate radii clamp to a half-turn.
fn max_angle_step(r_m: f64) -> f64 {
    let c = 1.0 - ARC_CHORD_TOL_M / r_m;
    if c <= -1.0 { PI } else { 2.0 * c.acos() }
}

/// A parsed `%AD…%` aperture definition (millimetre parameters
/// converted to metres). Templates outside `C`/`R` are bookkept as
/// [`Aperture::Other`] and rejected only if actually flashed.
#[derive(Debug, Clone, PartialEq)]
enum Aperture {
    /// Standard circle `C,<dia>` (no hole).
    Circle {
        /// Diameter in metres.
        d_m: f64,
    },
    /// Standard rectangle `R,<x>X<y>` (no hole).
    Rect {
        /// X extent in metres.
        x_m: f64,
        /// Y extent in metres.
        y_m: f64,
    },
    /// Anything else (obround, polygon, macro, holed C/R): kept with
    /// its raw definition text so a flash of it can name the offender.
    Other(String),
}

/// Linear vs circular interpolation mode (`G01`/`G02`/`G03`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Interp {
    /// `G01` — straight segments (the default).
    Linear,
    /// `G02` — clockwise circular interpolation.
    Cw,
    /// `G03` — counter-clockwise circular interpolation.
    Ccw,
}

/// Why a Gerber file could not be imported.
#[derive(Debug, Clone, PartialEq)]
pub enum GerberImportError {
    /// `%MOIN*%` — imperial units are outside the FS.3.0 subset.
    ImperialUnits,
    /// A word or extended command outside the supported subset, with the
    /// offending text.
    UnsupportedCommand(String),
    /// A `D01*` draw before any `D02*` move inside a region.
    DrawBeforeMove,
    /// The file ended inside an open `G36*` region.
    UnclosedRegion,
    /// A coordinate word failed to parse, with the offending text.
    BadCoordinate(String),
    /// The file parsed but contained no copper regions — a `Layout`
    /// needs at least one polygon (its bbox is undefined otherwise).
    NoCopper,
    /// A `D03*` flash inside an open `G36` region — forbidden by the
    /// Gerber specification.
    FlashInRegion,
    /// A `D03*` flash with no aperture selected, or with a D-code that
    /// was never defined by an `%AD…%`, with the offending text.
    UnknownAperture(String),
    /// A `D03*` flash of an aperture outside the supported `C`/`R`
    /// no-hole templates (obround, polygon, macro, holed, zero-size),
    /// with the offending definition.
    UnsupportedAperture(String),
    /// A circular-arc segment with broken geometry (zero radius, or
    /// start/end radii disagreeing beyond 10× the chord tolerance),
    /// with the offending text.
    BadArc(String),
}

impl std::fmt::Display for GerberImportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ImperialUnits => write!(f, "imperial units (%MOIN%) are not supported"),
            Self::UnsupportedCommand(w) => write!(f, "unsupported Gerber command: {w:?}"),
            Self::DrawBeforeMove => write!(f, "D01 draw before any D02 move in a region"),
            Self::UnclosedRegion => write!(f, "file ended inside an open G36 region"),
            Self::BadCoordinate(w) => write!(f, "malformed coordinate word: {w:?}"),
            Self::NoCopper => write!(f, "no copper regions found in the Gerber"),
            Self::FlashInRegion => write!(f, "D03 flash inside a G36 region is forbidden"),
            Self::UnknownAperture(w) => write!(f, "flash of unknown aperture: {w:?}"),
            Self::UnsupportedAperture(w) => write!(f, "flash of unsupported aperture: {w:?}"),
            Self::BadArc(w) => write!(f, "malformed arc segment: {w:?}"),
        }
    }
}

impl std::error::Error for GerberImportError {}

/// Parse the supported subset of an RS-274X Gerber into polygons
/// (metres): region fills (linear + circular-arc segments) and `C`/`R`
/// aperture flashes, in file order. Region contours drop the explicit
/// closing vertex, matching the [`Polygon`] convention (the writer
/// re-adds it on export, so the linear round-trip is byte-stable).
pub fn gerber_to_polygons(gerber: &str) -> Result<Vec<Polygon>, GerberImportError> {
    let mut polys = Vec::new();
    let mut in_region = false;
    let mut contour: Vec<Point2> = Vec::new();
    // Modal coordinate state (the subset always emits both X and Y, but
    // modal omission is core Gerber and cheap to honour).
    let (mut cur_x, mut cur_y) = (0.0_f64, 0.0_f64);
    // FS.3.2b state: aperture dictionary + selection, interpolation
    // mode, and whether multi-quadrant arc mode was declared.
    let mut apertures: HashMap<u32, Aperture> = HashMap::new();
    let mut current_aperture: Option<u32> = None;
    let mut interp = Interp::Linear;
    let mut multi_quadrant = false;

    for raw in gerber.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        // Extended commands: %...*%
        if let Some(ext) = line.strip_prefix('%').and_then(|l| l.strip_suffix("*%")) {
            if ext.starts_with("FSLA") {
                // Any AXmnYmn format is fine: we parse digits explicitly.
                continue;
            }
            if ext == "MOMM" {
                continue;
            }
            if ext == "MOIN" {
                return Err(GerberImportError::ImperialUnits);
            }
            if ext.starts_with("AD") {
                // Aperture definition: bookkept per D-code. Shapes we
                // cannot flash are stored as `Other` and rejected only
                // if a D03 actually uses them.
                if let Some((code, ap)) = parse_aperture_def(ext) {
                    apertures.insert(code, ap);
                }
                continue;
            }
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        }
        let Some(mut word) = line.strip_suffix('*') else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        // Inline interpolation-mode prefix on a coordinate word
        // (`G02X…I…J…D01`) — set the mode, keep parsing the rest.
        if word.len() > 3 {
            if let Some(rest) = word.strip_prefix("G01") {
                interp = Interp::Linear;
                word = rest;
            } else if let Some(rest) = word.strip_prefix("G02") {
                interp = Interp::Cw;
                word = rest;
            } else if let Some(rest) = word.strip_prefix("G03") {
                interp = Interp::Ccw;
                word = rest;
            }
        }
        match word {
            _ if word.starts_with("G04") => continue, // comment
            "G01" => {
                interp = Interp::Linear;
                continue;
            }
            "G02" => {
                interp = Interp::Cw;
                continue;
            }
            "G03" => {
                interp = Interp::Ccw;
                continue;
            }
            "G75" => {
                multi_quadrant = true;
                continue;
            }
            // Single-quadrant arc mode: legacy and ambiguous — named
            // rejection, never a mis-parse.
            "G74" => return Err(GerberImportError::UnsupportedCommand(line.into())),
            "G36" => {
                in_region = true;
                contour.clear();
                continue;
            }
            "G37" => {
                in_region = false;
                let mut verts = std::mem::take(&mut contour);
                // Drop the explicit closing vertex (Polygon convention).
                if verts.len() > 1 && verts[verts.len() - 1] == verts[0] {
                    verts.pop();
                }
                if verts.len() >= 3 {
                    polys.push(Polygon { verts });
                }
                continue;
            }
            "M02" => break,
            _ if word.starts_with('D')
                && word[1..].chars().all(|c| c.is_ascii_digit())
                && word[1..].parse::<u32>().is_ok_and(|c| c >= 10) =>
            {
                // Aperture select (D10 etc.): now bookkept for flashes.
                current_aperture = word[1..].parse().ok();
                continue;
            }
            _ => {}
        }
        // Coordinate words: [X<int>][Y<int>][I<int>][J<int>]D01|D02|D03
        // (bare D01/D02/D03 reuse the modal point).
        let (coords, op) = if let Some(c) = word.strip_suffix("D01") {
            (c, 1)
        } else if let Some(c) = word.strip_suffix("D02") {
            (c, 2)
        } else if let Some(c) = word.strip_suffix("D03") {
            (c, 3)
        } else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        let (start_x, start_y) = (cur_x, cur_y);
        let (i_off, j_off) = parse_coords_ij(coords, &mut cur_x, &mut cur_y)
            .map_err(|_| GerberImportError::BadCoordinate(line.into()))?;
        if op == 3 {
            // Flash. Forbidden inside a region by the Gerber spec.
            if in_region {
                return Err(GerberImportError::FlashInRegion);
            }
            let ap = current_aperture
                .and_then(|code| apertures.get(&code))
                .ok_or_else(|| GerberImportError::UnknownAperture(line.into()))?;
            polys.push(flash_polygon(ap, cur_x, cur_y)?);
            continue;
        }
        if !in_region {
            // Stroked draws (e.g. the outline layer) stay out of the
            // copper subset — linear or circular; a plain move outside
            // a region is harmless.
            if op == 1 {
                return Err(GerberImportError::UnsupportedCommand(line.into()));
            }
            continue;
        }
        match op {
            2 => {
                contour.clear();
                contour.push(Point2::new(cur_x, cur_y));
            }
            _ => {
                if contour.is_empty() {
                    return Err(GerberImportError::DrawBeforeMove);
                }
                match interp {
                    Interp::Linear => contour.push(Point2::new(cur_x, cur_y)),
                    Interp::Cw | Interp::Ccw => {
                        // Single-quadrant semantics would apply without
                        // G75 — reject rather than guess.
                        if !multi_quadrant {
                            return Err(GerberImportError::UnsupportedCommand(line.into()));
                        }
                        contour.extend(arc_vertices(
                            (start_x, start_y),
                            (cur_x, cur_y),
                            (start_x + i_off, start_y + j_off),
                            interp == Interp::Ccw,
                            line,
                        )?);
                    }
                }
            }
        }
    }
    if in_region {
        return Err(GerberImportError::UnclosedRegion);
    }
    Ok(polys)
}

/// Parse an `%AD…%` body (`ADD<code><template>,<params>`) into an
/// [`Aperture`] (mm → metres). Standard `C,<dia>` and `R,<x>X<y>`
/// without holes become flashable shapes; anything else — holed,
/// obround, polygon, macro, unparseable params — is stored as
/// [`Aperture::Other`] with the raw text so a flash of it can name the
/// offender. `None` only when no D-code can be extracted (the
/// definition is then simply ignored; flashing it is
/// [`GerberImportError::UnknownAperture`]).
fn parse_aperture_def(ext: &str) -> Option<(u32, Aperture)> {
    let body = ext.strip_prefix("ADD")?;
    let digits_end = body.find(|c: char| !c.is_ascii_digit())?;
    let code: u32 = body[..digits_end].parse().ok()?;
    if code < 10 {
        return None;
    }
    let rest = &body[digits_end..];
    let other = || Aperture::Other(ext.to_string());
    let ap = match rest.split_once(',') {
        Some(("C", params)) => {
            let p: Vec<Option<f64>> = params.split('X').map(|s| s.parse().ok()).collect();
            match p.as_slice() {
                [Some(d_mm)] if *d_mm > 0.0 => Aperture::Circle { d_m: d_mm * 1.0e-3 },
                _ => other(), // holed, zero-size, or malformed
            }
        }
        Some(("R", params)) => {
            let p: Vec<Option<f64>> = params.split('X').map(|s| s.parse().ok()).collect();
            match p.as_slice() {
                [Some(x_mm), Some(y_mm)] if *x_mm > 0.0 && *y_mm > 0.0 => Aperture::Rect {
                    x_m: x_mm * 1.0e-3,
                    y_m: y_mm * 1.0e-3,
                },
                _ => other(),
            }
        }
        _ => other(), // obround / polygon / macro / no params
    };
    Some((code, ap))
}

/// Convert a flash of `ap` at `(x, y)` metres into a [`Polygon`]:
/// rectangles exactly (4 corners, CCW from lower-left), circles
/// tessellated at [`ARC_CHORD_TOL_M`] (vertex 0 on the +x axis, CCW,
/// minimum 4 vertices).
fn flash_polygon(ap: &Aperture, x: f64, y: f64) -> Result<Polygon, GerberImportError> {
    match ap {
        Aperture::Circle { d_m } => {
            let r = 0.5 * d_m;
            let n = ((TAU / max_angle_step(r)).ceil() as usize).max(4);
            let verts = (0..n)
                .map(|k| {
                    let a = TAU * k as f64 / n as f64;
                    Point2::new(x + r * a.cos(), y + r * a.sin())
                })
                .collect();
            Ok(Polygon { verts })
        }
        Aperture::Rect { x_m, y_m } => {
            let (hx, hy) = (0.5 * x_m, 0.5 * y_m);
            Ok(Polygon {
                verts: vec![
                    Point2::new(x - hx, y - hy),
                    Point2::new(x + hx, y - hy),
                    Point2::new(x + hx, y + hy),
                    Point2::new(x - hx, y + hy),
                ],
            })
        }
        Aperture::Other(def) => Err(GerberImportError::UnsupportedAperture(def.clone())),
    }
}

/// Tessellate one multi-quadrant circular-arc segment from `start` to
/// `end` about `center` (all metres) into contour vertices at the
/// [`ARC_CHORD_TOL_M`] chord tolerance.
///
/// The start vertex is already in the contour (pushed by the previous
/// D02/D01), so only interior vertices plus the **exact** `end` point
/// are returned. Interior vertices sit at uniform angular steps with
/// the radius linearly interpolated start → end (the two radii may
/// disagree by up to the fixed-point quantum; beyond 10× the chord
/// tolerance it is a [`GerberImportError::BadArc`]). `start == end`
/// (bit-exact, as identical coordinate words decode identically) is a
/// full 360° circle per the Ucamco G75 rule.
fn arc_vertices(
    start: (f64, f64),
    end: (f64, f64),
    center: (f64, f64),
    ccw: bool,
    line: &str,
) -> Result<Vec<Point2>, GerberImportError> {
    let r0 = (start.0 - center.0).hypot(start.1 - center.1);
    let r1 = (end.0 - center.0).hypot(end.1 - center.1);
    if r0 <= ARC_CHORD_TOL_M || r1 <= ARC_CHORD_TOL_M {
        return Err(GerberImportError::BadArc(line.into()));
    }
    if (r0 - r1).abs() > 10.0 * ARC_CHORD_TOL_M {
        return Err(GerberImportError::BadArc(line.into()));
    }
    let a0 = (start.1 - center.1).atan2(start.0 - center.0);
    let a1 = (end.1 - center.1).atan2(end.0 - center.0);
    // Sweep magnitude in (0, 2π]; direction is applied via `sgn`.
    let sweep = if start == end {
        TAU
    } else {
        let mut d = if ccw { a1 - a0 } else { a0 - a1 };
        while d <= 0.0 {
            d += TAU;
        }
        d
    };
    let n = ((sweep / max_angle_step(r0.max(r1))).ceil() as usize).max(1);
    let sgn = if ccw { 1.0 } else { -1.0 };
    let mut verts = Vec::with_capacity(n);
    for k in 1..n {
        let t = k as f64 / n as f64;
        let a = a0 + sgn * sweep * t;
        let r = r0 + (r1 - r0) * t;
        verts.push(Point2::new(center.0 + r * a.cos(), center.1 + r * a.sin()));
    }
    // The endpoint is placed exactly as decoded, never recomputed.
    verts.push(Point2::new(end.0, end.1));
    Ok(verts)
}

/// Parse `[X<int>][Y<int>][I<int>][J<int>]` (4.6 fixed-point
/// millimetres) — `X`/`Y` update the modal state exactly (`n·1e-9`
/// metres); `I`/`J` are **non-modal** arc-centre offsets returned per
/// word, defaulting to 0 when omitted (the Ucamco G75 rule). An empty
/// string is a valid modal reuse (bare `D01*`/`D02*`/`D03*`).
fn parse_coords_ij(s: &str, x: &mut f64, y: &mut f64) -> Result<(f64, f64), ()> {
    let (mut i, mut j) = (0.0_f64, 0.0_f64);
    let mut rest = s;
    while !rest.is_empty() {
        let (axis, tail) = rest.split_at(1);
        let end = tail
            .find(|c: char| !(c.is_ascii_digit() || c == '-' || c == '+'))
            .unwrap_or(tail.len());
        let (num, next) = tail.split_at(end);
        let v: i64 = num.parse().map_err(|_| ())?;
        let metres = v as f64 * 1.0e-9;
        match axis {
            "X" => *x = metres,
            "Y" => *y = metres,
            "I" => i = metres,
            "J" => j = metres,
            _ => return Err(()),
        }
        rest = next;
    }
    Ok((i, j))
}

/// Parse `[X<int>][Y<int>]` (4.6 fixed-point millimetres) into the modal
/// state, exactly: `n·1e-9` metres.
fn parse_coords(s: &str, x: &mut f64, y: &mut f64) -> Result<(), ()> {
    let mut rest = s;
    let mut saw_any = false;
    while !rest.is_empty() {
        let (axis, tail) = rest.split_at(1);
        let end = tail
            .find(|c: char| !(c.is_ascii_digit() || c == '-' || c == '+'))
            .unwrap_or(tail.len());
        let (num, next) = tail.split_at(end);
        let v: i64 = num.parse().map_err(|_| ())?;
        let metres = v as f64 * 1.0e-9;
        match axis {
            "X" => *x = metres,
            "Y" => *y = metres,
            _ => return Err(()),
        }
        saw_any = true;
        rest = next;
    }
    if saw_any { Ok(()) } else { Err(()) }
}

/// Parse a **stroked outline** Gerber (FS.3.1a) — the Edge.Cuts dialect
/// [`crate::layout_to_gerber_outline`] emits: same header family, one
/// thin aperture, then a single `D02` move + `D01` draw chain OUTSIDE any
/// region. Returns the path vertices in order with the explicit closing
/// vertex dropped. Regions (`G36`) in an outline file are rejected — a
/// board profile is a cut path, not copper.
pub fn gerber_to_outline(gerber: &str) -> Result<Vec<Point2>, GerberImportError> {
    let mut path: Vec<Point2> = Vec::new();
    let (mut cur_x, mut cur_y) = (0.0_f64, 0.0_f64);
    for raw in gerber.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(ext) = line.strip_prefix('%').and_then(|l| l.strip_suffix("*%")) {
            if ext.starts_with("FSLA") || ext == "MOMM" || ext.starts_with("AD") {
                continue;
            }
            if ext == "MOIN" {
                return Err(GerberImportError::ImperialUnits);
            }
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        }
        let Some(word) = line.strip_suffix('*') else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        match word {
            _ if word.starts_with("G04") => continue,
            "M02" => break,
            "G36" | "G37" => {
                return Err(GerberImportError::UnsupportedCommand(line.into()));
            }
            _ if word.starts_with('D') && word[1..].chars().all(|c| c.is_ascii_digit()) => {
                continue;
            }
            _ => {}
        }
        let (coords, op) = if let Some(c) = word.strip_suffix("D01") {
            (c, 1)
        } else if let Some(c) = word.strip_suffix("D02") {
            (c, 2)
        } else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        parse_coords(coords, &mut cur_x, &mut cur_y)
            .map_err(|_| GerberImportError::BadCoordinate(line.into()))?;
        match op {
            2 => {
                path.clear();
                path.push(Point2::new(cur_x, cur_y));
            }
            _ => {
                if path.is_empty() {
                    return Err(GerberImportError::DrawBeforeMove);
                }
                path.push(Point2::new(cur_x, cur_y));
            }
        }
    }
    if path.len() > 1 && path[path.len() - 1] == path[0] {
        path.pop();
    }
    Ok(path)
}

/// Wrap imported polygons into a [`Layout`]. Gerber carries no stackup or
/// port information — the caller provides both (the studio's import flow
/// asks the user; gates pass the known originals).
pub fn gerber_to_layout(
    gerber: &str,
    substrate: Substrate,
    ports: Vec<PortRef>,
) -> Result<Layout, GerberImportError> {
    let traces = gerber_to_polygons(gerber)?;
    if traces.is_empty() {
        return Err(GerberImportError::NoCopper);
    }
    let bbox = yee_layout::BBox::from_polygons(&traces);
    Ok(Layout {
        substrate,
        traces,
        ports,
        bbox,
    })
}
