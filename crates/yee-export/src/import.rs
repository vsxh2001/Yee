//! RS-274X Gerber **import** (FS.3.0, ADR-0209): the region-fill subset
//! our own writer ([`crate::layout_to_gerber`]) emits, parsed back into
//! [`Polygon`]s — the walking skeleton of the "bring your own board"
//! door. Everything outside the subset is rejected with an explicit
//! error, never silently mis-parsed (see the spec's non-goals: arcs,
//! polarity, macro apertures, step-repeat, inches).
//!
//! Coordinates are exact: the writer's 4.6 fixed-point millimetre words
//! map back to metres as `n·1e-9`, so `export → import → export` is
//! byte-identical (gate `gerber-rt-001`).

use yee_layout::{Layout, Point2, Polygon, PortRef, Substrate};

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
        }
    }
}

impl std::error::Error for GerberImportError {}

/// Parse the region-fill subset of an RS-274X Gerber into polygons
/// (metres). Contours drop the explicit closing vertex, matching the
/// [`Polygon`] convention (the writer re-adds it on export, so the
/// round-trip is byte-stable).
pub fn gerber_to_polygons(gerber: &str) -> Result<Vec<Polygon>, GerberImportError> {
    let mut polys = Vec::new();
    let mut in_region = false;
    let mut contour: Vec<Point2> = Vec::new();
    // Modal coordinate state (the subset always emits both X and Y, but
    // modal omission is core Gerber and cheap to honour).
    let (mut cur_x, mut cur_y) = (0.0_f64, 0.0_f64);

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
                // Aperture definitions: bookkept by code only; regions
                // ignore aperture geometry entirely.
                continue;
            }
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        }
        let Some(word) = line.strip_suffix('*') else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        match word {
            _ if word.starts_with("G04") => continue, // comment
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
            _ if word.starts_with('D') && word[1..].chars().all(|c| c.is_ascii_digit()) => {
                // Aperture select (D10 etc.) — regions ignore it.
                continue;
            }
            _ => {}
        }
        // Coordinate words: [X<int>][Y<int>]D01|D02
        let (coords, op) = if let Some(c) = word.strip_suffix("D01") {
            (c, 1)
        } else if let Some(c) = word.strip_suffix("D02") {
            (c, 2)
        } else {
            return Err(GerberImportError::UnsupportedCommand(line.into()));
        };
        parse_coords(coords, &mut cur_x, &mut cur_y)
            .map_err(|_| GerberImportError::BadCoordinate(line.into()))?;
        if !in_region {
            // Stroked draws (e.g. the outline layer) are FS.3.1; a plain
            // move outside a region is harmless.
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
                contour.push(Point2::new(cur_x, cur_y));
            }
        }
    }
    if in_region {
        return Err(GerberImportError::UnclosedRegion);
    }
    Ok(polys)
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
