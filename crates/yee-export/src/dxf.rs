//! ASCII DXF (R12+ group-code format) **import** (FS.3.3, ADR-0230): the
//! same subset-plus-named-rejections discipline as the Gerber importer
//! ([`crate::import`]), targeted at trace-outline geometry instead of
//! copper regions.
//!
//! Supported subset: closed `LWPOLYLINE` (straight segments + bulge
//! arcs) and closed R12-style `POLYLINE`/`VERTEX`/`SEQEND` chains, on
//! any layer (with an optional [`DxfOptions::layer`] filter), under a
//! `$INSUNITS` of millimetres (`4`) or inches (`1`). Everything outside
//! that subset is a named, typed rejection — never a silent
//! mis-parse: open polylines, `CIRCLE`/`ARC`/`ELLIPSE`/`SPLINE`
//! entities, `TEXT`/`MTEXT`/`DIMENSION`, `INSERT` (block references),
//! nonzero `Z` / elevation, and any `$INSUNITS` other than mm/inch
//! (**including a missing header variable** — DXF's own default for an
//! absent `$INSUNITS` is "unitless", which this importer never guesses
//! at; every file must say what it means).
//!
//! Bulge (arc) segments tessellate at the same pinned
//! [`crate::import::ARC_CHORD_TOL_M`] chord tolerance as the Gerber
//! importer's circular arcs, reusing
//! [`crate::import::arc_vertices`]'s angle-stepping loop directly — only
//! the bulge → `(center, ccw)` conversion is DXF-specific (see
//! [`bulge_vertices`]).

use yee_layout::{Point2, Polygon};

use crate::import;

/// Options controlling the DXF importer.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DxfOptions {
    /// Restrict import to entities on this layer name (DXF group code
    /// `8`). `None` (the default) imports closed chains on every layer.
    pub layer: Option<String>,
}

/// Why a DXF file could not be imported.
#[derive(Debug, Clone, PartialEq)]
pub enum DxfImportError {
    /// `$INSUNITS` (header group `9`/`70`) is missing, or present with a
    /// value outside the supported mm (`4`) / inch (`1`) subset. Carries
    /// the raw value found, or the literal text `"missing"` when the
    /// variable was absent entirely.
    UnsupportedUnits(String),
    /// An entity type outside the supported `LWPOLYLINE` /
    /// `POLYLINE`+`VERTEX` subset (`CIRCLE`, `ARC`, `ELLIPSE`, `SPLINE`,
    /// `TEXT`, `MTEXT`, `DIMENSION`, `INSERT`, …), with the offending
    /// entity name.
    UnsupportedEntity(String),
    /// A `LWPOLYLINE`/`POLYLINE` whose closed flag (group `70` bit 0) is
    /// not set — an open path has no fillable outline area.
    OpenPolyline,
    /// Nonzero `Z` on a vertex, or a nonzero `LWPOLYLINE` constant
    /// elevation (group `38`) — the importer is 2-D only.
    NonzeroElevation,
    /// A coordinate, bulge, or flag value failed to parse, or a required
    /// group code was missing, with the offending text.
    BadValue(String),
    /// A `POLYLINE` entity's `VERTEX` chain never reached a `SEQEND`.
    UnclosedPolyline,
    /// A degenerate bulge segment (zero-length chord), with the
    /// offending segment description.
    BadBulge(String),
    /// The file parsed but contained no closed polylines — an outline
    /// import needs at least one polygon.
    NoOutline,
}

impl std::fmt::Display for DxfImportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedUnits(v) => write!(f, "unsupported or missing $INSUNITS: {v:?}"),
            Self::UnsupportedEntity(k) => write!(f, "unsupported DXF entity: {k:?}"),
            Self::OpenPolyline => write!(f, "polyline is not closed (group 70 bit 0 unset)"),
            Self::NonzeroElevation => write!(f, "nonzero Z / elevation — importer is 2-D only"),
            Self::BadValue(v) => write!(f, "malformed DXF value: {v:?}"),
            Self::UnclosedPolyline => write!(f, "POLYLINE vertex chain never reached SEQEND"),
            Self::BadBulge(v) => write!(f, "malformed bulge segment: {v:?}"),
            Self::NoOutline => write!(f, "no closed polylines found in the DXF"),
        }
    }
}

impl std::error::Error for DxfImportError {}

/// One `0`-delimited DXF entity: its type name plus every subsequent
/// `(group_code, value)` pair up to (not including) the next `0` group.
struct Entity<'a> {
    kind: &'a str,
    attrs: Vec<(i32, &'a str)>,
}

impl<'a> Entity<'a> {
    /// The value of the first occurrence of `code`, if any.
    fn get(&self, code: i32) -> Option<&'a str> {
        self.attrs.iter().find(|(c, _)| *c == code).map(|(_, v)| *v)
    }

    fn get_f64(&self, code: i32) -> Option<f64> {
        self.get(code)?.trim().parse().ok()
    }

    fn get_i32(&self, code: i32) -> Option<i32> {
        self.get(code)?.trim().parse().ok()
    }
}

/// Scan an ASCII DXF into `(group_code, value)` pairs: lines alternate a
/// group-code integer then its value, both trimmed. A group-code line
/// that fails to parse as an integer, or a dangling final line with no
/// value, is skipped rather than erroring — the scanner is deliberately
/// permissive about the boilerplate sections (`HEADER`/`TABLES`/
/// `BLOCKS`/`OBJECTS`) it never interprets; only entities actually used
/// inside `ENTITIES` are validated strictly.
fn group_pairs(dxf: &str) -> Vec<(i32, &str)> {
    let mut lines = dxf.lines().map(str::trim);
    let mut out = Vec::new();
    while let Some(code_line) = lines.next() {
        let Some(value) = lines.next() else { break };
        if let Ok(code) = code_line.parse::<i32>() {
            out.push((code, value));
        }
    }
    out
}

/// Find the `$INSUNITS` header variable's value (group `9` name
/// immediately followed by its value group, per the DXF header
/// convention), returning `None` if absent or unparseable.
fn header_insunits(dxf: &str) -> Option<i32> {
    let mut pending: Option<&str> = None;
    for (code, val) in group_pairs(dxf) {
        if code == 9 {
            pending = Some(val);
            continue;
        }
        if pending == Some("$INSUNITS") {
            return val.trim().parse().ok();
        }
        pending = None;
    }
    None
}

/// Map an `$INSUNITS` value to a metres-per-drawing-unit scale factor.
/// Only `4` (millimetres) and `1` (inches) are supported; anything else,
/// including a missing variable (`None`), is a named rejection — see the
/// module docs for why a missing `$INSUNITS` is not defaulted.
fn unit_scale_m(insunits: Option<i32>) -> Result<f64, DxfImportError> {
    match insunits {
        Some(4) => Ok(1.0e-3),
        Some(1) => Ok(0.0254),
        Some(v) => Err(DxfImportError::UnsupportedUnits(v.to_string())),
        None => Err(DxfImportError::UnsupportedUnits("missing".into())),
    }
}

/// Group the file's flat `(group_code, value)` stream into entities,
/// keeping only those inside the `SECTION 2 ENTITIES` block (so
/// `BLOCKS`-section entity definitions — never a supported subset item
/// per the `INSERT` rejection — are not mistaken for model-space
/// geometry).
fn scan_entities(dxf: &str) -> Vec<Entity<'_>> {
    let mut out = Vec::new();
    let mut section: Option<&str> = None;
    let mut expect_section_name = false;
    let mut cur: Option<Entity> = None;

    for (code, val) in group_pairs(dxf) {
        if code == 0 {
            if let Some(e) = cur.take()
                && section == Some("ENTITIES")
            {
                out.push(e);
            }
            match val {
                "SECTION" => expect_section_name = true,
                "ENDSEC" => section = None,
                "EOF" => break,
                _ => {
                    cur = Some(Entity {
                        kind: val,
                        attrs: Vec::new(),
                    });
                }
            }
            continue;
        }
        if expect_section_name && code == 2 {
            section = Some(val);
            expect_section_name = false;
            continue;
        }
        if let Some(e) = cur.as_mut() {
            e.attrs.push((code, val));
        }
    }
    if let Some(e) = cur.take()
        && section == Some("ENTITIES")
    {
        out.push(e);
    }
    out
}

/// Parse a `LWPOLYLINE`'s vertex list from its ordered attrs: group `10`
/// (X) starts a new vertex, `20` (Y) fills it, `42` (bulge) fills its
/// outgoing-segment bulge (default `0.0`), and a nonzero `30` (Z) or
/// `38` (constant elevation) is rejected. Order-dependent because `10`
/// repeats once per vertex — a `(code, value)` map would collapse them.
fn parse_lwpolyline_verts(attrs: &[(i32, &str)]) -> Result<Vec<(f64, f64, f64)>, DxfImportError> {
    let mut verts = Vec::new();
    let mut cur: Option<(f64, f64, f64)> = None;
    for &(code, val) in attrs {
        let f = || {
            val.trim()
                .parse::<f64>()
                .map_err(|_| DxfImportError::BadValue(val.into()))
        };
        match code {
            10 => {
                if let Some(v) = cur.take() {
                    verts.push(v);
                }
                cur = Some((f()?, 0.0, 0.0));
            }
            20 => {
                if let Some(v) = cur.as_mut() {
                    v.1 = f()?;
                }
            }
            42 => {
                if let Some(v) = cur.as_mut() {
                    v.2 = f()?;
                }
            }
            30 => {
                if f()? != 0.0 {
                    return Err(DxfImportError::NonzeroElevation);
                }
            }
            38 => {
                if f()? != 0.0 {
                    return Err(DxfImportError::NonzeroElevation);
                }
            }
            _ => {}
        }
    }
    if let Some(v) = cur.take() {
        verts.push(v);
    }
    Ok(verts)
}

/// Tessellate a DXF bulge segment from `start` to `end` (metres) into
/// interior arc vertices (the exact endpoint is dropped — the caller's
/// next chain vertex supplies it), reusing
/// [`crate::import::arc_vertices`]'s angle-stepping loop (same pinned
/// chord tolerance) — only the bulge → `(center, ccw)` conversion below
/// is DXF-specific.
///
/// DXF bulge (group `42`) is signed `tan(included_angle / 4)`: positive
/// traces the arc counter-clockwise from `start` to `end`, negative
/// clockwise (`bulge == 0` is a straight segment and must never reach
/// this function). The signed sagitta is `s = bulge · (chord_len / 2)`
/// measured along the chord's right-hand normal `n = (u.y, −u.x)` (`u`
/// the unit chord direction); the center then sits at `M + e·n` where
/// `e = (s² − h²) / (2s)` and `h` is the half-chord length (solved from
/// `r² = h² + e²` together with `|s − e| = r`). Verified against an
/// unambiguous, independently-checkable case: a CCW quarter-turn of the
/// unit circle from `(1, 0)` to `(0, 1)` with `bulge = tan(π/8)`
/// reproduces `center = (0, 0)` exactly.
fn bulge_vertices(
    start: Point2,
    end: Point2,
    bulge: f64,
    line: &str,
) -> Result<Vec<Point2>, DxfImportError> {
    let (dx, dy) = (end.x - start.x, end.y - start.y);
    let d = dx.hypot(dy);
    if d <= import::ARC_CHORD_TOL_M {
        return Err(DxfImportError::BadBulge(line.into()));
    }
    let (ux, uy) = (dx / d, dy / d);
    let (nx, ny) = (uy, -ux);
    let h = d / 2.0;
    let s = bulge * h;
    let e = (s * s - h * h) / (2.0 * s);
    let (mx, my) = ((start.x + end.x) / 2.0, (start.y + end.y) / 2.0);
    let center = (mx + e * nx, my + e * ny);
    let ccw = bulge > 0.0;
    import::arc_vertices((start.x, start.y), (end.x, end.y), center, ccw, line)
        .map_err(|err| DxfImportError::BadBulge(err.to_string()))
}

/// Build a closed [`Polygon`] from a raw `(x, y, bulge)` vertex chain
/// (drawing units, pre-scale) plus the `$INSUNITS` scale factor: each
/// vertex's outgoing segment is either straight (`bulge == 0.0`) or an
/// arc tessellated by [`bulge_vertices`], with the last segment wrapping
/// back to vertex 0 (the [`Polygon`] convention drops the explicit
/// closing vertex).
fn tessellate_chain(raw: &[(f64, f64, f64)], scale_m: f64) -> Result<Polygon, DxfImportError> {
    let pts: Vec<Point2> = raw
        .iter()
        .map(|&(x, y, _)| Point2::new(x * scale_m, y * scale_m))
        .collect();
    let n = pts.len();
    let mut verts = Vec::with_capacity(n);
    for i in 0..n {
        verts.push(pts[i]);
        let bulge = raw[i].2;
        if bulge != 0.0 {
            let end = pts[(i + 1) % n];
            let line = format!(
                "segment {i}: ({:.6},{:.6}) -> ({:.6},{:.6}) bulge {bulge}",
                pts[i].x, pts[i].y, end.x, end.y
            );
            let mut arc = bulge_vertices(pts[i], end, bulge, &line)?;
            // Drop the exact endpoint — either the next loop iteration's
            // `pts[i+1]` push, or (on the wrap-around segment) the
            // already-pushed `verts[0]`, supplies it identically.
            arc.pop();
            verts.extend(arc);
        }
    }
    Ok(Polygon { verts })
}

/// Parse the supported subset of an ASCII DXF into closed outline
/// polygons (metres): `LWPOLYLINE` and `POLYLINE`/`VERTEX` closed
/// chains, in file order, with bulge segments tessellated per the
/// module docs. `opts.layer` optionally restricts import to entities on
/// one layer (DXF group `8`); other layers are skipped, not rejected.
pub fn dxf_to_outline(dxf: &str, opts: &DxfOptions) -> Result<Vec<Polygon>, DxfImportError> {
    let scale_m = unit_scale_m(header_insunits(dxf))?;
    let entities = scan_entities(dxf);
    let mut polys = Vec::new();
    let mut it = entities.into_iter().peekable();

    while let Some(ent) = it.next() {
        match ent.kind {
            "LWPOLYLINE" => {
                if !layer_matches(&ent, opts) {
                    continue;
                }
                if ent.get_f64(38).unwrap_or(0.0) != 0.0 {
                    return Err(DxfImportError::NonzeroElevation);
                }
                let closed = ent.get_i32(70).unwrap_or(0) & 1 != 0;
                if !closed {
                    return Err(DxfImportError::OpenPolyline);
                }
                let raw = parse_lwpolyline_verts(&ent.attrs)?;
                if raw.len() >= 3 {
                    polys.push(tessellate_chain(&raw, scale_m)?);
                }
            }
            "POLYLINE" => {
                let closed = ent.get_i32(70).unwrap_or(0) & 1 != 0;
                let keep = layer_matches(&ent, opts);
                let mut raw = Vec::new();
                while matches!(it.peek(), Some(e) if e.kind == "VERTEX") {
                    let v = it.next().expect("peeked Some");
                    if keep {
                        if v.get_f64(30).unwrap_or(0.0) != 0.0 {
                            return Err(DxfImportError::NonzeroElevation);
                        }
                        let x = v
                            .get_f64(10)
                            .ok_or_else(|| DxfImportError::BadValue("VERTEX missing 10".into()))?;
                        let y = v
                            .get_f64(20)
                            .ok_or_else(|| DxfImportError::BadValue("VERTEX missing 20".into()))?;
                        raw.push((x, y, v.get_f64(42).unwrap_or(0.0)));
                    }
                }
                match it.next() {
                    Some(e) if e.kind == "SEQEND" => {}
                    _ => return Err(DxfImportError::UnclosedPolyline),
                }
                if !keep {
                    continue;
                }
                if !closed {
                    return Err(DxfImportError::OpenPolyline);
                }
                if raw.len() >= 3 {
                    polys.push(tessellate_chain(&raw, scale_m)?);
                }
            }
            "CIRCLE" | "ARC" | "ELLIPSE" | "SPLINE" | "TEXT" | "MTEXT" | "DIMENSION" | "INSERT"
            | "LINE" | "3DFACE" | "POINT" | "SOLID" | "HATCH" => {
                return Err(DxfImportError::UnsupportedEntity(ent.kind.into()));
            }
            other => return Err(DxfImportError::UnsupportedEntity(other.into())),
        }
    }
    if polys.is_empty() {
        return Err(DxfImportError::NoOutline);
    }
    Ok(polys)
}

/// Whether an entity's layer (group `8`) passes `opts.layer` (no filter
/// configured always passes).
fn layer_matches(ent: &Entity<'_>, opts: &DxfOptions) -> bool {
    match &opts.layer {
        None => true,
        Some(want) => ent.get(8) == Some(want.as_str()),
    }
}
