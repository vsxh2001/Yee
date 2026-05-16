//! Minimal KiCad `.kicad_pcb` importer.
//!
//! Phase 1.mesh.1 walking-skeleton scope (s-expression subset):
//! - Top-level `(kicad_pcb ...)` wrapper
//! - `(general (thickness ...))` → total board thickness
//! - `(layers (N "Name" type) ...)` → ordered layer list
//! - `(segment (start ...) (end ...) (width ...) (layer ...))` → traces
//! - `(zone ... (polygon (pts (xy ...) ...)))` → copper fills
//!
//! Everything else (footprints, vias, drills, silk, mask) is skipped. Full
//! parsing — including footprints with embedded copper, vias, arcs, and
//! masks — is a Phase 1.mesh.2 effort.
//!
//! No external s-expression crate is used; we hand-roll a tiny tokenizer +
//! recursive-descent walker over the resulting tree to keep `yee-mesh`'s
//! dependency footprint minimal.

use std::path::Path;

/// In-memory representation of a parsed `.kicad_pcb` file.
#[derive(Debug, Clone)]
pub struct KiCadBoard {
    /// Total board thickness in millimetres (from `(general (thickness ...))`).
    pub thickness_mm: f64,
    /// Ordered list of layers declared in the `(layers ...)` block.
    pub layers: Vec<LayerInfo>,
    /// Copper trace segments.
    pub segments: Vec<Segment>,
    /// Copper-fill zones with their polygon outline.
    pub zones: Vec<Zone>,
}

/// One entry from the top-level `(layers ...)` block.
#[derive(Debug, Clone)]
pub struct LayerInfo {
    /// Numeric ordinal (e.g. `0` for `F.Cu`, `31` for `B.Cu`).
    pub ordinal: u32,
    /// Human-readable name (e.g. `"F.Cu"`).
    pub name: String,
    /// Layer kind (e.g. `"signal"`, `"power"`, `"mixed"`, `"user"`).
    pub kind: String,
}

/// A copper trace segment.
#[derive(Debug, Clone)]
pub struct Segment {
    /// Start point in millimetres (KiCad-native units).
    pub start: (f64, f64),
    /// End point in millimetres.
    pub end: (f64, f64),
    /// Trace width in millimetres.
    pub width_mm: f64,
    /// Layer name the segment lives on.
    pub layer: String,
}

/// A copper-fill zone. We currently capture only the outline polygon; the
/// hatched fill pattern is reconstructed downstream.
#[derive(Debug, Clone)]
pub struct Zone {
    /// Layer the zone lives on.
    pub layer: String,
    /// Outline polygon as `(x, y)` vertices in millimetres.
    pub polygon: Vec<(f64, f64)>,
}

/// Errors produced while reading or parsing a `.kicad_pcb` file.
#[derive(Debug, thiserror::Error)]
pub enum KiCadError {
    /// Filesystem I/O failed.
    #[error("io error: {0}")]
    Io(String),
    /// Malformed s-expression or unexpected structure.
    #[error("parse error at byte {byte}: {msg}")]
    Parse {
        /// Byte offset within the input where the failure was detected.
        byte: usize,
        /// Human-readable explanation.
        msg: String,
    },
    /// The file parsed as s-expressions but does not look like a KiCad PCB.
    #[error("unsupported file: {0}")]
    Unsupported(String),
}

// -----------------------------------------------------------------------------
// Tokenizer
// -----------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Token {
    LParen,
    RParen,
    Atom(String),
}

struct Tokenizer<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Tokenizer<'a> {
    fn new(src: &'a str) -> Self {
        Self {
            src: src.as_bytes(),
            pos: 0,
        }
    }

    fn skip_ws(&mut self) {
        while self.pos < self.src.len() {
            let c = self.src[self.pos];
            if c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn next_token(&mut self) -> Result<Option<(usize, Token)>, KiCadError> {
        self.skip_ws();
        if self.pos >= self.src.len() {
            return Ok(None);
        }
        let start = self.pos;
        let c = self.src[self.pos];
        match c {
            b'(' => {
                self.pos += 1;
                Ok(Some((start, Token::LParen)))
            }
            b')' => {
                self.pos += 1;
                Ok(Some((start, Token::RParen)))
            }
            b'"' => {
                // Quoted string atom. Support `\"` and `\\` escapes; everything
                // else passes through verbatim.
                self.pos += 1;
                let mut buf = String::new();
                while self.pos < self.src.len() {
                    let ch = self.src[self.pos];
                    if ch == b'\\' {
                        if self.pos + 1 >= self.src.len() {
                            return Err(KiCadError::Parse {
                                byte: self.pos,
                                msg: "trailing escape in quoted atom".into(),
                            });
                        }
                        let esc = self.src[self.pos + 1];
                        buf.push(esc as char);
                        self.pos += 2;
                    } else if ch == b'"' {
                        self.pos += 1;
                        return Ok(Some((start, Token::Atom(buf))));
                    } else {
                        buf.push(ch as char);
                        self.pos += 1;
                    }
                }
                Err(KiCadError::Parse {
                    byte: start,
                    msg: "unterminated quoted atom".into(),
                })
            }
            _ => {
                // Bare atom: read until whitespace or paren.
                let begin = self.pos;
                while self.pos < self.src.len() {
                    let ch = self.src[self.pos];
                    if ch == b' '
                        || ch == b'\t'
                        || ch == b'\n'
                        || ch == b'\r'
                        || ch == b'('
                        || ch == b')'
                    {
                        break;
                    }
                    self.pos += 1;
                }
                let s = std::str::from_utf8(&self.src[begin..self.pos]).map_err(|_| {
                    KiCadError::Parse {
                        byte: begin,
                        msg: "non-utf8 atom".into(),
                    }
                })?;
                Ok(Some((start, Token::Atom(s.to_string()))))
            }
        }
    }
}

// -----------------------------------------------------------------------------
// S-expression tree
// -----------------------------------------------------------------------------

#[derive(Debug, Clone)]
enum Sexp {
    Atom(String),
    List(Vec<Sexp>),
}

impl Sexp {
    fn as_atom(&self) -> Option<&str> {
        match self {
            Sexp::Atom(s) => Some(s.as_str()),
            _ => None,
        }
    }

    fn as_list(&self) -> Option<&[Sexp]> {
        match self {
            Sexp::List(v) => Some(v.as_slice()),
            _ => None,
        }
    }

    /// Treat this node as a list whose first element is the head atom; return
    /// the head and the tail.
    fn head_tail(&self) -> Option<(&str, &[Sexp])> {
        let items = self.as_list()?;
        let head = items.first()?.as_atom()?;
        Some((head, &items[1..]))
    }
}

fn parse_tree(text: &str) -> Result<Sexp, KiCadError> {
    let mut tok = Tokenizer::new(text);
    let first = tok.next_token()?;
    let (start_byte, first_tok) = match first {
        Some(t) => t,
        None => {
            return Err(KiCadError::Parse {
                byte: 0,
                msg: "empty input".into(),
            });
        }
    };
    if first_tok != Token::LParen {
        return Err(KiCadError::Parse {
            byte: start_byte,
            msg: "expected '(' at start of file".into(),
        });
    }
    let root = parse_list(&mut tok, start_byte)?;
    // No trailing tokens (other than whitespace) are allowed.
    if let Some((b, _)) = tok.next_token()? {
        return Err(KiCadError::Parse {
            byte: b,
            msg: "unexpected trailing tokens after root list".into(),
        });
    }
    Ok(root)
}

/// Called after a `LParen` has been consumed; reads through the matching `RParen`.
fn parse_list(tok: &mut Tokenizer<'_>, open_byte: usize) -> Result<Sexp, KiCadError> {
    let mut items = Vec::new();
    loop {
        let next = tok.next_token()?;
        match next {
            None => {
                return Err(KiCadError::Parse {
                    byte: open_byte,
                    msg: "unterminated list".into(),
                });
            }
            Some((byte, Token::LParen)) => {
                items.push(parse_list(tok, byte)?);
            }
            Some((_, Token::RParen)) => return Ok(Sexp::List(items)),
            Some((_, Token::Atom(s))) => items.push(Sexp::Atom(s)),
        }
    }
}

// -----------------------------------------------------------------------------
// Domain extraction
// -----------------------------------------------------------------------------

fn atom_to_f64(s: &Sexp) -> Result<f64, KiCadError> {
    match s {
        Sexp::Atom(a) => a.parse::<f64>().map_err(|_| KiCadError::Parse {
            byte: 0,
            msg: format!("expected number, got `{a}`"),
        }),
        Sexp::List(_) => Err(KiCadError::Parse {
            byte: 0,
            msg: "expected number, got list".into(),
        }),
    }
}

fn atom_to_u32(s: &Sexp) -> Result<u32, KiCadError> {
    match s {
        Sexp::Atom(a) => a.parse::<u32>().map_err(|_| KiCadError::Parse {
            byte: 0,
            msg: format!("expected unsigned int, got `{a}`"),
        }),
        Sexp::List(_) => Err(KiCadError::Parse {
            byte: 0,
            msg: "expected unsigned int, got list".into(),
        }),
    }
}

fn atom_to_string(s: &Sexp) -> Result<String, KiCadError> {
    match s {
        Sexp::Atom(a) => Ok(a.clone()),
        Sexp::List(_) => Err(KiCadError::Parse {
            byte: 0,
            msg: "expected atom, got list".into(),
        }),
    }
}

/// Find the first child list whose head atom equals `key`. Returns `None` if
/// no such child exists.
fn find_child<'a>(items: &'a [Sexp], key: &str) -> Option<&'a [Sexp]> {
    for item in items {
        if let Sexp::List(inner) = item {
            if let Some(Sexp::Atom(head)) = inner.first() {
                if head == key {
                    return Some(&inner[1..]);
                }
            }
        }
    }
    None
}

fn extract_general_thickness(general_tail: &[Sexp]) -> Result<f64, KiCadError> {
    let thickness = find_child(general_tail, "thickness").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "missing (thickness ...) in (general ...)".into(),
    })?;
    let v = thickness.first().ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "empty (thickness ...) form".into(),
    })?;
    atom_to_f64(v)
}

fn extract_layers(layers_tail: &[Sexp]) -> Result<Vec<LayerInfo>, KiCadError> {
    let mut out = Vec::new();
    for item in layers_tail {
        let inner = match item.as_list() {
            Some(v) => v,
            None => continue,
        };
        // Each layer is `(N "Name" kind [maybe more])` — head is the ordinal,
        // not an atom keyword.
        if inner.len() < 3 {
            continue;
        }
        let ordinal = atom_to_u32(&inner[0])?;
        let name = atom_to_string(&inner[1])?;
        let kind = atom_to_string(&inner[2])?;
        out.push(LayerInfo {
            ordinal,
            name,
            kind,
        });
    }
    Ok(out)
}

fn extract_xy_pair(form: &[Sexp]) -> Result<(f64, f64), KiCadError> {
    if form.len() < 2 {
        return Err(KiCadError::Parse {
            byte: 0,
            msg: "expected at least 2 coordinates".into(),
        });
    }
    let x = atom_to_f64(&form[0])?;
    let y = atom_to_f64(&form[1])?;
    Ok((x, y))
}

fn extract_segment(tail: &[Sexp]) -> Result<Segment, KiCadError> {
    let start = find_child(tail, "start").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "segment missing (start ...)".into(),
    })?;
    let end = find_child(tail, "end").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "segment missing (end ...)".into(),
    })?;
    let width = find_child(tail, "width").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "segment missing (width ...)".into(),
    })?;
    let layer = find_child(tail, "layer").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "segment missing (layer ...)".into(),
    })?;

    let start_xy = extract_xy_pair(start)?;
    let end_xy = extract_xy_pair(end)?;
    let width_mm = atom_to_f64(width.first().ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "empty (width ...) form".into(),
    })?)?;
    let layer_name = atom_to_string(layer.first().ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "empty (layer ...) form".into(),
    })?)?;

    Ok(Segment {
        start: start_xy,
        end: end_xy,
        width_mm,
        layer: layer_name,
    })
}

/// Recursively search for the first `(polygon (pts (xy x y) ...))` form and
/// return its vertex list.
fn extract_zone_polygon(node: &Sexp) -> Option<Vec<(f64, f64)>> {
    let items = node.as_list()?;
    if let Some(Sexp::Atom(head)) = items.first() {
        if head == "polygon" {
            // Look for the (pts ...) child.
            let pts = find_child(&items[1..], "pts")?;
            let mut verts = Vec::new();
            for xy in pts {
                if let Sexp::List(inner) = xy {
                    if let Some(Sexp::Atom(h)) = inner.first() {
                        if h == "xy" && inner.len() >= 3 {
                            let x = inner[1].as_atom()?.parse::<f64>().ok()?;
                            let y = inner[2].as_atom()?.parse::<f64>().ok()?;
                            verts.push((x, y));
                        }
                    }
                }
            }
            return Some(verts);
        }
    }
    for child in &items[1..] {
        if let Some(v) = extract_zone_polygon(child) {
            return Some(v);
        }
    }
    None
}

fn extract_zone(tail: &[Sexp]) -> Result<Option<Zone>, KiCadError> {
    let layer = find_child(tail, "layer").ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "zone missing (layer ...)".into(),
    })?;
    let layer_name = atom_to_string(layer.first().ok_or_else(|| KiCadError::Parse {
        byte: 0,
        msg: "empty (layer ...) form".into(),
    })?)?;

    // Walk every child of the zone looking for the first (polygon ...).
    for item in tail {
        if let Some(verts) = extract_zone_polygon(item) {
            return Ok(Some(Zone {
                layer: layer_name,
                polygon: verts,
            }));
        }
    }
    // Zone without a polygon: silently skip — KiCad emits these for empty
    // unfilled zones.
    Ok(None)
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

impl KiCadBoard {
    /// Read and parse a `.kicad_pcb` file from disk.
    pub fn read(path: &Path) -> Result<Self, KiCadError> {
        let bytes = std::fs::read(path).map_err(|e| KiCadError::Io(e.to_string()))?;
        let text = std::str::from_utf8(&bytes)
            .map_err(|e| KiCadError::Io(format!("non-utf8 file: {e}")))?;
        Self::parse(text)
    }

    /// Parse a `.kicad_pcb` document already loaded into memory.
    pub fn parse(text: &str) -> Result<Self, KiCadError> {
        let tree = parse_tree(text)?;
        let items = tree
            .as_list()
            .ok_or_else(|| KiCadError::Unsupported("root sexp is not a list".into()))?;
        let head = items
            .first()
            .and_then(|s| s.as_atom())
            .ok_or_else(|| KiCadError::Unsupported("root list has no head atom".into()))?;
        if head != "kicad_pcb" {
            return Err(KiCadError::Unsupported(format!(
                "root form is `{head}`, expected `kicad_pcb`"
            )));
        }
        let tail = &items[1..];

        let mut thickness_mm = 0.0;
        let mut layers = Vec::new();
        let mut segments = Vec::new();
        let mut zones = Vec::new();

        for item in tail {
            let (head, body) = match item.head_tail() {
                Some(v) => v,
                None => continue,
            };
            match head {
                "general" => {
                    thickness_mm = extract_general_thickness(body)?;
                }
                "layers" => {
                    layers = extract_layers(body)?;
                }
                "segment" => {
                    segments.push(extract_segment(body)?);
                }
                "zone" => {
                    if let Some(z) = extract_zone(body)? {
                        zones.push(z);
                    }
                }
                _ => {
                    // Skip everything else (footprints, gr_*, vias, net, etc.)
                    // per Phase 1.mesh.1 scope.
                }
            }
        }

        Ok(Self {
            thickness_mm,
            layers,
            segments,
            zones,
        })
    }
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const MINI: &str = r#"(kicad_pcb (version 20221018) (generator pcbnew)
  (general (thickness 1.6))
  (layers
    (0 "F.Cu" signal)
    (31 "B.Cu" signal)
  )
  (segment (start 10 10) (end 50 10) (width 0.25) (layer "F.Cu") (net 1))
  (segment (start 50 10) (end 50 50) (width 0.25) (layer "F.Cu") (net 1))
  (zone (net 0) (net_name "GND") (layer "B.Cu")
    (polygon (pts (xy 0 0) (xy 100 0) (xy 100 100) (xy 0 100)))
  )
)"#;

    #[test]
    fn parse_thickness() {
        let board = KiCadBoard::parse(MINI).unwrap();
        assert!((board.thickness_mm - 1.6).abs() < 1e-9);
    }

    #[test]
    fn parse_two_layers() {
        let board = KiCadBoard::parse(MINI).unwrap();
        assert_eq!(board.layers.len(), 2);
        assert_eq!(board.layers[0].name, "F.Cu");
        assert_eq!(board.layers[1].name, "B.Cu");
    }

    #[test]
    fn parse_two_segments() {
        let board = KiCadBoard::parse(MINI).unwrap();
        assert_eq!(board.segments.len(), 2);
        assert!((board.segments[0].width_mm - 0.25).abs() < 1e-9);
    }

    #[test]
    fn parse_one_zone_with_four_corners() {
        let board = KiCadBoard::parse(MINI).unwrap();
        assert_eq!(board.zones.len(), 1);
        assert_eq!(board.zones[0].layer, "B.Cu");
        assert_eq!(board.zones[0].polygon.len(), 4);
    }

    #[test]
    fn malformed_input_rejected() {
        let bad = "(kicad_pcb (general";
        assert!(matches!(
            KiCadBoard::parse(bad),
            Err(KiCadError::Parse { .. })
        ));
    }
}
