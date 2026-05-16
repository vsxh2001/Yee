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

// The tokenizer + tree builder land before the domain extractor that consumes
// them; the `dead_code` allowance is dropped in the follow-up commit that
// wires `KiCadBoard::parse` to them.
#![allow(dead_code)]

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
