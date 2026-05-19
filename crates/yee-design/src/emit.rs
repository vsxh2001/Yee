//! Deterministic project-TOML emitter for the Phase 3.nl.0 design surface.
//!
//! Stage 5 of the spec §5 pipeline: consume a Stage-3 [`InitialEstimate`] plus
//! the [`DesignIntent`] that produced it and render two artefacts:
//!
//! - `ProjectFile::toml` — a `yee.toml` project file the existing `yee run`
//!   pipeline accepts (or will accept; see the cross-lane note below). The
//!   first two lines are metadata comments per spec §8:
//!     - `# nl-prompt: <verbatim source prompt>`
//!     - `# yee-design: <sha256(intent_json)> <provenance text>`
//! - `ProjectFile::intent_json` — the byte-for-byte `serde_json::to_string` of
//!   the [`DesignIntent`], used by the `<out>.intent.json` sidecar and by the
//!   spec §8 hash on line two.
//!
//! ## Determinism contract (spec §8)
//!
//! Two `emit(estimate, intent)` calls with `==` arguments produce
//! byte-identical `toml` and `intent_json` strings. This is enforced by:
//!
//! 1. The TOML body is hand-formatted; table order is fixed in code and
//!    keys *within* each table are emitted lexicographically via a
//!    [`BTreeMap`]. The `toml_edit` crate is consulted only as a parser in
//!    tests (round-trip), never as a writer — its float printer has
//!    historically tweaked across minor releases and its `Repr` setter is
//!    `pub(crate)` in 0.22.x.
//! 2. A single canonical float format — `format!("{:.6e}", x)` — is the
//!    one and only way an `f64` enters the emitted string. Integers
//!    (sweep-point counts, port indices) render as TOML integers.
//! 3. The two `#` comment lines sit at the very top, then exactly one
//!    blank line, then the body. No trailing whitespace ever.
//! 4. `serde_json::to_string` (compact, not `to_string_pretty`) for
//!    `intent_json`. `serde_json` preserves the `Serialize` impl's field
//!    order, which is the struct field declaration order, which is stable
//!    across compiler runs.
//!
//! ## Forward-compatible schema (cross-lane finding)
//!
//! At this base SHA, `yee-cli`'s `Run` subcommand is the Phase 0 stub
//! (`crates/yee-cli/src/main.rs:267` — `"yee run {} — Phase 0 stub."`), so
//! there is no canonical project-file schema to mirror. This emitter picks a
//! forward-compatible TOML shape that the R3 brief sketches:
//!
//! ```toml
//! # nl-prompt: ...
//! # yee-design: <hash> <provenance>
//!
//! [frequency]
//! center_hz = 2.400000e9
//! span_hz = 4.000000e8
//! sweep_points = 201
//!
//! [geometry]
//! feed_width_m = 3.060000e-3
//! inset_offset_m = 9.500000e-3
//! length_m = 2.943000e-2
//! type = "rectangular_inset_patch"
//! width_m = 3.804000e-2
//!
//! [substrate]
//! eps_r = 4.400000e0
//! h_m = 1.600000e-3
//! loss_tangent = 2.000000e-2
//!
//! [[ports]]
//! id = 1
//! inset_offset_m = 9.500000e-3
//! kind = "delta_gap"
//! location = "feed"
//! z0_ohm = 5.000000e1
//! ```
//!
//! When `yee run` lands its real parser, any field-name mismatch should be
//! reconciled by widening `yee-cli`'s lane, not this crate's. The shape
//! above is recorded as the forward-compatible target. // TBD: confirm with
//! yee-cli once Run parser lands.

use std::collections::BTreeMap;

use sha2::{Digest, Sha256};

use crate::estimate::InitialEstimate;
use crate::intent::{DesignIntent, GeometryFamily};

/// Emitted artefacts for one design pass.
///
/// Two fields, both `String`s: the TOML body of the project file and the
/// `serde_json` encoding of the [`DesignIntent`] that produced it. Both are
/// deterministic functions of `(estimate, intent)` per spec §8.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectFile {
    /// Project-file TOML body, with the two-line spec §8 metadata header.
    pub toml: String,
    /// `serde_json::to_string(&intent)` — the `<out>.intent.json` sidecar
    /// content. The same bytes appear hashed inside `toml`'s line-two header.
    pub intent_json: String,
}

/// A typed scalar suitable for a TOML key-value row.
///
/// We materialise the table contents as `BTreeMap<String, Scalar>` so keys
/// are sorted lexicographically by construction, and the serialiser walks
/// the map in order to produce the emitted TOML body.
#[derive(Debug, Clone, PartialEq)]
enum Scalar {
    /// A floating-point number, rendered via `format!("{:.6e}", x)`.
    F64(f64),
    /// A signed integer, rendered as a TOML integer literal.
    Int(i64),
    /// A bare string value, rendered as a double-quoted TOML string.
    Str(String),
}

impl Scalar {
    /// Render this scalar as the right-hand side of a TOML key-value row.
    fn render(&self) -> String {
        match self {
            Self::F64(x) => fmt_f64(*x),
            Self::Int(n) => n.to_string(),
            Self::Str(s) => format!("\"{}\"", escape_basic_string(s)),
        }
    }
}

/// Format an `f64` for the project file.
///
/// `format!("{:.6e}", x)` — one canonical scientific-notation form, six
/// fractional digits. This is the spec §8 numeric format; every float in
/// the emitted TOML uses it without exception.
fn fmt_f64(x: f64) -> String {
    format!("{:.6e}", x)
}

/// Escape a string for a TOML basic-string literal.
///
/// TOML basic strings (the double-quoted form) require escaping backslashes,
/// double quotes, and the standard control characters. Our inputs are model-
/// produced provenance strings + library names; the function is implemented
/// in full so the emitter cannot ever produce a malformed TOML.
fn escape_basic_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\x08' => out.push_str("\\b"),
            '\x0c' => out.push_str("\\f"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04X}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// Render a `BTreeMap<String, Scalar>` as a sequence of TOML key-value rows.
///
/// Keys come out in `BTreeMap`'s natural order (lexicographic on the `String`
/// keys). Caller is responsible for the leading `[table]` header line.
fn render_table_body(map: &BTreeMap<String, Scalar>) -> String {
    let mut out = String::new();
    for (k, v) in map {
        out.push_str(k);
        out.push_str(" = ");
        out.push_str(&v.render());
        out.push('\n');
    }
    out
}

/// Render Stage-3 outputs into a deterministic [`ProjectFile`].
///
/// Spec §8 byte-identity invariant: `emit(e, i) == emit(e, i)` for any
/// `(e, i)`. The function is total over the closed [`GeometryFamily`] enum
/// (`RectangularPatch` only at Phase 3.nl.0); future families are pattern-
/// matched explicitly so a new enum arm produces a compile error rather
/// than a silent wrong-shape TOML.
pub fn emit(estimate: &InitialEstimate, intent: &DesignIntent) -> ProjectFile {
    // 1. Serialize the intent first; its bytes feed both the sidecar and
    //    the line-two hash in the TOML header.
    let intent_json = serde_json::to_string(intent)
        .expect("yee-design::emit: DesignIntent is Serialize-infallible by construction");

    // 2. Build the TOML body. Top-level tables come out in a fixed order
    //    (`frequency`, `geometry`, `substrate`, `ports`); keys *within* each
    //    table sort lexicographically via the `BTreeMap`.
    let mut body = String::new();
    push_table(&mut body, "frequency", &build_frequency_table(intent));
    push_table(
        &mut body,
        "geometry",
        &build_geometry_table(estimate, intent),
    );
    push_table(&mut body, "substrate", &build_substrate_table(estimate));
    push_array_of_tables(&mut body, "ports", &build_ports_array(estimate));

    // 3. Prepend the spec §8 metadata header.
    let intent_hash = sha256_hex(&intent_json);
    let provenance = format!("source={}", intent.provenance.source);
    let header = format!(
        "# nl-prompt: {}\n# yee-design: {} {}\n\n",
        sanitise_comment(&intent.source_prompt),
        intent_hash,
        provenance,
    );
    let toml = format!("{header}{body}");

    ProjectFile { toml, intent_json }
}

/// Append `[name]\n<rows>\n` to `out`. One trailing blank line separates this
/// table from the next.
fn push_table(out: &mut String, name: &str, rows: &BTreeMap<String, Scalar>) {
    out.push('[');
    out.push_str(name);
    out.push_str("]\n");
    out.push_str(&render_table_body(rows));
    out.push('\n');
}

/// Append a `[[name]]` array-of-tables to `out`, one element at a time.
fn push_array_of_tables(out: &mut String, name: &str, rows_each: &[BTreeMap<String, Scalar>]) {
    for rows in rows_each {
        out.push_str("[[");
        out.push_str(name);
        out.push_str("]]\n");
        out.push_str(&render_table_body(rows));
        out.push('\n');
    }
}

/// SHA-256 hex digest of a string. Lower-case, no separators, 64 chars.
fn sha256_hex(s: &str) -> String {
    let mut h = Sha256::new();
    h.update(s.as_bytes());
    let digest = h.finalize();
    let mut out = String::with_capacity(64);
    for b in digest {
        out.push(nibble_to_hex(b >> 4));
        out.push(nibble_to_hex(b & 0x0f));
    }
    out
}

/// One-nibble (0..=15) → ASCII hex char ('0'..='9' | 'a'..='f').
fn nibble_to_hex(n: u8) -> char {
    match n {
        0..=9 => (b'0' + n) as char,
        10..=15 => (b'a' + (n - 10)) as char,
        _ => unreachable!("nibble out of range"),
    }
}

/// Strip newlines from a string so it can sit on a single TOML comment line.
///
/// The verbatim prompt is otherwise echoed unchanged; the prompt is preserved
/// byte-for-byte in `intent_json`, so this lossy single-line form in the
/// header is for human reading only.
fn sanitise_comment(s: &str) -> String {
    s.replace(['\r', '\n'], " ")
}

/// `[frequency]` table contents.
///
/// `center_hz` is the resonant frequency from the [`DesignIntent`]. `span_hz`
/// is a derived sweep span — 20 % of `center_hz` by default, or
/// `2 · bandwidth_target_mhz` if the intent expressed one (so the engineer's
/// stated bandwidth lands inside the swept band, not at its edge).
fn build_frequency_table(intent: &DesignIntent) -> BTreeMap<String, Scalar> {
    let f0 = intent.target_frequency_hz;
    let span = match intent.bandwidth_target_mhz {
        Some(bw_mhz) if bw_mhz > 0.0 => 2.0 * bw_mhz * 1.0e6,
        _ => 0.20 * f0,
    };
    let mut map = BTreeMap::new();
    map.insert("center_hz".to_string(), Scalar::F64(f0));
    map.insert("span_hz".to_string(), Scalar::F64(span));
    map.insert("sweep_points".to_string(), Scalar::Int(201));
    map
}

/// `[geometry]` table contents.
///
/// Field set is the [`InitialEstimate`] flat-struct projection plus the
/// geometry-family tag. Every dimension is in metres (SI; the spec §7
/// schema is metres-internal, millimetres only at the user-facing
/// substrate-library boundary).
fn build_geometry_table(
    estimate: &InitialEstimate,
    intent: &DesignIntent,
) -> BTreeMap<String, Scalar> {
    let kind = match intent.family {
        GeometryFamily::RectangularPatch => "rectangular_inset_patch",
    };
    let mut map = BTreeMap::new();
    map.insert("type".to_string(), Scalar::Str(kind.to_string()));
    map.insert("width_m".to_string(), Scalar::F64(estimate.width_m));
    map.insert("length_m".to_string(), Scalar::F64(estimate.length_m));
    map.insert(
        "inset_offset_m".to_string(),
        Scalar::F64(estimate.inset_offset_m),
    );
    map.insert(
        "feed_width_m".to_string(),
        Scalar::F64(estimate.feed_width_m),
    );
    map
}

/// `[substrate]` table contents.
///
/// Materialises the [`crate::ResolvedSubstrate`] (post library-lookup +
/// override merge) rather than the [`crate::Substrate`] enum, so the project
/// file is a flat set of named scalars. `roughness_m` is omitted at Phase
/// 3.nl.0 (R3 brief lists it as optional); it would land here if a later
/// phase surfaces it on the substrate library row.
fn build_substrate_table(estimate: &InitialEstimate) -> BTreeMap<String, Scalar> {
    let mut map = BTreeMap::new();
    map.insert("eps_r".to_string(), Scalar::F64(estimate.substrate.eps_r));
    map.insert("h_m".to_string(), Scalar::F64(estimate.substrate.h_m));
    map.insert(
        "loss_tangent".to_string(),
        Scalar::F64(estimate.substrate.loss_tangent),
    );
    map
}

/// `[[ports]]` array-of-tables contents.
///
/// Phase 3.nl.0 emits a single 50 Ω feed port. The shape is forward-
/// compatible with a multi-port geometry: each port is a TOML inline-style
/// table. Keys inside each port table are sorted via the same `BTreeMap`
/// mechanism as the top-level tables.
fn build_ports_array(estimate: &InitialEstimate) -> Vec<BTreeMap<String, Scalar>> {
    let mut port = BTreeMap::new();
    port.insert("id".to_string(), Scalar::Int(1));
    port.insert("kind".to_string(), Scalar::Str("delta_gap".to_string()));
    port.insert("z0_ohm".to_string(), Scalar::F64(50.0));
    port.insert("location".to_string(), Scalar::Str("feed".to_string()));
    // Symmetry with the geometry table: the port records the inset offset
    // it sits at, so a downstream consumer that ignores the geometry block
    // can still place the excitation correctly.
    port.insert(
        "inset_offset_m".to_string(),
        Scalar::F64(estimate.inset_offset_m),
    );
    vec![port]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intent::{NamedSubstrate, Provenance, Substrate, substrate_library};

    fn fr4_intent() -> DesignIntent {
        DesignIntent {
            family: GeometryFamily::RectangularPatch,
            target_frequency_hz: 2.4e9,
            substrate: Substrate::Named(NamedSubstrate {
                name: "FR4".to_string(),
                override_with: None,
            }),
            gain_target_dbi: None,
            bandwidth_target_mhz: Some(100.0),
            source_prompt: "2.4 GHz patch on FR4".to_string(),
            provenance: Provenance {
                source: "offline".to_string(),
                model: None,
                temperature: None,
                schema_version: "1".to_string(),
                substrate_library_version: substrate_library().version.clone(),
            },
        }
    }

    #[test]
    fn emit_header_has_two_comment_lines() {
        let intent = fr4_intent();
        let est = InitialEstimate::from_intent(&intent).expect("estimate");
        let out = emit(&est, &intent);
        let mut lines = out.toml.lines();
        let l0 = lines.next().expect("line 0");
        let l1 = lines.next().expect("line 1");
        assert!(l0.starts_with("# nl-prompt:"), "line 0: {l0}");
        assert!(l1.starts_with("# yee-design:"), "line 1: {l1}");
        // Hash on line 1 is 64 lowercase hex chars.
        let hash_field = l1
            .strip_prefix("# yee-design: ")
            .and_then(|rest| rest.split_whitespace().next())
            .expect("hash on line 1");
        assert_eq!(hash_field.len(), 64);
        assert!(hash_field.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn emit_floats_use_canonical_scientific_format() {
        let intent = fr4_intent();
        let est = InitialEstimate::from_intent(&intent).expect("estimate");
        let out = emit(&est, &intent);
        // Every f64 row contains `e` (scientific notation).
        // Spot-check one literal: eps_r for FR4 is 4.4, so the row should
        // be `eps_r = 4.400000e0`.
        assert!(
            out.toml.contains("eps_r = 4.400000e0"),
            "expected canonical eps_r literal; got:\n{}",
            out.toml
        );
    }

    #[test]
    fn emit_body_parses_as_toml() {
        let intent = fr4_intent();
        let est = InitialEstimate::from_intent(&intent).expect("estimate");
        let out = emit(&est, &intent);
        // Use `toml` (value-only) to confirm the emitted body is well-formed.
        let parsed: toml::Value = toml::from_str(&out.toml).expect("body must parse as TOML");
        assert!(parsed.get("geometry").is_some());
        assert!(parsed.get("substrate").is_some());
        assert!(parsed.get("frequency").is_some());
        assert!(parsed.get("ports").is_some());
    }
}
