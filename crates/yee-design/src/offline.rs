//! Deterministic offline parser for natural-language design prompts.
//!
//! Stage-1 fallback in the spec §5 pipeline: turns a free-form English prompt
//! such as `"2.4 GHz patch on FR4"` into a typed [`crate::DesignIntent`]
//! without consulting an LLM. The parser is intentionally forgiving — it scans
//! for keywords / regex matches and applies documented defaults for anything
//! it cannot extract — because Phase 3.nl.0's CI default and the offline
//! `--offline` flag both depend on a guaranteed-deterministic Stage-1 that
//! cannot reach the network.
//!
//! ## Grammar (informal)
//!
//! - **Frequency** — first match of `(\d+(\.\d+)?)\s*(GHz|MHz)`. Default
//!   `2.4 GHz` when absent.
//! - **Substrate** — first case-insensitive substring match against the
//!   canonical library names (`FR4`, `RO4003C`, `RO5880`, `AluminaTC`).
//!   Default `FR4` when absent. The match is whole-word-ish: we match against
//!   the lowercased prompt so `"ro4003c"` and `"RO4003C"` both work.
//! - **Geometry family** — Phase 3.nl.0 supports only
//!   [`crate::GeometryFamily::RectangularPatch`]; the parser always emits it.
//!
//! ## Determinism
//!
//! The parser is a pure function: same prompt in, bit-identical
//! [`crate::DesignIntent`] out. The substrate-library version is captured from
//! the process-wide [`crate::substrate_library`] singleton, which is itself
//! deterministic (baked-in `substrates.toml`). Provenance is always
//! `"offline"`; no LLM-specific metadata leaks in.
//!
//! See `crates/yee-design/validation/prompts.toml` for the 10 canonical
//! prompts that this parser must handle; they exercise every frequency-unit /
//! substrate-name combination the Phase 3.nl.0 walking skeleton supports.

use regex::Regex;
use std::sync::OnceLock;

use crate::intent::{
    DesignIntent, GeometryFamily, NamedSubstrate, Provenance, Substrate, substrate_library,
};

/// Schema version of the structured-output contract this parser emits.
///
/// Matches the value the LLM sidecar (R4) writes into its `Provenance`. Bumped
/// only when the `DesignIntent` shape changes incompatibly; substrate-library
/// drift is tracked separately via `substrate_library_version`.
const SCHEMA_VERSION: &str = "1";

/// Default operating frequency when the prompt does not specify one.
///
/// 2.4 GHz is the historical "default ISM" target on this repo (`mom-003`
/// patch resonance) and the most common implicit frequency in informal
/// antenna prompts.
const DEFAULT_FREQUENCY_HZ: f64 = 2.4e9;

/// Default substrate when the prompt does not name one. FR4 because it is
/// what every "I have a PCB" engineer assumes by default.
const DEFAULT_SUBSTRATE_NAME: &str = "FR4";

/// Canonical substrate library names recognised by the offline parser, in the
/// order we test them. Longer / more specific names come first so a prompt
/// like `"RO4003C"` is not partially-matched by a hypothetical shorter
/// `"RO40"` prefix.
const SUBSTRATE_NAMES: &[&str] = &["RO4003C", "RO5880", "AluminaTC", "FR4"];

/// Errors produced by [`parse`].
///
/// The parser is forgiving by design — most failure modes degrade to defaults
/// — so this enum is small. The one hard-fail is an empty prompt: returning
/// a "default" intent for an empty prompt would mask programmer error
/// (forgetting to pass the prompt at all).
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// The supplied prompt was empty or all-whitespace.
    #[error("yee-design::offline: prompt is empty")]
    EmptyPrompt,
}

/// Compiled regex matching a frequency literal with units.
///
/// Captures: 1 = numeric (`\d+(\.\d+)?`), 3 = unit (`GHz` / `MHz`, case-
/// insensitive via the `(?i)` flag at the head of the pattern). The middle
/// `\s*` allows `"2.4GHz"` and `"2.4 GHz"` equally.
fn frequency_regex() -> &'static Regex {
    static RE: OnceLock<Regex> = OnceLock::new();
    RE.get_or_init(|| {
        Regex::new(r"(?i)(\d+(\.\d+)?)\s*(GHz|MHz)")
            .expect("yee-design::offline: frequency regex must compile")
    })
}

/// Parse a free-form natural-language prompt into a [`DesignIntent`].
///
/// Forgiving keyword-based extraction; documented defaults fill in the gaps.
/// The geometry family is always [`GeometryFamily::RectangularPatch`] because
/// Phase 3.nl.0 ships no other family.
///
/// Returns [`Error::EmptyPrompt`] iff the prompt (after trimming) is empty.
pub fn parse(prompt: &str) -> Result<DesignIntent, Error> {
    if prompt.trim().is_empty() {
        return Err(Error::EmptyPrompt);
    }

    let target_frequency_hz = extract_frequency_hz(prompt).unwrap_or(DEFAULT_FREQUENCY_HZ);
    let substrate_name =
        extract_substrate_name(prompt).unwrap_or_else(|| DEFAULT_SUBSTRATE_NAME.to_string());

    let library_version = substrate_library().version.clone();

    Ok(DesignIntent {
        family: GeometryFamily::RectangularPatch,
        target_frequency_hz,
        substrate: Substrate::Named(NamedSubstrate {
            name: substrate_name,
            override_with: None,
        }),
        gain_target_dbi: None,
        bandwidth_target_mhz: None,
        source_prompt: prompt.to_string(),
        provenance: Provenance {
            source: "offline".to_string(),
            model: None,
            temperature: None,
            schema_version: SCHEMA_VERSION.to_string(),
            substrate_library_version: library_version,
        },
    })
}

/// Extract a frequency literal in hertz from the prompt, if any.
///
/// Returns the first match's value × {`GHz` → `1e9`, `MHz` → `1e6`}. Returns
/// `None` if no frequency literal is found, so the caller can pick a default.
fn extract_frequency_hz(prompt: &str) -> Option<f64> {
    let caps = frequency_regex().captures(prompt)?;
    let value: f64 = caps.get(1)?.as_str().parse().ok()?;
    let unit = caps.get(3)?.as_str().to_ascii_uppercase();
    let multiplier = match unit.as_str() {
        "GHZ" => 1.0e9,
        "MHZ" => 1.0e6,
        _ => return None,
    };
    Some(value * multiplier)
}

/// Find the first library-substrate name that appears in the prompt
/// (case-insensitive). Returns the canonical (capitalised) library key, not
/// the prompt's original capitalisation.
fn extract_substrate_name(prompt: &str) -> Option<String> {
    let needle = prompt.to_ascii_lowercase();
    for canonical in SUBSTRATE_NAMES {
        if needle.contains(&canonical.to_ascii_lowercase()) {
            return Some((*canonical).to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn named_substrate(intent: &DesignIntent) -> &str {
        match &intent.substrate {
            Substrate::Named(NamedSubstrate { name, .. }) => name.as_str(),
            Substrate::Explicit { .. } => panic!("offline parser must emit a Named substrate"),
        }
    }

    #[test]
    fn parses_canonical_2_4_ghz_fr4_prompt() {
        let intent = parse("2.4 GHz patch on FR4").expect("parse");
        assert_eq!(intent.target_frequency_hz, 2.4e9);
        assert_eq!(named_substrate(&intent), "FR4");
        assert_eq!(intent.family, GeometryFamily::RectangularPatch);
        assert_eq!(intent.source_prompt, "2.4 GHz patch on FR4");
        assert_eq!(intent.provenance.source, "offline");
    }

    #[test]
    fn parses_mhz_unit() {
        let intent = parse("915 MHz patch for IoT").expect("parse");
        assert_eq!(intent.target_frequency_hz, 915.0e6);
        // No substrate named → falls back to FR4 default.
        assert_eq!(named_substrate(&intent), "FR4");
    }

    #[test]
    fn parses_no_space_unit() {
        let intent = parse("5.8GHz patch on RO4003C").expect("parse");
        assert_eq!(intent.target_frequency_hz, 5.8e9);
        assert_eq!(named_substrate(&intent), "RO4003C");
    }

    #[test]
    fn defaults_when_no_frequency() {
        let intent = parse("patch antenna on RO5880").expect("parse");
        assert_eq!(intent.target_frequency_hz, DEFAULT_FREQUENCY_HZ);
        assert_eq!(named_substrate(&intent), "RO5880");
    }

    #[test]
    fn defaults_when_no_substrate() {
        let intent = parse("3.5 GHz patch").expect("parse");
        assert_eq!(intent.target_frequency_hz, 3.5e9);
        assert_eq!(named_substrate(&intent), "FR4");
    }

    #[test]
    fn matches_substrate_case_insensitive() {
        let intent = parse("10 GHz on ro4003c").expect("parse");
        // Canonical capitalisation, regardless of input casing.
        assert_eq!(named_substrate(&intent), "RO4003C");
    }

    #[test]
    fn empty_prompt_rejected() {
        assert!(matches!(parse(""), Err(Error::EmptyPrompt)));
        assert!(matches!(parse("   "), Err(Error::EmptyPrompt)));
    }

    #[test]
    fn deterministic_same_prompt_same_intent() {
        let a = parse("2.4 GHz patch on FR4").unwrap();
        let b = parse("2.4 GHz patch on FR4").unwrap();
        assert_eq!(a, b);
    }
}
