//! Typed surface for a parsed natural-language design intent.
//!
//! The shapes in this module mirror spec §6 (`DesignIntent`,
//! `GeometryFamily`, `Substrate`, `NamedSubstrate`, `Provenance`) and are the
//! sole boundary between Stage 1 (LLM / offline parser) and Stages 2–5 of the
//! Phase 3.nl.0 pipeline. Every type derives `serde::{Serialize,
//! Deserialize}` because spec §8 makes the `<out>.intent.json` artefact a
//! reproducibility-critical second-class citizen of every emitted project
//! file.
//!
//! ## Substrate library
//!
//! `substrates.toml` is `include_str!`-loaded at first access into a static
//! [`SubstrateLibrary`]; see [`substrate_library`]. The library has a version
//! string at the top of the file; [`Provenance::substrate_library_version`]
//! records the resolved value so a saved [`DesignIntent`] can detect drift
//! (spec §10 risk #2).

use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

/// Top-level structured representation of a parsed natural-language design
/// intent.
///
/// Spec §6. The verbatim source prompt is preserved in
/// [`DesignIntent::source_prompt`] so the spec §8 reproducibility invariant
/// (the project-file header reproduces the prompt) holds without a separate
/// round-trip from the LLM provider.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct DesignIntent {
    /// Closed-enum geometry family this intent targets.
    pub family: GeometryFamily,
    /// Target operating / resonant frequency in hertz.
    pub target_frequency_hz: f64,
    /// Substrate stack-up — either a library look-up or an explicit override.
    pub substrate: Substrate,
    /// Optional gain target (dBi). Not used by Phase 3.nl.0's initial
    /// estimator; carried for surrogate-refinement (Phase 3.nl.1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gain_target_dbi: Option<f64>,
    /// Optional fractional / absolute bandwidth target in megahertz. As with
    /// [`DesignIntent::gain_target_dbi`], Phase 3.nl.0 records but does not
    /// act on it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bandwidth_target_mhz: Option<f64>,
    /// Verbatim natural-language prompt that produced this intent. Echoed as
    /// the first-line comment of the emitted `yee.toml` (spec §8).
    pub source_prompt: String,
    /// Provenance metadata — model id, temperature, schema version,
    /// substrate-library version. Never carries secrets (spec §10 risk #3).
    pub provenance: Provenance,
}

/// Closed enumeration of supported geometry families.
///
/// Phase 3.nl.0 ships only `RectangularPatch`; spec §12 lists the families
/// that land in 3.nl.2.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum GeometryFamily {
    /// Rectangular inset-fed microstrip patch antenna (Balanis Ch. 14).
    RectangularPatch,
}

/// Substrate selection — either a name resolved against the canonical library
/// or an explicit `{eps_r, h_mm, loss_tangent}` override.
///
/// Mirrors the spec §7 `substrate` `oneOf` shape so the JSON-schema-validated
/// LLM tool-use response and the serde-derived shape stay byte-identical.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(untagged)]
pub enum Substrate {
    /// Named lookup against [`substrate_library`].
    Named(NamedSubstrate),
    /// Caller-supplied explicit substrate parameters.
    Explicit {
        /// Relative permittivity (unitless).
        eps_r: f64,
        /// Substrate height in millimetres.
        h_mm: f64,
        /// Loss tangent (tan δ, unitless).
        loss_tangent: f64,
    },
}

/// Name + optional explicit override pair.
///
/// The `override_with` field is present so a downstream stage can record a
/// substrate that was looked up but then over-ridden (e.g. the offline parser
/// recognised `"FR4"` but the prompt also specified `h = 0.8 mm`). Phase
/// 3.nl.0 may emit it with `None`; later stages may emit `Some` after
/// merging.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct NamedSubstrate {
    /// Library-key name — must match a row in `substrates.toml`.
    pub name: String,
    /// Optional override for one or more of the named substrate's fields.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub override_with: Option<SubstrateOverride>,
}

/// Partial override of a named substrate's parameters.
///
/// All fields are optional; only the ones the user explicitly specified are
/// populated.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SubstrateOverride {
    /// Override the relative permittivity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub eps_r: Option<f64>,
    /// Override the substrate height in millimetres.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub h_mm: Option<f64>,
    /// Override the loss tangent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loss_tangent: Option<f64>,
}

/// Provenance metadata about an intent's source.
///
/// Recorded so a re-run from `<out>.intent.json` is traceable (spec §8) and
/// so a substrate-library drift is detectable (spec §10 risk #2). Never
/// contains API keys or other secrets (spec §10 risk #3).
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Provenance {
    /// `"llm"` or `"offline"`; identifies the Stage 1 parser.
    pub source: String,
    /// LLM model id (e.g. `"claude-sonnet-4-5"`); `None` for the offline
    /// fallback.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub model: Option<String>,
    /// Sampling temperature; `None` when not applicable.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature: Option<f64>,
    /// Schema version of the structured-output contract this intent satisfies.
    pub schema_version: String,
    /// Version of `substrates.toml` resolved when this intent was constructed.
    pub substrate_library_version: String,
}

/// Canonical substrate library row.
///
/// Each row of `substrates.toml` deserialises into one of these. Field
/// semantics mirror the spec §7 explicit-substrate object.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SubstrateRecord {
    /// Library key (matches `NamedSubstrate::name`).
    pub name: String,
    /// Relative permittivity.
    pub eps_r: f64,
    /// Substrate height in millimetres.
    pub h_mm: f64,
    /// Loss tangent.
    pub loss_tangent: f64,
}

/// In-memory representation of `substrates.toml`.
///
/// Versioned at the top of the file; [`SubstrateLibrary::version`] is the
/// value recorded in [`Provenance::substrate_library_version`] when an intent
/// resolves a named substrate.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SubstrateLibrary {
    /// Library version string. Bump on any row change.
    pub version: String,
    /// All known substrates.
    #[serde(default, rename = "substrate")]
    pub substrates: Vec<SubstrateRecord>,
}

impl SubstrateLibrary {
    /// Look up a substrate by name (case-sensitive).
    ///
    /// Returns `None` for an unknown name; resolving the named-substrate /
    /// schema mismatch is Stage 2's job.
    pub fn get(&self, name: &str) -> Option<&SubstrateRecord> {
        self.substrates.iter().find(|s| s.name == name)
    }
}

/// Raw `substrates.toml` contents, baked in at build time.
const SUBSTRATES_TOML: &str = include_str!("../substrates.toml");

/// Access the canonical substrate library (lazy-initialised, process-wide).
///
/// Panics at startup only if the baked-in `substrates.toml` fails to parse —
/// which would be a build-time mistake, not a runtime input failure.
pub fn substrate_library() -> &'static SubstrateLibrary {
    static LIB: OnceLock<SubstrateLibrary> = OnceLock::new();
    LIB.get_or_init(|| {
        toml::from_str::<SubstrateLibrary>(SUBSTRATES_TOML)
            .expect("yee-design: substrates.toml failed to parse (baked-in asset is malformed)")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substrate_library_loads_and_contains_canonical_entries() {
        let lib = substrate_library();
        assert!(!lib.version.is_empty(), "library version must be non-empty");
        for name in ["FR4", "RO4003C", "RO5880", "AluminaTC"] {
            let row = lib
                .get(name)
                .unwrap_or_else(|| panic!("substrate '{name}' must be in library"));
            assert!(row.eps_r > 1.0, "{name} eps_r > 1");
            assert!(row.h_mm > 0.0, "{name} h_mm > 0");
            assert!(
                (0.0..=0.5).contains(&row.loss_tangent),
                "{name} loss_tangent in [0, 0.5]"
            );
        }
    }

    #[test]
    fn design_intent_named_substrate_round_trips() {
        let intent = DesignIntent {
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
        };
        let json = serde_json::to_string(&intent).expect("serialize");
        let back: DesignIntent = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(intent, back);
    }

    #[test]
    fn design_intent_explicit_substrate_round_trips() {
        let intent = DesignIntent {
            family: GeometryFamily::RectangularPatch,
            target_frequency_hz: 10.0e9,
            substrate: Substrate::Explicit {
                eps_r: 3.0,
                h_mm: 0.508,
                loss_tangent: 0.0027,
            },
            gain_target_dbi: Some(6.0),
            bandwidth_target_mhz: None,
            source_prompt: "10 GHz patch, explicit substrate".to_string(),
            provenance: Provenance {
                source: "llm".to_string(),
                model: Some("claude-sonnet-4-5".to_string()),
                temperature: Some(0.0),
                schema_version: "1".to_string(),
                substrate_library_version: substrate_library().version.clone(),
            },
        };
        let json = serde_json::to_string(&intent).expect("serialize");
        let back: DesignIntent = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(intent, back);
    }
}
