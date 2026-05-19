//! `nl-001` production validation gate — Phase 3.nl.0 R6.
//!
//! Implements the spec §9 four-sub-gate composition for the 10 canonical
//! prompts shipped in `crates/yee-design/validation/prompts.toml`:
//!
//! - **Sub-gate A (offline)** — every prompt parses to a [`DesignIntent`]
//!   via [`yee_design::parse_offline`].
//! - **Sub-gate B (schema)** — every parsed intent passes the spec §7
//!   structured-output schema, checked here against the baked-in
//!   [`yee_design::intent::INTENT_SCHEMA`] document via the small
//!   hand-rolled validator in [`validate_intent_against_schema`]. The R6
//!   brief endorses a hand-rolled validator as a documented escape hatch
//!   when a maintained pure-Rust `jsonschema` is not on the crate path;
//!   the schema is small enough (one closed `family` enum, one ranged
//!   `target_frequency_hz`, one `oneOf` substrate shape, two optional
//!   numeric ranges) that we cover every spec §7 required-field /
//!   enum / range invariant inline.
//! - **Sub-gate C (round-trip)** — `emit` once → re-parse the
//!   `intent.json` sidecar → `emit` again → byte-identical TOML and
//!   byte-identical `intent_json`. Mirrors the spec §8 determinism
//!   property.
//! - **Sub-gate D (solver, ±5 % f)** — the heavyweight `yee run`-style
//!   full solve on the emitted project TOML. Per the R6 brief, this is
//!   marked `#[ignore]` at the integration-test layer because the
//!   `mom-003` patch-resonance driver underneath is itself
//!   [`CaseStatus::Skipped`] pending Phase 1.1.1's real
//!   `MultilayerGreens`; the present module exposes the gate **as a
//!   structured `CaseResult`** for completeness so the slow-test entry
//!   point in `tests/nl_001_canonical_prompts.rs` can invoke it
//!   uniformly with sub-gates A–C.
//!
//! ## Default vs `#[ignore]`'d behaviour (CLAUDE.md §4 precedent)
//!
//! The fast smoke test [`run_nl_001_offline_schema_roundtrip`] exercises
//! sub-gates A + B + C for every prompt and runs in well under 30 s
//! wall-time in `--release`. The slow [`run_nl_001_solver_gate`] driver
//! is referenced from a `#[ignore]`'d integration test, matching the
//! `mom-001` precedent of separating the always-on lint-floor check
//! from the multi-minute solver invocation.

#![warn(missing_docs)]

use std::time::Instant;

use crate::{CaseResult, CaseStatus};
use yee_design::{DesignIntent, GeometryFamily, NamedSubstrate, Substrate, emit, parse_offline};

/// Raw `prompts.toml` shipped at `crates/yee-design/validation/prompts.toml`.
///
/// Baked in via `include_str!` so the production gate has no
/// filesystem dependency at test time. The TOML schema is a single
/// `prompts: [String]` array (spec §9, plan R5).
const PROMPTS_TOML: &str = include_str!("../../yee-design/validation/prompts.toml");

/// Container for the parsed `prompts.toml` manifest.
#[derive(serde::Deserialize)]
struct PromptManifest {
    prompts: Vec<String>,
}

/// Return the 10 canonical Phase 3.nl.0 prompts in declaration order.
///
/// Panics at startup only if the baked-in manifest fails to parse —
/// which would be a build-time mistake, not a runtime input failure.
pub fn canonical_prompts() -> Vec<String> {
    let manifest: PromptManifest = toml::from_str(PROMPTS_TOML)
        .expect("yee-validation::nl_001: prompts.toml failed to parse (baked-in asset malformed)");
    manifest.prompts
}

/// Errors produced by the hand-rolled schema validator.
///
/// Each variant carries enough context (field name, observed value) to
/// pinpoint the spec §7 invariant that fired. The full schema document
/// is available via [`yee_design::intent::INTENT_SCHEMA`] for cross-
/// reference.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum NlSchemaError {
    /// `target_frequency_hz` outside `[1 MHz, 1 THz]` (spec §7).
    #[error("nl-001 schema: target_frequency_hz = {0} Hz outside [1e6, 1e12]")]
    FrequencyOutOfRange(f64),
    /// `substrate.eps_r` outside `[1.0, 100.0]` (spec §7 explicit
    /// branch).
    #[error("nl-001 schema: substrate.eps_r = {0} outside [1.0, 100.0]")]
    EpsROutOfRange(f64),
    /// `substrate.h_mm` outside `[0.05, 10.0]` (spec §7 explicit
    /// branch).
    #[error("nl-001 schema: substrate.h_mm = {0} outside [0.05, 10.0]")]
    HMmOutOfRange(f64),
    /// `substrate.loss_tangent` outside `[0.0, 0.5]` (spec §7 explicit
    /// branch).
    #[error("nl-001 schema: substrate.loss_tangent = {0} outside [0.0, 0.5]")]
    LossTangentOutOfRange(f64),
    /// Named substrate's name is not in the spec §7 enum
    /// (`FR4`, `RO4003C`, `RO5880`, `AluminaTC`).
    #[error("nl-001 schema: substrate.name = '{0}' not in spec §7 enum")]
    NamedSubstrateNotInEnum(String),
    /// `gain_target_dbi` outside `[-10.0, 40.0]` (spec §7).
    #[error("nl-001 schema: gain_target_dbi = {0} outside [-10.0, 40.0]")]
    GainOutOfRange(f64),
    /// `bandwidth_target_mhz` negative (spec §7 `minimum = 0.0`).
    #[error("nl-001 schema: bandwidth_target_mhz = {0} < 0")]
    BandwidthNegative(f64),
}

/// Canonical names from the spec §7 substrate enum.
///
/// Matches `crates/yee-design/src/intent_schema.json`'s
/// `properties.substrate.oneOf[0].properties.name.enum`. Kept inline so
/// the validator does not have to re-parse the JSON schema at every
/// call; the unit test in this module cross-checks the inline list
/// against the schema document.
const SUBSTRATE_NAME_ENUM: &[&str] = &["FR4", "RO4003C", "RO5880", "AluminaTC"];

/// Hand-rolled validator for a [`DesignIntent`] against the spec §7
/// schema.
///
/// Returns `Ok(())` iff every required field is present (the type
/// system already enforces this — `family`, `target_frequency_hz`,
/// `substrate` are non-`Option`), every closed enum value is in range,
/// and every ranged number is within its spec §7 bounds. The function
/// **does not** re-check that `family == "rectangular_patch"` —
/// [`GeometryFamily`] is a closed enum in the Phase 3.nl.0 type
/// surface, so any other variant is a compile-time impossibility.
///
/// Sub-gate B of the spec §9 gate composition. The schema document
/// itself is available via [`yee_design::intent::INTENT_SCHEMA`].
pub fn validate_intent_against_schema(intent: &DesignIntent) -> Result<(), NlSchemaError> {
    // 1. target_frequency_hz ∈ [1e6, 1e12]
    let f = intent.target_frequency_hz;
    if !(1.0e6..=1.0e12).contains(&f) || !f.is_finite() {
        return Err(NlSchemaError::FrequencyOutOfRange(f));
    }

    // 2. substrate (oneOf: named with optional override, or explicit)
    match &intent.substrate {
        Substrate::Named(NamedSubstrate {
            name,
            override_with,
        }) => {
            if !SUBSTRATE_NAME_ENUM.contains(&name.as_str()) {
                return Err(NlSchemaError::NamedSubstrateNotInEnum(name.clone()));
            }
            if let Some(ov) = override_with {
                if let Some(v) = ov.eps_r
                    && !(1.0..=100.0).contains(&v)
                {
                    return Err(NlSchemaError::EpsROutOfRange(v));
                }
                if let Some(v) = ov.h_mm
                    && !(0.05..=10.0).contains(&v)
                {
                    return Err(NlSchemaError::HMmOutOfRange(v));
                }
                if let Some(v) = ov.loss_tangent
                    && !(0.0..=0.5).contains(&v)
                {
                    return Err(NlSchemaError::LossTangentOutOfRange(v));
                }
            }
        }
        Substrate::Explicit {
            eps_r,
            h_mm,
            loss_tangent,
        } => {
            if !(1.0..=100.0).contains(eps_r) {
                return Err(NlSchemaError::EpsROutOfRange(*eps_r));
            }
            if !(0.05..=10.0).contains(h_mm) {
                return Err(NlSchemaError::HMmOutOfRange(*h_mm));
            }
            if !(0.0..=0.5).contains(loss_tangent) {
                return Err(NlSchemaError::LossTangentOutOfRange(*loss_tangent));
            }
        }
    }

    // 3. Optional gain_target_dbi ∈ [-10.0, 40.0]
    if let Some(g) = intent.gain_target_dbi
        && !(-10.0..=40.0).contains(&g)
    {
        return Err(NlSchemaError::GainOutOfRange(g));
    }

    // 4. Optional bandwidth_target_mhz ≥ 0
    if let Some(bw) = intent.bandwidth_target_mhz
        && bw < 0.0
    {
        return Err(NlSchemaError::BandwidthNegative(bw));
    }

    // family is closed-enum in Rust; the schema's `family ∈
    // {"rectangular_patch"}` constraint is structurally enforced by
    // [`GeometryFamily`]'s single variant.
    let _ = GeometryFamily::RectangularPatch;

    Ok(())
}

/// Per-prompt sub-gate results.
///
/// The fast-test entry point materialises one of these per canonical
/// prompt. `solver` is `None` when the slow solver gate is not
/// exercised (default `cargo test`); see [`run_nl_001_solver_gate`] for
/// the heavyweight invocation.
#[derive(Debug, Clone)]
pub struct NlSubGateResults {
    /// The verbatim prompt this row corresponds to.
    pub prompt: String,
    /// Sub-gate A: `parse_offline(prompt)` succeeded.
    pub offline_passed: bool,
    /// Sub-gate B: the parsed intent satisfies the spec §7 schema.
    pub schema_passed: bool,
    /// Sub-gate C: `emit` → re-parse `intent.json` → `emit` produces
    /// byte-identical output (TOML body **and** `intent_json` sidecar).
    pub roundtrip_passed: bool,
    /// Diagnostic message (one of "ok", or the first failing sub-gate).
    pub notes: String,
}

impl NlSubGateResults {
    /// `true` iff every sub-gate this row records ran and passed.
    pub fn all_passed(&self) -> bool {
        self.offline_passed && self.schema_passed && self.roundtrip_passed
    }
}

/// Run sub-gates A + B + C for one prompt.
///
/// Pure function over the offline pipeline; no network, no
/// filesystem, sub-second per prompt. Used by both the smoke-test entry
/// point and the slow solver gate (which threads its solve through
/// after this returns `all_passed()`).
pub fn run_nl_001_offline_schema_roundtrip(prompt: &str) -> NlSubGateResults {
    // Sub-gate A: parse via the offline parser.
    let intent = match parse_offline(prompt) {
        Ok(intent) => intent,
        Err(e) => {
            return NlSubGateResults {
                prompt: prompt.to_string(),
                offline_passed: false,
                schema_passed: false,
                roundtrip_passed: false,
                notes: format!("sub-gate A (offline) failed: {e}"),
            };
        }
    };

    // Sub-gate B: schema validation.
    if let Err(e) = validate_intent_against_schema(&intent) {
        return NlSubGateResults {
            prompt: prompt.to_string(),
            offline_passed: true,
            schema_passed: false,
            roundtrip_passed: false,
            notes: format!("sub-gate B (schema) failed: {e}"),
        };
    }

    // Sub-gate C: emit → re-parse intent.json → emit → byte-identical.
    let estimate = match yee_design::InitialEstimate::from_intent(&intent) {
        Ok(est) => est,
        Err(e) => {
            return NlSubGateResults {
                prompt: prompt.to_string(),
                offline_passed: true,
                schema_passed: true,
                roundtrip_passed: false,
                notes: format!("sub-gate C precondition (estimate) failed: {e}"),
            };
        }
    };
    let first = emit(&estimate, &intent);
    let reparsed_intent: DesignIntent = match serde_json::from_str(&first.intent_json) {
        Ok(v) => v,
        Err(e) => {
            return NlSubGateResults {
                prompt: prompt.to_string(),
                offline_passed: true,
                schema_passed: true,
                roundtrip_passed: false,
                notes: format!("sub-gate C (round-trip JSON parse) failed: {e}"),
            };
        }
    };
    let reparsed_estimate = match yee_design::InitialEstimate::from_intent(&reparsed_intent) {
        Ok(est) => est,
        Err(e) => {
            return NlSubGateResults {
                prompt: prompt.to_string(),
                offline_passed: true,
                schema_passed: true,
                roundtrip_passed: false,
                notes: format!("sub-gate C (round-trip estimate) failed: {e}"),
            };
        }
    };
    let second = emit(&reparsed_estimate, &reparsed_intent);

    if first != second {
        return NlSubGateResults {
            prompt: prompt.to_string(),
            offline_passed: true,
            schema_passed: true,
            roundtrip_passed: false,
            notes: "sub-gate C (round-trip): emit not byte-identical".to_string(),
        };
    }

    NlSubGateResults {
        prompt: prompt.to_string(),
        offline_passed: true,
        schema_passed: true,
        roundtrip_passed: true,
        notes: "ok".to_string(),
    }
}

/// Driver for sub-gate D (`±5 %` frequency gate) per the spec §9
/// solver-gate composition.
///
/// The R6 brief notes:
///
/// > For R6 scope, the "solver gate" is the full mom-002 / mom-003-
/// > style solve which is SLOW (~minutes per prompt × 10 = ~30+ min).
/// > Implement two test entry points: ... Slow production gate
/// > (`#[ignore]`'d by default per CLAUDE.md §4 mom-001 precedent).
///
/// Per CLAUDE.md §10, `mom-003` (the 2.4 GHz patch resonance — the
/// closest existing patch driver to what `nl-001` would invoke) is
/// itself [`CaseStatus::Skipped`] pending Phase 1.1.1's real
/// `MultilayerGreens`. The Phase 1.1.0 placeholder cannot produce a
/// meaningful `|S11|` minimum within ±5 % of `f_target` for arbitrary
/// patch dimensions. Wiring the gate to the existing
/// `MultilayerGreens` would assert against numerics that the project
/// already declared too loose to trust.
///
/// This driver therefore returns a [`CaseStatus::Skipped`]
/// [`CaseResult`] with a notes string spelling out the upstream
/// dependency. When Phase 1.1.1 lands, this body switches to invoke
/// the real patch-resonance solve (the I/O shape is already in place
/// via the emitted project TOML at `<out>.toml`, parsed in
/// step R5's CLI); the test entry point in
/// `tests/nl_001_canonical_prompts.rs` stays `#[ignore]`'d.
pub fn run_nl_001_solver_gate(prompt: &str) -> CaseResult {
    let t0 = Instant::now();

    // Sub-gates A + B + C first; if any of those fail the solver gate
    // is meaningless, surface that immediately.
    let sub_abc = run_nl_001_offline_schema_roundtrip(prompt);
    if !sub_abc.all_passed() {
        return CaseResult {
            id: format!("nl-001/{prompt}"),
            description: "Phase 3.nl.0 production gate — sub-gates A–C precondition failed".into(),
            status: CaseStatus::Failed,
            notes: sub_abc.notes,
            wall_time_seconds: t0.elapsed().as_secs_f64(),
            plot_paths: Vec::new(),
        };
    }

    // Sub-gate D: would invoke the patch-resonance solve here. Per
    // CLAUDE.md §10 the only available patch driver is `mom-003`
    // (Skipped); see the docstring for the full rationale.
    CaseResult {
        id: format!("nl-001/{prompt}"),
        description: "Phase 3.nl.0 production gate — solver sub-gate (±5 % f)".into(),
        status: CaseStatus::Skipped,
        notes: "sub-gate D (solver, ±5 % f): the existing `mom-003` patch driver \
             is itself `CaseStatus::Skipped` per CLAUDE.md §10 (`MultilayerGreens` \
             one-image DCIM placeholder; awaiting Phase 1.1.1 Sommerfeld-integral / \
             multi-image DCIM extraction). The R6 production gate inherits that \
             deferral; the sub-gate is documented `#[ignore]`'d at the integration- \
             test layer and re-runs through this driver once Phase 1.1.1 ships."
            .into(),
        wall_time_seconds: t0.elapsed().as_secs_f64(),
        plot_paths: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_prompts_loads_ten() {
        let prompts = canonical_prompts();
        assert_eq!(
            prompts.len(),
            10,
            "spec §9 ships exactly 10 canonical prompts; got {}",
            prompts.len()
        );
    }

    #[test]
    fn schema_validator_accepts_offline_default_intent() {
        let intent = parse_offline("2.4 GHz patch on FR4").unwrap();
        validate_intent_against_schema(&intent).expect("default offline intent must pass schema");
    }

    #[test]
    fn schema_validator_rejects_low_frequency() {
        let mut intent = parse_offline("2.4 GHz patch on FR4").unwrap();
        intent.target_frequency_hz = 1.0; // 1 Hz — below the 1 MHz floor.
        assert!(matches!(
            validate_intent_against_schema(&intent),
            Err(NlSchemaError::FrequencyOutOfRange(_))
        ));
    }

    #[test]
    fn schema_validator_rejects_unknown_named_substrate() {
        let mut intent = parse_offline("2.4 GHz patch on FR4").unwrap();
        intent.substrate = Substrate::Named(NamedSubstrate {
            name: "NotASubstrate".to_string(),
            override_with: None,
        });
        assert!(matches!(
            validate_intent_against_schema(&intent),
            Err(NlSchemaError::NamedSubstrateNotInEnum(_))
        ));
    }

    #[test]
    fn schema_validator_substrate_name_enum_matches_schema_doc() {
        // Cross-check that the inline SUBSTRATE_NAME_ENUM here matches
        // the spec §7 schema document baked into yee-design.
        let v: serde_json::Value = serde_json::from_str(yee_design::intent::INTENT_SCHEMA)
            .expect("INTENT_SCHEMA is valid JSON");
        let names = v
            .pointer("/properties/substrate/oneOf/0/properties/name/enum")
            .and_then(|x| x.as_array())
            .expect("schema /properties/substrate/oneOf/0/.../enum present");
        let from_schema: Vec<&str> = names.iter().filter_map(|x| x.as_str()).collect();
        assert_eq!(
            from_schema, SUBSTRATE_NAME_ENUM,
            "inline SUBSTRATE_NAME_ENUM must mirror the schema document"
        );
    }

    #[test]
    fn run_nl_001_all_three_sub_gates_pass_on_canonical_prompt() {
        let r = run_nl_001_offline_schema_roundtrip("2.4 GHz patch on FR4");
        assert!(r.offline_passed, "offline: {}", r.notes);
        assert!(r.schema_passed, "schema: {}", r.notes);
        assert!(r.roundtrip_passed, "round-trip: {}", r.notes);
        assert!(r.all_passed(), "all sub-gates: {}", r.notes);
    }
}
