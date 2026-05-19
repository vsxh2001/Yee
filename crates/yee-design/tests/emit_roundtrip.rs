//! Phase 3.nl.0 R3 — integration tests for the deterministic project-TOML
//! emitter.
//!
//! The emitter is the spec §8 reproducibility gatekeeper: two calls with the
//! same `(estimate, intent)` MUST produce byte-identical `toml` and
//! `intent_json` strings. The four tests in this module pin that contract:
//!
//! 1. `emit_roundtrip_2g4_fr4_patch` — 2.4 GHz FR-4 patch round-trips through
//!    the `toml` parser back to a `toml::Value` with all geometry / substrate
//!    / frequency fields preserved bit-equal to the original
//!    [`InitialEstimate`] + [`DesignIntent`].
//! 2. `emit_roundtrip_5g8_ro4003c_patch` — same for 5.8 GHz on RO4003C
//!    (different substrate, different frequency decade).
//! 3. `emit_is_byte_identical_across_two_calls` — the spec §8 determinism
//!    gate: `emit(e, i) == emit(e, i)` byte-for-byte on both fields.
//! 4. `emit_intent_json_serde_matches_field_by_field` — the
//!    `ProjectFile::intent_json` field is exactly
//!    `serde_json::to_string(intent)` byte-for-byte. This is the contract
//!    the `<out>.intent.json` sidecar artefact relies on.
//!
//! References:
//! - `docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md` §8.
//! - `docs/superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md` R3.

use yee_design::{
    DesignIntent, GeometryFamily, InitialEstimate, NamedSubstrate, ProjectFile, Provenance,
    Substrate, emit, substrate_library,
};

/// Build a Phase 3.nl.0 [`DesignIntent`] for the given named substrate +
/// target frequency. The bandwidth / gain targets are populated with stable
/// non-trivial values so the deterministic-rerender test exercises every
/// optional field.
fn make_intent(name: &str, f_hz: f64, prompt: &str) -> DesignIntent {
    DesignIntent {
        family: GeometryFamily::RectangularPatch,
        target_frequency_hz: f_hz,
        substrate: Substrate::Named(NamedSubstrate {
            name: name.to_string(),
            override_with: None,
        }),
        gain_target_dbi: Some(6.0),
        bandwidth_target_mhz: Some(100.0),
        source_prompt: prompt.to_string(),
        provenance: Provenance {
            source: "offline".to_string(),
            model: None,
            temperature: None,
            schema_version: "1".to_string(),
            substrate_library_version: substrate_library().version.clone(),
        },
    }
}

/// Pull the float at `path[..]` (sequence of TOML keys) out of a parsed
/// `toml::Value`. Panics with a descriptive message if any step misses; this
/// is a test helper, not production code.
fn get_f64(root: &toml::Value, path: &[&str]) -> f64 {
    let mut cur = root;
    for (i, key) in path.iter().enumerate() {
        cur = cur
            .get(key)
            .unwrap_or_else(|| panic!("path {path:?} missing key '{key}' at step {i}"));
    }
    cur.as_float()
        .unwrap_or_else(|| panic!("path {path:?} is not a float; got {cur:?}"))
}

/// Pull a TOML string at `path[..]`.
fn get_str<'a>(root: &'a toml::Value, path: &[&str]) -> &'a str {
    let mut cur = root;
    for (i, key) in path.iter().enumerate() {
        cur = cur
            .get(key)
            .unwrap_or_else(|| panic!("path {path:?} missing key '{key}' at step {i}"));
    }
    cur.as_str()
        .unwrap_or_else(|| panic!("path {path:?} is not a string; got {cur:?}"))
}

/// Compare two floats that should be bit-identical *after* the `{:.6e}`
/// round trip — i.e. equal to within the 6-decimal-digit canonical format.
///
/// The emitter formats every f64 via `format!("{:.6e}", x)`, then a TOML
/// parser reads it back. We assert relative error ≤ 1e-6 (one ULP of the
/// 6-significant-digit literal), which is the natural fidelity ceiling of
/// the canonical format.
fn assert_close(actual: f64, expected: f64, name: &str) {
    let denom = expected.abs().max(1.0e-12);
    let rel = (actual - expected).abs() / denom;
    assert!(
        rel <= 1.0e-6,
        "{name}: actual = {actual:.10e}, expected = {expected:.10e}, rel err = {rel:.3e}"
    );
}

/// Emit a 2.4 GHz FR-4 patch project file, re-parse the TOML, and confirm
/// every geometry / substrate / frequency scalar survives the round trip
/// within the `{:.6e}` canonical-format fidelity ceiling. The two metadata
/// header lines are stripped before parsing so the round-trip uses the
/// same `toml` value-parser as `yee run` would.
#[test]
fn emit_roundtrip_2g4_fr4_patch() {
    let intent = make_intent("FR4", 2.4e9, "2.4 GHz patch on FR4 with 100 MHz bandwidth");
    let estimate = InitialEstimate::from_intent(&intent).expect("estimate");
    let out = emit(&estimate, &intent);

    let parsed: toml::Value = toml::from_str(&out.toml).expect("emitted TOML must parse");

    // Geometry: every field of `InitialEstimate` is reachable.
    assert_eq!(
        get_str(&parsed, &["geometry", "type"]),
        "rectangular_inset_patch"
    );
    assert_close(
        get_f64(&parsed, &["geometry", "width_m"]),
        estimate.width_m,
        "geometry.width_m",
    );
    assert_close(
        get_f64(&parsed, &["geometry", "length_m"]),
        estimate.length_m,
        "geometry.length_m",
    );
    assert_close(
        get_f64(&parsed, &["geometry", "inset_offset_m"]),
        estimate.inset_offset_m,
        "geometry.inset_offset_m",
    );
    assert_close(
        get_f64(&parsed, &["geometry", "feed_width_m"]),
        estimate.feed_width_m,
        "geometry.feed_width_m",
    );

    // Substrate: matches the resolved (post library-lookup) values.
    assert_close(
        get_f64(&parsed, &["substrate", "eps_r"]),
        estimate.substrate.eps_r,
        "substrate.eps_r",
    );
    assert_close(
        get_f64(&parsed, &["substrate", "h_m"]),
        estimate.substrate.h_m,
        "substrate.h_m",
    );
    assert_close(
        get_f64(&parsed, &["substrate", "loss_tangent"]),
        estimate.substrate.loss_tangent,
        "substrate.loss_tangent",
    );

    // Frequency: centre + span + sweep-point count.
    assert_close(
        get_f64(&parsed, &["frequency", "center_hz"]),
        intent.target_frequency_hz,
        "frequency.center_hz",
    );
    // Span = 2 * bandwidth (intent has 100 MHz).
    let expected_span = 2.0 * 100.0e6;
    assert_close(
        get_f64(&parsed, &["frequency", "span_hz"]),
        expected_span,
        "frequency.span_hz",
    );
    let sweep_points = parsed
        .get("frequency")
        .and_then(|t| t.get("sweep_points"))
        .and_then(|v| v.as_integer())
        .expect("frequency.sweep_points must be int");
    assert_eq!(sweep_points, 201);

    // Ports: at least one entry, and Z₀ = 50 Ω.
    let ports = parsed
        .get("ports")
        .and_then(|p| p.as_array())
        .expect("ports must be an array");
    assert!(!ports.is_empty());
    let port0 = &ports[0];
    assert_eq!(
        port0.get("kind").and_then(|v| v.as_str()),
        Some("delta_gap")
    );
    assert_close(
        port0
            .get("z0_ohm")
            .and_then(|v| v.as_float())
            .expect("z0_ohm"),
        50.0,
        "ports[0].z0_ohm",
    );
}

/// Same shape as [`emit_roundtrip_2g4_fr4_patch`] but on a different
/// substrate (RO4003C) and frequency (5.8 GHz). Confirms the emitter does
/// not bake in any FR4-specific assumptions.
#[test]
fn emit_roundtrip_5g8_ro4003c_patch() {
    let intent = make_intent(
        "RO4003C",
        5.8e9,
        "5.8 GHz patch on RO4003C with 200 MHz bandwidth",
    );
    let estimate = InitialEstimate::from_intent(&intent).expect("estimate");
    let out = emit(&estimate, &intent);

    let parsed: toml::Value = toml::from_str(&out.toml).expect("emitted TOML must parse");

    // Geometry: resolved RO4003C dimensions round-trip.
    assert_eq!(
        get_str(&parsed, &["geometry", "type"]),
        "rectangular_inset_patch"
    );
    assert_close(
        get_f64(&parsed, &["geometry", "width_m"]),
        estimate.width_m,
        "geometry.width_m",
    );
    assert_close(
        get_f64(&parsed, &["geometry", "length_m"]),
        estimate.length_m,
        "geometry.length_m",
    );

    // Substrate: resolved eps_r matches the library row.
    let lib = substrate_library();
    let ro = lib.get("RO4003C").expect("RO4003C in library");
    assert_close(get_f64(&parsed, &["substrate", "eps_r"]), ro.eps_r, "eps_r");
    assert_close(
        get_f64(&parsed, &["substrate", "h_m"]),
        ro.h_mm * 1.0e-3,
        "h_m",
    );

    // Frequency: target_frequency_hz lives on the line after parsing.
    assert_close(
        get_f64(&parsed, &["frequency", "center_hz"]),
        5.8e9,
        "frequency.center_hz",
    );
}

/// Spec §8 determinism gate: a second `emit(e, i)` call must produce a
/// `ProjectFile` that is byte-for-byte identical to the first on both
/// `toml` and `intent_json` fields. This is the property the spec calls
/// out as the distinguishing line between "convenience over a script"
/// and "magic".
#[test]
fn emit_is_byte_identical_across_two_calls() {
    let intent = make_intent("FR4", 2.4e9, "2.4 GHz patch on FR4 with 100 MHz bandwidth");
    let estimate = InitialEstimate::from_intent(&intent).expect("estimate");
    let a: ProjectFile = emit(&estimate, &intent);
    let b: ProjectFile = emit(&estimate, &intent);
    assert_eq!(a.toml, b.toml, "ProjectFile::toml not deterministic");
    assert_eq!(
        a.intent_json, b.intent_json,
        "ProjectFile::intent_json not deterministic"
    );
    // For total safety, also check the structurally-derived `PartialEq`.
    assert_eq!(a, b);
}

/// `ProjectFile::intent_json` must be exactly `serde_json::to_string(intent)`
/// — no extra whitespace, no field reorder, no pretty-printing. This is the
/// contract the `<out>.intent.json` sidecar relies on so the same hash that
/// goes into the TOML header lines up with what re-reading the JSON
/// produces.
#[test]
fn emit_intent_json_serde_matches_field_by_field() {
    let intent = make_intent("FR4", 2.4e9, "2.4 GHz patch on FR4 with 100 MHz bandwidth");
    let estimate = InitialEstimate::from_intent(&intent).expect("estimate");
    let out = emit(&estimate, &intent);
    let direct = serde_json::to_string(&intent).expect("serialize intent");
    assert_eq!(out.intent_json, direct);
    // And a JSON-decoded round-trip yields the original DesignIntent.
    let back: DesignIntent = serde_json::from_str(&out.intent_json).expect("deserialize");
    assert_eq!(back, intent);
}
