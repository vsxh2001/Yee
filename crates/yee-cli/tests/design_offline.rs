//! Integration test for `yee design --offline` over the 10 canonical prompts.
//!
//! Phase 3.nl.0 plan R5 DoD. Loops over
//! `crates/yee-design/validation/prompts.toml`, invokes the `yee` binary
//! once per prompt, and asserts four sub-gates per prompt:
//!
//! 1. Exit code 0.
//! 2. The project-TOML output file exists and is non-empty.
//! 3. The `<output>.intent.json` sidecar exists and is non-empty.
//! 4. The emitted TOML parses, and `[frequency].center_hz` matches the
//!    frequency the offline parser extracted from the prompt within ±0.1 %
//!    (e.g. `"2.4 GHz"` → `2.4e9 ± 2.4 MHz`).
//!
//! Output paths land under the per-test `CARGO_TARGET_TMPDIR` (set by cargo
//! for integration tests) so parallel test runs do not collide and no
//! `/tmp/...` cleanup is required.

use std::path::PathBuf;
use std::process::Command;

use serde::Deserialize;

/// Path to the canonical-prompts manifest, baked at compile time via
/// `env!("CARGO_MANIFEST_DIR")` so the test does not depend on the cwd at
/// runtime. The path resolves to
/// `<repo>/crates/yee-design/validation/prompts.toml`.
fn prompts_manifest_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("yee-design")
        .join("validation")
        .join("prompts.toml")
}

#[derive(Deserialize)]
struct PromptsManifest {
    prompts: Vec<String>,
}

#[derive(Deserialize)]
struct EmittedFrequency {
    center_hz: f64,
}

#[derive(Deserialize)]
struct EmittedProject {
    frequency: EmittedFrequency,
}

/// Re-implement the offline parser's frequency-extraction rule so this test
/// asserts the *prompt's* expected frequency, not the parser's output. The
/// regex is deliberately the same shape as `yee_design::offline` — keep them
/// in sync if either changes.
fn expected_frequency_hz(prompt: &str) -> f64 {
    let re = regex::Regex::new(r"(?i)(\d+(\.\d+)?)\s*(GHz|MHz)").unwrap();
    if let Some(caps) = re.captures(prompt) {
        let value: f64 = caps.get(1).unwrap().as_str().parse().unwrap();
        let unit = caps.get(3).unwrap().as_str().to_ascii_uppercase();
        let mult = match unit.as_str() {
            "GHZ" => 1.0e9,
            "MHZ" => 1.0e6,
            _ => 2.4e9,
        };
        return value * mult;
    }
    2.4e9
}

#[test]
fn ten_canonical_prompts_round_trip_offline() {
    let manifest_path = prompts_manifest_path();
    let manifest_text = std::fs::read_to_string(&manifest_path)
        .unwrap_or_else(|e| panic!("read {}: {e}", manifest_path.display()));
    let manifest: PromptsManifest = toml::from_str(&manifest_text).expect("parse prompts.toml");
    assert_eq!(
        manifest.prompts.len(),
        10,
        "spec §9 / plan R5 require exactly 10 canonical prompts, got {}",
        manifest.prompts.len()
    );

    let tmp_root = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(&tmp_root).expect("create tmp dir");

    for (i, prompt) in manifest.prompts.iter().enumerate() {
        let out_path = tmp_root.join(format!("yee-cli-design-test-{i}.toml"));
        let intent_path = {
            let mut s = out_path.as_os_str().to_owned();
            s.push(".intent.json");
            PathBuf::from(s)
        };
        // Force a clean slate so a stale artefact from a previous run cannot
        // accidentally satisfy the existence assertion below.
        let _ = std::fs::remove_file(&out_path);
        let _ = std::fs::remove_file(&intent_path);

        // 1. Exit code 0.
        let output = Command::new(env!("CARGO_BIN_EXE_yee"))
            .args(["design", prompt, "-o"])
            .arg(&out_path)
            .arg("--offline")
            // Scrub ANTHROPIC_API_KEY so a developer with the env var set
            // does not accidentally route through the LLM-not-wired stub.
            .env_remove("ANTHROPIC_API_KEY")
            .output()
            .expect("invoke yee");

        assert!(
            output.status.success(),
            "prompt #{i} `{prompt}`: yee exited non-zero. stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );

        // 2. Project TOML exists and is non-empty.
        let toml_bytes = std::fs::read(&out_path).unwrap_or_else(|e| {
            panic!(
                "prompt #{i} `{prompt}`: project TOML {} missing: {e}",
                out_path.display()
            )
        });
        assert!(
            !toml_bytes.is_empty(),
            "prompt #{i} `{prompt}`: project TOML is empty"
        );

        // 3. Intent JSON sidecar exists and is non-empty.
        let intent_bytes = std::fs::read(&intent_path).unwrap_or_else(|e| {
            panic!(
                "prompt #{i} `{prompt}`: intent.json sidecar {} missing: {e}",
                intent_path.display()
            )
        });
        assert!(
            !intent_bytes.is_empty(),
            "prompt #{i} `{prompt}`: intent.json sidecar is empty"
        );

        // 4. Emitted TOML parses and [frequency].center_hz matches the prompt
        //    within ±0.1 %.
        let toml_text = std::str::from_utf8(&toml_bytes)
            .unwrap_or_else(|_| panic!("prompt #{i} `{prompt}`: TOML is not UTF-8"));
        let emitted: EmittedProject = toml::from_str(toml_text).unwrap_or_else(|e| {
            panic!("prompt #{i} `{prompt}`: emitted TOML failed to parse: {e}\n---\n{toml_text}")
        });

        let expected = expected_frequency_hz(prompt);
        let actual = emitted.frequency.center_hz;
        let rel = (actual - expected).abs() / expected;
        assert!(
            rel < 1.0e-3,
            "prompt #{i} `{prompt}`: center_hz = {actual} differs from expected {expected} by {rel:.4} (>0.1%)"
        );
    }
}
