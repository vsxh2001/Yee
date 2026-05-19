//! `nl-001` production validation gate — Phase 3.nl.0 R6 integration
//! tests.
//!
//! Two entry points, matching the R6 brief's "Fast smoke test
//! (always-on)" / "Slow production gate (`#[ignore]`'d)" split:
//!
//! - [`nl_001_offline_schema_roundtrip_passes_all_ten_prompts`] —
//!   exercises spec §9 sub-gates A (offline), B (schema), and C
//!   (round-trip) for every canonical prompt shipped in
//!   `crates/yee-design/validation/prompts.toml`. Pure-Rust, no
//!   network, no filesystem; finishes well under 30 s wall in
//!   `--release`. Always-on in default CI.
//! - [`nl_001_solver_gate_runs_under_thirty_minutes_all_ten_prompts`]
//!   — exercises sub-gate D (full solve, `±5 %` on `f_min`). Marked
//!   `#[ignore]` per the R6 brief's "wall-time < 30 min `--release`;
//!   hardware-gate behind `#[ignore]` if overrun" and per CLAUDE.md §4
//!   `mom-001` precedent (the only existing patch driver, `mom-003`,
//!   is itself [`yee_validation::Status::Skipped`] pending Phase 1.1.1
//!   per CLAUDE.md §10).

use yee_validation::nl_001::{
    canonical_prompts, run_nl_001_offline_schema_roundtrip, run_nl_001_solver_gate,
};

/// Fast smoke test — always-on in default CI.
///
/// Asserts sub-gates A + B + C pass for every prompt in spec §9's
/// canonical 10. Runs in sub-second wall-time in `--release` (pure-
/// Rust offline parse + serde round-trip + emit). The slow solver
/// sub-gate D is the responsibility of the separate `#[ignore]`'d test
/// below.
#[test]
fn nl_001_offline_schema_roundtrip_passes_all_ten_prompts() {
    let prompts = canonical_prompts();
    assert_eq!(
        prompts.len(),
        10,
        "spec §9 ships exactly 10 canonical prompts; manifest now has {}",
        prompts.len()
    );

    let mut failed: Vec<String> = Vec::new();
    for prompt in &prompts {
        let r = run_nl_001_offline_schema_roundtrip(prompt);
        if !r.all_passed() {
            failed.push(format!("'{}' — {}", prompt, r.notes));
        }
    }

    assert!(
        failed.is_empty(),
        "nl-001 sub-gates A+B+C failed on {} prompt(s):\n  - {}",
        failed.len(),
        failed.join("\n  - ")
    );
}

/// Per-prompt visibility: a separate assertion per sub-gate so a
/// failure in the offline parser is distinguishable from a failure in
/// the schema validator or in the round-trip emit.
///
/// This complements the bulk pass/fail above; both can run in default
/// CI because the cost is sub-second.
#[test]
fn nl_001_each_prompt_passes_each_of_three_sub_gates_individually() {
    for prompt in canonical_prompts() {
        let r = run_nl_001_offline_schema_roundtrip(&prompt);
        assert!(
            r.offline_passed,
            "sub-gate A (offline) failed on '{}': {}",
            prompt, r.notes
        );
        assert!(
            r.schema_passed,
            "sub-gate B (schema) failed on '{}': {}",
            prompt, r.notes
        );
        assert!(
            r.roundtrip_passed,
            "sub-gate C (round-trip) failed on '{}': {}",
            prompt, r.notes
        );
    }
}

/// Slow production gate — sub-gate D (`±5 %` on `f_min`).
///
/// `#[ignore]`'d per:
///
/// 1. R6 brief — "wall-time < 30 min `--release`; hardware-gate behind
///    `#[ignore]` if overrun." The patch-resonance solve underneath is
///    measured at ~minutes per prompt × 10 prompts.
/// 2. CLAUDE.md §4 — `mom-001` precedent for separating multi-minute
///    solver tests from the always-on lint floor.
/// 3. CLAUDE.md §10 — `MultilayerGreens` is a Phase 1.1.0 placeholder;
///    the closest existing patch driver (`mom-003`) is itself
///    [`yee_validation::Status::Skipped`] pending Phase 1.1.1's real
///    Sommerfeld-integral / multi-image DCIM extraction. The current
///    body of [`run_nl_001_solver_gate`] therefore returns
///    [`yee_validation::Status::Skipped`] with a notes string spelling
///    out the upstream dependency; this test does not assert
///    `Passed` (which would require a working patch driver) — it
///    asserts the slow-gate plumbing returns a structured
///    `CaseResult` and that sub-gates A+B+C are met as a
///    precondition.
///
/// Run explicitly with:
/// `cargo test -p yee-validation --release --test nl_001_canonical_prompts \
///     -- --include-ignored`.
#[test]
#[ignore = "slow: nl-001 solver gate is `#[ignore]`'d pending Phase 1.1.1 MultilayerGreens"]
fn nl_001_solver_gate_runs_under_thirty_minutes_all_ten_prompts() {
    let prompts = canonical_prompts();
    let start = std::time::Instant::now();
    let thirty_minutes = std::time::Duration::from_secs(30 * 60);

    for prompt in &prompts {
        let r = run_nl_001_solver_gate(prompt);
        // Per the R6 brief and CLAUDE.md §10, `Skipped` is the
        // expected status until Phase 1.1.1 lands; `Failed` here
        // means a sub-gate A/B/C precondition fired, which is a real
        // bug.
        assert!(
            !matches!(r.status, yee_validation::Status::Failed),
            "nl-001 solver gate FAILED on '{}': {} ({})",
            prompt,
            r.notes,
            r.id
        );
    }

    let elapsed = start.elapsed();
    assert!(
        elapsed < thirty_minutes,
        "nl-001 solver gate exceeded 30 min wall-time budget (R6 DoD): \
            ran {} prompts in {:?}",
        prompts.len(),
        elapsed
    );
}
