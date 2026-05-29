//! Smoke tests for `yee validate`.
//!
//! The fast (`help`) test runs every CI invocation; the
//! `yee_validate_mom_json_runs` test is `#[ignore]` because it invokes
//! the real `yee_validation::Report::run_all`, which executes the
//! `mom-001` 24x176 cylinder dipole solve (~7-8 min wall time in
//! `--release` — see CLAUDE.md §4 / §10). Opt in with
//! `cargo test -p yee-cli -- --include-ignored`.

use std::process::Command;

#[test]
#[ignore = "slow: invokes the real aggregator (mom-001 ~8 min)"]
fn yee_validate_mom_json_runs() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["validate", "mom", "--json"])
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("\"mom-001\""));
}

/// `yee validate --list` must print the registered-case inventory and
/// exit 0 **without running any solver**, so this test is fast and NOT
/// `#[ignore]`'d. Asserts the canonical `mom-001` and `fem-eig-006`
/// ids appear in the table.
#[test]
fn yee_validate_list_runs() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["validate", "--list"])
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee validate --list exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("mom-001"),
        "--list output missing mom-001; got: {stdout}"
    );
    assert!(
        stdout.contains("fem-eig-006"),
        "--list output missing fem-eig-006; got: {stdout}"
    );
}

/// `yee validate --list --json` emits the registered-case inventory as a JSON
/// array, exit 0, running no solver (ADR-0083). Fast, NOT `#[ignore]`'d.
/// Substring assertions avoid needing a JSON parser in the test deps.
#[test]
fn yee_validate_list_json_runs() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["validate", "--list", "--json"])
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee validate --list --json exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let trimmed = stdout.trim();
    assert!(
        trimmed.starts_with('[') && trimmed.ends_with(']'),
        "--list --json output is not a JSON array; got: {stdout}"
    );
    for needle in [
        "\"mom-001\"",
        "\"fem-eig-006\"",
        "\"SkippedGateOpen\"",
        "\"Run\"",
    ] {
        assert!(
            stdout.contains(needle),
            "--list --json output missing {needle}; got: {stdout}"
        );
    }
}

#[test]
fn yee_validate_help_lists_target() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["validate", "--help"])
        .output()
        .expect("invoke yee");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.to_lowercase().contains("mom"));
    assert!(stdout.to_lowercase().contains("fdtd"));
    assert!(stdout.to_lowercase().contains("fem"));
    assert!(stdout.to_lowercase().contains("all"));
}

/// `yee validate fem --json` must exit 0 and emit a JSON array that
/// contains a `fem-eig-001` entry. The FEM drivers that run (001/002/
/// 004/005) complete in seconds in `--release`; 003 and 006 are
/// registered Skipped (wall-time / open gate) so no `Failed` case is
/// present and the command exits 0.
///
/// Gated `#[ignore]` because the FEM sparse-LU solves are slow in
/// debug mode. Run with:
/// `cargo test -p yee-cli --release -- --include-ignored`.
#[test]
#[ignore = "slow in debug: FEM solves (fem-eig-001 ~7 s release, longer debug); run with --release --include-ignored"]
fn yee_validate_fem_json_exits_0_and_contains_fem_eig_001() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["validate", "fem", "--json"])
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee validate fem --json exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("\"fem-eig-001\""),
        "JSON output does not contain fem-eig-001 entry; got: {stdout}"
    );
}
