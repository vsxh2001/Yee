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
    assert!(stdout.to_lowercase().contains("all"));
}
