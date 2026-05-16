//! Integration tests for the `yee` binary.
//!
//! Every test asserts both an exit code AND meaningful stdout/stderr content.

use assert_cmd::Command;
use predicates::str::contains;

/// `yee --help` exits 0 and lists every subcommand.
#[test]
fn help_lists_every_subcommand() {
    Command::cargo_bin("yee")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(contains("validate"))
        .stdout(contains("mesh"))
        .stdout(contains("export"))
        .stdout(contains("run"));
}

/// `yee --version` exits 0 and prints the workspace version.
#[test]
fn version_prints_workspace_version() {
    Command::cargo_bin("yee")
        .unwrap()
        .arg("--version")
        .assert()
        .success()
        .stdout(contains("0.0.0"));
}

/// `yee validate mom` exits 0 and prints the three planned mom cases.
#[test]
fn validate_mom_prints_planned_cases() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["validate", "mom"])
        .assert()
        .success()
        .stdout(contains("yee validate mom"))
        .stdout(contains("planned cases"))
        .stdout(contains("mom-001"))
        .stdout(contains("mom-002"))
        .stdout(contains("mom-003"));
}

/// `yee validate fdtd` exits 0 and prints the Phase 2 deliverable notice.
#[test]
fn validate_fdtd_prints_phase2_notice() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["validate", "fdtd"])
        .assert()
        .success()
        .stdout(contains("Phase 2 deliverable"))
        .stdout(contains("yee-fdtd"));
}

/// `yee validate all` exits 0, runs both, and stdout includes the Phase 2 line.
#[test]
fn validate_all_runs_both() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["validate", "all"])
        .assert()
        .success()
        .stdout(contains("mom-001"))
        .stdout(contains("Phase 2 deliverable"));
}

/// `yee mesh` without the `gmsh` feature exits with code 2 and mentions gmsh.
#[test]
fn mesh_without_gmsh_feature_exits_2() {
    let output = Command::cargo_bin("yee")
        .unwrap()
        .args(["mesh", "/tmp/nonexistent.step"])
        .assert()
        .failure()
        .code(2)
        .get_output()
        .clone();
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    let combined = format!("{stdout}{stderr}");
    assert!(
        combined.contains("gmsh"),
        "expected output to mention gmsh, got: stdout={stdout:?} stderr={stderr:?}"
    );
}

/// Unknown subcommand exits non-zero with a clap error/suggestion.
#[test]
fn unknown_subcommand_suggests() {
    let output = Command::cargo_bin("yee")
        .unwrap()
        .arg("garbage-subcmd")
        .assert()
        .failure()
        .get_output()
        .clone();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unrecognized") || stderr.contains("similar") || stderr.contains("error"),
        "expected an error/suggestion from clap, got: {stderr}"
    );
}
