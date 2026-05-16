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
        .stdout(contains("completions"))
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

/// `yee export --format hdf5` exits with code 2 and mentions hdf5 on stderr.
#[test]
fn export_hdf5_exits_2() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["export", "/tmp/in.s1p", "--format", "hdf5", "/tmp/out.h5"])
        .assert()
        .failure()
        .code(2)
        .stderr(contains("hdf5"));
}

/// `yee export --format touchstone` reads a real Touchstone file and writes a new one.
#[test]
fn export_touchstone_roundtrip() {
    let tmp = TempDir::new();
    let input = tmp.path().join("in.s1p");
    let output = tmp.path().join("out.s1p");
    std::fs::write(
        &input,
        "! example one-port\n# GHz S RI R 50\n1.0 0.5 -0.25\n2.0 0.4 -0.20\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "export",
            input.to_str().unwrap(),
            "--format",
            "touchstone",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let written = std::fs::read_to_string(&output).unwrap();
    assert!(
        written.contains("# GHz S RI"),
        "expected option line in output, got: {written}"
    );
    assert!(written.contains("0.5"), "expected data, got: {written}");
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

/// `yee completions bash` writes a bash completion script to stdout.
#[test]
fn completions_bash_emits_script() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["completions", "bash"])
        .assert()
        .success()
        .stdout(contains("yee"));
}

/// `yee completions zsh` writes a zsh completion script.
#[test]
fn completions_zsh_emits_script() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["completions", "zsh"])
        .assert()
        .success()
        .stdout(contains("yee"));
}

/// `yee completions fish` writes a fish completion script.
#[test]
fn completions_fish_emits_script() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["completions", "fish"])
        .assert()
        .success()
        .stdout(contains("yee"));
}

/// RAII wrapper around a unique scratch directory under `std::env::temp_dir()`.
///
/// The directory and all of its contents are removed when the `TempDir` is
/// dropped (test exit, scope exit, or panic-unwind) so repeated CI runs
/// don't accumulate `yee-cli-test-<pid>-<n>` directories in `$TMPDIR`.
/// `remove_dir_all` errors are ignored on drop — a leak in a failure path is
/// preferable to a destructor that panics inside an unwinding test.
struct TempDir {
    path: std::path::PathBuf,
}

impl TempDir {
    fn new() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let pid = std::process::id();
        let path = std::env::temp_dir().join(format!("yee-cli-test-{pid}-{n}"));
        std::fs::create_dir_all(&path).unwrap();
        Self { path }
    }

    fn path(&self) -> &std::path::Path {
        &self.path
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}
