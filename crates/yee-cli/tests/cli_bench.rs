//! `yee bench` help-smoke integration tests.
//!
//! These tests deliberately only exercise clap's help-rendering path: they
//! confirm that the `Bench` subcommand exists and that its `BenchTarget`
//! value-enum advertises every supported target. Running an actual
//! benchmark would spend minutes per CI invocation and is therefore left
//! to developers / nightly runs.

use std::process::Command;

/// `yee bench --help` exits 0 and lists every `BenchTarget` variant.
///
/// The help output renders value-enum variants in lowercase
/// (clap's default for `ValueEnum`), so we lowercase the captured stdout
/// before substring-matching to stay robust to clap's formatting.
#[test]
fn yee_bench_help_lists_targets() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["bench", "--help"])
        .output()
        .expect("invoke yee");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lower = stdout.to_lowercase();
    for target in &["mom", "fdtd", "gmres", "gp", "bo", "all"] {
        assert!(
            lower.contains(target),
            "expected target '{target}' in `yee bench --help`, got:\n{stdout}"
        );
    }
}
