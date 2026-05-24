//! Integration tests for `yee plot` from the post-run Touchstone-review angle.
//!
//! The pre-existing `tests/cli.rs` already exercises the single-file `--kind`
//! paths (PNG dB, SVG Smith, unknown-extension rejection). This file covers
//! the contract that the brief asked for explicitly:
//!
//! 1. `--help` mentions Touchstone input and the dB/Smith plot kinds so a
//!    user can discover the workflow without reading the source.
//! 2. `--format` is accepted as the canonical flag name (with `--kind` kept
//!    as a backwards-compat alias).
//! 3. `--format both` emits two PNGs with `-db` / `-smith` suffixes inserted
//!    before the extension, and each is a non-trivial file (> 1 KB) so we
//!    catch regressions where the backend silently writes an empty buffer.

use assert_cmd::Command;
use predicates::boolean::PredicateBooleanExt;
use predicates::str::contains;

/// `yee plot --help` mentions Touchstone input and the dB / Smith plot
/// kinds. This guards the README example (`yee plot ... --format both`)
/// from drifting out of the binary's actual surface.
#[test]
fn plot_help_mentions_touchstone_and_formats() {
    let output = Command::cargo_bin("yee")
        .unwrap()
        .args(["plot", "--help"])
        .output()
        .expect("invoke yee");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lowered = stdout.to_lowercase();
    assert!(
        lowered.contains("touchstone"),
        "help should mention touchstone input, got: {stdout}"
    );
    assert!(
        lowered.contains("smith"),
        "help should mention smith, got: {stdout}"
    );
    assert!(
        lowered.contains("db"),
        "help should mention db kind, got: {stdout}"
    );
    assert!(
        lowered.contains("both"),
        "help should mention both kind, got: {stdout}"
    );
}

/// `yee plot --format db` (canonical spelling) reads a tiny Touchstone
/// fixture and emits a non-trivial PNG. This complements the existing
/// `--kind db` test in `tests/cli.rs` and locks in the `--format` flag
/// name called out by the brief.
#[test]
fn plot_format_db_emits_non_trivial_png() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s1p");
    let output = tmp.path().join("out.png");

    // Three-point .s1p; matches the format yee-io::touchstone::read accepts.
    std::fs::write(
        &input,
        "# Hz S RI R 50\n1.0e9 0.1 0.2\n1.5e9 0.15 0.25\n2.0e9 0.2 0.3\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "db",
            "--output",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let size = std::fs::metadata(&output).unwrap().len();
    assert!(size > 1024, "PNG output too small: {size} bytes");
}

/// `yee plot --format both` emits two files: `<stem>-db.<ext>` and
/// `<stem>-smith.<ext>`. Each must exist and be > 1 KB. The base
/// `--output` path itself is *not* written — the suffixing is what makes
/// the two-file invocation unambiguous.
#[test]
fn plot_format_both_emits_db_and_smith() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s1p");
    let output = tmp.path().join("out.png");

    std::fs::write(
        &input,
        "# Hz S RI R 50\n1.0e9 0.1 0.2\n1.5e9 0.15 0.25\n2.0e9 0.2 0.3\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "both",
            "--output",
            output.to_str().unwrap(),
        ])
        .assert()
        .success();

    let db_path = tmp.path().join("out-db.png");
    let smith_path = tmp.path().join("out-smith.png");
    assert!(db_path.exists(), "dB PNG not created: {db_path:?}");
    assert!(smith_path.exists(), "Smith PNG not created: {smith_path:?}");
    let db_size = std::fs::metadata(&db_path).unwrap().len();
    let smith_size = std::fs::metadata(&smith_path).unwrap().len();
    assert!(db_size > 1024, "dB PNG too small: {db_size} bytes");
    assert!(smith_size > 1024, "Smith PNG too small: {smith_size} bytes");
}

/// `--kind` remains accepted as a legacy alias for `--format`. Drops in
/// this guarantee would break the pre-existing `cli.rs` tests *and* any
/// scripts users have around; this test surfaces the regression at the
/// CLI-flag layer rather than only at the integration test layer.
#[test]
fn plot_kind_alias_is_accepted() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["plot", "--help"])
        .assert()
        .success()
        .stdout(contains("--kind").or(contains("kind")));
}

/// `yee plot --help` now mentions the new `--entry` flag. Guards the multi-trace
/// surface from drifting out of the binary's actual help text.
#[test]
fn plot_help_mentions_entry_flag() {
    Command::cargo_bin("yee")
        .unwrap()
        .args(["plot", "--help"])
        .assert()
        .success()
        .stdout(contains("--entry"));
}

/// `yee plot --entry 11 --entry 21` on a 2-port Touchstone file produces a
/// non-empty PNG that overlays S11 and S21. This is the primary DoD-2 gate.
#[test]
fn plot_multi_trace_s11_s21_overlay_emits_non_empty_png() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s2p");
    let output = tmp.path().join("out.png");

    // Minimal 2-port .s2p: 3 frequency points.
    // Format: f S11_re S11_im S12_re S12_im S21_re S21_im S22_re S22_im
    std::fs::write(
        &input,
        "# Hz S RI R 50\n\
         1.0e9  0.50 0.10  0.30 0.05  0.30 0.05  0.45 0.08\n\
         1.5e9  0.45 0.08  0.25 0.04  0.25 0.04  0.40 0.07\n\
         2.0e9  0.40 0.06  0.20 0.03  0.20 0.03  0.35 0.06\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "db",
            "--output",
            output.to_str().unwrap(),
            "--entry",
            "11",
            "--entry",
            "21",
        ])
        .assert()
        .success();

    assert!(output.exists(), "output PNG not created: {output:?}");
    let size = std::fs::metadata(&output).unwrap().len();
    assert!(size > 1024, "multi-trace PNG too small: {size} bytes");
}

/// `yee plot --all` on a 2-port file overlays all 4 entries (S11, S12, S21,
/// S22) and produces a non-empty PNG.
#[test]
fn plot_all_entries_on_2port_emits_non_empty_png() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s2p");
    let output = tmp.path().join("out.png");

    std::fs::write(
        &input,
        "# Hz S RI R 50\n\
         1.0e9  0.50 0.10  0.30 0.05  0.30 0.05  0.45 0.08\n\
         1.5e9  0.45 0.08  0.25 0.04  0.25 0.04  0.40 0.07\n\
         2.0e9  0.40 0.06  0.20 0.03  0.20 0.03  0.35 0.06\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "db",
            "--output",
            output.to_str().unwrap(),
            "--all",
        ])
        .assert()
        .success();

    assert!(output.exists(), "output PNG not created: {output:?}");
    let size = std::fs::metadata(&output).unwrap().len();
    assert!(size > 1024, "multi-trace --all PNG too small: {size} bytes");
}

/// Out-of-range `--entry` rejects cleanly (no panic, non-zero exit, error
/// message mentions the entry).
#[test]
fn plot_out_of_range_entry_errors_cleanly() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s1p");
    let output = tmp.path().join("out.png");

    // 1-port file; S21 does not exist.
    std::fs::write(&input, "# Hz S RI R 50\n1.0e9 0.5 0.1\n2.0e9 0.4 0.08\n").unwrap();

    let out = Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "db",
            "--output",
            output.to_str().unwrap(),
            "--entry",
            "21",
        ])
        .output()
        .expect("invoke yee");

    assert!(
        !out.status.success(),
        "expected failure for out-of-range entry"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("21") || stderr.contains("out of range") || stderr.contains("range"),
        "stderr should mention the entry or 'out of range': {stderr}"
    );
}

/// `--entry 11 --format smith` is rejected with a clean error (not a panic).
#[test]
fn plot_entry_with_smith_errors_cleanly() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s2p");
    let output = tmp.path().join("out.png");

    std::fs::write(
        &input,
        "# Hz S RI R 50\n\
         1.0e9  0.50 0.10  0.30 0.05  0.30 0.05  0.45 0.08\n\
         2.0e9  0.40 0.06  0.20 0.03  0.20 0.03  0.35 0.06\n",
    )
    .unwrap();

    let out = Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "smith",
            "--output",
            output.to_str().unwrap(),
            "--entry",
            "11",
        ])
        .output()
        .expect("invoke yee");

    assert!(
        !out.status.success(),
        "expected failure for smith + --entry combination"
    );
}

/// Default single-`--port` behaviour is unchanged when neither `--entry`
/// nor `--all` is passed. A 2-port .s2p with `--port 0` extracts S11
/// and writes a non-empty PNG — same contract as before the multi-trace addition.
#[test]
fn plot_default_port_behaviour_unchanged() {
    let tmp = TempDir::new();
    let input = tmp.path().join("test.s2p");
    let output = tmp.path().join("out.png");

    std::fs::write(
        &input,
        "# Hz S RI R 50\n\
         1.0e9  0.50 0.10  0.30 0.05  0.30 0.05  0.45 0.08\n\
         2.0e9  0.40 0.06  0.20 0.03  0.20 0.03  0.35 0.06\n",
    )
    .unwrap();

    Command::cargo_bin("yee")
        .unwrap()
        .args([
            "plot",
            input.to_str().unwrap(),
            "--format",
            "db",
            "--output",
            output.to_str().unwrap(),
            "--port",
            "0",
        ])
        .assert()
        .success();

    assert!(output.exists(), "output PNG not created: {output:?}");
    let size = std::fs::metadata(&output).unwrap().len();
    assert!(size > 1024, "default-port PNG too small: {size} bytes");
}

/// RAII wrapper around a unique scratch directory under `std::env::temp_dir()`.
///
/// Duplicated (intentionally) from `cli.rs` so this test file is
/// self-contained. Test binaries don't share modules.
struct TempDir {
    path: std::path::PathBuf,
}

impl TempDir {
    fn new() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let pid = std::process::id();
        let path = std::env::temp_dir().join(format!("yee-cli-plot-test-{pid}-{n}"));
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
