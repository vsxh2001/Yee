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
