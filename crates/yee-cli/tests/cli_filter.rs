//! Smoke test for `yee filter synth` (Filter Phase F0).
//!
//! Runs the CLI against the committed satisfiable Chebyshev bandpass fixture
//! and asserts: exit 0, stdout reports the coupling matrix and a `PASS`
//! verdict, and a Touchstone `.s2p` file is written. Fast (pure math, no EM),
//! so NOT `#[ignore]`'d.

use std::path::PathBuf;
use std::process::Command;

/// Absolute path to the committed fixture spec.
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cheb_bpf.toml")
}

#[test]
fn yee_filter_synth_passes_and_writes_touchstone() {
    // Write the Touchstone into the cargo target dir so the test is hermetic
    // and leaves no artifact in `tests/fixtures/`.
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_filter_cheb_bpf.s2p");
    let _ = std::fs::remove_file(&out);

    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(fixture())
        .arg("--output")
        .arg(&out)
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee filter synth exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("coupling matrix"),
        "stdout missing the coupling matrix; got:\n{stdout}"
    );
    assert!(
        stdout.contains("PASS"),
        "stdout missing PASS verdict; got:\n{stdout}"
    );

    assert!(
        out.exists(),
        "Touchstone output {} was not written",
        out.display()
    );
    let s2p = std::fs::read_to_string(&out).expect("read written Touchstone");
    assert!(
        s2p.contains("# Hz S RI R 50"),
        "written file is not a Z0=50 Touchstone option line; got:\n{}",
        &s2p[..s2p.len().min(200)]
    );

    let _ = std::fs::remove_file(&out);
}

#[test]
fn yee_filter_help_lists_synth() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "--help"])
        .output()
        .expect("invoke yee");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.to_lowercase().contains("synth"));
}
