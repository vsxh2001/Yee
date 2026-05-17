//! Integration tests for `yee fdtd-run`.
//!
//! The fast `--help` smoke test asserts every documented flag is wired
//! through clap; the slow end-to-end test (ignored by default) drives the
//! real [`yee_fdtd::FdtdDriver`] on a 60³ grid for a few hundred steps and
//! checks the JSON it writes has the expected θ-sweep shape.

use std::process::Command;

#[test]
fn fdtd_run_help_lists_args() {
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["fdtd-run", "--help"])
        .output()
        .expect("invoke yee");
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let lowered = stdout.to_lowercase();
    for arg in &[
        "--grid",
        "--dx",
        "--steps",
        "--source",
        "--dipole-length",
        "--freq",
        "--ntff-pad",
        "--cpml",
        "--output",
    ] {
        assert!(lowered.contains(arg), "help missing {arg}");
    }
}

#[test]
#[ignore = "slow: ~5s for 60³ × 800 steps in release"]
fn fdtd_run_emits_radiation_pattern_json() {
    let tmp = std::env::temp_dir();
    let out = tmp.join(format!("yee-cli-fdtd-{}.json", std::process::id()));
    let status = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args([
            "fdtd-run",
            "--steps",
            "400",
            "--output",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("invoke yee");
    assert!(status.success());
    let text = std::fs::read_to_string(&out).unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert!(v["theta_deg"].as_array().unwrap().len() >= 19);
    assert!(v["e_theta_phi0"].as_array().unwrap().len() >= 19);
    std::fs::remove_file(out).ok();
}
