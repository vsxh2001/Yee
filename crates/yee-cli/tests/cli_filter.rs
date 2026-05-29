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
fn yee_filter_synth_plot_writes_png() {
    // `--plot` renders the |S21| response with the spec mask overlaid (no EM).
    let out = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_filter_cheb_bpf.png");
    let s2p = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_filter_plot.s2p");
    let _ = std::fs::remove_file(&out);

    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(fixture())
        .arg("--output")
        .arg(&s2p)
        .arg("--plot")
        .arg(&out)
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee filter synth --plot exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stdout).contains("wrote plot"),
        "stdout missing 'wrote plot'"
    );
    assert!(out.exists(), "plot {} was not written", out.display());
    let len = std::fs::metadata(&out).expect("plot metadata").len();
    assert!(len > 1024, "plot PNG must be non-trivial, got {len} bytes");

    let _ = std::fs::remove_file(&out);
    let _ = std::fs::remove_file(&s2p);
}

/// F1.2.0: `yee filter synth` emits physical edge-coupled dimensions for the
/// committed Chebyshev 0.5 dB N=5 fixture on the default FR-4 substrate, and
/// `--layout-svg` writes a well-formed SVG. Pure closed-form math (no EM/FDTD),
/// so NOT `#[ignore]`'d.
#[test]
fn cli_dims() {
    let s2p = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_dims.s2p");
    let svg = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_dims.svg");
    let _ = std::fs::remove_file(&s2p);
    let _ = std::fs::remove_file(&svg);

    // FR-4 defaults (eps_r=4.4, h=1.6 mm) are supplied by the CLI — exercise
    // them implicitly by omitting --eps-r/--h-mm.
    let output = Command::new(env!("CARGO_BIN_EXE_yee"))
        .args(["filter", "synth"])
        .arg(fixture())
        .arg("--output")
        .arg(&s2p)
        .arg("--layout-svg")
        .arg(&svg)
        .output()
        .expect("invoke yee");

    assert!(
        output.status.success(),
        "yee filter synth --layout-svg exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("physical dimensions"),
        "stdout missing the dimensions block; got:\n{stdout}"
    );

    // Parse the line width / resonator length / gaps back out of the SI metres
    // column (the value right before the `m  (` mm-annotation) and assert all
    // are strictly positive.
    let line_width = parse_si_after(&stdout, "line width");
    let res_len = parse_si_after(&stdout, "resonator length");
    assert!(line_width > 0.0, "line width must be > 0; got {line_width}");
    assert!(res_len > 0.0, "resonator length must be > 0; got {res_len}");

    let n_gaps = stdout.matches("gap[").count();
    assert_eq!(n_gaps, 4, "N=5 filter has N-1=4 inter-resonator gaps");
    for i in 0..n_gaps {
        let gap = parse_si_after(&stdout, &format!("gap[{i}]"));
        assert!(gap > 0.0, "gap[{i}] must be > 0; got {gap}");
    }

    // The layout SVG must be a well-formed, non-empty SVG document.
    assert!(svg.exists(), "layout SVG {} was not written", svg.display());
    let svg_text = std::fs::read_to_string(&svg).expect("read layout SVG");
    assert!(
        svg_text.contains("<svg"),
        "layout SVG missing opening <svg tag; got:\n{}",
        &svg_text[..svg_text.len().min(200)]
    );
    assert!(
        svg_text.contains("</svg>"),
        "layout SVG missing closing </svg> tag"
    );

    let _ = std::fs::remove_file(&s2p);
    let _ = std::fs::remove_file(&svg);
}

/// Pull the first SI-metres value out of a dimensions line whose text starts
/// with `label` and has the shape `... = <value> m  (<mm> mm) ...`.
fn parse_si_after(stdout: &str, label: &str) -> f64 {
    let line = stdout
        .lines()
        .find(|l| l.trim_start().starts_with(label))
        .unwrap_or_else(|| panic!("no line starting with {label:?} in:\n{stdout}"));
    // Take the token immediately preceding the literal `m` units marker.
    let (lhs, _) = line
        .split_once(" m  (")
        .unwrap_or_else(|| panic!("line {line:?} not in expected `= <v> m  (` shape"));
    let value = lhs
        .rsplit(['=', ' '])
        .find(|t| !t.is_empty())
        .unwrap_or_else(|| panic!("no value token in {line:?}"));
    value
        .parse::<f64>()
        .unwrap_or_else(|_| panic!("value {value:?} (from {line:?}) is not f64"))
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
