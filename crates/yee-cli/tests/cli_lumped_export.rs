//! Gate for `yee filter synth --lumped` (Filter Phase F2.2-cli, ADR-0158).
//!
//! Runs the CLI against the committed satisfiable Chebyshev N=5 bandpass fixture
//! with `--lumped --footprint 0603 --gerber <tmp>` and asserts the written
//! Gerber is a structurally-valid RS-274X file (`%FS` / `G36` / `G37` / `M02`)
//! whose filled-region count covers every L/C component pad plus the signal
//! line + ground rail, AND that it **differs** from the planar (`--lumped`
//! omitted) Gerber for the same spec — proving the `--lumped` branch is taken
//! rather than silently falling through to the distributed edge-coupled layout.
//!
//! Pure closed-form geometry (no EM/FDTD), so NOT `#[ignore]`'d.

use std::path::PathBuf;
use std::process::Command;

/// Absolute path to the committed Chebyshev N=5 bandpass fixture (shared with
/// the planar `cli_gerber` gate).
fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cheb_bpf.toml")
}

/// Run `yee filter synth <fixture> --output <s2p> --gerber <gbr> [--lumped
/// --footprint 0603]` and return the written Gerber text. Asserts a zero exit.
fn synth_gerber(tag: &str, lumped: bool) -> String {
    let s2p = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("cli_lumped_{tag}.s2p"));
    let gbr = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("cli_lumped_{tag}.gbr"));
    let _ = std::fs::remove_file(&s2p);
    let _ = std::fs::remove_file(&gbr);

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_yee"));
    cmd.args(["filter", "synth"])
        .arg(fixture())
        .arg("--output")
        .arg(&s2p)
        .arg("--gerber")
        .arg(&gbr);
    if lumped {
        cmd.arg("--lumped").args(["--footprint", "0603"]);
    }
    let output = cmd.output().expect("invoke yee");

    assert!(
        output.status.success(),
        "yee filter synth (lumped={lumped}) exited non-zero; stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(
        gbr.exists(),
        "layout Gerber {} was not written",
        gbr.display()
    );
    let gerber = std::fs::read_to_string(&gbr).expect("read layout Gerber");
    let _ = std::fs::remove_file(&s2p);
    let _ = std::fs::remove_file(&gbr);
    gerber
}

#[test]
fn cli_lumped_gerber() {
    let lumped = synth_gerber("lumped", true);

    // Structurally-valid RS-274X: coordinate-format header, region open/close,
    // and the end-of-file marker.
    assert!(!lumped.is_empty(), "lumped Gerber is empty");
    assert!(
        lumped.contains("%FSLAX46Y46*%"),
        "lumped Gerber missing %FS coordinate-format header; got:\n{}",
        &lumped[..lumped.len().min(200)]
    );
    assert!(
        lumped.contains("G36*"),
        "lumped Gerber missing G36 region-open"
    );
    assert!(
        lumped.contains("G37*"),
        "lumped Gerber missing G37 region-close"
    );
    assert!(
        lumped.contains("M02*"),
        "lumped Gerber missing M02 end-of-file marker"
    );

    // The N=5 fixture synthesizes 5 LC resonators → 2·5 = 10 components (an L
    // and a C per resonator), each a two-pad SMD footprint → ≥ 2 copper regions
    // per component, plus the signal line and ground rail. So the filled-region
    // (`G36*`) count must be ≥ 2·N_components = 2·(2·resonators) = 20 for the
    // pads alone, with the line + rail on top.
    const RESONATORS: usize = 5; // order-5 Chebyshev fixture
    const N_COMPONENTS: usize = 2 * RESONATORS; // an L + a C per resonator
    let n_regions = lumped.matches("G36*").count();
    assert!(
        n_regions >= 2 * N_COMPONENTS,
        "lumped Gerber has {n_regions} filled regions; expected ≥ {} \
         (2 pads × {N_COMPONENTS} components) plus signal line + ground rail",
        2 * N_COMPONENTS
    );

    // The `--lumped` branch must actually be taken: the lumped board Gerber must
    // differ from the planar edge-coupled Gerber for the SAME spec/substrate. A
    // silent fall-through to the distributed layout would make these identical.
    let planar = synth_gerber("planar", false);
    let n_planar = planar.matches("G36*").count();
    assert_ne!(
        lumped, planar,
        "lumped Gerber is byte-identical to the planar Gerber — the --lumped \
         branch was not taken (lumped regions={n_regions}, planar regions={n_planar})"
    );
}
