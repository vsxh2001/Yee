//! Gate `cli-finite-q-s2p` for `yee filter synth --q-unloaded` (Filter Phase
//! F2-Q, ADR-0161).
//!
//! `--q-unloaded <Q>` makes the exported Touchstone `.s2p` carry the realistic
//! **finite-Q lumped-LC** response (`ladder_s21_lossy`) instead of the ideal
//! lossless one. This gate runs the CLI against the standard 3-pole 0.5 dB
//! Chebyshev bandpass (f0 = 2 GHz, FBW = 0.10, Z0 = 50), **reads the written
//! `.s2p` back** via `yee_io::touchstone::read`, and asserts:
//!
//!   1. the midband (`f0`) insertion loss `−20·log10|S21|` matches **Cohn's**
//!      dissipation-loss formula `4.343·Σg/(Q_u·FBW)` (≈ 1.86 dB) to ≤ 15 %.
//!      This is non-circular: Cohn is computed independently from the prototype
//!      g-sum (`yee_synth::prototype`, g[1..=3]) / Q_u / FBW, while the `.s2p`
//!      value comes from the CLI's own finite-Q sweep + the Touchstone
//!      round-trip.
//!   2. the finite-Q `.s2p` **byte-differs** from the default (ideal) `.s2p`
//!      for the same spec — proving the `--q-unloaded` branch is taken (the
//!      ADR-0158 byte-diff pattern), AND the default file's midband IL ≈ 0 dB
//!      (the lossless response).
//!   3. a non-positive / non-finite `--q-unloaded` is rejected (non-zero exit).
//!
//! Pure closed-form (no EM/FDTD), so NOT `#[ignore]`'d.

use std::path::{Path, PathBuf};
use std::process::Command;

use num_complex::Complex64;
use yee_synth::{Approximation, prototype};

/// Standard 3-pole 0.5 dB Chebyshev bandpass shared by every case here.
const F0_HZ: f64 = 2.0e9;
const FBW: f64 = 0.10;
const Z0_OHM: f64 = 50.0;
const ORDER: usize = 3;
const RIPPLE_DB: f64 = 0.5;
/// Per-resonator unloaded Q exercised by the finite-Q case.
const Q_UNLOADED: f64 = 100.0;

/// Write the standard 3-pole spec TOML into the cargo target tmp dir and return
/// its path. Self-contained (not the committed N=5 fixture) so the Cohn Σg —
/// computed below over g[1..=3] — is the like-for-like reference for THIS spec.
fn write_spec() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_finite_q_spec.toml");
    let toml = format!(
        "response = \"Bandpass\"\n\
         f0_hz = {F0_HZ:e}\n\
         fbw = {FBW}\n\
         order = {ORDER}\n\
         z0_ohm = {Z0_OHM}\n\
         \n\
         [approximation.Chebyshev]\n\
         ripple_db = {RIPPLE_DB}\n\
         \n\
         [mask]\n\
         passband_ripple_db = {RIPPLE_DB}\n\
         return_loss_db = 9.0\n\
         stopband = [[2.4e9, 20.0]]\n"
    );
    std::fs::write(&path, toml).expect("write 3-pole spec TOML");
    path
}

/// Run `yee filter synth <spec> --output <out> [--q-unloaded <q>]` and return
/// `(success, stdout)`. Does not assert the exit status (so the rejection case
/// can inspect a failure).
fn run_synth(out: &Path, q_unloaded: Option<f64>) -> (bool, String) {
    let spec = write_spec();
    let _ = std::fs::remove_file(out);
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_yee"));
    cmd.args(["filter", "synth"])
        .arg(&spec)
        .arg("--output")
        .arg(out);
    if let Some(q) = q_unloaded {
        cmd.args(["--q-unloaded", &q.to_string()]);
    }
    let output = cmd.output().expect("invoke yee");
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

/// Read the `.s2p` back and return the midband insertion loss `−20·log10|S21|`
/// at the swept frequency nearest `f0`. S21 of a 2-port row-major S-matrix is
/// index `[1*2 + 0] = 2`.
fn midband_il_db(s2p: &Path) -> f64 {
    let file = yee_io::touchstone::read(s2p).expect("read written .s2p");
    assert_eq!(file.n_ports, 2, ".s2p must be a 2-port file");
    assert!(!file.freq_hz.is_empty(), ".s2p has no frequency points");
    let idx = file
        .freq_hz
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (*a - F0_HZ).abs().partial_cmp(&(*b - F0_HZ).abs()).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let s21: Complex64 = file.data[idx][2];
    -20.0 * s21.norm().max(1e-12).log10()
}

/// Cohn dissipation-loss reference `4.343·Σg/(Q_u·FBW)` (dB), with Σg the sum
/// of the reactive prototype elements g[1..=N] — computed independently of the
/// CLI from `yee_synth::prototype` so the gate is non-circular.
fn cohn_il_db() -> f64 {
    let proto = prototype(
        Approximation::Chebyshev {
            ripple_db: RIPPLE_DB,
        },
        ORDER,
    );
    let sum_g: f64 = proto.g[1..=ORDER].iter().sum();
    4.343 * sum_g / (Q_UNLOADED * FBW)
}

#[test]
fn cli_finite_q_s2p() {
    // ---- (1) finite-Q .s2p: midband IL ≈ Cohn ---------------------------
    let out_q = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_finite_q_out_q.s2p");
    let (ok_q, stdout_q) = run_synth(&out_q, Some(Q_UNLOADED));
    assert!(ok_q, "finite-Q synth exited non-zero; stdout:\n{stdout_q}");
    assert!(
        out_q.exists(),
        "finite-Q Touchstone {} was not written",
        out_q.display()
    );
    assert!(
        stdout_q.contains("finite-Q response"),
        "stdout missing the finite-Q response notice; got:\n{stdout_q}"
    );

    let il_q = midband_il_db(&out_q);
    let il_cohn = cohn_il_db();

    // ---- (2) ideal .s2p: byte-differs + IL ≈ 0 dB -----------------------
    let out_ideal = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_finite_q_out_ideal.s2p");
    let (ok_i, stdout_i) = run_synth(&out_ideal, None);
    assert!(ok_i, "ideal synth exited non-zero; stdout:\n{stdout_i}");
    let il_ideal = midband_il_db(&out_ideal);

    let bytes_q = std::fs::read(&out_q).expect("read finite-Q .s2p bytes");
    let bytes_ideal = std::fs::read(&out_ideal).expect("read ideal .s2p bytes");

    println!(
        "cli-finite-q-s2p: IL_q = {il_q:.4} dB | IL_cohn = {il_cohn:.4} dB \
         (rel err {:.2} %) | IL_ideal = {il_ideal:.4} dB | byte-differ = {}",
        100.0 * (il_q - il_cohn).abs() / il_cohn,
        bytes_q != bytes_ideal,
    );

    // Assertion (1): finite-Q midband IL matches Cohn within 15 %.
    let rel_err = (il_q - il_cohn).abs() / il_cohn;
    assert!(
        rel_err <= 0.15,
        "finite-Q midband IL {il_q:.4} dB disagrees with Cohn {il_cohn:.4} dB by \
         {:.2} % (> 15 %)",
        100.0 * rel_err
    );

    // Assertion (2a): the finite-Q file byte-differs from the ideal file —
    // proves the `--q-unloaded` branch is actually taken (ADR-0158 pattern).
    assert_ne!(
        bytes_q, bytes_ideal,
        "finite-Q .s2p is byte-identical to the ideal .s2p — the --q-unloaded \
         branch was not taken (IL_q={il_q:.4} dB, IL_ideal={il_ideal:.4} dB)"
    );
    // Assertion (2b): the default (ideal) response is ~lossless at midband.
    assert!(
        il_ideal <= 0.2,
        "ideal midband IL {il_ideal:.4} dB is not ~0 dB (lossless response expected)"
    );
    // The finite-Q loss must be materially above the lossless floor (sanity).
    assert!(
        il_q > il_ideal + 0.5,
        "finite-Q IL {il_q:.4} dB is not materially above ideal {il_ideal:.4} dB"
    );

    // ---- (3) non-positive / non-finite Q is rejected --------------------
    let out_bad = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("cli_finite_q_out_bad.s2p");
    for bad_q in [-1.0_f64, 0.0_f64] {
        let (ok_bad, _) = run_synth(&out_bad, Some(bad_q));
        assert!(
            !ok_bad,
            "--q-unloaded {bad_q} should be rejected with a non-zero exit, but synth succeeded"
        );
    }

    let _ = std::fs::remove_file(&out_q);
    let _ = std::fs::remove_file(&out_ideal);
    let _ = std::fs::remove_file(&out_bad);
}
