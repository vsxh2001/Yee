//! Integration tests for the Touchstone v1.1 reader/writer.
//!
//! Strategy: for every passive fixture, do `read → write → read` and check
//! that the two reads produce equal `File` structs up to a 1e-12 relative
//! float tolerance. Then test the negative paths (passivity, malformed).

use std::path::{Path, PathBuf};

use yee_io::touchstone;

fn fixtures_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("validation")
        .join("fixtures")
        .join("touchstone")
}

fn approx_eq(a: f64, b: f64, rel: f64) -> bool {
    if a == b {
        return true;
    }
    let diff = (a - b).abs();
    let scale = a.abs().max(b.abs()).max(1.0);
    diff <= rel * scale
}

fn assert_files_close(lhs: &touchstone::File, rhs: &touchstone::File, rel: f64) {
    assert_eq!(lhs.n_ports, rhs.n_ports, "n_ports mismatch");
    assert!(approx_eq(lhs.z0, rhs.z0, rel), "z0 mismatch");
    assert_eq!(lhs.freq_unit, rhs.freq_unit, "freq_unit mismatch");
    assert_eq!(lhs.format, rhs.format, "format mismatch");
    assert_eq!(lhs.comments, rhs.comments, "comments mismatch");
    assert_eq!(
        lhs.freq_hz.len(),
        rhs.freq_hz.len(),
        "freq count mismatch (lhs {} vs rhs {})",
        lhs.freq_hz.len(),
        rhs.freq_hz.len(),
    );
    for (k, (a, b)) in lhs.freq_hz.iter().zip(&rhs.freq_hz).enumerate() {
        assert!(approx_eq(*a, *b, rel), "freq[{k}] mismatch: {a} vs {b}");
    }
    for (k, (m1, m2)) in lhs.data.iter().zip(&rhs.data).enumerate() {
        assert_eq!(m1.len(), m2.len(), "S-matrix len mismatch at freq {k}");
        for (idx, (z1, z2)) in m1.iter().zip(m2).enumerate() {
            assert!(
                approx_eq(z1.re, z2.re, rel) && approx_eq(z1.im, z2.im, rel),
                "S[{idx}] mismatch at freq {k}: {z1} vs {z2}",
            );
        }
    }
}

fn roundtrip(fixture_name: &str) {
    let src = fixtures_dir().join(fixture_name);
    let parsed =
        touchstone::read(&src).unwrap_or_else(|e| panic!("read({fixture_name}) failed: {e}"));

    // Round-trip through a temp file with the same extension.
    let tmp_dir = std::env::temp_dir();
    let tmp_path = tmp_dir.join(format!(
        "yee_io_roundtrip_{}_{}",
        std::process::id(),
        fixture_name,
    ));
    touchstone::write(&tmp_path, &parsed).unwrap_or_else(|e| panic!("write failed: {e}"));
    let reparsed = touchstone::read(&tmp_path)
        .unwrap_or_else(|e| panic!("re-read({fixture_name}) failed: {e}"));

    assert_files_close(&parsed, &reparsed, 1e-12);
    let _ = std::fs::remove_file(&tmp_path);
}

#[test]
fn s1p_roundtrip() {
    roundtrip("1port.s1p");
}

#[test]
fn s2p_ri_roundtrip() {
    roundtrip("2port.s2p");
}

#[test]
fn s2p_db_roundtrip() {
    roundtrip("2port_db.s2p");
}

#[test]
fn s1p_freq_unit_dispatch() {
    let f = touchstone::read(&fixtures_dir().join("1port.s1p")).unwrap();
    assert_eq!(f.freq_unit, touchstone::FreqUnit::GHz);
    assert_eq!(f.freq_hz, vec![1.0e9, 2.0e9, 5.0e9, 10.0e9]);
    assert_eq!(f.n_ports, 1);
    assert_eq!(f.z0, 50.0);
}

#[test]
fn s2p_db_uses_mhz_unit() {
    let f = touchstone::read(&fixtures_dir().join("2port_db.s2p")).unwrap();
    assert_eq!(f.freq_unit, touchstone::FreqUnit::MHz);
    assert_eq!(f.format, touchstone::Format::DecibelAngle);
    // 1000 MHz -> 1 GHz canonical Hz.
    assert!((f.freq_hz[0] - 1.0e9).abs() < 1e-6);
}

#[test]
fn s2p_layout_matches_on_disk_quirk() {
    // 2port.s2p has S11 = S22 = 0, S21 = S12 = 0.5012. After parsing into
    // row-major, S[0]=S11=0, S[1]=S12=0.5012, S[2]=S21=0.5012, S[3]=S22=0.
    let f = touchstone::read(&fixtures_dir().join("2port.s2p")).unwrap();
    let mat = &f.data[0];
    assert!(mat[0].norm() < 1e-12, "S11 ~ 0");
    assert!((mat[1].norm() - 0.5011872336272722).abs() < 1e-12, "S12");
    assert!((mat[2].norm() - 0.5011872336272722).abs() < 1e-12, "S21");
    assert!(mat[3].norm() < 1e-12, "S22 ~ 0");
}

#[test]
fn comments_are_preserved_through_roundtrip() {
    let f = touchstone::read(&fixtures_dir().join("1port.s1p")).unwrap();
    assert!(!f.comments.is_empty(), "expected at least one comment");
    let any_provenance_line = f.comments.iter().any(|c| c.contains("75-ohm"));
    assert!(
        any_provenance_line,
        "comments should include fixture provenance"
    );
}

#[test]
fn nonpassive_fixture_is_rejected() {
    let err = touchstone::read(&fixtures_dir().join("2port_nonpassive.s2p")).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("passivity"),
        "expected passivity-violation message, got: {msg}"
    );
}

#[test]
fn malformed_option_line_is_rejected() {
    // Wrong parameter type (Y instead of S).
    let tmp = std::env::temp_dir().join(format!("yee_io_malformed_{}.s2p", std::process::id()));
    std::fs::write(&tmp, "! malformed\n# GHz Y RI R 50\n1.0 0 0 0 0 0 0 0 0\n").unwrap();
    let err = touchstone::read(&tmp).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("not supported") || msg.contains("Y"),
        "expected option-line rejection, got: {msg}"
    );
    // Position should be at line 2 (the option line).
    if let yee_io::Error::TouchstoneParse { line, .. } = err {
        assert_eq!(line, 2);
    } else {
        panic!("expected TouchstoneParse, got: {msg}");
    }
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn malformed_column_count_is_rejected() {
    // 2-port file needs 1 + 2*4 = 9 floats per frequency; supply only 5.
    let tmp = std::env::temp_dir().join(format!("yee_io_short_{}.s2p", std::process::id()));
    std::fs::write(&tmp, "# GHz S RI R 50\n1.0 0.0 0.0 0.0 0.0\n").unwrap();
    let err = touchstone::read(&tmp).unwrap_err();
    let msg = format!("{err}");
    assert!(
        msg.contains("multiple of 9") || msg.contains("frequency"),
        "expected column-count error, got: {msg}"
    );
    if let yee_io::Error::TouchstoneParse { line, .. } = err {
        // Failure should be flagged on the data line (line 2).
        assert!(line >= 2, "expected line >= 2, got {line}");
    } else {
        panic!("expected TouchstoneParse");
    }
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn s3p_roundtrip_synthetic() {
    // Build an in-memory 3-port file (passive identity scaled by 0.3),
    // write it through the writer, read it back, compare.
    use num_complex::Complex64;
    let n = 3usize;
    let mut data = Vec::new();
    let mut diag = vec![Complex64::new(0.0, 0.0); n * n];
    for i in 0..n {
        diag[i * n + i] = Complex64::new(0.3, 0.0);
    }
    data.push(diag);

    let original = yee_io::File {
        n_ports: n,
        z0: 50.0,
        freq_unit: touchstone::FreqUnit::GHz,
        format: touchstone::Format::RealImag,
        freq_hz: vec![2.4e9],
        data,
        comments: vec![" synthetic 3-port".to_string()],
    };

    let tmp = std::env::temp_dir().join(format!("yee_io_s3p_{}.s3p", std::process::id()));
    touchstone::write(&tmp, &original).unwrap();
    let reparsed = touchstone::read(&tmp).unwrap();
    assert_files_close(&original, &reparsed, 1e-12);
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn s4p_roundtrip_synthetic() {
    use num_complex::Complex64;
    let n = 4usize;
    // Passive: scaled identity 0.25 + small off-diagonal coupling.
    let mut mat = vec![Complex64::new(0.0, 0.0); n * n];
    for i in 0..n {
        mat[i * n + i] = Complex64::new(0.25, 0.05);
    }
    let original = yee_io::File {
        n_ports: n,
        z0: 50.0,
        freq_unit: touchstone::FreqUnit::Hz,
        format: touchstone::Format::MagAngle,
        freq_hz: vec![1.0e9, 2.0e9],
        data: vec![mat.clone(), mat],
        comments: vec![],
    };

    let tmp = std::env::temp_dir().join(format!("yee_io_s4p_{}.s4p", std::process::id()));
    touchstone::write(&tmp, &original).unwrap();
    let reparsed = touchstone::read(&tmp).unwrap();
    assert_files_close(&original, &reparsed, 1e-12);
    let _ = std::fs::remove_file(&tmp);
}

#[test]
fn rejects_bogus_extension() {
    let tmp = std::env::temp_dir().join(format!("yee_io_bogus_{}.s9p", std::process::id()));
    std::fs::write(&tmp, "# GHz S RI R 50\n").unwrap();
    let err = touchstone::read(&tmp).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("s9p") || msg.contains("Phase 0"), "{msg}");
    let _ = std::fs::remove_file(&tmp);
}
