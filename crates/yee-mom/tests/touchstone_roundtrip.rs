//! Integration test: `SParameters` ↔ `yee_io::touchstone::File` round-trip.
//!
//! Builds a small passive 2-port S-parameter dataset, writes it via
//! `SParameters::write_touchstone`, reads the file back through
//! `yee_io::touchstone::read`, converts back into `SParameters` and asserts
//! struct equality to 1e-12 relative tolerance.

use num_complex::Complex64;
use tempfile::TempDir;
use yee_mom::SParameters;

/// Build a tiny passive 2-port S dataset at three frequencies. All entries
/// have |S_ij| well below 1.0, so passivity (|σ_max(S)| ≤ 1) is satisfied.
fn synthetic_2port() -> SParameters {
    // Row-major per frequency: [S11, S12, S21, S22].
    let f0 = vec![
        Complex64::new(0.10, 0.02),
        Complex64::new(0.30, -0.05),
        Complex64::new(0.30, -0.05),
        Complex64::new(0.12, 0.01),
    ];
    let f1 = vec![
        Complex64::new(0.15, -0.04),
        Complex64::new(0.28, 0.06),
        Complex64::new(0.28, 0.06),
        Complex64::new(0.14, -0.02),
    ];
    let f2 = vec![
        Complex64::new(0.20, 0.03),
        Complex64::new(0.25, -0.07),
        Complex64::new(0.25, -0.07),
        Complex64::new(0.18, 0.04),
    ];
    SParameters {
        freq_hz: vec![1.0e9, 1.5e9, 2.0e9],
        data: vec![f0, f1, f2],
        n_ports: 2,
    }
}

fn relative_close(a: Complex64, b: Complex64, tol: f64) -> bool {
    let diff = (a - b).norm();
    let scale = a.norm().max(b.norm()).max(1.0);
    diff <= tol * scale
}

#[test]
fn write_read_back_roundtrip_2port() {
    let original = synthetic_2port();

    let dir = TempDir::new().expect("tempdir");
    let path = dir.path().join("network.s2p");

    original
        .write_touchstone(&path, 50.0)
        .expect("write_touchstone");

    let file = yee_io::touchstone::read(&path).expect("read");

    // Sanity-check the on-disk metadata before converting back.
    assert_eq!(file.n_ports, original.n_ports);
    assert_eq!(file.freq_hz.len(), original.freq_hz.len());
    assert_eq!(file.z0, 50.0);

    let round_tripped = SParameters::from_touchstone(&file);

    assert_eq!(round_tripped.n_ports, original.n_ports);
    assert_eq!(
        round_tripped.freq_hz.len(),
        original.freq_hz.len(),
        "frequency count must match"
    );
    for (k, (f_round, f_orig)) in round_tripped
        .freq_hz
        .iter()
        .zip(original.freq_hz.iter())
        .enumerate()
    {
        let denom = f_orig.abs().max(1.0);
        let rel = (f_round - f_orig).abs() / denom;
        assert!(rel <= 1e-12, "freq[{k}] mismatch: {f_round} vs {f_orig}");
    }
    assert_eq!(round_tripped.data.len(), original.data.len());
    for (k, (row_round, row_orig)) in round_tripped
        .data
        .iter()
        .zip(original.data.iter())
        .enumerate()
    {
        assert_eq!(row_round.len(), row_orig.len(), "row {k} length");
        for (idx, (a, b)) in row_round.iter().zip(row_orig.iter()).enumerate() {
            assert!(
                relative_close(*a, *b, 1e-12),
                "S-matrix mismatch at freq {k}, slot {idx}: {a} vs {b}",
            );
        }
    }
}
