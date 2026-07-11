//! Gate `net-001` (FS.6.0, ADR-0212): the 2-port network algebra
//! reproduces textbook identities exactly:
//!
//! 1. S↔T round-trip;
//! 2. thru is the cascade identity, both sides;
//! 3. matched attenuators compose in dB with phases summed;
//! 4. cascade is associative;
//! 5. de-embedding a fixture recovers the DUT;
//! 6. mismatch physics: two series impedances cascade to their sum
//!    (the ABCD identity [[1,Z1],[0,1]]·[[1,Z2],[0,1]] = [[1,Z1+Z2],[0,1]]
//!    expressed in S-parameters — reflections are non-zero throughout);
//! 7. `cascade_files` happy path + explicit rejections (ports, z0, grid).
//!
//! All closed-form, instant, non-ignored.

use num_complex::Complex64;
use yee_io::network::TwoPort;
use yee_io::touchstone::{File, Format, FreqUnit};
use yee_io::{
    Error, cascade, cascade_files, deembed_left, deembed_right, renormalize, renormalize_file,
    s_to_t, t_to_s,
};

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn assert_close(a: &TwoPort, b: &TwoPort, tol: f64, what: &str) {
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        assert!(
            (x - y).norm() < tol,
            "{what}: element {i}: {x} vs {y} (|diff| {})",
            (x - y).norm()
        );
    }
}

/// A deliberately non-symmetric, non-reciprocal, lossy 2-port.
fn messy() -> TwoPort {
    [
        c(0.31, -0.12),
        c(0.55, 0.20),
        c(0.72, -0.31),
        c(-0.18, 0.09),
    ]
}

fn thru() -> TwoPort {
    [c(0.0, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)]
}

/// Matched attenuator: |s21| = 10^(−dB/20) with phase θ, s11 = s22 = 0.
fn attenuator(db: f64, theta: f64) -> TwoPort {
    let a = Complex64::from_polar(10f64.powf(-db / 20.0), -theta);
    [c(0.0, 0.0), a, a, c(0.0, 0.0)]
}

/// Series impedance Z in a z0 system: s11 = Z/(Z+2z0), s21 = 2z0/(Z+2z0).
fn series(z: Complex64, z0: f64) -> TwoPort {
    let d = z + 2.0 * z0;
    [z / d, 2.0 * z0 / d, 2.0 * z0 / d, z / d]
}

#[test]
fn net_001_s_t_roundtrip_is_exact() {
    let s = messy();
    let back = t_to_s(&s_to_t(&s).unwrap()).unwrap();
    assert_close(&s, &back, 1e-15, "S -> T -> S");
    // Isolation has no chain representation — explicit rejection.
    let iso = [c(0.9, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.9, 0.0)];
    assert!(matches!(s_to_t(&iso), Err(Error::Network(_))));
}

#[test]
fn net_001_thru_is_cascade_identity() {
    let x = messy();
    assert_close(&cascade(&thru(), &x).unwrap(), &x, 1e-14, "thru . X");
    assert_close(&cascade(&x, &thru()).unwrap(), &x, 1e-14, "X . thru");
}

#[test]
fn net_001_attenuators_compose_in_db_and_phase() {
    let out = cascade(&attenuator(3.0, 0.4), &attenuator(3.0, 0.7)).unwrap();
    let s21 = out[2];
    let db = -20.0 * s21.norm().log10();
    assert!((db - 6.0).abs() < 1e-12, "3 dB + 3 dB = {db} dB");
    assert!(
        (-s21.arg() - 1.1).abs() < 1e-12,
        "phases must sum: {}",
        -s21.arg()
    );
    assert!(
        out[0].norm() < 1e-15 && out[3].norm() < 1e-15,
        "stays matched"
    );
}

#[test]
fn net_001_cascade_is_associative() {
    let (a, b) = (messy(), attenuator(2.0, 0.3));
    let d = series(c(30.0, 45.0), 50.0);
    let left = cascade(&cascade(&a, &b).unwrap(), &d).unwrap();
    let right = cascade(&a, &cascade(&b, &d).unwrap()).unwrap();
    assert_close(&left, &right, 1e-12, "(A.B).D vs A.(B.D)");
}

#[test]
fn net_001_deembed_recovers_dut() {
    let (fixture, dut) = (series(c(12.0, -20.0), 50.0), messy());
    let measured = cascade(&fixture, &dut).unwrap();
    let recovered = deembed_left(&fixture, &measured).unwrap();
    assert_close(&recovered, &dut, 1e-12, "deembed(F, F.D)");
}

#[test]
fn net_001_series_impedances_cascade_to_their_sum() {
    let z0 = 50.0;
    let (z1, z2) = (c(20.0, 35.0), c(5.0, -60.0));
    let out = cascade(&series(z1, z0), &series(z2, z0)).unwrap();
    assert_close(&out, &series(z1 + z2, z0), 1e-13, "series(Z1).series(Z2)");
}

#[test]
fn net_002_deembed_right_recovers_dut() {
    let (fixture, dut) = (series(c(12.0, -20.0), 50.0), messy());
    let measured = cascade(&dut, &fixture).unwrap();
    let recovered = deembed_right(&measured, &fixture).unwrap();
    assert_close(&recovered, &dut, 1e-12, "deembed_right(D.F, F)");
}

#[test]
fn net_002_renormalize_matches_closed_form() {
    // Same-impedance renormalization is the identity, bit-exact.
    let s = messy();
    assert_eq!(renormalize(&s, 50.0, 50.0).unwrap(), s);

    // A series impedance has a closed form at ANY reference impedance:
    // renormalizing the 50-ohm construction to 75 ohm must equal the
    // direct 75-ohm construction.
    let z = c(20.0, 35.0);
    let re75 = renormalize(&series(z, 50.0), 50.0, 75.0).unwrap();
    assert_close(&re75, &series(z, 75.0), 1e-13, "series Z at 75 ohm");

    // Round-trip 50 -> 75 -> 50 is the identity.
    let back = renormalize(&renormalize(&s, 50.0, 75.0).unwrap(), 75.0, 50.0).unwrap();
    assert_close(&back, &s, 1e-14, "50 -> 75 -> 50 round-trip");

    // Non-positive impedance is an explicit rejection.
    assert!(matches!(renormalize(&s, 50.0, 0.0), Err(Error::Network(_))));
}

#[test]
fn net_002_renormalize_file_unblocks_strict_cascade() {
    let grid: Vec<f64> = (1..=5).map(|k| k as f64 * 1.0e9).collect();
    let a = two_port_file(0.9, 50.0, grid.clone());
    let b75 = two_port_file(0.8, 75.0, grid);
    // Strict path rejects the z0 mismatch...
    assert!(matches!(cascade_files(&a, &b75), Err(Error::Network(_))));
    // ...explicit renormalization unblocks it.
    let b50 = renormalize_file(&b75, 50.0).unwrap();
    assert_eq!(b50.z0, 50.0);
    let out = cascade_files(&a, &b50).unwrap();
    assert_eq!(out.freq_hz, a.freq_hz);
    // And the renormalization itself round-trips through the File layer.
    let b75_again = renormalize_file(&b50, 75.0).unwrap();
    for (sa, sb) in b75.data.iter().zip(&b75_again.data) {
        assert_close(
            &[sa[0], sa[1], sa[2], sa[3]],
            &[sb[0], sb[1], sb[2], sb[3]],
            1e-14,
            "file renormalize round-trip",
        );
    }
}

fn two_port_file(scale: f64, z0: f64, freq_hz: Vec<f64>) -> File {
    let data = freq_hz
        .iter()
        .map(|&f| {
            let phase = -2.0e-10 * f;
            let s21 = Complex64::from_polar(scale, phase);
            let s11 = c(0.1 * scale, 0.02);
            vec![s11, s21, s21, s11]
        })
        .collect();
    File {
        n_ports: 2,
        z0,
        freq_unit: FreqUnit::GHz,
        format: Format::RealImag,
        freq_hz,
        data,
        comments: vec![],
    }
}

#[test]
fn net_001_cascade_files_happy_and_rejections() {
    let grid: Vec<f64> = (1..=5).map(|k| k as f64 * 1.0e9).collect();
    let a = two_port_file(0.9, 50.0, grid.clone());
    let b = two_port_file(0.8, 50.0, grid.clone());
    let out = cascade_files(&a, &b).unwrap();
    assert_eq!(out.freq_hz, a.freq_hz);
    for (k, s) in out.data.iter().enumerate() {
        let ta: TwoPort = [a.data[k][0], a.data[k][1], a.data[k][2], a.data[k][3]];
        let tb: TwoPort = [b.data[k][0], b.data[k][1], b.data[k][2], b.data[k][3]];
        let expect = cascade(&ta, &tb).unwrap();
        assert_close(&[s[0], s[1], s[2], s[3]], &expect, 1e-14, "file point");
    }

    // Rejections: port count, z0, frequency grid.
    let mut one_port = a.clone();
    one_port.n_ports = 1;
    assert!(matches!(
        cascade_files(&one_port, &b),
        Err(Error::Network(_))
    ));
    let z75 = two_port_file(0.8, 75.0, grid.clone());
    assert!(matches!(cascade_files(&a, &z75), Err(Error::Network(_))));
    let shifted = two_port_file(0.8, 50.0, grid.iter().map(|f| f + 1.0e6).collect());
    assert!(matches!(
        cascade_files(&a, &shifted),
        Err(Error::Network(_))
    ));
}
