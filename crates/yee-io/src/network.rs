//! Two-port network algebra (FS.6.0, ADR-0212): S↔T conversion, cascade,
//! and fixture de-embedding — the composition layer every commercial
//! suite has and Touchstone I/O alone cannot provide.
//!
//! ## Convention
//!
//! A 2-port S-matrix is `[s11, s12, s21, s22]` row-major, matching the
//! flattening of [`crate::touchstone::File::data`]. The transfer (chain
//! scattering) matrix `T` is defined by
//!
//! ```text
//! [b1]       [a2]            1  [ -det S   s11 ]
//! [a1]  =  T [b2],    T  =  ---  [                ]
//!                           s21  [  -s22     1  ]
//! ```
//!
//! Derivation: from `b2 = s21·a1 + s22·a2` solve `a1 = (b2 − s22·a2)/s21`;
//! substitute into `b1 = s11·a1 + s12·a2` to get
//! `b1 = (−det S/s21)·a2 + (s11/s21)·b2`. Connecting port 2 of A to
//! port 1 of B identifies `(a2, b2)` of A with `(b1, a1)` of B, so the
//! chain is a plain matrix product: **`T_cascade = T_A · T_B`** — the
//! reason this convention (and not one of its transposes) is used here.
//!
//! A network with `s21 = 0` transmits nothing and has no transfer matrix;
//! conversions reject it with [`Error::Network`] instead of emitting NaN.

use crate::{Error, Result, touchstone::File};
use num_complex::Complex64;

/// 2-port matrix, row-major `[m11, m12, m21, m22]`.
pub type TwoPort = [Complex64; 4];

/// Relative tolerance for "same frequency" in [`cascade_files`] and the
/// exact-zero test for `s21` singularity.
const EPS: f64 = 1e-12;

/// Convert a 2-port S-matrix to its transfer (chain) matrix.
///
/// Rejects `|s21| ≈ 0` (no transmission ⇒ no chain representation).
pub fn s_to_t(s: &TwoPort) -> Result<TwoPort> {
    let [s11, s12, s21, s22] = *s;
    if s21.norm() < EPS {
        return Err(Error::Network(
            "s21 = 0: an isolating network has no transfer matrix".into(),
        ));
    }
    let det = s11 * s22 - s12 * s21;
    Ok([
        -det / s21,
        s11 / s21,
        -s22 / s21,
        Complex64::new(1.0, 0.0) / s21,
    ])
}

/// Convert a transfer (chain) matrix back to its S-matrix.
///
/// Rejects `|t22| ≈ 0` (the image of `s21 → ∞`, not a physical network).
pub fn t_to_s(t: &TwoPort) -> Result<TwoPort> {
    let [t11, t12, t21, t22] = *t;
    if t22.norm() < EPS {
        return Err(Error::Network(
            "t22 = 0: transfer matrix has no S-parameter image".into(),
        ));
    }
    let det = t11 * t22 - t12 * t21;
    Ok([
        t12 / t22,
        det / t22,
        Complex64::new(1.0, 0.0) / t22,
        -t21 / t22,
    ])
}

fn mat_mul(a: &TwoPort, b: &TwoPort) -> TwoPort {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

fn mat_inv(m: &TwoPort, what: &str) -> Result<TwoPort> {
    let det = m[0] * m[3] - m[1] * m[2];
    if det.norm() < EPS {
        return Err(Error::Network(format!("{what} is singular")));
    }
    Ok([m[3] / det, -m[1] / det, -m[2] / det, m[0] / det])
}

/// Cascade two 2-ports (port 2 of `a` into port 1 of `b`, same reference
/// impedance): `S(a·b) = t_to_s(T_a · T_b)`.
pub fn cascade(a: &TwoPort, b: &TwoPort) -> Result<TwoPort> {
    t_to_s(&mat_mul(&s_to_t(a)?, &s_to_t(b)?))
}

/// Remove a known fixture from the **left** (input side) of a measured
/// cascade: given `measured = fixture · dut`, recover `dut` as
/// `t_to_s(T_fixture⁻¹ · T_measured)`.
pub fn deembed_left(fixture: &TwoPort, measured: &TwoPort) -> Result<TwoPort> {
    let t_f = mat_inv(&s_to_t(fixture)?, "fixture transfer matrix")?;
    t_to_s(&mat_mul(&t_f, &s_to_t(measured)?))
}

/// Cascade two 2-port Touchstone [`File`]s frequency-by-frequency.
///
/// FS.6.0 requires **identical** frequency grids (relative tolerance
/// 1e-12; no interpolation) and identical reference impedance — anything
/// else is an explicit [`Error::Network`], never a silent resample.
/// Comments and formatting metadata are taken from `a`.
pub fn cascade_files(a: &File, b: &File) -> Result<File> {
    if a.n_ports != 2 || b.n_ports != 2 {
        return Err(Error::Network(format!(
            "cascade needs two 2-ports, got {}-port and {}-port",
            a.n_ports, b.n_ports
        )));
    }
    if (a.z0 - b.z0).abs() > EPS * a.z0.abs().max(1.0) {
        return Err(Error::Network(format!(
            "reference impedances differ: {} vs {} ohm (renormalization is FS.6.1)",
            a.z0, b.z0
        )));
    }
    if a.freq_hz.len() != b.freq_hz.len()
        || a.freq_hz
            .iter()
            .zip(&b.freq_hz)
            .any(|(fa, fb)| (fa - fb).abs() > EPS * fa.abs().max(1.0))
    {
        return Err(Error::Network(
            "frequency grids differ (FS.6.0 does not interpolate)".into(),
        ));
    }
    let mut data = Vec::with_capacity(a.data.len());
    for (k, (sa, sb)) in a.data.iter().zip(&b.data).enumerate() {
        let ta: TwoPort = [sa[0], sa[1], sa[2], sa[3]];
        let tb: TwoPort = [sb[0], sb[1], sb[2], sb[3]];
        let s = cascade(&ta, &tb)
            .map_err(|e| Error::Network(format!("at {} Hz (point {k}): {e}", a.freq_hz[k])))?;
        data.push(s.to_vec());
    }
    Ok(File {
        n_ports: 2,
        z0: a.z0,
        freq_unit: a.freq_unit,
        format: a.format,
        freq_hz: a.freq_hz.clone(),
        data,
        comments: a.comments.clone(),
    })
}
