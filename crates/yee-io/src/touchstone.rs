//! Touchstone v1.1 reader/writer.
//!
//! Implements the option-line grammar `# <freq_unit> <param_type> <format> R <Z0>`
//! with these constraints (Phase 0):
//!
//! - `freq_unit` ∈ {Hz, kHz, MHz, GHz} (case-insensitive)
//! - `param_type` = `S` only; Y/Z/H/G are rejected with a descriptive error
//! - `format` ∈ {RI, MA, DB} (case-insensitive)
//! - `R <Z0>` is optional; default `Z0 = 50.0`
//!
//! Comments begin with `!` and may appear anywhere; they are preserved on
//! [`File`] in their original order so a `read → write` round-trip is faithful.
//!
//! Port count is determined from the file extension (`.s1p` … `.s4p`).
//!
//! S-parameter data layout per spec:
//!
//! - **n = 1**: each frequency line is `f S11.a S11.b`
//! - **n = 2**: each frequency line is `f S11 S21 S12 S22` (note the swapped
//!   off-diagonal order — this is a Touchstone v1 oddity)
//! - **n ≥ 3**: row-major `S11..S1n / S21..S2n / ...`, the first frequency on
//!   the leading row; per spec, four S-parameters per source line maximum,
//!   continued on subsequent lines. We parse permissively (any whitespace
//!   layout) and emit one full frequency record per output line.
//!
//! Passivity (|σ_max(S)| ≤ 1 + 1e-9) is checked on every read.

use crate::{Error, Result};
use num_complex::Complex64;
use std::path::Path;

/// One sample of S-parameter data: an `n_ports × n_ports` matrix stored
/// row-major.
pub type SMatrix = Vec<Complex64>;

/// On-disk numeric format for each complex datum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Real / Imaginary pair.
    RealImag,
    /// Magnitude / phase (degrees).
    MagAngle,
    /// 20·log10(|S|) / phase (degrees).
    DecibelAngle,
}

/// Frequency units recognised in the option line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FreqUnit {
    /// Hz.
    Hz,
    /// kHz.
    KHz,
    /// MHz.
    MHz,
    /// GHz.
    GHz,
}

impl FreqUnit {
    /// Multiplier to convert a number expressed in this unit into Hz.
    pub fn to_hz_multiplier(self) -> f64 {
        match self {
            FreqUnit::Hz => 1.0,
            FreqUnit::KHz => 1.0e3,
            FreqUnit::MHz => 1.0e6,
            FreqUnit::GHz => 1.0e9,
        }
    }

    /// Canonical string for the option-line emitter.
    pub fn as_str(self) -> &'static str {
        match self {
            FreqUnit::Hz => "Hz",
            FreqUnit::KHz => "kHz",
            FreqUnit::MHz => "MHz",
            FreqUnit::GHz => "GHz",
        }
    }

    fn parse(token: &str) -> Option<Self> {
        match token.to_ascii_lowercase().as_str() {
            "hz" => Some(FreqUnit::Hz),
            "khz" => Some(FreqUnit::KHz),
            "mhz" => Some(FreqUnit::MHz),
            "ghz" => Some(FreqUnit::GHz),
            _ => None,
        }
    }
}

impl Format {
    /// Canonical string for the option-line emitter.
    pub fn as_str(self) -> &'static str {
        match self {
            Format::RealImag => "RI",
            Format::MagAngle => "MA",
            Format::DecibelAngle => "DB",
        }
    }

    fn parse(token: &str) -> Option<Self> {
        match token.to_ascii_lowercase().as_str() {
            "ri" => Some(Format::RealImag),
            "ma" => Some(Format::MagAngle),
            "db" => Some(Format::DecibelAngle),
            _ => None,
        }
    }
}

/// Parsed Touchstone file content.
#[derive(Debug, Clone, PartialEq)]
pub struct File {
    /// Number of ports (1..=4 in Phase 0).
    pub n_ports: usize,
    /// Reference impedance (Ω). Default 50.0.
    pub z0: f64,
    /// Original frequency unit from the option line. Preserved so a
    /// `read → write` round-trip emits the same unit string.
    pub freq_unit: FreqUnit,
    /// Original on-disk numeric format. Preserved for round-trip fidelity.
    pub format: Format,
    /// Frequencies in Hz (canonical SI).
    pub freq_hz: Vec<f64>,
    /// `data[k]` is the `n_ports × n_ports` S-matrix at `freq_hz[k]`,
    /// row-major in physical (mathematical) order: index `[i*n + j] = S_{i,j}`
    /// with `i,j` 0-based.
    pub data: Vec<SMatrix>,
    /// Comments (lines beginning with `!`) preserved in source order,
    /// stripped of the leading `!` but with internal whitespace retained.
    pub comments: Vec<String>,
}

/// Read a Touchstone file from `path`. Port count is inferred from the
/// extension (`.s1p` … `.s4p`); any other extension is rejected.
pub fn read(path: &Path) -> Result<File> {
    let n_ports = port_count_from_extension(path)?;
    let bytes = std::fs::read(path).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
    let text =
        std::str::from_utf8(&bytes).map_err(|e| Error::Io(format!("{}: {e}", path.display())))?;
    parse(text, n_ports)
}

/// Write `file` to `path`. Output is deterministic: preserved comments, then
/// the option line, then one row per frequency.
pub fn write(path: &Path, file: &File) -> Result<()> {
    let text = render(file)?;
    std::fs::write(path, text).map_err(|e| Error::Io(format!("{}: {e}", path.display())))
}

// ----------------------------------------------------------------------------
// Internals
// ----------------------------------------------------------------------------

fn port_count_from_extension(path: &Path) -> Result<usize> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .ok_or_else(|| Error::TouchstoneParse {
            line: 0,
            col: 0,
            msg: format!(
                "missing or non-UTF8 file extension on {}; expected .s1p..=.s4p",
                path.display()
            ),
        })?
        .to_ascii_lowercase();
    if ext.len() != 3 || !ext.starts_with('s') || !ext.ends_with('p') {
        return Err(Error::TouchstoneParse {
            line: 0,
            col: 0,
            msg: format!("unrecognised Touchstone extension `.{ext}`; expected .s1p..=.s4p"),
        });
    }
    let digit = ext.as_bytes()[1] as char;
    match digit {
        '1'..='4' => Ok((digit as u8 - b'0') as usize),
        _ => Err(Error::TouchstoneParse {
            line: 0,
            col: 0,
            msg: format!(
                "Touchstone extension `.{ext}` not supported in Phase 0 (only .s1p..=.s4p)"
            ),
        }),
    }
}

/// Strip an inline `!`-comment from a non-comment line, returning
/// `(payload, comment_after_bang_or_None)`. The bang itself is not included
/// in the returned comment.
fn split_inline_comment(line: &str) -> (&str, Option<&str>) {
    match line.find('!') {
        Some(idx) => (&line[..idx], Some(&line[idx + 1..])),
        None => (line, None),
    }
}

fn parse(text: &str, n_ports: usize) -> Result<File> {
    let mut comments: Vec<String> = Vec::new();
    let mut option_line: Option<(usize, String)> = None;
    let mut data_lines: Vec<(usize, String)> = Vec::new();

    for (idx, raw_line) in text.lines().enumerate() {
        let line_no = idx + 1;
        let trimmed = raw_line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('!') {
            comments.push(rest.to_string());
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('#') {
            if option_line.is_some() {
                return Err(Error::TouchstoneParse {
                    line: line_no,
                    col: 1,
                    msg: "multiple option lines (`#`) found; Touchstone v1.1 allows only one"
                        .into(),
                });
            }
            option_line = Some((line_no, rest.trim().to_string()));
            continue;
        }
        // Strip any trailing inline comment, then keep the payload.
        let (payload, inline) = split_inline_comment(trimmed);
        if let Some(c) = inline {
            comments.push(c.to_string());
        }
        let payload = payload.trim();
        if payload.is_empty() {
            continue;
        }
        data_lines.push((line_no, payload.to_string()));
    }

    let (opt_line_no, opt_text) = option_line.ok_or_else(|| Error::TouchstoneParse {
        line: 0,
        col: 0,
        msg: "no option line (`# ...`) found".into(),
    })?;
    let (freq_unit, format, z0) = parse_option_line(opt_line_no, &opt_text)?;

    // Flatten data tokens.
    let mut tokens: Vec<(usize, usize, f64)> = Vec::new();
    for (line_no, payload) in &data_lines {
        let mut col = 1usize;
        for tok in payload.split_whitespace() {
            // Column is approximate when the source uses multi-character
            // whitespace between tokens: we advance `col` past the previous
            // token and `find` the next non-whitespace match in the
            // remaining slice, so a tab- or multi-space-separated layout
            // yields the column of the first byte of the token rather than
            // a precise visual column.
            let pos = payload[col - 1..].find(tok).map(|p| col + p).unwrap_or(col);
            let val: f64 = tok.parse().map_err(|_| Error::TouchstoneParse {
                line: *line_no,
                col: pos,
                msg: format!("expected float, got `{tok}`"),
            })?;
            tokens.push((*line_no, pos, val));
            col = pos + tok.len();
        }
    }

    // Each frequency record is 1 + 2*n^2 floats.
    let per_record = 1 + 2 * n_ports * n_ports;
    if tokens.is_empty() {
        return Err(Error::TouchstoneParse {
            line: 0,
            col: 0,
            msg: "no data lines found".into(),
        });
    }
    if !tokens.len().is_multiple_of(per_record) {
        // Safe: we already returned above when `tokens.is_empty()`.
        let (line, col, _) = tokens.last().copied().unwrap_or((0, 0, 0.0));
        return Err(Error::TouchstoneParse {
            line,
            col,
            msg: format!(
                "expected a multiple of {per_record} floats per frequency record \
                 ({n_ports}-port file: 1 frequency + 2·{n_ports}² S-parameter scalars), \
                 found {} total floats",
                tokens.len()
            ),
        });
    }

    let n_records = tokens.len() / per_record;
    let mut freq_hz = Vec::with_capacity(n_records);
    let mut data: Vec<SMatrix> = Vec::with_capacity(n_records);
    let mult = freq_unit.to_hz_multiplier();
    for record_idx in 0..n_records {
        let base = record_idx * per_record;
        let f_raw = tokens[base].2;
        freq_hz.push(f_raw * mult);
        let mut mat: SMatrix = vec![Complex64::new(0.0, 0.0); n_ports * n_ports];
        // Read 2·n² scalars, decode to n² complex values in source order.
        for k in 0..(n_ports * n_ports) {
            let a = tokens[base + 1 + 2 * k].2;
            let b = tokens[base + 1 + 2 * k + 1].2;
            mat[k] = decode_complex(format, a, b);
        }
        // Reorder Touchstone-on-disk layout into mathematical row-major.
        let mat = on_disk_to_row_major(n_ports, &mat);
        data.push(mat);
    }

    let file = File {
        n_ports,
        z0,
        freq_unit,
        format,
        freq_hz,
        data,
        comments,
    };

    check_passivity(&file)?;
    Ok(file)
}

/// Convert the on-disk per-frequency S-matrix layout into mathematical
/// row-major. For n = 2 the on-disk order is `S11 S21 S12 S22`; otherwise it
/// is already row-major.
fn on_disk_to_row_major(n: usize, on_disk: &[Complex64]) -> SMatrix {
    if n == 2 {
        // On-disk: [S11, S21, S12, S22] → row-major [S11, S12, S21, S22]
        vec![on_disk[0], on_disk[2], on_disk[1], on_disk[3]]
    } else {
        on_disk.to_vec()
    }
}

/// Inverse of [`on_disk_to_row_major`].
fn row_major_to_on_disk(n: usize, row_major: &[Complex64]) -> SMatrix {
    if n == 2 {
        // Row-major: [S11, S12, S21, S22] → on-disk [S11, S21, S12, S22]
        vec![row_major[0], row_major[2], row_major[1], row_major[3]]
    } else {
        row_major.to_vec()
    }
}

fn parse_option_line(line_no: usize, body: &str) -> Result<(FreqUnit, Format, f64)> {
    let mut iter = body.split_whitespace();
    let unit_tok = iter.next().ok_or_else(|| Error::TouchstoneParse {
        line: line_no,
        col: 1,
        msg: "option line missing frequency unit".into(),
    })?;
    let freq_unit = FreqUnit::parse(unit_tok).ok_or_else(|| Error::TouchstoneParse {
        line: line_no,
        col: 1,
        msg: format!("unknown frequency unit `{unit_tok}`; expected Hz/kHz/MHz/GHz"),
    })?;

    let param_tok = iter.next().ok_or_else(|| Error::TouchstoneParse {
        line: line_no,
        col: 1,
        msg: "option line missing parameter type".into(),
    })?;
    match param_tok.to_ascii_uppercase().as_str() {
        "S" => {}
        "Y" | "Z" | "H" | "G" => {
            return Err(Error::TouchstoneParse {
                line: line_no,
                col: 1,
                msg: format!("parameter type `{param_tok}` not supported in Phase 0 (S only)"),
            });
        }
        other => {
            return Err(Error::TouchstoneParse {
                line: line_no,
                col: 1,
                msg: format!("unknown parameter type `{other}`; expected S"),
            });
        }
    }

    let fmt_tok = iter.next().ok_or_else(|| Error::TouchstoneParse {
        line: line_no,
        col: 1,
        msg: "option line missing data format".into(),
    })?;
    let format = Format::parse(fmt_tok).ok_or_else(|| Error::TouchstoneParse {
        line: line_no,
        col: 1,
        msg: format!("unknown data format `{fmt_tok}`; expected RI/MA/DB"),
    })?;

    // Optional "R <Z0>"
    let z0 = match iter.next() {
        None => 50.0_f64,
        Some(tok) => {
            if !tok.eq_ignore_ascii_case("R") {
                return Err(Error::TouchstoneParse {
                    line: line_no,
                    col: 1,
                    msg: format!("expected `R` token after format, got `{tok}`"),
                });
            }
            let val_tok = iter.next().ok_or_else(|| Error::TouchstoneParse {
                line: line_no,
                col: 1,
                msg: "`R` token not followed by an impedance value".into(),
            })?;
            val_tok.parse::<f64>().map_err(|_| Error::TouchstoneParse {
                line: line_no,
                col: 1,
                msg: format!("reference impedance `{val_tok}` is not a valid float"),
            })?
        }
    };

    if iter.next().is_some() {
        return Err(Error::TouchstoneParse {
            line: line_no,
            col: 1,
            msg: "option line has trailing tokens after `R <Z0>`".into(),
        });
    }

    Ok((freq_unit, format, z0))
}

fn decode_complex(format: Format, a: f64, b: f64) -> Complex64 {
    match format {
        Format::RealImag => Complex64::new(a, b),
        Format::MagAngle => {
            let theta = b.to_radians();
            Complex64::from_polar(a, theta)
        }
        Format::DecibelAngle => {
            let mag = 10.0_f64.powf(a / 20.0);
            let theta = b.to_radians();
            Complex64::from_polar(mag, theta)
        }
    }
}

fn encode_complex(format: Format, z: Complex64) -> (f64, f64) {
    match format {
        Format::RealImag => (z.re, z.im),
        Format::MagAngle => {
            let r = z.norm();
            let theta = z.arg().to_degrees();
            (r, theta)
        }
        Format::DecibelAngle => {
            let r = z.norm();
            let db = if r > 0.0 {
                20.0 * r.log10()
            } else {
                f64::NEG_INFINITY
            };
            let theta = z.arg().to_degrees();
            (db, theta)
        }
    }
}

fn render(file: &File) -> Result<String> {
    if !(1..=4).contains(&file.n_ports) {
        return Err(Error::InvalidFile(format!(
            "n_ports = {} is outside the Phase 0 range (1..=4)",
            file.n_ports
        )));
    }
    if file.freq_hz.len() != file.data.len() {
        return Err(Error::InvalidFile(format!(
            "freq_hz.len() = {} != data.len() = {}",
            file.freq_hz.len(),
            file.data.len()
        )));
    }
    let n = file.n_ports;
    for (i, mat) in file.data.iter().enumerate() {
        if mat.len() != n * n {
            return Err(Error::InvalidFile(format!(
                "S-matrix at index {i} has length {} but expected {} for an {}-port file",
                mat.len(),
                n * n,
                n
            )));
        }
    }
    if !file.z0.is_finite() {
        return Err(Error::InvalidFile(format!(
            "reference impedance Z0 = {} is not finite",
            file.z0
        )));
    }
    for (k, f) in file.freq_hz.iter().enumerate() {
        if !f.is_finite() {
            return Err(Error::InvalidFile(format!(
                "frequency at index {k} = {f} is not finite"
            )));
        }
    }

    let mut out = String::new();
    for c in &file.comments {
        out.push('!');
        out.push_str(c);
        out.push('\n');
    }
    // Option line.
    out.push_str(&format!(
        "# {} S {} R {}\n",
        file.freq_unit.as_str(),
        file.format.as_str(),
        format_g(file.z0),
    ));

    let mult = file.freq_unit.to_hz_multiplier();
    for (k, freq) in file.freq_hz.iter().enumerate() {
        let f_in_unit = freq / mult;
        let on_disk = row_major_to_on_disk(n, &file.data[k]);
        let mut line = String::new();
        line.push_str(&format_g(f_in_unit));
        // S-parameter cell indices are 0-based in row-major math order;
        // for n=2 we already swapped to on-disk order above, but the
        // diagnostic below refers to the on-disk slot order so a user can
        // locate the offending value in the emitted file.
        for (slot, z) in on_disk.iter().enumerate() {
            let (a, b) = encode_complex(file.format, *z);
            if !a.is_finite() || !b.is_finite() {
                return Err(Error::InvalidFile(format!(
                    "S-matrix entry at freq index {k} (slot {slot}) produced \
                     a non-finite value ({a}, {b}) under format {:?}; use a \
                     finite dB floor (e.g. -200 dB) for zero-magnitude entries \
                     when emitting `DB` files",
                    file.format,
                )));
            }
            line.push(' ');
            line.push_str(&format_g(a));
            line.push(' ');
            line.push_str(&format_g(b));
        }
        out.push_str(&line);
        out.push('\n');
    }

    Ok(out)
}

/// Render an `f64` in Rust's shortest-decimal form, which losslessly
/// round-trips any finite value since Rust 1.55. Analogous to C's `%g`
/// for finite inputs. Returns `"inf"` / `"-inf"` / `"NaN"` for non-finite
/// inputs — callers are responsible for validating finiteness before
/// passing values to this function (see the finite-value checks at the
/// top of [`render`]).
fn format_g(x: f64) -> String {
    format!("{x}")
}

/// Check passivity of every frequency sample. We use power iteration on
/// `S† S` to estimate `λ_max`. For n ≤ 4 this converges in a handful of
/// iterations.
fn check_passivity(file: &File) -> Result<()> {
    let tol = 1.0 + 1e-9;
    let n = file.n_ports;
    for (k, mat) in file.data.iter().enumerate() {
        let lambda_max = max_eig_s_dagger_s(n, mat);
        if lambda_max > tol {
            return Err(Error::TouchstoneParse {
                line: 0,
                col: 0,
                msg: format!(
                    "passivity violation at frequency index {k} \
                     (f = {} Hz): max eigenvalue of S†S = {lambda_max:.6e} > 1 + 1e-9",
                    file.freq_hz[k],
                ),
            });
        }
    }
    Ok(())
}

/// Largest eigenvalue of the Hermitian PSD matrix `S† S` via power
/// iteration. `mat` is row-major `n×n`. Sufficient for n ≤ 4.
///
/// Starting-vector caveat: we use the all-ones vector, which has non-zero
/// projection onto the dominant eigenvector for every physical S-matrix
/// that arises in Phase 0 (`n ≤ 4`, finite reciprocal networks). A
/// pathological input perfectly orthogonal to that eigenvector (for
/// example a strictly antisymmetric matrix whose dominant eigenvector
/// sits in the antisymmetric subspace) would leave `lambda` at 0.0 and
/// silently miss the passivity violation. Phase 1 should swap this for a
/// proper Hermitian eigensolver (e.g. `nalgebra::SymmetricEigen`) once
/// the workspace tolerates a heavier linear-algebra dep here.
fn max_eig_s_dagger_s(n: usize, mat: &[Complex64]) -> f64 {
    // Compute M = S† S (n×n Hermitian PSD), row-major real-imag.
    let m = compute_s_dagger_s(n, mat);
    // Start with all-ones vector. See the function doc-comment for the
    // caveat on antisymmetric edge cases.
    let mut v: Vec<Complex64> = vec![Complex64::new(1.0, 0.0); n];
    // Normalise.
    normalise(&mut v);
    let mut lambda = 0.0f64;
    for _ in 0..200 {
        // w = M v
        let mut w: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n];
        for i in 0..n {
            for j in 0..n {
                w[i] += m[i * n + j] * v[j];
            }
        }
        // Rayleigh quotient (M is Hermitian → real).
        let mut num = 0.0f64;
        for i in 0..n {
            num += (w[i].conj() * v[i]).re;
        }
        let den = norm_sq(&v);
        let new_lambda = if den > 0.0 { num / den } else { 0.0 };
        normalise(&mut w);
        v = w;
        if (new_lambda - lambda).abs() < 1e-14 * (1.0 + new_lambda.abs()) {
            lambda = new_lambda;
            break;
        }
        lambda = new_lambda;
    }
    lambda
}

fn compute_s_dagger_s(n: usize, mat: &[Complex64]) -> Vec<Complex64> {
    let mut out = vec![Complex64::new(0.0, 0.0); n * n];
    // (S† S)_{i,j} = Σ_k conj(S_{k,i}) * S_{k,j}
    for i in 0..n {
        for j in 0..n {
            let mut acc = Complex64::new(0.0, 0.0);
            for k in 0..n {
                acc += mat[k * n + i].conj() * mat[k * n + j];
            }
            out[i * n + j] = acc;
        }
    }
    out
}

fn norm_sq(v: &[Complex64]) -> f64 {
    v.iter().map(|z| z.norm_sqr()).sum()
}

fn normalise(v: &mut [Complex64]) {
    let n = norm_sq(v).sqrt();
    if n > 0.0 {
        for z in v.iter_mut() {
            *z /= n;
        }
    }
}

// ----------------------------------------------------------------------------
// Unit tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_option_line_minimal() {
        let (u, f, z) = parse_option_line(1, "GHz S RI").unwrap();
        assert_eq!(u, FreqUnit::GHz);
        assert_eq!(f, Format::RealImag);
        assert_eq!(z, 50.0);
    }

    #[test]
    fn parse_option_line_with_r() {
        let (u, f, z) = parse_option_line(1, "MHz S MA R 75").unwrap();
        assert_eq!(u, FreqUnit::MHz);
        assert_eq!(f, Format::MagAngle);
        assert_eq!(z, 75.0);
    }

    #[test]
    fn parse_option_line_rejects_y_param() {
        let err = parse_option_line(1, "GHz Y RI R 50").unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not supported"), "{msg}");
    }

    #[test]
    fn parse_option_line_case_insensitive() {
        let (u, f, _) = parse_option_line(1, "ghz s db").unwrap();
        assert_eq!(u, FreqUnit::GHz);
        assert_eq!(f, Format::DecibelAngle);
    }

    #[test]
    fn freq_unit_multipliers() {
        assert_eq!(FreqUnit::Hz.to_hz_multiplier(), 1.0);
        assert_eq!(FreqUnit::KHz.to_hz_multiplier(), 1.0e3);
        assert_eq!(FreqUnit::MHz.to_hz_multiplier(), 1.0e6);
        assert_eq!(FreqUnit::GHz.to_hz_multiplier(), 1.0e9);
    }

    #[test]
    fn complex_codec_ri_roundtrip() {
        let z = Complex64::new(0.5, -0.25);
        let (a, b) = encode_complex(Format::RealImag, z);
        let z2 = decode_complex(Format::RealImag, a, b);
        assert!((z - z2).norm() < 1e-15);
    }

    #[test]
    fn complex_codec_ma_roundtrip() {
        let z = Complex64::from_polar(0.6, 1.2);
        let (a, b) = encode_complex(Format::MagAngle, z);
        let z2 = decode_complex(Format::MagAngle, a, b);
        assert!((z - z2).norm() < 1e-12);
    }

    #[test]
    fn complex_codec_db_roundtrip() {
        let z = Complex64::from_polar(0.6, 1.2);
        let (a, b) = encode_complex(Format::DecibelAngle, z);
        let z2 = decode_complex(Format::DecibelAngle, a, b);
        assert!((z - z2).norm() < 1e-12);
    }

    #[test]
    fn on_disk_layout_n2_swaps_off_diagonals() {
        // Row-major: S11=1, S12=2, S21=3, S22=4
        let row_major = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(3.0, 0.0),
            Complex64::new(4.0, 0.0),
        ];
        let disk = row_major_to_on_disk(2, &row_major);
        // On-disk for n=2 is S11, S21, S12, S22
        assert_eq!(disk[0], Complex64::new(1.0, 0.0)); // S11
        assert_eq!(disk[1], Complex64::new(3.0, 0.0)); // S21
        assert_eq!(disk[2], Complex64::new(2.0, 0.0)); // S12
        assert_eq!(disk[3], Complex64::new(4.0, 0.0)); // S22
        let back = on_disk_to_row_major(2, &disk);
        assert_eq!(back, row_major);
    }

    #[test]
    fn passivity_passes_for_identity_like() {
        let file = File {
            n_ports: 2,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz: vec![1e9],
            data: vec![vec![
                Complex64::new(0.5, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.5, 0.0),
            ]],
            comments: vec![],
        };
        check_passivity(&file).unwrap();
    }

    #[test]
    fn passivity_rejects_gain() {
        let file = File {
            n_ports: 2,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz: vec![1e9],
            data: vec![vec![
                Complex64::new(2.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(2.0, 0.0),
            ]],
            comments: vec![],
        };
        let err = check_passivity(&file).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("passivity"), "{msg}");
    }

    #[test]
    fn render_rejects_zero_magnitude_under_db_format() {
        // DB encoding of a zero S-parameter would emit -inf dB, which is
        // not representable in Touchstone. render() must surface this as
        // Error::InvalidFile so the file is never written.
        let file = File {
            n_ports: 1,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::DecibelAngle,
            freq_hz: vec![1.0e9],
            data: vec![vec![Complex64::new(0.0, 0.0)]],
            comments: vec![],
        };
        let err = render(&file).unwrap_err();
        match &err {
            Error::InvalidFile(msg) => {
                assert!(
                    msg.contains("non-finite") && msg.contains("-200 dB"),
                    "expected dB-floor guidance, got: {msg}"
                );
            }
            other => panic!("expected InvalidFile, got: {other:?}"),
        }
    }

    #[test]
    fn render_accepts_zero_magnitude_under_ri_format() {
        // The same zero S-parameter under RI must succeed — only DB fails.
        let file = File {
            n_ports: 1,
            z0: 50.0,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz: vec![1.0e9],
            data: vec![vec![Complex64::new(0.0, 0.0)]],
            comments: vec![],
        };
        let s = render(&file).expect("render should succeed for RI zero entry");
        assert!(s.contains("# GHz S RI R 50"));
        assert!(s.contains("0 0"), "expected `0 0` data row, got: {s}");
    }

    #[test]
    fn render_rejects_non_finite_z0() {
        let file = File {
            n_ports: 1,
            z0: f64::INFINITY,
            freq_unit: FreqUnit::GHz,
            format: Format::RealImag,
            freq_hz: vec![1.0e9],
            data: vec![vec![Complex64::new(0.0, 0.0)]],
            comments: vec![],
        };
        let err = render(&file).unwrap_err();
        match err {
            Error::InvalidFile(msg) => {
                assert!(msg.contains("Z0"), "{msg}");
            }
            other => panic!("expected InvalidFile, got: {other:?}"),
        }
    }
}
