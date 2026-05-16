//! # yee-core
//!
//! Shared types, traits, errors, and units for the Yee electromagnetic simulation studio.
//!
//! This crate intentionally has no CUDA, no GUI, and no I/O dependencies. It is the stable
//! foundation that every other Yee crate depends on. Keep it small and well-documented.
//!
//! Phase 0 scope:
//! - Physical units and constants (`units` module)
//! - Frequency/range and time-step types
//! - The `Solver` trait skeleton
//! - The crate-wide [`Error`] type

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Crate-wide result alias.
///
/// Equivalent to [`core::result::Result<T, Error>`].
pub type Result<T> = core::result::Result<T, Error>;

/// Crate-wide error type.
///
/// Specific solver and I/O errors compose into one of these variants so callers
/// can match on a single, stable surface.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Invalid input from the caller (out-of-range frequencies, malformed geometry, etc.).
    #[error("invalid input: {0}")]
    Invalid(String),

    /// Underlying numerical failure (singular matrix, divergence, NaN).
    #[error("numerical failure: {0}")]
    Numerical(String),

    /// Feature not yet implemented in this Phase.
    #[error("unimplemented: {0}")]
    Unimplemented(&'static str),

    /// Input/output failure surfaced by an upstream crate (file, network, etc.).
    #[error("I/O error: {0}")]
    Io(String),
}

/// Physical units and CODATA 2018 reference constants.
///
/// All values are populated from the CODATA 2018 recommended set; consult
/// `tests/units.rs` for the tolerance bounds each constant is verified against.
pub mod units {
    /// Speed of light in vacuum (m/s). Defined exactly.
    pub const C0: f64 = 299_792_458.0;
    /// Vacuum permittivity (F/m), CODATA 2018.
    pub const EPS0: f64 = 8.854_187_812_8e-12;
    /// Vacuum permeability (H/m), CODATA 2018.
    pub const MU0: f64 = 1.256_637_062_12e-6;
    /// Free-space wave impedance (Ω), CODATA 2018.
    pub const ETA0: f64 = 376.730_313_668;
}

/// A linear frequency sweep specified by `start_hz`, `stop_hz`, and `n_points`.
///
/// Construct with [`FreqRange::new`] for input validation, or build the struct
/// directly when validation has already been performed. Iterate the sample
/// points with [`FreqRange::iter`].
///
/// # Examples
///
/// ```
/// use yee_core::FreqRange;
///
/// let band = FreqRange::new(1.0e9, 2.0e9, 5).unwrap();
/// assert_eq!(band.iter().count(), 5);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreqRange {
    /// Start frequency (Hz).
    pub start_hz: f64,
    /// Stop frequency (Hz).
    pub stop_hz: f64,
    /// Number of points (inclusive of endpoints).
    pub n_points: usize,
}

impl FreqRange {
    /// Construct a new [`FreqRange`], validating the inputs.
    ///
    /// Returns [`Error::Invalid`] when `start_hz >= stop_hz` or `n_points == 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use yee_core::FreqRange;
    ///
    /// let band = FreqRange::new(2.4e9, 2.5e9, 11).unwrap();
    /// assert_eq!(band.n_points, 11);
    ///
    /// assert!(FreqRange::new(2.5e9, 2.4e9, 11).is_err());
    /// assert!(FreqRange::new(1.0e9, 2.0e9, 0).is_err());
    /// ```
    pub fn new(start_hz: f64, stop_hz: f64, n_points: usize) -> Result<Self> {
        if start_hz >= stop_hz || start_hz.is_nan() || stop_hz.is_nan() {
            return Err(Error::Invalid(format!(
                "FreqRange requires start_hz < stop_hz, got start_hz = {start_hz}, stop_hz = {stop_hz}"
            )));
        }
        if n_points == 0 {
            return Err(Error::Invalid(
                "FreqRange requires n_points >= 1, got 0".to_string(),
            ));
        }
        Ok(Self {
            start_hz,
            stop_hz,
            n_points,
        })
    }

    /// Iterate the `n_points` linearly spaced frequencies in Hz.
    ///
    /// - `n_points == 1` yields `[start_hz]`.
    /// - `n_points == 2` yields `[start_hz, stop_hz]` (endpoints exact).
    /// - `n_points >= 3` yields `n_points` evenly spaced samples with both
    ///   endpoints reproduced exactly.
    ///
    /// # Examples
    ///
    /// ```
    /// use yee_core::FreqRange;
    ///
    /// let band = FreqRange::new(1.0e9, 2.0e9, 3).unwrap();
    /// let pts: Vec<f64> = band.iter().collect();
    /// assert_eq!(pts[0], 1.0e9);
    /// assert_eq!(pts[2], 2.0e9);
    /// assert_eq!(pts.len(), 3);
    /// ```
    pub fn iter(&self) -> FreqRangeIter {
        FreqRangeIter {
            start_hz: self.start_hz,
            stop_hz: self.stop_hz,
            n_points: self.n_points,
            index: 0,
        }
    }
}

/// Iterator over the sample points produced by [`FreqRange::iter`].
///
/// Pins both endpoints exactly: `pts[0] == start_hz` and, when `n_points >= 2`,
/// the final yielded value equals `stop_hz` bit-for-bit.
#[derive(Debug, Clone)]
pub struct FreqRangeIter {
    start_hz: f64,
    stop_hz: f64,
    n_points: usize,
    index: usize,
}

impl Iterator for FreqRangeIter {
    type Item = f64;

    fn next(&mut self) -> Option<f64> {
        if self.index >= self.n_points {
            return None;
        }
        let i = self.index;
        self.index += 1;

        // Pin endpoints exactly to avoid floating-point rounding at the bounds.
        if i == 0 {
            return Some(self.start_hz);
        }
        if i + 1 == self.n_points {
            return Some(self.stop_hz);
        }
        // Interior samples for n_points >= 3.
        let denom = (self.n_points - 1) as f64;
        let t = (i as f64) / denom;
        Some(self.start_hz + t * (self.stop_hz - self.start_hz))
    }
}

/// Solver-agnostic skeleton. Concrete solvers (planar MoM, 3D FDTD) implement this.
///
/// # Examples
///
/// ```
/// use yee_core::{FreqRange, Result, Solver};
///
/// struct NullSolver;
/// impl Solver for NullSolver {
///     type Geometry = ();
///     type Output = usize;
///     fn run(&self, _geometry: &(), freq: FreqRange) -> Result<usize> {
///         Ok(freq.n_points)
///     }
/// }
///
/// let band = FreqRange::new(1.0e9, 2.0e9, 7).unwrap();
/// assert_eq!(NullSolver.run(&(), band).unwrap(), 7);
/// ```
pub trait Solver {
    /// Geometry type accepted by this solver.
    type Geometry;
    /// Output type produced by [`Solver::run`].
    type Output;

    /// Run the solver to completion. May be GPU-bound, CPU-bound, or both.
    fn run(&self, geometry: &Self::Geometry, freq: FreqRange) -> Result<Self::Output>;
}
