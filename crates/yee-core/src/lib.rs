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
pub type Result<T> = core::result::Result<T, Error>;

/// Crate-wide error type. Specific solver/IO errors compose into this.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Invalid input from the caller (out-of-range frequencies, malformed geometry, etc.).
    #[error("invalid input: {0}")]
    Invalid(String),

    /// Underlying numerical failure (singular matrix, divergence, NaN).
    #[error("numerical failure: {0}")]
    Numerical(String),

    /// Generic placeholder while the crate is under construction.
    #[error("unimplemented: {0}")]
    Unimplemented(&'static str),
}

/// Physical units and constants.
pub mod units {
    /// Speed of light in vacuum (m/s).
    pub const C0: f64 = 299_792_458.0;
    /// Vacuum permittivity (F/m).
    pub const EPS0: f64 = 8.854_187_812_8e-12;
    /// Vacuum permeability (H/m), CODATA 2018.
    pub const MU0: f64 = 1.256_637_062_12e-6;
    /// Free-space impedance (Ω).
    pub const ETA0: f64 = 376.730_313_668;
}

/// A linear frequency range: `start`, `stop`, `n_points`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreqRange {
    /// Start frequency (Hz).
    pub start_hz: f64,
    /// Stop frequency (Hz).
    pub stop_hz: f64,
    /// Number of points (inclusive of endpoints).
    pub n_points: usize,
}

/// Solver-agnostic skeleton. Concrete solvers (planar MoM, 3D FDTD) implement this.
pub trait Solver {
    /// Geometry type accepted by this solver.
    type Geometry;
    /// Output type produced by [`Solver::run`].
    type Output;

    /// Run the solver to completion. May be GPU-bound, CPU-bound, or both.
    fn run(&self, geometry: &Self::Geometry, freq: FreqRange) -> Result<Self::Output>;
}
