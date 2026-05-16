//! # yee-fdtd
//!
//! 3D FDTD on the Yee staggered grid. This crate is being built up as the
//! **Phase 2 walking skeleton**: a CPU-only, single-threaded, scalar (FP64)
//! Yee solver.
//!
//! This commit lands the grid + Courant stability machinery. Update kernels,
//! sources, and the solver wrapper follow in subsequent commits.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod grid;

pub use grid::YeeGrid;

/// FDTD-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Invalid input from the caller (out-of-range size, bad time step, etc.).
    #[error("invalid input: {0}")]
    Invalid(String),

    /// Numerical failure (NaN, divergence, instability).
    #[error("numerical failure: {0}")]
    Numerical(String),
}

/// FDTD-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;
