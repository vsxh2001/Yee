//! # yee-fdtd
//!
//! 3D FDTD on the Yee staggered grid. This crate is being built up as the
//! **Phase 2 walking skeleton**: a CPU-only, single-threaded, scalar (FP64)
//! Yee solver.
//!
//! This commit adds the Gaussian point source and the hard PEC outer
//! boundary. CPML (the *real* absorbing boundary) is **Phase 2.1+ work** —
//! see [`boundary`] for the limitation. The solver wrapper follows in the
//! next commit.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod boundary;
pub mod grid;
pub mod sources;
pub mod update;

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
