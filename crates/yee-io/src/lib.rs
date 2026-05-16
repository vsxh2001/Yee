//! # yee-io
//!
//! File format I/O for Yee:
//! - Touchstone v1.1 reader/writer (`.s1p` through `.s4p` in Phase 0)
//! - CAD ingestion (STEP / IGES / KiCad) via `opencascade-rs` (feature `opencascade`)
//! - Solver output containers (HDF5 / Arrow IPC arrive Phase 1)
//!
//! Phase 0 ships Touchstone only. Everything else is feature-gated and stubs return
//! [`Error::NotEnabled`] until the corresponding feature is toggled.
//!
//! ## Touchstone
//!
//! See [`touchstone`] for the parser, writer, and [`touchstone::File`] container.
//! The IBIS Touchstone v1.1 specification is the reference grammar:
//! <https://ibis.org/connector/touchstone_spec11.pdf>.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod touchstone;

pub use touchstone::File;

/// I/O-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Underlying OS / file-system error.
    #[error("io error: {0}")]
    Io(String),
    /// Touchstone file failed parsing or read-time validation (including
    /// passivity).
    ///
    /// `line` and `col` are 1-based; `col = 0` means "column unknown / N/A".
    #[error("touchstone parse error at line {line}, col {col}: {msg}")]
    TouchstoneParse {
        /// 1-based line number in the source file.
        line: usize,
        /// 1-based column number, or `0` when not meaningful.
        col: usize,
        /// Human-readable explanation of the failure.
        msg: String,
    },
    /// In-memory file structure failed validation before writing.
    ///
    /// Used by write-side checks where there is no source line to point at —
    /// e.g. `n_ports` / S-matrix length mismatch, or a numeric encoding
    /// failure such as a `-inf` dB magnitude that Touchstone cannot
    /// represent.
    #[error("invalid touchstone file: {0}")]
    InvalidFile(String),
    /// Feature flag not enabled in this build.
    #[error("io feature `{0}` not enabled; rebuild with that feature on")]
    NotEnabled(&'static str),
}

/// I/O-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;
