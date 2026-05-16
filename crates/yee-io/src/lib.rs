//! # yee-io
//!
//! File format I/O for Yee:
//! - Touchstone v1.1 reader/writer (`.s1p` through `.sNp`)
//! - CAD ingestion (STEP / IGES / KiCad) via `opencascade-rs` (feature `opencascade`)
//! - Solver output containers (HDF5 / Arrow IPC arrive Phase 1)
//!
//! Phase 0 ships Touchstone only. Everything else is feature-gated and stubs return
//! `Error::NotEnabled` until the corresponding feature is toggled.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use num_complex::Complex64;

/// I/O-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Underlying OS / file-system error.
    #[error("io error: {0}")]
    Io(String),
    /// Touchstone file failed parsing.
    #[error("touchstone parse: {0}")]
    TouchstoneParse(String),
    /// Feature flag not enabled in this build.
    #[error("io feature `{0}` not enabled; rebuild with that feature on")]
    NotEnabled(&'static str),
}

/// I/O-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;

/// Touchstone v1.1 reader/writer skeleton.
pub mod touchstone {
    use super::*;

    /// Parsed Touchstone file content.
    #[derive(Debug, Clone)]
    pub struct File {
        /// Reference impedance (Ω). Default 50.0.
        pub z0: f64,
        /// Number of ports.
        pub n_ports: usize,
        /// Frequencies (Hz).
        pub freq_hz: Vec<f64>,
        /// `data[k]` is the n×n S-matrix at `freq_hz[k]`, row-major flat.
        pub data: Vec<Vec<Complex64>>,
    }

    /// Read a Touchstone file. Phase 0 stub.
    pub fn read(_path: &std::path::Path) -> Result<File> {
        Err(Error::TouchstoneParse(
            "Touchstone reader not yet implemented".into(),
        ))
    }

    /// Write a Touchstone file. Phase 0 stub.
    pub fn write(_path: &std::path::Path, _file: &File) -> Result<()> {
        Err(Error::TouchstoneParse(
            "Touchstone writer not yet implemented".into(),
        ))
    }
}
