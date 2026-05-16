//! # yee-fdtd
//!
//! 3D FDTD on the Yee staggered grid. **Phase 2 deliverable**; this crate is a stub
//! during Phase 0 / 1 so the workspace builds end-to-end.
//!
//! See `README.md` and `ROADMAP.md` in this crate for full scope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// FDTD-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Functionality not yet available in the current phase.
    #[error("yee-fdtd is a Phase-2 deliverable; not available yet")]
    NotYet,
}

/// FDTD-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;
