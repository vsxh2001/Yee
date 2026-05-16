//! # yee-mom
//!
//! Planar Method of Moments solver — the **Yee v1 beachhead**.
//!
//! Phase 0 ships a lossless, single-layer, PEC-only solver with a CPU dense LU via
//! `faer` and a GPU port via cuSOLVER hidden behind the `cuda` feature. Phase 1 adds
//! multilayer dielectric stack-ups, RWG/rooftop basis functions, lumped + wave ports,
//! TRL/SOLT de-embedding, and the production GPU path.
//!
//! See `README.md` and `ROADMAP.md` in this crate for full scope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use num_complex::Complex64;
use yee_core::{FreqRange, Solver};
use yee_mesh::TriMesh;

/// MoM-layer errors.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Mesh did not pass validation.
    #[error("mesh invalid: {0}")]
    BadMesh(String),
    /// Numerical failure during fill or solve.
    #[error("numerical: {0}")]
    Numerical(String),
    /// Phase 0 placeholder.
    #[error("not yet implemented in this phase: {0}")]
    Unimplemented(&'static str),
}

/// MoM-layer result alias.
pub type Result<T> = core::result::Result<T, Error>;

/// Multi-port S-parameter container — Phase 0 placeholder.
#[derive(Debug, Clone)]
pub struct SParameters {
    /// Frequencies (Hz) corresponding to each S-matrix row in `data`.
    pub freq_hz: Vec<f64>,
    /// `data[k]` is the n×n S-matrix at `freq_hz[k]`, row-major flat.
    pub data: Vec<Vec<Complex64>>,
    /// Number of ports (n).
    pub n_ports: usize,
}

/// The planar MoM solver. Phase 0: empty shell.
#[derive(Debug, Default)]
pub struct PlanarMoM {
    // TODO(phase-0): mesh, ports, Green's function evaluator, GPU context.
}

impl Solver for PlanarMoM {
    type Geometry = TriMesh;
    type Output = SParameters;

    fn run(&self, _geometry: &Self::Geometry, _freq: FreqRange) -> yee_core::Result<Self::Output> {
        Err(yee_core::Error::Unimplemented(
            "PlanarMoM::run not implemented in phase 0",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_constructs() {
        // Phase 0 sanity: the empty-shell solver must be default-constructible.
        let _solver = PlanarMoM::default();
    }

    #[test]
    fn run_returns_unimplemented_with_exact_message() {
        // The Phase 0 contract is that `run` returns the variant
        // `yee_core::Error::Unimplemented` with this exact static message.
        let solver = PlanarMoM::default();
        let mesh = TriMesh::default();
        let freq = FreqRange::new(1.0e9, 2.0e9, 3).expect("valid FreqRange");
        let err = solver.run(&mesh, freq).expect_err("run must return Err in Phase 0");
        match err {
            yee_core::Error::Unimplemented(msg) => {
                assert_eq!(msg, "PlanarMoM::run not implemented in phase 0");
            }
            other => panic!("expected Unimplemented, got: {other:?}"),
        }
    }
}
