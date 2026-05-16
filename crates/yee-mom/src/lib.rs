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

/// Boundary mapping: any failure surfaced by `yee_io` while writing a
/// Touchstone file is rendered into `yee_core::Error::Io` so callers higher
/// in the stack (the CLI, solver drivers, etc.) only need to match a single
/// crate-wide error surface. The full `yee_io::Error` message text — line
/// and column hints included — is preserved verbatim inside the wrapped
/// string.
fn io_to_core(e: yee_io::Error) -> yee_core::Error {
    yee_core::Error::Io(e.to_string())
}

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

impl SParameters {
    /// Build an [`SParameters`] from a parsed [`yee_io::touchstone::File`].
    ///
    /// `yee_io` already canonicalises frequencies to Hz and reorders the
    /// S-matrix into mathematical row-major (including the n = 2 off-diagonal
    /// swap), so this is a structural copy — no numeric transformation.
    pub fn from_touchstone(file: &yee_io::touchstone::File) -> Self {
        Self {
            freq_hz: file.freq_hz.clone(),
            data: file.data.clone(),
            n_ports: file.n_ports,
        }
    }

    /// Build a [`yee_io::touchstone::File`] from `self` using the Phase 0
    /// defaults: `Format::RealImag` numeric encoding and `FreqUnit::Hz` for
    /// frequencies.
    ///
    /// `FreqUnit::Hz` is hard-coded because [`SParameters::freq_hz`] is the
    /// canonical SI Hz representation — writing under any other unit would
    /// silently misinterpret the values (e.g. emitting 1e9 Hz as 1 GHz numerically
    /// is fine, but as "1e9 GHz" in the option line is a unit-mismatch bug).
    /// Callers that need a non-Hz on-disk unit or a non-RI numeric format
    /// must use [`SParameters::to_touchstone_with`] explicitly.
    ///
    /// Comments are intentionally left empty — this constructor exists for
    /// the simulation → file path where there is no source commentary to
    /// preserve.
    pub fn to_touchstone(&self, z0: f64) -> yee_io::touchstone::File {
        self.to_touchstone_with(
            z0,
            yee_io::touchstone::Format::RealImag,
            yee_io::touchstone::FreqUnit::Hz,
        )
    }

    /// Advanced-caller form of [`SParameters::to_touchstone`] that exposes
    /// the on-disk numeric format and frequency unit. Most callers want
    /// the spec-default [`SParameters::to_touchstone`] instead; reach for
    /// this only when emitting a file targeting a specific consumer's
    /// expectations (e.g. a GHz-MA legacy tool).
    ///
    /// Note: the in-memory `freq_hz` is always Hz; choosing `freq_unit`
    /// here only affects how those numbers are rendered on disk — the
    /// writer divides by the unit's multiplier when emitting.
    pub fn to_touchstone_with(
        &self,
        z0: f64,
        format: yee_io::touchstone::Format,
        freq_unit: yee_io::touchstone::FreqUnit,
    ) -> yee_io::touchstone::File {
        yee_io::touchstone::File {
            n_ports: self.n_ports,
            z0,
            freq_unit,
            format,
            freq_hz: self.freq_hz.clone(),
            data: self.data.clone(),
            comments: Vec::new(),
        }
    }

    /// Write `self` to `path` as a Touchstone v1.1 file using the same
    /// defaults as [`SParameters::to_touchstone`]: `Format::RealImag` and
    /// `FreqUnit::Hz`. Errors from `yee_io` are mapped to
    /// [`yee_core::Error::Io`] via the boundary helper documented at module
    /// level.
    pub fn write_touchstone(&self, path: &std::path::Path, z0: f64) -> yee_core::Result<()> {
        let file = self.to_touchstone(z0);
        yee_io::touchstone::write(path, &file).map_err(io_to_core)
    }
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
