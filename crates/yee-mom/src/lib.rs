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

pub(crate) mod basis;
pub(crate) mod fill;
pub(crate) mod greens;
pub(crate) mod quadrature;
pub(crate) mod solve;

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

    fn run(&self, geometry: &Self::Geometry, freq: FreqRange) -> yee_core::Result<Self::Output> {
        let basis = basis::RwgBasis::from_mesh(geometry.clone())?;
        let file = solve::s_parameters_sweep(&basis, 1, freq, 50.0)?;
        Ok(SParameters::from_touchstone(&file))
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
    fn run_without_port_tags_returns_numerical_error() {
        use nalgebra::Vector3;
        use yee_mesh::TriMesh;

        let mesh = TriMesh::new(
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.1, 0.0, 0.0),
                Vector3::new(0.1, 0.1, 0.0),
                Vector3::new(0.0, 0.1, 0.0),
            ],
            vec![[0u32, 1, 2], [0u32, 2, 3]],
            vec![0u32, 0u32], // no port tags → port edges empty → port current vanishes
        )
        .unwrap();
        let freq = FreqRange::new(1.0e9, 2.0e9, 2).unwrap();
        let result = PlanarMoM::default().run(&mesh, freq);
        match result {
            Err(yee_core::Error::Numerical(msg)) => {
                assert!(msg.contains("port current"), "got: {msg}");
            }
            other => panic!("expected Numerical error, got {other:?}"),
        }
    }
}

#[doc(hidden)]
pub mod __internal {
    //! Test-helper surface. Not stable API; do not depend on it.

    use crate::basis::RwgBasis;
    use crate::fill::impedance_matrix;
    use crate::greens::FreeSpaceGreen;
    use yee_core::Error;
    use yee_mesh::TriMesh;

    /// Build the impedance matrix and return its condition number via
    /// `cond = sigma_max / sigma_min`. Helper for the condition-number
    /// regression test; not a public API.
    ///
    /// The `_port_tag` argument is reserved for future per-port conditioning
    /// diagnostics; the matrix itself depends only on the mesh and the
    /// excitation frequency, so it is intentionally unused today.
    pub fn condition_number_at_freq(
        mesh: &TriMesh,
        _port_tag: u32,
        freq_hz: f64,
    ) -> Result<f64, Error> {
        let basis = RwgBasis::from_mesh(mesh.clone())?;
        let green = FreeSpaceGreen::new(freq_hz);
        let z = impedance_matrix(&basis, &green);

        // faer 0.23 ships a `MatRef::singular_values()` shortcut that
        // computes the SVD and returns the singular values as a plain
        // `Vec<f64>` (real, nonnegative, descending). This avoids juggling
        // the lower-level `Svd::new(...).S()` / `DiagRef::column_vector()`
        // chain — see
        // https://docs.rs/faer/0.23/faer/struct.MatRef.html#method.singular_values.
        let s = z
            .as_ref()
            .singular_values()
            .map_err(|e| Error::Numerical(format!("SVD failed: {e:?}")))?;

        let mut max_s: f64 = 0.0;
        let mut min_s: f64 = f64::INFINITY;
        for sv in s.iter().copied() {
            if sv > max_s {
                max_s = sv;
            }
            if sv > 0.0 && sv < min_s {
                min_s = sv;
            }
        }
        if min_s <= 0.0 || !min_s.is_finite() {
            return Err(Error::Numerical("Z is singular".into()));
        }
        Ok(max_s / min_s)
    }
}
