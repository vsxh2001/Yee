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
pub(crate) mod gpof;
pub(crate) mod greens;
pub mod iterative;
pub(crate) mod multilayer;
pub mod ports;
pub(crate) mod quadrature;
pub mod roughness;
pub(crate) mod solve;

pub use iterative::{GmresParams, GmresResult, gmres_jacobi};
pub use roughness::{RoughnessModel, SIGMA_COPPER};

/// Public re-export of the [`Greens`] trait — the abstraction over
/// MPIE Green's-function kernels. Phase 1.1 promoted this from
/// `__internal` to the stable surface so callers (including
/// `yee-validation`) can name it in [`GreensSpec`] documentation
/// without reaching through the test-only module.
pub use crate::greens::{FreeSpaceGreen, Greens};

/// Public re-export of the multilayer placeholder Green's function.
/// Promoted from `__internal` alongside [`Greens`] so that
/// [`GreensSpec::Microstrip`] can be documented in terms of the
/// concrete kernel it builds.
pub use crate::multilayer::MultilayerGreens;

use num_complex::Complex64;
use yee_core::{FreqRange, Solver};
use yee_mesh::TriMesh;

/// Frequency-agnostic specification of the MPIE Green's-function kernel
/// to use in [`PlanarMoM::run`].
///
/// The MPIE assembly in [`crate::fill::impedance_matrix`] is generic over
/// the [`crate::greens::Greens`] trait, but every concrete `Greens`
/// instance is **constructed at a specific frequency** (the wave number
/// `k₀ = ω / c` is baked in). [`PlanarMoM::run`] sweeps across an
/// entire [`FreqRange`], so what it stores cannot be a single `Greens`
/// — it must be a freq-agnostic *spec* that the sweep can re-evaluate
/// at each point.
///
/// This is the spec type. It is intentionally a small enum rather than
/// a `Box<dyn GreensFactory>` trait object so callers can construct it
/// inline without naming a separate factory type, and so the variant
/// list captures the supported kernels at the type-system level.
///
/// New kernels (Phase 1.1.1 multi-image / Sommerfeld DCIM, future
/// stratified-media spec) extend this enum; do not introduce a parallel
/// trait-object hierarchy.
#[derive(Debug, Clone, Copy, Default)]
pub enum GreensSpec {
    /// Free-space scalar Green's function — `exp(-j k₀ R) / (4 π R)`,
    /// identical for vector and scalar potentials. The Phase 1.0
    /// default, used by mom-001 (free-space half-wave dipole).
    #[default]
    FreeSpace,
    /// Phase 1.1.0 multilayer placeholder: substrate slab of relative
    /// permittivity `eps_r` and thickness `h_m` over a PEC ground
    /// plane, evaluated with the one-image DCIM approximation in
    /// [`crate::multilayer::MultilayerGreens`]. The TE / TM split is
    /// collapsed; the N-image fit lives under
    /// [`GreensSpec::MicrostripDcim`].
    Microstrip {
        /// Relative permittivity of the substrate slab.
        eps_r: f64,
        /// Substrate thickness in metres; PEC ground at `z = -h_m`.
        h_m: f64,
    },
    /// Phase 1.1.1.0 multi-image DCIM: same substrate geometry as
    /// [`GreensSpec::Microstrip`] but the [`MultilayerGreens`] kernel
    /// is built with `n_images` complex image pairs fitted via GPOF
    /// against the slab's TE / TM spectral reflection coefficients
    /// ([`MultilayerGreens::new_microstrip_with_n_images`]). The
    /// TE/TM split that Phase 1.1.0 collapsed is resolved here, so
    /// the vector and scalar potentials use independent image trains.
    ///
    /// `n_images = 1` reduces to the [`GreensSpec::Microstrip`] path
    /// bit-for-bit; the recommended value for FR-4 microstrip is
    /// `n_images = 5` (Aksun 1996). Phase 1.1.1.1 (real Sommerfeld
    /// extraction with surface-wave pole subtraction) supersedes
    /// this; until then this is the preferred multilayer spec.
    MicrostripDcim {
        /// Relative permittivity of the substrate slab.
        eps_r: f64,
        /// Substrate thickness in metres; PEC ground at `z = -h_m`.
        h_m: f64,
        /// Number of DCIM image pairs to fit. Typical: 5.
        n_images: usize,
    },
}

impl GreensSpec {
    /// Convenience constructor for the microstrip placeholder. Mirrors
    /// [`MultilayerGreens::new_microstrip`]'s parameter order minus the
    /// frequency, which the sweep supplies at evaluation time.
    pub fn microstrip(eps_r: f64, h_m: f64) -> Self {
        Self::Microstrip { eps_r, h_m }
    }

    /// Convenience constructor for the N-image DCIM microstrip kernel.
    /// Routes through
    /// [`MultilayerGreens::new_microstrip_with_n_images`] at sweep time.
    /// `n_images = 1` is functionally identical to
    /// [`Self::microstrip`]; values in `2..=10` exercise the GPOF fit.
    pub fn microstrip_dcim(eps_r: f64, h_m: f64, n_images: usize) -> Self {
        Self::MicrostripDcim {
            eps_r,
            h_m,
            n_images,
        }
    }

    /// Build a concrete [`Greens`] kernel at `freq_hz`. Used by the
    /// per-frequency hot loop inside `s_parameters_sweep`.
    pub(crate) fn build(&self, freq_hz: f64) -> Box<dyn Greens + Send + Sync> {
        match *self {
            Self::FreeSpace => Box::new(FreeSpaceGreen::new(freq_hz)),
            Self::Microstrip { eps_r, h_m } => {
                Box::new(MultilayerGreens::new_microstrip(freq_hz, eps_r, h_m))
            }
            Self::MicrostripDcim {
                eps_r,
                h_m,
                n_images,
            } => Box::new(MultilayerGreens::new_microstrip_with_n_images(
                eps_r, h_m, freq_hz, n_images,
            )),
        }
    }
}

// Boxed Greens satisfies the Greens trait by forwarding through the
// pointer. `impedance_matrix` requires `G: Greens + Sync`, and
// `Box<dyn Greens + Send + Sync>` is both `Sync` (because the inner trait
// object is `Sync`) and `Greens` (via this blanket impl). The blanket is
// gated on `?Sized` so it covers both sized boxes and trait-object
// boxes.
impl<T: Greens + ?Sized> Greens for Box<T> {
    fn k0(&self) -> Complex64 {
        (**self).k0()
    }
    fn eta0(&self) -> f64 {
        (**self).eta0()
    }
    fn scalar_vector(&self, r1: nalgebra::Vector3<f64>, r2: nalgebra::Vector3<f64>) -> Complex64 {
        (**self).scalar_vector(r1, r2)
    }
    fn scalar_scalar(&self, r1: nalgebra::Vector3<f64>, r2: nalgebra::Vector3<f64>) -> Complex64 {
        (**self).scalar_scalar(r1, r2)
    }
    fn scalar_vector_smooth(
        &self,
        r1: nalgebra::Vector3<f64>,
        r2: nalgebra::Vector3<f64>,
    ) -> Complex64 {
        (**self).scalar_vector_smooth(r1, r2)
    }
    fn scalar_scalar_smooth(
        &self,
        r1: nalgebra::Vector3<f64>,
        r2: nalgebra::Vector3<f64>,
    ) -> Complex64 {
        (**self).scalar_scalar_smooth(r1, r2)
    }
}

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

/// The planar MoM solver.
///
/// Holds the Green's-function spec used during impedance-matrix
/// assembly. The default is [`GreensSpec::FreeSpace`], which preserves
/// the mom-001 (free-space half-wave dipole) numerics bit-for-bit.
/// Switch to a multilayer kernel via [`PlanarMoM::with_greens`].
#[derive(Debug, Default)]
pub struct PlanarMoM {
    /// The Green's-function spec to use when filling the impedance
    /// matrix at each frequency. Defaults to [`GreensSpec::FreeSpace`].
    greens: GreensSpec,
    // TODO(phase-0): mesh, ports, GPU context.
}

impl PlanarMoM {
    /// Replace the default [`GreensSpec::FreeSpace`] kernel with the
    /// supplied spec. The Greens kernel is rebuilt per frequency inside
    /// the sweep, so the spec itself remains frequency-agnostic — see
    /// [`GreensSpec`] for the rationale.
    ///
    /// This is the entry point Phase 1.1.1 (real Sommerfeld extraction)
    /// will exercise once the production multilayer kernel lands; in
    /// Phase 1.1.0 it routes through the one-image DCIM placeholder
    /// documented at the [`MultilayerGreens`](crate::multilayer::MultilayerGreens)
    /// level.
    pub fn with_greens(mut self, greens: GreensSpec) -> Self {
        self.greens = greens;
        self
    }
}

impl Solver for PlanarMoM {
    type Geometry = TriMesh;
    type Output = SParameters;

    fn run(&self, geometry: &Self::Geometry, freq: FreqRange) -> yee_core::Result<Self::Output> {
        let basis = basis::RwgBasis::from_mesh(geometry.clone())?;
        // mom-001 / Phase 1.0 default excitation: 1 V delta-gap on port_tag 1.
        // Phase 1.3 routes this through the `Port` trait without changing
        // numerics — `DeltaGapPort { tag: 1, voltage: 1+0i }` reproduces the
        // legacy `delta_gap_rhs(..., 1)` bit-for-bit.
        let port = ports::DeltaGapPort {
            tag: 1,
            voltage: Complex64::new(1.0, 0.0),
        };
        let file = solve::s_parameters_sweep(&basis, &port, freq, 50.0, None, &self.greens)?;
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

    use crate::fill::impedance_matrix;
    use crate::solve::delta_gap_rhs;
    use faer::linalg::solvers::{PartialPivLu, Solve};
    use num_complex::Complex64;
    use yee_core::Error;
    use yee_mesh::TriMesh;

    /// Public re-export of the crate-private RWG basis type so integration
    /// tests can inspect port edges, lengths, and counts without forcing the
    /// real `basis` module to be `pub`.
    pub use crate::basis::RwgBasis;

    /// Public re-export of [`MultilayerGreens`] for the Phase 1.1
    /// integration tests. Not part of the stable API.
    pub use crate::multilayer::MultilayerGreens;

    /// Public re-export of the `Greens` trait and `FreeSpaceGreen` struct
    /// so `__internal` callers (the integration tests) can name them in
    /// generic signatures. Not part of the stable API.
    pub use crate::greens::{FreeSpaceGreen, Greens};

    /// Test-only constructor for [`RwgBasis`] — wraps the crate-private
    /// `from_mesh` so integration tests can build a basis without making
    /// `basis::RwgBasis::from_mesh` itself public.
    pub fn build_basis(mesh: &TriMesh) -> Result<RwgBasis, Error> {
        RwgBasis::from_mesh(mesh.clone())
    }

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

    /// Diagnostic helper: solve `Z·i = b` at `freq_hz` and return both
    /// `Z_in = V_port / I_port` and the relative LU residual
    /// `||Z·i - b||_2 / ||b||_2`. A clean LU should produce a residual
    /// well below `1e-10` on this geometry; anything larger indicates the
    /// solve itself is broken rather than the formulation.
    ///
    /// This mirrors `solve::s_parameters_at_freq` but exposes `Z_in`
    /// directly (instead of `S11`) and the LU residual instead of
    /// swallowing it. It is intentionally a separate helper rather than a
    /// public surface change on `solve` because the residual diagnostic
    /// is not something callers should depend on long term.
    pub fn z_in_and_residual_at_freq(
        mesh: &TriMesh,
        port_tag: u32,
        freq_hz: f64,
        _z0_ref: f64,
    ) -> Result<(Complex64, f64), Error> {
        let basis = RwgBasis::from_mesh(mesh.clone())?;
        let green = FreeSpaceGreen::new(freq_hz);
        let z = impedance_matrix(&basis, &green);
        let b = delta_gap_rhs(&basis, port_tag);

        let lu = PartialPivLu::new(z.as_ref());
        let i = lu.solve(b.as_ref());

        // Residual norm: ||Z·i - b||_2 / ||b||_2. faer's Mat does not have a
        // direct `*` operator producing a Mat (it uses MatRef × MatRef on the
        // generic level), so we hand-roll the matvec to keep this self-
        // contained — n is small enough (~1k) that the cost is irrelevant.
        let n = z.nrows();
        let mut residual_sq = 0.0_f64;
        let mut b_norm_sq = 0.0_f64;
        for m in 0..n {
            let mut zi_m = Complex64::new(0.0, 0.0);
            for k in 0..n {
                zi_m += z[(m, k)] * i[(k, 0)];
            }
            let diff = zi_m - b[(m, 0)];
            residual_sq += diff.norm_sqr();
            b_norm_sq += b[(m, 0)].norm_sqr();
        }
        let rel_residual = (residual_sq / b_norm_sq.max(f64::MIN_POSITIVE)).sqrt();

        let mut i_port = Complex64::new(0.0, 0.0);
        for k in basis.port_basis_indices(port_tag) {
            i_port += b[(k, 0)] * i[(k, 0)];
        }
        if i_port.norm() < 1e-30 {
            return Err(Error::Numerical(
                "port current vanished; check port tagging".into(),
            ));
        }
        let v_port = Complex64::new(1.0, 0.0);
        let z_in = v_port / i_port;
        Ok((z_in, rel_residual))
    }

    /// Diagnostic: return per-port-edge `(length_k, i_k)` pairs at `freq_hz`.
    /// Used to check whether the RWG +/- orientation around the cylinder
    /// port ring is consistent (all `i_k` same sign and magnitude on a
    /// symmetric mesh) or whether some port edges are flipped and partially
    /// cancelling in the `Σ b_k · i_k` sum.
    pub fn port_edge_currents(
        mesh: &TriMesh,
        port_tag: u32,
        freq_hz: f64,
    ) -> Result<Vec<(f64, Complex64)>, Error> {
        let basis = RwgBasis::from_mesh(mesh.clone())?;
        let green = FreeSpaceGreen::new(freq_hz);
        let z = impedance_matrix(&basis, &green);
        let b = delta_gap_rhs(&basis, port_tag);
        let lu = PartialPivLu::new(z.as_ref());
        let i = lu.solve(b.as_ref());
        Ok(basis
            .port_basis_indices(port_tag)
            .map(|k| (basis.edges[k].length, i[(k, 0)]))
            .collect())
    }

    /// Generic Phase 1.1 helper: assemble the MPIE impedance matrix using
    /// the supplied [`Greens`] implementation, solve a delta-gap excitation
    /// at `port_tag`, and return `Z_in = V_port / I_port`. Identical to the
    /// `z_in_and_residual_at_freq` helper but parameterised over the
    /// Green's-function kernel so multilayer integration tests can compare
    /// free-space and multilayer evaluations on the same mesh.
    pub fn z_in_with_greens<G: Greens + Sync>(
        mesh: &TriMesh,
        port_tag: u32,
        green: &G,
    ) -> Result<Complex64, Error> {
        let basis = RwgBasis::from_mesh(mesh.clone())?;
        let z = impedance_matrix(&basis, green);
        let b = delta_gap_rhs(&basis, port_tag);

        let lu = PartialPivLu::new(z.as_ref());
        let i = lu.solve(b.as_ref());

        let mut i_port = Complex64::new(0.0, 0.0);
        for k in basis.port_basis_indices(port_tag) {
            i_port += b[(k, 0)] * i[(k, 0)];
        }
        if i_port.norm() < 1e-30 {
            return Err(Error::Numerical(
                "port current vanished; check port tagging".into(),
            ));
        }
        let v_port = Complex64::new(1.0, 0.0);
        Ok(v_port / i_port)
    }

    /// Free-space convenience wrapper around [`z_in_with_greens`]: build
    /// the basis, instantiate [`FreeSpaceGreen`] at `freq_hz`, and solve.
    /// Mirrors the multilayer entry point so call sites in tests stay
    /// uniform.
    pub fn z_in_free_space(
        mesh: &TriMesh,
        port_tag: u32,
        freq_hz: f64,
    ) -> Result<Complex64, Error> {
        let green = FreeSpaceGreen::new(freq_hz);
        z_in_with_greens(mesh, port_tag, &green)
    }
}
