//! Phase 4.fem.eig.1 dispersive eigenmode tracker — Newton-Raphson
//! frequency tracking of `K(ω) e = (ω/c)² M(ω) e` with complex
//! `ε(ω)`, `μ(ω)` per tet.
//!
//! ## D4 scope
//!
//! [`DispersiveSolver::solve_at_frequency`] solves the **linearised**
//! complex generalised eigenproblem `K(ω₀) e = θ M(ω₀) e` at a single
//! trial angular frequency ω₀ via the complex sparse inverse-power
//! eigensolver from `crates/yee-fem/src/solve.rs` (D2). The returned
//! eigenvalues are linearised wavenumber-squareds at the trial ω; they
//! are **not yet** self-consistent dispersive eigenmodes.
//!
//! ## D5 scope (this file)
//!
//! [`DispersiveSolver::solve_with_newton`] wraps
//! [`DispersiveSolver::solve_at_frequency`] in an outer fixed-point
//! iteration that closes the self-consistency relation
//! `ω² ε(ω) μ_0 ε_0 = k²` for a single physical mode. At each
//! iteration the linearised solver returns the FEM generalised
//! eigenvalue `λ(ω_n)` from `K(ω_n) e = λ M(ω_n) e`, where
//! `K ∋ (1/μ)·curl·curl` is the stiffness matrix and
//! `M ∋ ε(ω_n)·basis·basis` is the mass matrix — i.e. `ε(ω)` is
//! **already** baked into `M` at assembly time. At a self-consistent
//! dispersive eigenmode `λ(ω*) = (ω*/c)²`, so the fixed-point update
//! is
//!
//! ```text
//!     ω_{n+1} = c · sqrt( λ(ω_n) )       (equivalently ω² = c² · λ).
//! ```
//!
//! This is **not** `ω² = k²/(μ₀ ε₀ ε(ω))` — applying that form would
//! divide by `ε(ω)` a *second* time after the FEM `M` matrix already
//! accounts for it, collapsing the converged `Re(f_FEM)` to
//! `Re(f_analytic) / √ε_∞`. See `crates/yee-fem/validation/README.md`
//! "fem-eig-002 D6 finding 1 — `solve_with_newton` ε double-divide"
//! for the QQQQQQQQ-track bug history and the 1/√3.78 measured ratio
//! that pinned the diagnosis. When the relative step
//! `|ω_{n+1} − ω_n| / |ω_n|` falls below
//! [`DispersiveSolver::newton_tol`], the iteration returns
//! [`DispersiveEigenpair`] with the converged complex `k = ω · √(μ₀ ε₀ ε(ω))`.
//! The two convergence knobs (inner-solver
//! [`DispersiveSolver::tol`] and outer-Newton
//! [`DispersiveSolver::newton_tol`]) are intentionally separate — the
//! inner Rayleigh-quotient tolerance is far tighter than the outer
//! fixed-point tolerance under any sane configuration.
//!
//! The convergence behaviour is mathematically a fixed-point iteration
//! rather than a full Hellmann–Feynman Newton step (spec §4.2 ships the
//! latter as a later milestone). For the v1 validation regime —
//! lightly-dispersive cavities warm-started from a free-space resonance
//! — the fixed-point form is monotone-convergent and avoids the
//! complex-symmetric derivative book-keeping that the full Newton step
//! requires. Step labels: this file delivers the D5 "Newton ω-tracker"
//! gate per the implementation plan; future tracks may upgrade the
//! inner update rule without changing the public surface.
//!
//! ## Pipeline at a trial ω
//!
//! 1. Build [`crate::FemEigenAssembly`] over the mesh with placeholder
//!    real-valued `ε_r = μ_r = 1` arrays (the assembler stores them
//!    but the complex path below ignores them).
//! 2. Call [`crate::FemEigenAssembly::assemble_complex`] passing
//!    `omega` and the stored [`MaterialDatabase`] — this evaluates the
//!    per-tet `ε(ω)`, `μ(ω)` via the material database and emits
//!    [`crate::AssembledMatricesComplex`].
//! 3. Invoke
//!    [`crate::ComplexInverseIterEigen::solve`](crate::solve::ComplexInverseIterEigen::solve)
//!    on the assembled `(K, M)` at the caller-supplied complex shift
//!    `σ` for the requested number of eigenvalues.
//!
//! Bit-for-bit, when the supplied [`MaterialDatabase`] returns purely
//! real `(ε, μ)` at `omega`, the linearised eigenvalues match the v0
//! free-space [`crate::FemEigenAssembly::assemble`] +
//! [`crate::InverseIterEigen`] path to solver tolerance.  This is the
//! load-bearing backward-compatibility invariant per ADR-0039 §4.
//!
//! ## References
//!
//! * `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
//!   §4 (theory anchor), §7 (public API surface).
//! * `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`
//!   step D4.
//! * `docs/src/decisions/0039-phase-4-fem-eig-1-dispersive-scope.md`.

use nalgebra::DVector;
use num_complex::Complex64;
use yee_core::units::{C0, EPS0, MU0};
use yee_mesh::{MaterialTag, TetMesh3D};

use crate::assembly::FemEigenAssembly;
use crate::material::MaterialDatabase;
use crate::solve::{ComplexInverseIterEigen, EigenpairListComplex, SparseEigenComplex};

/// Default bulk-material tag used by [`DispersiveSolver::solve_with_newton`]
/// when evaluating `ε(ω)` for the fixed-point self-consistency update.
///
/// Per ADR-0039 §4 the v1 Newton tracker drives the dispersive
/// self-consistency relation against the permittivity of a single
/// "bulk" material — i.e. the dispersive filler of the cavity. By
/// project convention every test fixture, validation driver, and
/// example tags the bulk filler as [`MaterialTag`] `0`; the constant
/// here pins that convention in the source so tag drift in a future
/// fixture is loud rather than silent. A more flexible
/// `solve_with_newton_with_bulk_tag(_, bulk_tag)` overload is left as
/// a Phase 4.fem.eig.1.1 extension if mixed-bulk Newton tracking ever
/// shows up as a validation requirement.
pub const BULK_TAG: MaterialTag = 0;

/// Number of inner eigenvalues requested per Newton step.
///
/// Mirrors the D4 fixture convention in `tests/dispersive_solve.rs`,
/// which uses `solve_at_frequency(_, _, 10, _)` at
/// `sigma_factor = 2.5` so the smallest-Re(k²) physical mode (TE_{101}
/// for a rectangular cavity) is reliably the first entry of the
/// ascending-Re-sorted output. Asking for fewer eigenvalues at this
/// shift lets the inner shift-invert lock onto a higher physical mode
/// whose `|k² − σ|` is smaller, which would feed the wrong
/// wavenumber-squared into the dispersive self-consistency update.
const NEWTON_INNER_NUM_EIGS: usize = 10;

/// Phase 4.fem.eig.1 dispersive eigenmode solver — wraps a
/// [`MaterialDatabase`] and a configured inner
/// [`ComplexInverseIterEigen`] for the linearised solve at each trial
/// frequency.
///
/// D4 ships [`Self::solve_at_frequency`] (one-trial linearised solve).
/// The outer Newton loop closing the fixed-point
/// `θ = (ω/c)²` condition lands as `Self::track_mode` in plan step D5;
/// D4 deliberately stops short of the Newton step so the linearised
/// path can be validated against fem-eig-001 in isolation.
///
/// ## Tuning
///
/// The default constructor [`Self::new`] uses `max_iter = 1000` and
/// `tol = 1e-8` — the same defaults as the v0
/// [`crate::InverseIterEigen`]. These are sufficient for the fem-eig-001
/// scale (~2 k interior DoFs) and for the unit-test fixtures in
/// `crates/yee-fem/tests/dispersive_solve.rs`.
#[derive(Debug, Clone)]
pub struct DispersiveSolver {
    /// Per-tet material database keyed by [`yee_mesh::MaterialTag`].
    pub material_db: MaterialDatabase,
    /// Per-mode iteration cap forwarded to the inner
    /// [`ComplexInverseIterEigen`].
    pub max_iter: usize,
    /// Relative Rayleigh-quotient convergence tolerance forwarded to
    /// the inner [`ComplexInverseIterEigen`].
    pub tol: f64,
    /// Outer Newton-tracker iteration cap used by
    /// [`DispersiveSolver::solve_with_newton`]. Defaults to `20` — well
    /// above the 3–5-iteration budget the spec §4.2 reports as typical
    /// for warm-started lightly-dispersive cavities, leaving headroom
    /// for cold-started or strongly-Drude-loaded fixtures.
    pub newton_max_iter: usize,
    /// Relative-step convergence tolerance for
    /// [`DispersiveSolver::solve_with_newton`]: the iteration returns
    /// when `|ω_{n+1} − ω_n| / |ω_n| < newton_tol`. Defaults to `1e-6`
    /// — same order as the spec §9 ±0.5 % / ±5 % Re/Im gates with
    /// three decades of headroom.
    pub newton_tol: f64,
}

impl DispersiveSolver {
    /// Build a dispersive solver from a [`MaterialDatabase`] with
    /// default inner-solver tuning (`max_iter = 1000`, `tol = 1e-8`).
    pub fn new(material_db: MaterialDatabase) -> Self {
        Self {
            material_db,
            max_iter: 1000,
            tol: 1e-8,
            newton_max_iter: 20,
            newton_tol: 1e-6,
        }
    }

    /// Build a dispersive solver with explicit inner-solver tuning.
    /// See [`crate::ComplexInverseIterEigen::new`] for the meaning of
    /// `max_iter` and `tol`. The outer-Newton parameters
    /// (`newton_max_iter`, `newton_tol`) remain at their defaults; set
    /// them on the returned value if a tighter / looser Newton budget
    /// is required.
    pub fn with_tuning(material_db: MaterialDatabase, max_iter: usize, tol: f64) -> Self {
        Self {
            material_db,
            max_iter,
            tol,
            newton_max_iter: 20,
            newton_tol: 1e-6,
        }
    }

    /// Solve the **linearised** complex generalised eigenproblem
    /// `K(ω) e = θ M(ω) e` at a single trial real angular frequency
    /// `omega` (rad/s), for the `num_eigs` eigenvalues nearest the
    /// complex shift `sigma`.
    ///
    /// The returned eigenvalues are linearised wavenumber-squareds
    /// `k² = θ`; the dispersive Newton fixed-point condition
    /// `θ = (ω/c)²` is **not** enforced here. That is the job of the
    /// outer Newton loop in plan step D5.
    ///
    /// # Arguments
    ///
    /// * `mesh` — tet mesh with per-tet
    ///   [`yee_mesh::TetMesh3D::tetrahedron_material`] tags consumed by
    ///   the stored [`MaterialDatabase`].
    /// * `omega` — trial angular frequency (rad/s). Real-valued; see
    ///   [`crate::FemEigenAssembly::assemble_complex`] for why a
    ///   complex `omega` is not supported in v1.
    /// * `num_eigs` — number of eigenvalues to request from the inner
    ///   solver. Must be ≥ 1 and ≤ interior DoF count.
    /// * `sigma` — complex shift `σ` for the shift-invert pencil
    ///   `(K − σM)`. Typical choice: `σ = (ω/c)²` (the analytic
    ///   free-space wavenumber-squared at the trial ω, which is also
    ///   the fixed-point target for D5's Newton step).
    ///
    /// # Errors
    ///
    /// Propagates [`yee_core::Error`] from
    /// [`crate::FemEigenAssembly::new`],
    /// [`crate::FemEigenAssembly::assemble_complex`], and
    /// [`ComplexInverseIterEigen::solve`]. The first surfaces mesh /
    /// material-database shape mismatches; the second guards against
    /// an empty mesh; the third reports inner sparse-LU failures or
    /// per-mode non-convergence in the iteration budget.
    pub fn solve_at_frequency(
        &self,
        mesh: &TetMesh3D,
        omega: f64,
        num_eigs: usize,
        sigma: Complex64,
    ) -> Result<EigenpairListComplex, yee_core::Error> {
        // The real `eps_r`, `mu_r` arrays are placeholders for the
        // complex path — `assemble_complex` ignores them and looks up
        // ε(ω), μ(ω) per-tet via the supplied MaterialDatabase. We
        // pass `1.0` everywhere so the assembler's length-check
        // succeeds; the values themselves are unused on this code
        // path.
        let n_tets = mesh.tetrahedra.len();
        let assembly = FemEigenAssembly::new(mesh, vec![1.0; n_tets], vec![1.0; n_tets])?;
        let assembled = assembly.assemble_complex(omega, &self.material_db)?;

        let solver = ComplexInverseIterEigen::new(self.max_iter, self.tol);
        solver.solve(&assembled.k, &assembled.m, num_eigs, sigma)
    }

    /// Outer fixed-point ("Newton") tracker — converge a single physical
    /// dispersive eigenmode from a real warm-start angular frequency.
    ///
    /// Iterates the dispersive self-consistency relation
    ///
    /// ```text
    ///     ω² ε_r(ω) μ_0 ε_0 = k_phys²(ω),
    /// ```
    ///
    /// where the FEM generalised eigenvalue
    /// `λ(ω) := K(ω) e / (M(ω) e)` already bakes `ε(ω)` into `M`, so at
    /// a self-consistent dispersive eigenmode `λ(ω*) = (ω*/c)²` —
    /// `ε(ω)` does **not** appear in the update rule. The fixed-point
    /// form is therefore
    ///
    /// ```text
    ///     ω_{n+1} = c · sqrt( λ(ω_n) ).
    /// ```
    ///
    /// This corrects a previous form `ω_{n+1}² = λ / (μ₀ε₀ε(ω))` that
    /// divided by `ε(ω)` a *second* time after `M` already accounted
    /// for it; that bug collapsed the converged `Re(f_FEM)` to
    /// `Re(f_analytic) / √ε_∞` on the fem-eig-002 cavity (Track
    /// QQQQQQQQ D6 finding 1 — see
    /// `crates/yee-fem/validation/README.md`). The bulk-material `ε(ω)`
    /// is still consumed by [`DispersiveEigenpair::k_complex`] to
    /// compose the converged physical wavenumber, and is looked up at
    /// tag [`BULK_TAG`].
    ///
    /// The returned [`DispersiveEigenpair`] carries the converged
    /// complex `k = ω · √(μ₀ ε₀ ε(ω))` and the M-orthonormalised
    /// complex eigenvector from the final inner solve.
    ///
    /// # Arguments
    ///
    /// * `mesh` — tet mesh consumed by the inner
    ///   [`Self::solve_at_frequency`] call.
    /// * `omega_0` — initial real-valued warm-start angular frequency
    ///   (rad/s). Typical sources: a Phase 4.0 free-space resonance
    ///   from [`crate::InverseIterEigen`] or the Pozar §3.1 analytic
    ///   cavity mode. The convergence basin is small for strongly
    ///   dispersive media; for the v1 validation regime (lightly-loaded
    ///   cavities) the v0 air resonance is comfortably inside the
    ///   monotone-convergence range.
    /// * `sigma_factor` — shift multiplier for the inner shift-invert
    ///   pencil. The inner solver targets eigenvalues nearest
    ///   `σ = sigma_factor · Re(ω_n)² / c²` (real-valued — keeps the
    ///   numerics conditioned even when `Im(ω) > 0`). Typical value
    ///   `2.5` — matches the
    ///   `dispersive_solve.rs` D4 fixture convention and sits between
    ///   the 8th and 9th physical modes on the canonical WR-90 mesh.
    ///
    /// # Errors
    ///
    /// * [`DispersiveError::Underlying`] wraps any
    ///   [`yee_core::Error`] surfaced by the inner
    ///   [`Self::solve_at_frequency`] (mesh / material-database shape
    ///   mismatch; inner sparse-LU failure; inner inverse-power
    ///   non-convergence).
    /// * [`DispersiveError::NewtonDidNotConverge`] when the outer
    ///   fixed-point iteration reaches [`Self::newton_max_iter`]
    ///   without the relative step `|Δω/ω|` falling below
    ///   [`Self::newton_tol`]. Carries the last iterate and the last
    ///   residual for caller diagnosis.
    pub fn solve_with_newton(
        &self,
        mesh: &TetMesh3D,
        omega_0: Complex64,
        sigma_factor: f64,
    ) -> Result<DispersiveEigenpair, DispersiveError> {
        let c0 = yee_core::units::C0;
        let mut omega = omega_0;
        let mut last_residual = f64::INFINITY;
        let mut last_k_sq = Complex64::new(f64::NAN, f64::NAN);

        for _iter in 0..self.newton_max_iter {
            let omega_re = omega.re;
            // Real shift for numerical conditioning per ADR-0039
            // §"Hellmann–Feynman transposed vs. Hermitian" — keeps the
            // shift-invert factor real-positive even when the running
            // ω develops a small imaginary part from a lossy material.
            let k0_sq_re = (omega_re / c0).powi(2);
            let sigma = Complex64::new(sigma_factor * k0_sq_re, 0.0);

            // Inner shift-invert returns eigenvalues nearest σ, sorted
            // ascending by Re(k²). The smallest physical mode (TE_{101}
            // for a rectangular cavity) is the lowest of the first ten
            // — this matches the D4 fixture convention in
            // `tests/dispersive_solve.rs` which uses
            // `solve_at_frequency(_, _, 10, _)` and takes `pairs.k[0]`.
            // Asking for fewer than ten eigenvalues at `sigma_factor =
            // 2.5` lets the inner solver lock onto a higher physical
            // mode whose `|k² − σ|` is smaller, which is the wrong
            // eigenvalue for the dispersive self-consistency relation.
            let num_eigs_inner = NEWTON_INNER_NUM_EIGS;
            let pairs = self
                .solve_at_frequency(mesh, omega_re, num_eigs_inner, sigma)
                .map_err(DispersiveError::Underlying)?;

            // Pick the smallest **physical** mode — `Re(k²) > 0`. The
            // inner shift-invert can return spurious near-null-space
            // gradient modes with `Re(k²) ≤ 0` ahead of the lowest
            // physical TE / TM mode when the trial ω sits well below
            // the cavity resonance and the shift-invert pencil is far
            // from the physical spectrum. Match the D4 fixture
            // convention from `tests/dispersive_solve.rs` —
            // `pairs.k[0]` is fine when the warm-start ω is close to
            // the physical resonance, but a cold-started Newton step
            // (warm-start an order of magnitude off) needs the
            // gradient-mode filter. Threshold `Re(k²) > 0` is
            // permissive: the physical modes carry `Re(k²) ≫ 0` at
            // any reasonable mesh.
            let (k_sq_lin, e_col) = pairs
                .k
                .iter()
                .zip(pairs.e.column_iter())
                .find(|(k_sq, _)| k_sq.re > 0.0)
                .ok_or_else(|| {
                    DispersiveError::Underlying(yee_core::Error::Numerical(format!(
                        "DispersiveSolver::solve_with_newton: inner solver \
                             returned no eigenvalue with Re(k²) > 0 at trial \
                             ω = {omega_re} rad/s (sigma = {sigma}); all returned \
                             k² are non-physical / spurious gradient modes"
                    )))
                })?;
            let k_sq_lin = *k_sq_lin;
            let e_vec = DVector::from_column_slice(e_col.as_slice());

            // Fixed-point update: ω_{n+1} = c · √λ, equivalently
            // ω_{n+1}² = c² · λ. The FEM generalised eigenvalue
            // `λ = k_sq_lin` returned by `solve_at_frequency` is
            // `(ω_phys/c)²` at a self-consistent dispersive eigenmode
            // because `ε(ω)` is already baked into the `M` matrix at
            // assembly time. The earlier `ω² = λ/(μ₀ε₀ε(ω))` form
            // divided by `ε(ω)` a *second* time and collapsed the
            // converged `Re(f_FEM)` to `Re(f_analytic)/√ε_∞` — see
            // `crates/yee-fem/validation/README.md` "fem-eig-002 D6
            // finding 1" for the bug history (Track QQQQQQQQ).
            let c_sq = Complex64::new(C0 * C0, 0.0);
            let omega_sq_new = c_sq * k_sq_lin;
            let omega_new = omega_sq_new.sqrt();

            // Convergence: relative step on |ω|.
            let denom_norm = omega.norm();
            let residual = if denom_norm > 0.0 {
                (omega_new - omega).norm() / denom_norm
            } else {
                (omega_new - omega).norm()
            };
            last_residual = residual;
            last_k_sq = k_sq_lin;

            if residual < self.newton_tol {
                // Compose the physical complex k for the converged ω.
                // k = ω · √(μ₀ ε₀ ε(ω_new.re))
                let eps_final = self.material_db.eps_at(BULK_TAG, omega_new.re);
                let k_complex = omega_new * (Complex64::new(MU0 * EPS0, 0.0) * eps_final).sqrt();
                return Ok(DispersiveEigenpair {
                    omega: omega_new,
                    k_complex,
                    e_vec,
                });
            }

            omega = omega_new;
        }

        Err(DispersiveError::NewtonDidNotConverge {
            last_omega: omega,
            last_k_sq,
            last_residual,
        })
    }
}

/// Output of [`DispersiveSolver::solve_with_newton`] — a single
/// converged dispersive eigenmode.
///
/// All fields are populated only on a successful return; the failure
/// path carries its own last-iterate / last-residual diagnostics inside
/// the [`DispersiveError::NewtonDidNotConverge`] variant.
#[derive(Debug, Clone)]
pub struct DispersiveEigenpair {
    /// Converged complex angular frequency (rad/s). `Re(omega)` is the
    /// physical resonance; `Im(omega) < 0` encodes the lossy-mode decay
    /// rate under the engineering `exp(+jωt)` convention.
    pub omega: Complex64,
    /// Physical complex wavenumber at the converged frequency:
    /// `k = ω · √(μ₀ ε₀ ε_r(ω))`. The complex `k²` satisfies the
    /// dispersive self-consistency relation
    /// `k² = ω² ε(ω) μ₀ ε₀` to within
    /// [`DispersiveSolver::newton_tol`].
    pub k_complex: Complex64,
    /// M-orthonormalised complex eigenvector on the interior-DoF basis
    /// at the converged frequency.
    pub e_vec: DVector<Complex64>,
}

/// Error type for [`DispersiveSolver::solve_with_newton`] — wraps the
/// underlying [`yee_core::Error`] surface and adds a
/// fixed-point-non-convergence variant.
///
/// A dedicated `yee-fem`-local error type is used here rather than a
/// new variant on `yee_core::Error` to keep the
/// `crates/yee-fem/src/dispersive.rs` lane self-contained (per the D5
/// agent brief escape-hatch — see ADR-0039 §"Material relocation" for
/// the cross-lane lift policy). Lifting this surface into
/// `yee_core::Error::NewtonDidNotConverge` is a follow-up task once a
/// second consumer of the variant appears in the workspace.
#[derive(Debug, thiserror::Error)]
pub enum DispersiveError {
    /// An error from the inner [`DispersiveSolver::solve_at_frequency`]
    /// call (mesh / material-database shape mismatch, inner sparse-LU
    /// failure, inner inverse-power non-convergence).
    #[error("dispersive inner solve failed: {0}")]
    Underlying(#[from] yee_core::Error),

    /// The outer fixed-point ("Newton") iteration reached
    /// [`DispersiveSolver::newton_max_iter`] without satisfying the
    /// relative-step convergence test `|Δω/ω| < newton_tol`.
    ///
    /// Carries the last iterate, the last inner linearised
    /// wavenumber-squared, and the last relative residual for caller
    /// diagnosis. The standard remediation when this fires on a
    /// production workload is to (a) increase
    /// [`DispersiveSolver::newton_max_iter`], (b) widen
    /// [`DispersiveSolver::newton_tol`], or (c) supply a closer
    /// warm-start `omega_0` to
    /// [`DispersiveSolver::solve_with_newton`].
    #[error("Newton ω-tracker did not converge within {last_residual:e} (last ω = {last_omega})")]
    NewtonDidNotConverge {
        /// Last trial angular frequency in the iteration (rad/s).
        last_omega: Complex64,
        /// Last linearised inner-solver eigenvalue at `last_omega`.
        last_k_sq: Complex64,
        /// Last relative step `|ω_{n+1} − ω_n| / |ω_n|`.
        last_residual: f64,
    },
}
