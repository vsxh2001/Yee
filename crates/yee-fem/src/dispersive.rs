//! Phase 4.fem.eig.1 dispersive eigenmode tracker — Newton-Raphson
//! frequency tracking of `K(ω) e = (ω/c)² M(ω) e` with complex
//! `ε(ω)`, `μ(ω)` per tet.
//!
//! ## D4 scope (this file)
//!
//! [`DispersiveSolver::solve_at_frequency`] solves the **linearised**
//! complex generalised eigenproblem `K(ω₀) e = θ M(ω₀) e` at a single
//! trial angular frequency ω₀ via the complex sparse inverse-power
//! eigensolver from `crates/yee-fem/src/solve.rs` (D2). The returned
//! eigenvalues are linearised wavenumber-squareds at the trial ω; they
//! are **not yet** self-consistent dispersive eigenmodes — the Newton
//! outer loop in plan step D5 closes the fixed-point condition
//! `θ = (ω/c)²`.
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

use num_complex::Complex64;
use yee_mesh::TetMesh3D;

use crate::assembly::FemEigenAssembly;
use crate::material::MaterialDatabase;
use crate::solve::{ComplexInverseIterEigen, EigenpairListComplex, SparseEigenComplex};

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
}

impl DispersiveSolver {
    /// Build a dispersive solver from a [`MaterialDatabase`] with
    /// default inner-solver tuning (`max_iter = 1000`, `tol = 1e-8`).
    pub fn new(material_db: MaterialDatabase) -> Self {
        Self {
            material_db,
            max_iter: 1000,
            tol: 1e-8,
        }
    }

    /// Build a dispersive solver with explicit inner-solver tuning.
    /// See [`crate::ComplexInverseIterEigen::new`] for the meaning of
    /// `max_iter` and `tol`.
    pub fn with_tuning(material_db: MaterialDatabase, max_iter: usize, tol: f64) -> Self {
        Self {
            material_db,
            max_iter,
            tol,
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
}
