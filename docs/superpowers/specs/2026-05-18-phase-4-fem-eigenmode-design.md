# Phase 4 — 3-D FEM Eigenmode Solver for Resonant Cavities

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.0 (walking skeleton); 4.fem.eig.1+ deferred (see §13)
**Depends on:** Phase 1.3.1.1 (2-D Nedelec cross-section eigensolver, shipped), Phase 1.mesh.0/1 (Gmsh FFI + KiCad import, shipped)
**Blocks:** Phase 4.fem.eig.2 (lossy / Q-factor extraction), Phase 4.fem.eig.3 (dielectric-resonator antenna validation), Phase 4 driven-FEM solver (separate sub-project)

## 1. Motivation

The shipped Yee solver portfolio leaves a real gap for resonant 3-D problems:

- **FDTD** can extract resonant modes of an arbitrary closed cavity by exciting with a wideband pulse and Fourier-analyzing the ringdown — but that workflow needs many time-steps to resolve closely-spaced modes, fails outright on lossless cavities (the ringdown never decays), and gives only a noisy estimate of Q-factor for lightly-loaded cavities. Mode-by-mode extraction beyond the lowest few resonances is impractical.
- **Planar MoM** is the wrong tool entirely. It is a surface integral-equation method on conducting boundaries embedded in a stratified open background — it has no concept of a volumetric resonant interior.
- The 2-D Nedelec cross-section eigensolver from Phase 1.3.1.1 finds *guided-wave* modes on a translationally-invariant cross-section (`E ∝ e^{-jβz}`), not the *resonant* modes of a 3-D enclosure.

A 3-D FEM eigenmode solver fills this gap directly. The use cases the project has on its roadmap that require it:

- Resonant-cavity filter design (loaded waveguide cavities, combline / iris-coupled bandpass).
- Q-factor extraction for lossy or dielectric-loaded cavities via complex eigenvalues (`Q = ω' / (2 ω'')`).
- Multi-mode dielectric-resonator antenna (DRA) analysis, where the dielectric puck and the surrounding air share a single closed FEM domain.
- Particle-accelerator cavity design (RF gun / pillbox / elliptical cell), called out as an "application pack" in `ROADMAP.md` Phase 4.

The mathematics — Maxwell eigenproblem in weak form on curl-conforming Nedelec edge elements — is textbook (Jin 3rd ed. Ch. 9–10, Boffi-Brezzi-Demkowicz 2013). The engineering risk is sparse generalized-eigensolve plumbing in pure Rust at production scale, not physics.

## 2. Non-goals (Phase 4.fem.eig.0 walking skeleton)

Explicitly out of scope for v0:

- **Driven 3-D FEM** (right-hand-side from wave-port, current density, or incident plane wave). That is a separate Phase 4 sub-project sharing the same `K` assembly but a different solve path.
- **Open-region radiators.** FEM eigenmode is for closed cavities; open antennas remain MoM's job. Coupling an FEM interior to an MoM exterior (FEM-BEM hybrid) is Phase 4+ and not promised here.
- **Adaptive remeshing / mesh refinement.** The mesh is an input contract from `yee-mesh`.
- **Moving meshes / time-domain FEM.**
- **GPU acceleration.** v0 is CPU/FP64 scalar; a `Backend`-style abstraction is left for v1 if/when sparse eigensolve becomes the bottleneck.
- **Higher-order Nedelec elements.** v0 ships first-order (lowest-order Whitney-1) edge elements only. Hierarchical p-refinement is 4.fem.eig.1+.
- **Periodic / Floquet boundary conditions.** v0 ships PEC and PMC (natural Neumann) only. Phase-shift periodic BCs are 4.fem.eig.4+ if user demand materializes.
- **Tensor / anisotropic / dispersive permittivity inside the cavity.** Scalar isotropic complex `ε_r`, `μ_r` per tetrahedron only.

## 3. Scope decision — Phase 4.fem.eig.0 walking skeleton

Walking-skeleton-first per `CLAUDE.md` §3. The v0 deliverable:

- Tetrahedral mesh `TetMesh3D` of a closed 3-D domain, supplied by the caller.
- First-order Nedelec curl-conforming edge elements (one DoF per tet edge, six DoFs per tet).
- Assemble the generalized eigenproblem `K · e = k² · M · e`, where
  - `K_ij = ∫_Ω (1/μ_r) (∇×N_i) · (∇×N_j) dV` — curl-curl stiffness,
  - `M_ij = ∫_Ω ε_r N_i · N_j dV` — vector mass.
- PEC tangential-`E`-zero Dirichlet on boundary edges (row/column elimination).
- Sparse generalized eigensolve via shift-invert Arnoldi targeting smallest physical eigenvalues, ten lowest eigenpairs extracted.
- One end-to-end gate: rectangular metallic air-filled cavity, lowest mode matches Pozar TE_{101} within ±0.3% (§9, fem-eig-001).
- CPU, FP64, scalar Rust. No GPU. No higher-order. No dielectric loss.

What v0 does **not** ship, deferred to 4.fem.eig.1+ (see §13):

- Production-scale sparse eigensolve > 100k DoFs.
- Lossy / complex-`ε_r` Q-factor extraction.
- Higher-order Nedelec.
- Periodic BCs.
- GPU backend.

## 4. Theory anchor

The source-free Maxwell curl equations in a closed lossless domain `Ω` with PEC boundary `∂Ω` reduce, after eliminating `H`, to the vector wave equation

```
∇ × ( (1/μ_r) ∇ × E )  =  k₀² ε_r E      in Ω
n × E  =  0                              on ∂Ω
```

with `k₀² = ω² ε₀ μ₀ = (ω/c)²`. The variational (weak) form, multiplying by a test field `v ∈ H₀(curl; Ω)` and integrating by parts, is:

> Find `(k₀², E) ∈ (ℝ₊, H₀(curl; Ω))` such that for all `v ∈ H₀(curl; Ω)`,
>
> ```
> ∫_Ω (1/μ_r) (∇×E) · (∇×v) dV  =  k₀² ∫_Ω ε_r E · v dV.
> ```

Discretizing `E ≈ Σ_j e_j N_j(x)` with `{N_j}` a basis of `H₀(curl; Ω)` yields the matrix eigenproblem `K e = k₀² M e`.

**Why Nedelec edge elements** — and not nodal Lagrange. Nodal `H¹` discretizations of this problem are well-known to admit *spurious modes* whose `∇×` vanishes but which are not gradients of `H¹` scalars (i.e. the discrete kernel of the curl operator is too large). Nedelec curl-conforming elements (Whitney 1, lowest order) are the canonical fix: their kernel is exactly `∇ H¹`, so the discrete spurious modes coincide with the continuous gauge freedom and cluster harmlessly at `k₀² = 0`. Shift-invert targeting `σ > 0` skips them automatically. This is the same argument that motivated the 2-D Nedelec basis in Phase 1.3.1.1 (see `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md` §Approach); the 3-D Whitney-1 basis is its direct generalization.

References for §4:

- Jin, *The Finite Element Method in Electromagnetics*, 3rd ed., Wiley 2014, Ch. 9 (Nedelec elements on tetrahedra) and Ch. 10 (eigenvalue problems).
- Boffi, Brezzi, Fortin, *Mixed Finite Element Methods and Applications*, Springer 2013, Ch. 7 (curl-conforming discretizations and the discrete compactness property).
- Demkowicz, *Computing with hp-Adaptive Finite Elements*, vol. 2, CRC 2007, Ch. 4 (de Rham complex on tetrahedra).
- Bossavit, *Computational Electromagnetism*, Academic Press 1998 (Whitney forms, geometric origin of edge elements).

## 5. Mesh requirements — `TetMesh3D`

New type in `yee-mesh`, paralleling the existing `TriMesh2D` from Phase 1.3.1.1:

- `vertices: Vec<Point3<f64>>` — node coordinates in metres.
- `tetrahedra: Vec<[u32; 4]>` — four vertex indices per tet, ordered so the signed volume `(1/6) (v₁−v₀) · ((v₂−v₀) × (v₃−v₀))` is positive (consistent CCW orientation per the right-hand rule). Construction validates orientation and rejects degenerate / inverted tets.
- `vertex_material: Option<Vec<MaterialTag>>` — optional per-vertex tag (PEC marker for boundary classification).
- `tet_material: Option<Vec<MaterialTag>>` — optional per-tet tag (dielectric region).
- Methods: `signed_volume(tet) -> f64`, `centroid(tet) -> Point3<f64>`, `edges() -> impl Iterator<Item = [u32; 2]>` returning the global edge list with a canonical orientation (lower-index endpoint first), and `boundary_edges() -> &[EdgeId]` returning edges lying on `∂Ω` for Dirichlet application.
- Construction: `TetMesh3D::new(vertices, tetrahedra, vertex_material, tet_material) -> Result<Self, yee_core::Error>` validates non-degeneracy, orientation, and tag-length consistency.

Mesh sources:

- **Hand-rolled fixtures** for the v0 gate (rectangular cavity decomposed into ≤ 6 tets per cell on a regular grid — well-conditioned and easy to refine).
- **Gmsh-generated `.msh` import** via the existing `gmsh` feature on `yee-mesh`. Gmsh's `tetgen`-backed mesher produces orientation-consistent tets and second-order surface elements are flattened to first-order for v0.

`TetMesh3D` is purely a geometric/topological container; no FEM logic lives in `yee-mesh`. The eigensolver consumes it via reference.

## 6. Assembly + solver pipeline

End-to-end pipeline, mirroring the structure of `crates/yee-mom/src/eigensolver/{mesh,assembly,solve}.rs` from Phase 1.3.1.1:

1. **Edge enumeration.** Walk `TetMesh3D::tetrahedra`, build a global edge list with canonical orientation. Each tet contributes six local edges `(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)`. Local→global edge map stored once per assembly.
2. **Local element matrices.** For each tet:
   - Compute the 4×4 matrix of gradients `∇λ_i` of the barycentric (linear nodal) basis functions.
   - The Nedelec basis on edge `(i,j)` is `N_{ij} = λ_i ∇λ_j − λ_j ∇λ_i`; its curl is the constant vector `∇λ_i × ∇λ_j − ∇λ_j × ∇λ_i = 2 ∇λ_i × ∇λ_j` (Jin §9.4, eq. 9.43).
   - Local stiffness `K^e_{αβ} = (1/μ_r,e) · V_e · (∇×N_α) · (∇×N_β)` (the curls are constant per tet, so the integral is exact).
   - Local mass `M^e_{αβ} = ε_r,e · ∫_{tet} N_α · N_β dV`, computed with 4-point Gauss quadrature on the reference tet (exact for the integrand's polynomial degree).
3. **Global assembly.** Scatter local 6×6 blocks into a sparse `CsrMatrix<f64>` (or `faer::sparse::SparseColMat`, decision in §8). Signs are tracked via the canonical edge orientation: when a tet's local edge runs against the global orientation, the corresponding row/column of the local matrix is negated.
4. **Boundary conditions.** For each edge in `boundary_edges()`, eliminate the row and column from both `K` and `M`, *or* apply the standard penalty / row-replacement technique (TBD: penalty avoids resizing the index map at the cost of one large eigenvalue per Dirichlet edge that shift-invert skips trivially; row elimination keeps spectra clean). v0 uses row-elimination.
5. **Generalized sparse eigensolve.** Target the ten smallest *positive* eigenvalues of `K e = k² M e` via shift-invert Arnoldi at shift `σ = (k₀_lowest_expected)² · 0.5` (i.e. below the smallest physical mode but above the gradient-kernel cluster at `0`). Concretely solve `(K − σ M)^{-1} M e = θ e` with `k² = σ + 1/θ`. Sparse LU of `(K − σ M)` is provided by `faer::sparse::FaerLuSolver` or equivalent.
6. **Post-processing.** Convert eigenvalues to resonant frequencies via `f_n = c · √(k_n²) / (2π)`. Eigenvectors `e_n` are real per-edge coefficients; the modal `E`-field at an arbitrary point is reconstructed via the Nedelec interpolation `E(x) = Σ_j e_j N_j(x)` and is what callers visualise.

Spurious-mode handling: the gradient kernel of `∇×` is non-empty (any `E = ∇φ` for `φ ∈ H₀¹(Ω)` has `∇×E = 0`). After PEC elimination it manifests as roughly `(N_vertices − 1)` eigenpairs clustered at `k² ≈ 0`. Shift-invert at `σ > 0` simply does not see them; they are correctly identified as null-space modes.

## 7. API sketch

```rust
//! crates/yee-mesh/src/lib.rs — extend with TetMesh3D

/// 3-D tetrahedral mesh for FEM eigenmode + driven FEM solvers.
///
/// Mirrors [`TriMesh2D`] in invariants and tag conventions; sources include
/// hand-rolled fixtures and Gmsh `.msh` import via the `gmsh` feature.
pub struct TetMesh3D {
    pub vertices: Vec<Point3<f64>>,
    pub tetrahedra: Vec<[u32; 4]>,
    pub vertex_material: Option<Vec<MaterialTag>>,
    pub tet_material: Option<Vec<MaterialTag>>,
}

impl TetMesh3D {
    pub fn new(
        vertices: Vec<Point3<f64>>,
        tetrahedra: Vec<[u32; 4]>,
        vertex_material: Option<Vec<MaterialTag>>,
        tet_material: Option<Vec<MaterialTag>>,
    ) -> Result<Self, yee_core::Error>;

    pub fn signed_volume(&self, tet_idx: usize) -> f64;
    pub fn centroid(&self, tet_idx: usize) -> Point3<f64>;
    pub fn edges(&self) -> &[EdgeKey];
    pub fn boundary_edges(&self) -> &[usize]; // indices into edges()
}

//! New crate or module — TBD whether this lives in `yee-mom::eigensolver`
//! (alongside the 2-D solver) or in a new `yee-fem` crate. See §11 and §14.

pub struct FemEigenAssembly<'m> {
    mesh: &'m TetMesh3D,
    eps_r: HashMap<TetId, Complex64>, // v0: imag part must be zero
    mu_r:  HashMap<TetId, Complex64>, // v0: imag part must be zero
}

impl<'m> FemEigenAssembly<'m> {
    pub fn new(
        mesh: &'m TetMesh3D,
        eps_r: HashMap<TetId, Complex64>,
        mu_r: HashMap<TetId, Complex64>,
    ) -> Result<Self, yee_core::Error>;

    /// Assemble the global K (curl-curl stiffness) and M (vector mass)
    /// sparse matrices, with PEC tangential-E-zero Dirichlet applied.
    pub fn assemble(&self) -> Result<(SparseCsr64, SparseCsr64), yee_core::Error>;
}

pub struct FemEigenSolver { /* shift, num_eigs, max_iter, tol */ }

impl FemEigenSolver {
    pub fn new(num_eigs: usize) -> Self;
    pub fn with_shift(self, sigma_k2: f64) -> Self;

    /// Solve K e = k² M e for the `num_eigs` smallest positive eigenvalues.
    pub fn solve(
        &self,
        k: &SparseCsr64,
        m: &SparseCsr64,
    ) -> Result<EigenpairList, yee_core::Error>;
}

pub struct EigenpairList {
    /// Resonant wavenumbers (k = ω/c), one per mode, sorted ascending.
    pub k: Vec<f64>, // v0: real. Complex in 4.fem.eig.2 for lossy cavities.
    /// Mode-coefficient vectors stacked column-wise: e[:, n] is mode n.
    pub e: DMatrix<f64>,
}
```

End-to-end example (informally):

```rust
let mesh = TetMesh3D::rectangular_cavity(a, b, d, nx, ny, nz)?;
let eps = HashMap::from_iter(mesh.tets().map(|t| (t, Complex64::new(1.0, 0.0))));
let mu  = HashMap::from_iter(mesh.tets().map(|t| (t, Complex64::new(1.0, 0.0))));
let asm = FemEigenAssembly::new(&mesh, eps, mu)?;
let (k_mat, m_mat) = asm.assemble()?;
let pairs = FemEigenSolver::new(10).solve(&k_mat, &m_mat)?;
let f1 = SPEED_OF_LIGHT * pairs.k[0] / (2.0 * PI);
// f1 should match Pozar TE_{101} = (c/2)·sqrt((1/a)² + (1/d)²) to ±0.3%.
```

## 8. Sparse-eigen library decision

The 2-D Phase 1.3.1.1 solver escape-hatched to a **dense `nalgebra::SymmetricEigen`** at ≤ 500 DoFs because `arpack-rs` is unmaintained and pulls a system LAPACK that the CI lint-clean policy rejects. The 3-D problem has dramatically more DoFs (≈ 7 edges per tet × 10³–10⁶ tets), so dense fallback is **not** viable in general. Survey of pure-Rust / Rust-binding options as of 2026-05:

| Option | Pros | Cons |
|---|---|---|
| `arpack-rs` (FORTRAN binding to ARPACK) | Mature, battle-tested shift-invert Arnoldi, gold-standard for sparse generalized eigenproblems. | Unmaintained crate; pulls system LAPACK; CI-lint-policy conflict per the Phase 1.3.1.1 escape hatch. |
| `lobpcg` (pure-Rust LOBPCG) | Pure-Rust, MIT, no system deps, works for symmetric generalized problems with a preconditioner. | LOBPCG converges to a few extreme eigenvalues; shift-invert requires a sparse LU preconditioner we'd build on top. Less proven on indefinite shift-invert spectra. |
| `faer` sparse eigensolve | Same ecosystem as our dense LU. | As of 2026-05 `faer` provides sparse LU and symmetric *dense* eigen but **no native generalized sparse eigen** (TBD: re-check at implementation time; `faer` is moving quickly). |
| Hand-rolled inverse-power iteration on `faer` sparse LU | Pure-Rust, no new dep, swap-point trait exposed. | One mode at a time + deflation; convergence rate sensitive to spectrum clustering near the shift. |
| `slepc-rs` / PETSc binding | Production-grade Krylov-Schur. | Heavy native dependency; out of character with the rest of the workspace; opens an MPI surface we don't otherwise want at v0. |

**Recommendation for v0:** ship a `SparseEigen` trait that abstracts the solve, with two concrete impls behind feature flags:

- `LobpcgEigen` (default, pure-Rust, no system deps) — `lobpcg` crate plus a `faer` sparse LU as the inner preconditioner for shift-invert.
- `ArpackEigen` (feature `eig-arpack`, off by default) — if a future maintainer revives an `arpack-rs`-equivalent binding or PETSc/SLEPc becomes acceptable.

The trait is the load-bearing decision; the library is the swap point. This mirrors the `yee_cuda::backend::Backend` trait that exists precisely so the cudarc choice is not load-bearing (`TECH_STACK.md`).

`TBD: confirm at implementation time whether faer has added native sparse generalized eigen — if so, prefer it for ecosystem coherence.`

## 9. Validation gate — fem-eig-001 rectangular metallic cavity

Canonical Pozar §6.3 example. Lossless air-filled rectangular metallic cavity with dimensions

- `a = 22.86 mm` (WR-90 broad wall),
- `b = 10.16 mm` (WR-90 narrow wall),
- `d = 30 mm` (cavity length).

The analytic TE_{mnp} resonant frequencies are

```
f_{mnp} = (c/2) · √( (m/a)² + (n/b)² + (p/d)² ).
```

The fundamental TE_{101} mode resonates at **`f_{101} ≈ 9.660 GHz`** for this geometry (Pozar 4th ed. eq. 6.42).

**Gate criteria:**

1. Lowest extracted eigenvalue `k_1` satisfies `|f_FEM − 9.660 GHz| / 9.660 GHz ≤ 0.3%`. The tolerance is tighter than the Phase 1.3.1.1 TE10 cross-section gate (0.1%) because the 3-D problem has more discretization error per DoF; matching Pozar's 4-significant-digit tabulation is the target.
2. The lowest ten eigenvalues, sorted ascending, match the analytic Pozar TE/TM table for this geometry within ±1% pairwise (mode-by-mode RMS error; ordering must agree).
3. No spurious mode appears below the analytic TE_{101} (i.e. shift-invert correctly skips the gradient-kernel cluster at `k² ≈ 0`).
4. The eigenmode-001 binary runs end-to-end in `< 60 s` on a ~30k-edge mesh in `--release`. (Informational; not a CI gate.)
5. Standard verification chain green: `cargo build`, `cargo clippy -- -D warnings`, `cargo test --release`, `cargo fmt --check` on the touched crates.

## 10. Higher-applications roadmap

Beyond fem-eig-001, the eigensolver targets two further validation cases before Phase 4.fem.eig is called done:

- **fem-eig-002 — lossy-cavity Q-factor.** Same rectangular cavity, walls modeled with a finite surface resistance (alternatively: bulk complex `ε_r = ε_r' − j ε_r''` in a thin lossy dielectric liner). Extract Q via `Q = Re(ω) / (2 Im(ω))` from the complex eigenvalues, compare to Pozar §6.3 closed-form Q for the dominant mode (TE_{101} loaded Q from wall losses). Tolerance ±5% on Q. Requires complex `ε_r` end-to-end and a complex generalized eigensolve (which by linearity reduces to the same sparse Arnoldi machinery applied to a complex `(K − σM)`).
- **fem-eig-003 — cylindrical dielectric-resonator antenna in air-filled outer cavity.** Validation case from Petosa, *Dielectric Resonator Antenna Handbook*, ch. 3: cylindrical ε_r = 9.8 puck inside a metallic outer cylinder; extract the lowest HEM_{11δ} and TE_{01δ} modes and compare to Petosa's tabulated frequencies. Tolerance ±2% on resonance; exercises the full piecewise-`ε_r` path.

These are deferred to Phase 4.fem.eig.2 and 4.fem.eig.3 respectively (§13).

## 11. Risks and open questions

- **Spurious modes from the gradient kernel of `∇×`.** Mitigated by Nedelec basis (kernel is exactly the discrete gradient subspace), shift-invert targeting `σ > 0`, and a sanity-check post-filter that rejects any returned mode with `‖∇·(ε_r E)‖ > tol`. Risk: shift `σ` chosen too small accidentally captures gradient modes; mitigation is to set `σ ≥ 0.1 · k_lowest_expected²`.
- **Pure-Rust sparse generalized eigensolve maturity.** As of 2026-05 no single crate is the obvious choice (§8). The `SparseEigen` trait isolates this so a library swap is one-PR. Worst-case fallback is hand-rolled inverse-power iteration on a `faer` sparse LU; documented as the escape hatch.
- **Mesh-quality sensitivity.** Cavity-mode accuracy is known to be sensitive to tet aspect ratio; Pozar TE_{101} on a sliver-tet mesh can show 5%+ error even with first-order Nedelec correct. Gmsh's default mesher with `MeshSizeFactor ≤ 0.05·λ` produces acceptable meshes for the gate. Document in `validation/README.md`.
- **Memory cost.** At ~1 M edges (a realistic DRA mesh) the sparse `K` and `M` each hold ~14·10⁶ nnz (six edges per tet, ~7 nonzeros per edge column), which is ≈ 0.5 GB per matrix in CSR FP64. Sparse LU at that scale is the bottleneck; LOBPCG with an algebraic-multigrid preconditioner would scale better but is out of scope for v0. `TBD: confirm whether `lobpcg` crate supports a custom preconditioner trait or expects a static matrix.`
- **FFI surface for Gmsh tet meshing.** `yee-mesh`'s Gmsh FFI is feature-gated and well-tested for triangular surface meshes. The 3-D path uses Gmsh's `tetgen` algorithm, which is already in the bundled SDK — no new FFI symbols expected, but quality flags (`Mesh.Algorithm3D`, `Mesh.Optimize`) need to be exposed.
- **First-order Nedelec accuracy floor.** First-order edge elements converge at `O(h²)` for eigenvalues but the constant is large; matching Pozar's 4-significant-digit table at the fem-eig-001 mesh size is tight. If the gate proves marginal, the right next step is hierarchical `p`-refinement (4.fem.eig.1), not `h`-refinement to absurd mesh counts.

## 12. Dependencies

- **`yee-mesh` extension** — new `TetMesh3D` type (parallel to `TriMesh2D`), Gmsh tet-meshing exposure behind the existing `gmsh` feature flag. No new external dep.
- **Sparse linear algebra** — `faer` sparse LU (already in the workspace) for the inner shift-invert solve.
- **Sparse eigensolve** — `lobpcg` crate (new direct dep, pure-Rust, MIT) for the v0 outer Krylov iteration. `arpack-rs` remains a non-default alternative gated behind feature `eig-arpack`.
- **`yee-core`** — error type, units, no API changes expected.
- **Gmsh SDK** — only when the caller uses `.msh` import; the hand-rolled-cavity fixture path is feature-independent.

No strict ordering constraint relative to other Phase 4 sub-projects: 4.fem.eig.0 stands alone behind its own walking-skeleton gate.

## 13. Phase numbering ladder

- **Phase 4.fem.eig.0** — walking skeleton (this spec): rectangular metallic cavity, first-order Nedelec, lossless, fem-eig-001 passes.
- **Phase 4.fem.eig.1** — production-scale sparse eigensolve: higher-order Nedelec (hierarchical `p ≤ 2`), preconditioned LOBPCG or revived ARPACK binding, ≥ 100k-edge meshes, performance budget published.
- **Phase 4.fem.eig.2** — lossy-cavity Q-factor (fem-eig-002): complex `ε_r`, complex eigenvalues, validated against Pozar wall-loss Q.
- **Phase 4.fem.eig.3** — dielectric-resonator antenna (fem-eig-003): piecewise dielectric, Petosa DRA Handbook validation.
- **Phase 4.fem.eig.4+** — periodic / Floquet BCs, GPU sparse eigensolve, FEM-BEM hybrid for open radiators. Open-ended.

## 14. Lane

Spec file:

```
docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md
```

Implementation lane (declared here for the follow-up plan, not edited by this spec):

- `crates/yee-mesh/**` — new `TetMesh3D` type, optional Gmsh `.msh` tet-mesh ingestion.
- `crates/yee-fem/**` *(new crate)* — `FemEigenAssembly`, `FemEigenSolver`, `EigenpairList`, `SparseEigen` trait. Rationale for a new crate (rather than extending `yee-mom::eigensolver`): the FEM solver is a peer of MoM and FDTD, not a child of MoM. A `yee-fem` crate matches the workspace shape used by `yee-mom` and `yee-fdtd` and avoids loading non-MoM users with FEM compile cost.
- `crates/yee-core/**` — possibly one new error variant (`Eigenproblem`); no API breaks expected.
- Out-of-lane (do not touch in the implementation PR): `yee-py` Python binding for the eigensolver is a follow-up yee-py-lane PR; `yee-cli` subcommand likewise.

## 15. References

- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd ed., Wiley 2014. Ch. 9 (Nedelec tetrahedral elements), Ch. 10 (eigenvalue problems and cavity resonators).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012. §6.3 (rectangular cavity resonator, analytic TE_{mnp} and TM_{mnp} frequencies, wall-loss Q).
- Boffi, D., Brezzi, F., Fortin, M., *Mixed Finite Element Methods and Applications*, Springer 2013. Ch. 7 (curl-conforming spaces, discrete compactness).
- Demkowicz, L., *Computing with hp-Adaptive Finite Elements*, vol. 2, CRC 2007. Ch. 4 (de Rham complex, hierarchical Nedelec).
- Bossavit, A., *Computational Electromagnetism*, Academic Press 1998. Ch. 5 (Whitney forms, geometric origin of edge elements).
- Petosa, A., *Dielectric Resonator Antenna Handbook*, Artech House 2007. Ch. 3 (cylindrical DRA modal frequencies — fem-eig-003 reference).
- Lee, J.-F., Sun, D.-K., Cendes, Z. J., "Full-wave analysis of dielectric waveguides using tangential vector finite elements", IEEE Trans. Microwave Theory Tech., 39(8), 1991, pp. 1262–1271 (the 2-D analog cited in Phase 1.3.1.1).
- `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md` — Phase 1.3.1.1 2-D Nedelec cross-section eigensolver spec; direct analog for the 3-D Whitney-1 generalization in §4 and §6.
