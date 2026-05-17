# Phase 4.fem.eig.0 — 3-D FEM Eigenmode Solver (walking skeleton) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan track-by-track.

**Companion spec:** `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`
**Base SHA:** `3364b0c` (CLAUDE.md §2 — yee-surrogate status reconciliation).
**Target phase:** 4.fem.eig.0 only. 4.fem.eig.1–4 are explicitly deferred — see §"Out of scope".
**Tech-stack additions:** new workspace member `crates/yee-fem/`; new direct deps `nalgebra-sparse` (already in lock from Phase 1.3.1.1), `lobpcg` (pure-Rust LOBPCG, new); `faer` sparse LU (already in workspace) for the shift-invert inner solve. No system / FFI additions.

---

## Goal

Phase 4.fem.eig.0 ships a single end-to-end pipeline: a hand-rolled tetrahedral mesh of a closed rectangular metallic air-filled cavity (a = 22.86 mm, b = 10.16 mm, d = 30 mm) is consumed by a new `crates/yee-fem` crate which assembles the first-order Nedelec curl-curl stiffness `K` and vector mass `M` sparse matrices, applies PEC tangential-`E`-zero Dirichlet via row/column elimination, and solves the generalized eigenproblem `K e = k² M e` for the ten smallest positive eigenvalues via shift-invert LOBPCG. Validation gate `fem-eig-001` enforces TE_{101} resonance at 9.660 GHz within ±0.3% and the lowest ten eigenvalues match the Pozar §6.3 analytic table within ±1% pairwise. CPU-only, single-threaded, scalar FP64, scalar isotropic real `ε_r` / `μ_r` per tet, no losses, no higher-order elements, no GPU — same execution model as the 2-D Phase 1.3.1.1 cross-section eigensolver, generalised one dimension up.

## Pre-flight — sparse-linalg library decision review

Spec §8 recommends `lobpcg` (pure-Rust) for v0 with a `SparseEigen` trait swap-point. Before Step T5 starts, confirm at the implementation base SHA:

1. `lobpcg` crate is published, MIT-licensed, no system deps, and its public API accepts a user-supplied matrix-vector product closure (we need `(K − σM)^{-1} M · x` for shift-invert). If the crate has regressed to a static matrix-only API, fall back to a hand-rolled inverse-power iteration on `faer` sparse LU as documented in spec §8 (escape hatch: still one-PR behind the `SparseEigen` trait).
2. `nalgebra-sparse::CsrMatrix<f64>` is interop-compatible with `lobpcg`'s expected BLAS / dense block format (LOBPCG inherently iterates on a small dense block of vectors). The trait abstraction lets T4 emit `nalgebra-sparse` while T5 translates if needed.
3. `faer::sparse::SparseColMat<f64>` + `FaerLuSolver` (or equivalent) is the inner shift-invert preconditioner. If `faer` has added a native sparse generalized eigen between spec-write and implementation, prefer it for ecosystem coherence (spec §8 `TBD`).
4. Fallback to `arpack-rs` is **not** pursued — same CI-lint-policy escape hatch as Phase 1.3.1.1 (`crates/yee-mom/src/eigensolver/solve.rs` ships dense-only). Worst-case fallback is the hand-rolled inverse-power iteration documented in the spec.

If any of (1)–(3) blocks, escape-hatch per the standard >15-min rule (CLAUDE.md §5) and surface as a Phase 4.fem.eig.0.1 finding; do **not** weaken the fem-eig-001 gate to compensate.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `Cargo.toml` (workspace root) | Modify | Add `yee-fem` workspace member; promote `lobpcg` to `[workspace.dependencies]`. |
| `crates/yee-fem/Cargo.toml` | Create | New crate; deps: `yee-core`, `yee-mesh`, `nalgebra`, `nalgebra-sparse`, `faer`, `lobpcg`, `thiserror`, `num-complex`. |
| `crates/yee-fem/src/lib.rs` | Create | Crate root, `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`, module declarations + re-exports. |
| `crates/yee-fem/src/element.rs` | Create | Local 6-edge Nedelec tet element matrices (`K_local`, `M_local`) on a reference / physical tet. |
| `crates/yee-fem/src/assembly.rs` | Create | Global edge enumeration, scatter into sparse `K`, `M`; PEC row/column elimination. |
| `crates/yee-fem/src/solve.rs` | Create | `SparseEigen` trait + `LobpcgEigen` impl with shift-invert via `faer` sparse LU. |
| `crates/yee-fem/src/lib.rs` (cont.) | Create | `FemEigenAssembly`, `FemEigenSolver`, `EigenpairList` public surface. |
| `crates/yee-mesh/src/lib.rs` | Modify | Add `TetMesh3D` struct + constructor + `signed_volume` / `centroid` / `edges` / `boundary_edges`. |
| `crates/yee-mesh/src/cavity.rs` | Create | `TetMesh3D::cavity_uniform(a, b, d, nx, ny, nz)` — axis-aligned box tetrahedralization. |
| `crates/yee-fem/tests/element_local_matrices.rs` | Create | T3 unit test against Jin Ch. 9 published reference. |
| `crates/yee-fem/tests/assembly_dimensions.rs` | Create | T4 unit test — global matrix shape matches interior-edge count. |
| `crates/yee-fem/tests/lobpcg_smoke.rs` | Create | T5 unit test — recover smallest eigenvalue of a known symmetric pencil. |
| `crates/yee-mesh/tests/tet_mesh_3d.rs` | Create | T2 unit tests — orientation, volume, boundary classification. |
| `crates/yee-mesh/tests/cavity_uniform.rs` | Create | T6 unit test — cavity mesh sanity (cell count, boundary edges, total volume). |
| `crates/yee-validation/src/lib.rs` | Modify | T7 — `run_fem_eig_001_rectangular_cavity` driver. |
| `crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs` | Create | T7 production-gate test. |
| `crates/yee-fem/validation/README.md` | Create | `fem-eig-001 (TE_{101})` and `fem-eig-001 (mode-10 ordering)` rows. |
| `crates/yee-py/src/fem.rs` | Create (T8, optional) | Python binding `yee.fem.solve_cavity(...)`. |
| `crates/yee-py/src/lib.rs` | Modify (T8, optional) | Register `yee.fem` submodule. |
| `docs/src/tutorials/04-fem-cavity-eigenmode.md` | Create (T9, optional) | mdBook tutorial. |
| `docs/src/SUMMARY.md` | Modify (T9, optional) | Add tutorial entry. |

No changes to `yee-fdtd`, `yee-mom`, `yee-cuda`, `yee-gui`, `yee-plotters`, `yee-io`, `yee-cli`, `yee-surrogate`. The `#![forbid(unsafe_code)]` floor is preserved across every touched crate.

## Step ladder

### Step T1 — `crates/yee-fem/` scaffold (workspace member, lib.rs, Cargo.toml)

- **Brief:** Create `crates/yee-fem/` workspace member with `Cargo.toml`, `src/lib.rs`. Crate header sets `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`. Declare `pub mod element;`, `pub mod assembly;`, `pub mod solve;` as empty stubs that compile. Public surface placeholders: `FemEigenAssembly<'m>`, `FemEigenSolver`, `EigenpairList` matching the spec §7 API sketch — methods return `yee_core::Error::Unimplemented("Phase 4.fem.eig.0 step T1 stub")` for now. **Cross-lane:** root `Cargo.toml` (members array only). Call this out explicitly in the agent's report — CLAUDE.md §6 cross-lane convention.
- **Lane:** `Cargo.toml` (root, members + `[workspace.dependencies]` for `lobpcg`), `crates/yee-fem/**`.
- **Base SHA dep:** none — branches off `3364b0c` directly.
- **DoD:** `cargo check -p yee-fem` exits 0; `cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0; `cargo doc -p yee-fem --no-deps` is `missing_docs`-clean; public types `FemEigenAssembly`, `FemEigenSolver`, `EigenpairList` are documented and exported from the crate root.
- **Verification:** `cargo check -p yee-fem && cargo clippy -p yee-fem --all-targets -- -D warnings && cargo doc -p yee-fem --no-deps` exits 0.
- **Escape hatch:** blocked > 15 min on workspace `Cargo.lock` resolver churn → `git checkout --theirs Cargo.lock && cargo check --workspace && git add Cargo.lock` per CLAUDE.md §5. Do not hand-merge the lock.
- **LOC:** ~140.

### Step T2 — `TetMesh3D` extension to `yee-mesh`

- **Brief:** Add `TetMesh3D` to `crates/yee-mesh/src/lib.rs`, parallelling the existing `TriMesh2D` (spec §5, §7). Fields: `vertices: Vec<Point3<f64>>`, `tetrahedra: Vec<[u32; 4]>`, `vertex_material: Option<Vec<MaterialTag>>`, `tet_material: Option<Vec<MaterialTag>>`. Constructor `TetMesh3D::new(...)` validates: ≥ 4 vertices, ≥ 1 tet, every tet's signed volume `(1/6)·(v₁−v₀)·((v₂−v₀)×(v₃−v₀)) > ε`, tag-length consistency. Auto-reorient inverted tets by swapping `v₂ ↔ v₃` (vs. rejecting outright — easier for callers; document the silent re-ordering). Methods: `signed_volume(tet)`, `centroid(tet)`, `edges()` returning canonical-orientation edge list (`lower_endpoint` first), `boundary_edges()` returning indices into `edges()` for edges lying on `∂Ω` (an edge is boundary iff it bounds exactly two face-coincident triangles on `< 2` tets, i.e. the face is shared by `< 2` tets). The eigensolver consumes this purely by reference; no FEM logic lives in `yee-mesh`.
- **Lane:** `crates/yee-mesh/src/lib.rs`, `crates/yee-mesh/tests/tet_mesh_3d.rs` (create). No `cavity.rs` yet — that is T6.
- **Base SHA dep:** none — branches off `3364b0c` directly; parallel with T1.
- **DoD:** unit tests pass — `single_tet_signed_volume_positive` (reference tet `[(0,0,0), (1,0,0), (0,1,0), (0,0,1)]` → `V = 1/6` ± ε); `inverted_tet_is_reoriented_silently`; `boundary_edges_of_single_tet_has_6` (every edge of a free tet is boundary); `two_tet_shared_face_has_no_interior_boundary_edges` (two tets sharing a triangular face → 9 edges total, 8 boundary, 1 interior); `edges_canonical_orientation_lower_index_first`. `cargo doc -p yee-mesh --no-deps` is `missing_docs`-clean.
- **Verification:** `cargo test -p yee-mesh --release tet_mesh_3d && cargo clippy -p yee-mesh --all-targets -- -D warnings && cargo doc -p yee-mesh --no-deps` exits 0.
- **Escape hatch:** blocked > 15 min on the boundary-edge classifier (the "edge belongs to N tets via M shared faces" count is fiddly) → first implement the face-incidence map `HashMap<[u32;3], Vec<usize>>` (sorted triplet → tets touching), then boundary-edge is any edge whose two endpoint vertices co-appear in a face whose tet-count is 1. Validate against the two single-cube hand fixtures before generalising.
- **LOC:** ~360.

### Step T3 — local 6-edge tet element matrices (Nedelec K_local + M_local)

- **Brief:** Create `crates/yee-fem/src/element.rs`. For a physical tet with vertices `v₀..v₃` and per-tet `ε_r`, `μ_r`: compute the 4 barycentric gradients `∇λ_i` (closed form, `∇λ_i = (face-normal opposite i) / (3 · V)`), build the six Nedelec edge bases `N_{ij} = λ_i ∇λ_j − λ_j ∇λ_i` (Jin §9.4 eq. 9.43), and emit the 6×6 local stiffness `K^e_{αβ} = (1/μ_r,e) · V_e · (∇×N_α)·(∇×N_β)` (curls are constant per tet, integral is exact) and 6×6 local mass `M^e_{αβ} = ε_r,e · ∫ N_α · N_β dV` via 4-point Gauss quadrature on the reference tet (exact for the integrand's polynomial degree per Jin §9.4 quadrature table). Sign convention follows the canonical edge orientation (lower-endpoint-first); when the tet's local edge runs against global orientation, the row/column gets negated at assembly time, not here.
- **Lane:** `crates/yee-fem/src/element.rs`, `crates/yee-fem/tests/element_local_matrices.rs` (create).
- **Base SHA dep:** T1 merged.
- **DoD:** unit tests pass — `reference_tet_K_local_matches_jin_table_9_x` (hand-tabulated 6×6 from Jin 3rd ed. Table 9.x for the unit reference tet at `μ_r = 1`; entries agree to 1e-10); `M_local_symmetric_positive_definite` (smallest eigenvalue > 0 by `nalgebra::SymmetricEigen` on the 6×6 dense block); `K_local_kernel_is_dimension_3` (the kernel of `K^e` is the discrete gradient subspace, dimension `n_verts_tet − 1 = 3`); `scaling_with_tet_volume` (uniform scale of vertices by `α` → `K_local` scales by `α`, `M_local` by `α³`).
- **Verification:** `cargo test -p yee-fem --release element && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the Jin Table 9.x reference values (the book's notation conflates global / reference / barycentric coordinates) → first verify on a tet whose curls have analytic form (the reference tet at `(0,0,0),(1,0,0),(0,1,0),(0,0,1)` has `∇λ_0 = (-1,-1,-1)`, etc.); back out the 6×6 by hand; cross-check with the 2-D analog in `crates/yee-mom/src/eigensolver/assembly.rs` for sign conventions.
- **LOC:** ~280.

### Step T4 — global sparse `K`, `M` assembly + PEC Dirichlet elimination

- **Brief:** Create `crates/yee-fem/src/assembly.rs`. Walk `TetMesh3D::tetrahedra`; for each tet build the local→global edge map (six entries) using `TetMesh3D::edges()`; track local-vs-global orientation per edge (sign flip in scatter). Scatter local 6×6 blocks into `nalgebra_sparse::coo::CooMatrix<f64>`, convert to `CsrMatrix<f64>` at the end. Apply PEC Dirichlet: collect `boundary_edges()`, drop those rows and columns from both `K` and `M` (row elimination per spec §6 — keeps spectra clean; penalty alternative is explicitly out of scope for v0). Build the interior-DoF index map so eigenvectors can be lifted back to the full edge basis. Expose `FemEigenAssembly::assemble(&self) -> Result<AssembledMatrices, Error>` where `AssembledMatrices { k: CsrMatrix<f64>, m: CsrMatrix<f64>, interior_edges: Vec<usize> }`.
- **Lane:** `crates/yee-fem/src/assembly.rs`, `crates/yee-fem/src/lib.rs` (re-exports), `crates/yee-fem/tests/assembly_dimensions.rs` (create).
- **Base SHA dep:** T1 + T2 + T3 all merged.
- **DoD:** unit tests pass — `assembled_K_dimensions_match_interior_edge_count` (single-tet fixture → 6 edges, 0 interior, K is 0×0; two-tet fixture sharing one face → 9 edges, 1 interior, K is 1×1); `assembled_K_is_symmetric` (`(K - K.transpose()).max_norm() < 1e-12`); `assembled_M_is_symmetric_positive_definite` (smallest eigenvalue > 0 on a 4-tet cube split); `assembled_K_kernel_dimension_equals_interior_vertex_count_minus_one` (the discrete gradient subspace has dimension `N_int_verts − 1`; verify by counting null-space dimension via dense eigen on a small fixture).
- **Verification:** `cargo test -p yee-fem --release assembly && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the local-to-global edge orientation sign → write the four-cube-split fixture out by hand on paper (vertex coordinates, expected edge list, expected signs); validate the orientation-aware scatter on it before any larger mesh. The 2-D analog in `crates/yee-mom/src/eigensolver/assembly.rs` handles the same problem one dimension down — match its sign convention exactly.
- **LOC:** ~420.

### Step T5 — `SparseEigen` trait + `LobpcgEigen` shift-invert impl

- **Brief:** Create `crates/yee-fem/src/solve.rs`. Define `pub trait SparseEigen { fn solve(&self, k: &CsrMatrix<f64>, m: &CsrMatrix<f64>, num_eigs: usize, sigma: f64) -> Result<EigenpairList, Error>; }`. Implement `pub struct LobpcgEigen { max_iter: usize, tol: f64 }` against this trait: build `(K − σM)` as a `faer` sparse matrix, factor once via `faer::sparse::FaerLuSolver` (or `nalgebra-sparse` LU equivalent if `faer` doesn't expose the API at base SHA — see Pre-flight (3)). LOBPCG iterates on the operator `(K − σM)^{-1} M`, returning the `num_eigs` largest eigenvalues `θ` of that operator; physical eigenvalues are `k² = σ + 1/θ`. Sort ascending, return `EigenpairList { k: Vec<f64>, e: DMatrix<f64> }` with eigenvectors lifted to interior-DoF indexing (caller lifts further to full edges). The trait is the load-bearing decision; `LobpcgEigen` is the swap point — spec §8.
- **Lane:** `crates/yee-fem/src/solve.rs`, `crates/yee-fem/tests/lobpcg_smoke.rs` (create).
- **Base SHA dep:** T1 merged (uses `EigenpairList` type from T1's stubs). Dependency-sibling of T3/T4 but no shared file.
- **DoD:** unit tests pass — `lobpcg_recovers_smallest_eigenvalue_on_known_dense_pencil` (hand-built 4×4 symmetric pencil with known eigenvalues `[0.5, 1.2, 3.4, 7.8]`; LOBPCG at shift `σ = 0.1` returns the lowest three within 1e-6); `lobpcg_on_scaled_identity_pencil` (`K = αI`, `M = βI` → all eigenvalues are `α/β`, shift-invert converges in one iteration); `lobpcg_converges_within_max_iter_for_3d_laplacian` (small 3-D scalar Laplacian as a sanity case — not the FEM eigenproblem, just the solver wiring). `cargo doc -p yee-fem --no-deps` clean.
- **Verification:** `cargo test -p yee-fem --release solve && cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the `lobpcg` crate API mismatch (pre-flight (1)/(2) failed) → implement `InverseIterEigen` as the `SparseEigen` impl instead: deflated inverse-power iteration on `(K − σM)^{-1} M` using `faer` sparse LU, one mode at a time with explicit deflation against converged eigenvectors. Document in `validation/README.md`. The trait keeps the solve abstract; gate downstream is unaffected.
- **LOC:** ~320.

### Step T6 — rectangular-cavity mesh constructor

- **Brief:** Add `crates/yee-mesh/src/cavity.rs`: `pub fn cavity_uniform(a: f64, b: f64, d: f64, nx: usize, ny: usize, nz: usize) -> Result<TetMesh3D, Error>`. Build a regular `(nx+1) × (ny+1) × (nz+1)` grid of vertices inside `[0, a] × [0, b] × [0, d]`, decompose each axis-aligned brick into six tetrahedra via the canonical 6-tet Kuhn decomposition (preserves orientation, no slivers, well-conditioned for first-order Nedelec — spec §5 / spec §11 "mesh-quality sensitivity"). Tag all `(nx × ny × nz × 6)` tets with `MaterialTag::AIR` (or eventually a caller-supplied tag map; v0 uses uniform tagging). Return a `TetMesh3D` that passes `TetMesh3D::new` validation (positive signed volumes, consistent orientation).
- **Lane:** `crates/yee-mesh/src/cavity.rs`, `crates/yee-mesh/src/lib.rs` (one `pub mod cavity;`), `crates/yee-mesh/tests/cavity_uniform.rs` (create).
- **Base SHA dep:** T2 merged.
- **DoD:** unit tests pass — `cavity_uniform_cell_count` (`cavity_uniform(a,b,d, 2,2,2)` → 8 bricks × 6 tets = 48 tets, 27 vertices); `cavity_uniform_total_volume` (sum of `signed_volume` over all tets equals `a·b·d` within 1e-12); `cavity_uniform_boundary_edges_form_closed_surface` (every boundary edge connects two boundary vertices; boundary vertex count = `2·((nx+1)(ny+1) + (nx+1)(nz+1) + (ny+1)(nz+1)) − 4·(nx+ny+nz+1) − ...`-style Euler-formula identity, or — pragmatically — assert boundary-vertex set equals the geometric box surface). `cargo doc -p yee-mesh --no-deps` clean.
- **Verification:** `cargo test -p yee-mesh --release cavity_uniform && cargo clippy -p yee-mesh --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the Kuhn 6-tet decomposition orientation table → use the published table (e.g. <https://en.wikipedia.org/wiki/Tetrahedron#Subdivision_into_orthoschemes>, Kuhn 1960) verbatim and validate orientation via `TetMesh3D::new`'s built-in check (T2 auto-reorients inverted tets silently, so even a sign-flipped table is recoverable).
- **LOC:** ~280.

### Step T7 — fem-eig-001 validation gate (rectangular metallic cavity)

- **Brief:** Implement spec §9 validation. Construct `TetMesh3D::cavity_uniform(0.02286, 0.01016, 0.030, 8, 6, 10)` (~480 bricks × 6 = 2880 tets, ~12k edges, ~2k DoFs after PEC elimination — well within v0 sparse-eigen budget). Build `FemEigenAssembly` with `ε_r = μ_r = 1.0` (air) on every tet. Assemble `K`, `M`. Solve `FemEigenSolver::new(10).with_shift(σ).solve(&k, &m)` where `σ = 0.5 · k₀_TE101² = 0.5 · (2π · 9.66e9 / c)²`. Convert eigenvalues to frequencies via `f_n = c · sqrt(k_n²) / (2π)`. Sort ascending. Assert: (1) `|f_1 − 9.660 GHz| / 9.660 GHz ≤ 0.003` (±0.3%); (2) the lowest ten `f_n` match the Pozar TE/TM analytic table for `(a, b, d)` within ±1% pairwise (mode-by-mode RMS); (3) no spurious mode below `f_1` (the lowest returned eigenvalue is `> 0.5 · k_0_TE101²`); (4) wall-time `< 60 s` `--release` (informational; not a hard gate). Register the case in `crates/yee-validation/src/lib.rs` as `run_fem_eig_001_rectangular_cavity`, mirroring the existing per-validation drivers.
- **Lane:** `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs` (create), `crates/yee-fem/validation/README.md` (create).
- **Base SHA dep:** T4 + T5 + T6 all merged.
- **DoD:** test passes within the ±0.3% TE_{101} bound and the ±1% mode-10 bound; validation README has `fem-eig-001 (TE_{101})` and `fem-eig-001 (mode-10 ordering)` rows; CI workflow `ci.yml` test step picks up the new validation crate test automatically (`cargo test --workspace --release`).
- **Verification:** `cargo test -p yee-validation --release fem_eig_001_rectangular_cavity` exits 0.
- **Escape hatch:** blocked > 15 min with TE_{101} error between 0.3% and 1% on the (8,6,10) mesh → refine to (12,9,15) and retry; if still > 0.3% the failure is not mesh resolution — investigate (a) shift-invert finding gradient-kernel modes (lower `σ` to verify, then raise above the cluster), (b) sign error in T4 orientation scatter, (c) Jin Table 9.x value used in T3 was wrong reference. Do **not** weaken the ±0.3% bound — that is the spec §9 gate. Per CLAUDE.md §4 "no solver feature ships without a published-benchmark validation case"; if the gate cannot be met, the feature does not ship and the failure is a spec-level decision.
- **LOC:** ~320.

### Step T8 (optional) — Python binding `yee.fem.solve_cavity(...)`

- **Brief:** Add `crates/yee-py/src/fem.rs` exposing `yee.fem.solve_cavity(a: float, b: float, d: float, nx: int, ny: int, nz: int, num_eigs: int = 10) -> tuple[list[float], np.ndarray]` (resonant frequencies in Hz + mode-coefficient block). Register `yee.fem` submodule in `crates/yee-py/src/lib.rs`. Mirror the existing `yee.eigensolver` PyO3 binding pattern from Phase 1.3.1.1. Pytest case `crates/yee-py/tests/test_fem.py` re-runs the fem-eig-001 frequency check from Python with the same ±0.3% tolerance.
- **Lane:** `crates/yee-py/src/fem.rs`, `crates/yee-py/src/lib.rs`, `crates/yee-py/tests/test_fem.py`.
- **Base SHA dep:** T7 merged.
- **DoD:** `maturin develop -p yee-py --release` succeeds; `pytest crates/yee-py/tests/test_fem.py` exits 0; `yee.fem.solve_cavity(0.02286, 0.01016, 0.030, 8, 6, 10)[0][0]` is `9.66e9 ± 0.3%`.
- **Verification:** `cd crates/yee-py && maturin develop --release && pytest tests/test_fem.py` exits 0.
- **Escape hatch:** blocked > 15 min on PyO3 0.28 ABI / abi3-py310 mismatch with the new `yee-fem` dep → check the existing `yee.eigensolver` binding for the canonical pattern (Phase 1.3.1.1 follow-up). If the LOBPCG dep transitively pulls a non-abi3-compatible crate, defer the binding to Phase 4.fem.eig.0.1 and surface as a finding.
- **LOC:** ~180.

### Step T9 (optional) — mdBook tutorial `04-fem-cavity-eigenmode.md`

- **Brief:** Add `docs/src/tutorials/04-fem-cavity-eigenmode.md` walking through the fem-eig-001 problem end-to-end: build the WR-90-based rectangular cavity, assemble, solve, compare against Pozar §6.3, plot the mode profile (use the existing `yee-plotters` static PNG path). Register the chapter in `docs/src/SUMMARY.md`. Pure docs lane; no Rust touched.
- **Lane:** `docs/src/tutorials/04-fem-cavity-eigenmode.md`, `docs/src/SUMMARY.md`.
- **Base SHA dep:** T7 merged (tutorial references the published gate).
- **DoD:** `mdbook build docs/` exits 0; new chapter renders; cross-link from `docs/src/theory/` (TBD: add a theory chapter on FEM eigenmode in Phase 4.fem.eig.1) is left as a `// TBD: link once theory chapter lands` comment.
- **Verification:** `mdbook build docs/` exits 0; SUMMARY.md lists the new chapter under the Tutorials section.
- **Escape hatch:** blocked > 15 min on the mode-profile plot (3-D vector field rendering is non-trivial) → ship the chapter with a 2-D mid-plane slice plot (sample `E(x, b/2, z)` on a 2-D grid, plot `|E|` heatmap via `yee-plotters`); leave 3-D quiver as a Phase 4.fem.eig.0.1 finding.
- **LOC:** ~220.

## Track sequencing

Critical path: `T1 → T3 → T4 → T7` and `T2 → T6 → T7` and `T1 → T5 → T7`.

```
                    T1 ──┬── T3 ──┬── T4 ──┐
                         │        │        │
                         └── T5 ──┼────────┤
                                  │        │
T2 ──── T6 ──────────────────────────────── T7 ──┬── T8 (optional)
                                                 │
                                                 └── T9 (optional)
```

- **T1 and T2 run in parallel** (no shared file, no base-SHA dep). Both branch off `3364b0c`.
- **T3 depends on T1** (uses `yee-fem` crate scaffold).
- **T4 depends on T2 + T3** (consumes `TetMesh3D` + local element matrices).
- **T5 depends on T1** (lives in `yee-fem`, uses the trait surface). Dependency-sibling of T3/T4 — runs in parallel with either.
- **T6 depends on T2** (extends `yee-mesh` with the cavity constructor).
- **T7 depends on T4 + T5 + T6** (production gate consumes the full pipeline).
- **T8 and T9 depend on T7** and are independent of each other; both optional, can run in parallel.

Within CLAUDE.md §5's "up to 5 parallel agents" envelope: peak parallelism is at the start (T1 ‖ T2), then T3 ‖ T5 once T1 lands, then T4 ‖ T5 ‖ T6 once T2 + T3 land. Serial bottleneck is the T4 → T7 edge; T7 cannot start until the production-scale solve is wired.

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| **fem-eig-001 TE_{101}** — lowest eigenvalue vs Pozar | T7 | `|f_FEM − 9.660 GHz| / 9.660 GHz ≤ 0.3%` | `< 60 s` `--release` |
| **fem-eig-001 mode-10 ordering** — ten lowest modes vs Pozar table | T7 | ±1% pairwise RMS; no spurious mode below TE_{101} | covered by same test |
| No-spurious-modes sanity | T7 (sub-assertion) | smallest returned `k² > 0.5 · k₀_TE101²` (shift skipped the gradient cluster) | covered by same test |

Both rows land in `crates/yee-fem/validation/README.md`. Per CLAUDE.md §4 "no solver feature ships without a published-benchmark validation case" — fem-eig-001 is the published benchmark (Pozar §6.3 eq. 6.42). Higher-application gates are scoped to later phases:

- **fem-eig-002** — lossy-cavity Q-factor validation. Phase 4.fem.eig.2. Out of scope here.
- **fem-eig-003** — cylindrical DRA validation. Phase 4.fem.eig.3. Out of scope here.

## Lane / file inventory

| Step | Files |
|------|-------|
| T1 | `Cargo.toml` (root, cross-lane), `crates/yee-fem/Cargo.toml` (create), `crates/yee-fem/src/lib.rs` (create) |
| T2 | `crates/yee-mesh/src/lib.rs`, `crates/yee-mesh/tests/tet_mesh_3d.rs` (create) |
| T3 | `crates/yee-fem/src/element.rs` (create), `crates/yee-fem/tests/element_local_matrices.rs` (create) |
| T4 | `crates/yee-fem/src/assembly.rs` (create), `crates/yee-fem/src/lib.rs`, `crates/yee-fem/tests/assembly_dimensions.rs` (create) |
| T5 | `crates/yee-fem/src/solve.rs` (create), `crates/yee-fem/tests/lobpcg_smoke.rs` (create) |
| T6 | `crates/yee-mesh/src/cavity.rs` (create), `crates/yee-mesh/src/lib.rs`, `crates/yee-mesh/tests/cavity_uniform.rs` (create) |
| T7 | `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fem_eig_001_rectangular_cavity.rs` (create), `crates/yee-fem/validation/README.md` (create) |
| T8 (opt) | `crates/yee-py/src/fem.rs` (create), `crates/yee-py/src/lib.rs`, `crates/yee-py/tests/test_fem.py` (create) |
| T9 (opt) | `docs/src/tutorials/04-fem-cavity-eigenmode.md` (create), `docs/src/SUMMARY.md` |

Cross-lane consumers (`yee-cli`, `yee-gui`) are not touched in 4.fem.eig.0. CLI / GUI exposure lands as follow-up Phase 4.fem.eig.0.1 once the gate is green.

## Risk register

Spec §11 risks mapped to steps:

1. **Spurious modes from the gradient kernel of `∇×`** (spec §11). Mitigated by Nedelec basis (kernel is exactly the discrete gradient subspace) + shift-invert targeting `σ > 0`. **Materialises in Step T5** as the shift selection (`σ = 0.5 · k₀_TE101²` for fem-eig-001); the no-spurious-modes assertion in T7 is the canary. If T7 catches a gradient mode below TE_{101}, raise `σ` past the kernel cluster (spec §11 mitigation: `σ ≥ 0.1 · k_lowest_expected²`).
2. **Pure-Rust sparse generalized eigensolve maturity** (spec §11). `lobpcg` is published but not battle-tested at FEM scale. **Materialises in Step T5** as the `SparseEigen` trait swap-point; escape hatch is hand-rolled inverse-power iteration on `faer` sparse LU. The trait is the load-bearing decision; the library is replaceable in one PR.
3. **Mesh-quality dependence** (spec §11). First-order Nedelec on sliver tets shows 5%+ error even on TE_{101}. **Materialises in Step T6** — the Kuhn 6-tet uniform decomposition is well-conditioned by construction (no slivers), so fem-eig-001 on the uniform mesh hits ±0.3%. Arbitrary Gmsh-imported meshes are deferred to Phase 4.fem.eig.1 along with a `MeshSizeFactor ≤ 0.05·λ` recommendation; v0 uses the hand-rolled cavity only.
4. **Memory cost at ~1M edges** (spec §11). The fem-eig-001 mesh is ~12k edges (~2k DoFs after PEC), well within v0 budget. Production scale (≥100k edges) is **deferred to Phase 4.fem.eig.1**; v0 does not promise it. No mitigation needed here; called out so the gate run-time stays inside `< 60 s`.
5. **First-order Nedelec accuracy floor** (spec §11). The (8,6,10) mesh hits Pozar's 4-significant-digit table tightly. If the gate proves marginal in practice the right next step is hierarchical `p`-refinement (Phase 4.fem.eig.1), not `h`-refinement to absurd mesh counts. **Materialises in Step T7** escape hatch.
6. **Gmsh FFI surface for tet meshing** (spec §11). Not touched in v0 — the hand-rolled `cavity_uniform` constructor is the only mesh source. Gmsh tet ingestion is deferred. **No risk on this step ladder.**

## Out of scope

Explicit non-goals for this plan, per spec §2 and §13:

- **No driven 3-D FEM** (right-hand-side from wave-port / current / plane wave). Separate Phase 4 sub-project sharing `K` assembly with a different solve path.
- **No open-region radiators** (FEM-BEM hybrid). Phase 4.fem.eig.4+.
- **No higher-order Nedelec elements.** First-order Whitney-1 only. Phase 4.fem.eig.1.
- **No lossy / complex `ε_r` Q-factor extraction.** Real eigenvalues only. Phase 4.fem.eig.2 (`fem-eig-002` lossy-cavity gate).
- **No dielectric-resonator-antenna validation.** Phase 4.fem.eig.3 (`fem-eig-003` Petosa DRA gate).
- **No periodic / Floquet boundary conditions.** PEC + PMC natural Neumann only. Phase 4.fem.eig.4+.
- **No GPU.** CPU-only, single-threaded, scalar FP64 throughout — same execution model as Phase 1.3.1.1.
- **No adaptive remeshing / mesh refinement.** Mesh is an input contract from `yee-mesh`.
- **No MLFMA / ACA / hierarchical compression.** Direct sparse formats only at v0.
- **No tensor / anisotropic / dispersive permittivity.** Scalar isotropic real `ε_r`, `μ_r` per tet only.
- **No Gmsh `.msh` import path in v0.** Hand-rolled `cavity_uniform` only; Gmsh tet ingestion is a Phase 4.fem.eig.0.2 follow-up.
- **No CLI / GUI exposure** beyond the optional Python binding (T8). Direct Rust API + Python only in 4.fem.eig.0.

## Final verification

```bash
cargo build  -p yee-fem -p yee-mesh -p yee-validation
cargo clippy -p yee-fem -p yee-mesh -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-fem --release
cargo test   -p yee-mesh --release
cargo test   -p yee-validation --release fem_eig_001
cargo fmt    --check --all
cargo doc    --no-deps -p yee-fem -p yee-mesh
```

All seven must exit 0. Every existing `crates/yee-mom/`, `crates/yee-fdtd/`, `crates/yee-mesh/` test stays green — `yee-fem` is a new peer crate, not a refactor of any shipped solver. The mom-001 dipole gate (`dipole_z_at_resonance`) is untouched.

## Estimated total

- LOC: ~2 400 core (T1 ~140, T2 ~360, T3 ~280, T4 ~420, T5 ~320, T6 ~280, T7 ~320, T8 ~180, T9 ~220).
- Wall-time per agent: 5–7 days end-to-end at one-engineer pace. Critical path T1 → T3 → T4 → T7 is ~4 days; T2 → T6 → T7 runs in parallel (~3 days). T8 + T9 add ~1 day.
- Risk concentration: Step T5 (sparse generalized eigensolve in pure Rust) is the load-bearing engineering risk per spec §8; the `SparseEigen` trait isolates the library choice so a swap is one PR. Step T4 (orientation-aware scatter) is the load-bearing correctness risk; the 2-D analog in `crates/yee-mom/src/eigensolver/assembly.rs` is the pattern file.
