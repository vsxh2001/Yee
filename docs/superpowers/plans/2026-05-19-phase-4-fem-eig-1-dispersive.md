# Phase 4.fem.eig.1 — Dispersive `ε_r(ω)` on the Tet-Mesh FEM Eigensolver — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan track-by-track.

**Companion spec:** `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
**Base SHA:** `5602609` (CLAUDE.md §1 — Phase 4.fem.eig.0 walking skeleton shipped end-to-end against `fem-eig-001` at 0.09 % rel.err on TE_{101}).
**Target phase:** 4.fem.eig.1 only. 4.fem.eig.1.1 / 1.2 / 1.5 / 2 / 3 are explicitly deferred — see §"Out of scope".
**Tech-stack additions:** none new. `faer 0.x` already in workspace; verify `Complex64` sparse LU surface area in pre-flight (1). No new direct dep.

---

## Goal

Phase 4.fem.eig.1 extends the shipped Phase 4.fem.eig.0 free-space FEM
eigensolver to **single-pole dispersive lossy media**, reusing the
Phase 2.fdtd.3 ADE `Material` enum (Drude / Lorentz / Debye) verbatim. The
delivered pipeline: a hand-rolled lossy-SiO₂-filled rectangular metallic
cavity (a = 10 mm, b = 5 mm, d = 20 mm) is consumed by a new
`yee-fem::dispersive` module which assembles complex `K(ω)` and `M(ω)` at
each Newton trial frequency, solves the linearised generalised eigenproblem
`K(ω₀) e = θ M(ω₀) e` via a complex peer of the v0 `InverseIterEigen`, and
takes an analytic Hellmann–Feynman Newton step in ω until the fixed-point
condition `θ = (ω/c)²` holds. Validation gate `fem-eig-002` enforces complex
TE_{101} resonance within ±0.5 % on Re(f) and ±5 % on Im(f) (Q-factor)
against a hand-derived Drude-model reference. CPU-only, single-threaded,
scalar FP64 complex, no GPU, single-pole dispersion only, scalar isotropic
complex `ε_r(ω)` and real `μ_r` per tet, PEC closed cavity only — same
execution model as Phase 4.fem.eig.0, lifted one axis (real → complex) and
wrapped in an outer Newton loop.

## Pre-flight — sparse complex linalg availability

Spec §8 sketches the complex inner solver as the load-bearing risk. Before
Step D2 starts, confirm at the implementation base SHA `5602609`:

1. `faer::sparse::SparseColMat<Complex64>` is constructible and
   `faer::sparse::FaerLuSolver<Complex64>` (or whatever the equivalent is in
   the workspace pin) factors and back-substitutes. If the API regressed,
   fall back to dense `nalgebra::DMatrix<Complex64>` LU on a per-mode
   shift-invert block — fem-eig-002 is ~2 k DoFs and dense complex is
   tractable.
2. `num-complex::Complex64` is already pulled transitively; no workspace
   `Cargo.toml` dep change needed (verify with `cargo tree -p yee-fem`).
3. The v0 `InverseIterEigen<f64>` algorithm in
   `crates/yee-fem/src/solve.rs` is a near-textbook deflated inverse-power
   iteration — its lift to `Complex64` is a search-and-replace plus
   complex-norm and `e^T M e` (transposed, not Hermitian) inner-product
   substitutions. If the v0 code has been written in a tightly `f64`-coupled
   way that breaks the lift, factor out the common iteration loop into a
   trait method first; the trait was introduced precisely for this purpose
   (spec §8).
4. `Material::permittivity` returns `Complex64` already (Phase 2.fdtd.3
   `crates/yee-fdtd/src/material.rs`). The plan moves the type to
   `yee-core::material::Material` and adds `permittivity_derivative`; the
   move is mechanical (one `pub use` re-export keeps `yee-fdtd` callers
   green).

If any of (1)–(3) blocks, escape-hatch per the standard >15-min rule
(CLAUDE.md §5) and surface as a Phase 4.fem.eig.1.0.1 finding; do **not**
weaken the fem-eig-002 gate to compensate.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-core/src/material.rs` | Create | Move `Material` enum from `yee-fdtd`; add `permittivity_derivative(omega)`. |
| `crates/yee-core/src/lib.rs` | Modify | `pub mod material;` |
| `crates/yee-fdtd/src/material.rs` | Modify | `pub use yee_core::material::{Material, MaterialMap as _LegacyMap};` — preserve the historical local `MaterialMap` (grid-tag struct; D3 keeps it in-tree). |
| `crates/yee-fem/Cargo.toml` | Modify | Add `num-complex` as a direct dep (already transitively present; promote for clarity). |
| `crates/yee-fem/src/element.rs` | Modify | Complex-coefficient lift: `assemble_tet_element(verts, eps_omega: Complex64, mu_omega: Complex64) -> NedelecTetElement<Complex64>`. Real callers wrap with `Complex64::from(real)`. |
| `crates/yee-fem/src/assembly.rs` | Modify | Complex matrix path; preserve real `FemEigenAssembly::new_free_space` via the literal lift. |
| `crates/yee-fem/src/solve.rs` | Modify | Add `ComplexInverseIterEigen` + `SparseEigenComplex` trait peer. |
| `crates/yee-fem/src/dispersive.rs` | Create | `MaterialDatabase`, `DispersiveSolver::{solve_at_frequency, track_mode}`, `DispersiveEigenpair(s)`. Newton + bisection-fallback. |
| `crates/yee-fem/src/lib.rs` | Modify | `pub mod dispersive;`, re-export the new public types. |
| `crates/yee-fem/tests/element_complex.rs` | Create | D1 unit test — complex-coefficient local matrices reduce to v0 real matrices when `eps_omega` / `mu_omega` are pure real. |
| `crates/yee-fem/tests/complex_inverse_iter_smoke.rs` | Create | D2 unit test — `ComplexInverseIterEigen` recovers known eigenvalues of a small hand-built complex pencil. |
| `crates/yee-fem/tests/material_database.rs` | Create | D3 unit tests — `permittivity_at`, `permittivity_derivative_at` match closed-form Drude / Lorentz / Debye. |
| `crates/yee-fem/tests/solve_at_frequency.rs` | Create | D4 unit test — single-frequency linearised solve on a lossy-cavity fixture. |
| `crates/yee-fem/tests/newton_tracker.rs` | Create | D5 unit test — `track_mode` converges on a hand-fixture with a known analytic complex root. |
| `crates/yee-validation/src/lib.rs` | Modify | D6 — `run_fem_eig_002_lossy_sio2_cavity` driver. |
| `crates/yee-validation/tests/fem_eig_002_lossy_sio2_cavity.rs` | Create | D6 production-gate test. |
| `crates/yee-fem/validation/README.md` | Modify | `fem-eig-002 (lossy-SiO₂)` row. |
| `crates/yee-py/src/fem.rs` | Modify (D7, optional) | Python binding `yee.fem.solve_cavity_dispersive(...)`. |
| `crates/yee-py/tests/test_fem_dispersive.py` | Create (D7, optional) | Python pytest re-running fem-eig-002 from Python. |

No changes to `yee-mom`, `yee-cuda`, `yee-gui`, `yee-plotters`, `yee-io`,
`yee-cli`, `yee-surrogate`, `yee-mesh`. The `#![forbid(unsafe_code)]` floor
is preserved across every touched crate.

## Step ladder

### Step D1 — complex-coefficient lift of `assemble_tet_element`

- **Brief:** Re-type `crates/yee-fem/src/element.rs`'s
  `assemble_tet_element` and `NedelecTetElement` to carry `Complex64` matrix
  entries (spec §6). Barycentric gradients (real), Nedelec basis curls
  (real), 4-point Gauss quadrature weights (real) are unchanged; only the
  scalar pre-multiplier `(1/mu_omega) * V` on `K_local` and `eps_omega * V`
  on `M_local` become complex. The v0 real signature stays available as a
  thin wrapper `assemble_tet_element_real(verts, eps_r: f64, mu_r: f64)`
  that calls the complex path with `Complex64::from`. Pattern file:
  `crates/yee-fem/src/element.rs` itself (lines 83 onward — preserve the
  reference-tet ordering doc).
- **Lane:** `crates/yee-fem/src/element.rs`, `crates/yee-fem/Cargo.toml`,
  `crates/yee-fem/tests/element_complex.rs` (create).
- **Base SHA dep:** none — branches off `5602609` directly.
- **DoD:** unit tests pass — `complex_matches_real_for_pure_real_eps_mu`
  (`assemble_tet_element` with `Complex64::new(2.5, 0.0)` reproduces the v0
  real path to `1e-12`); `imaginary_eps_produces_imaginary_M_diagonal`
  (`assemble_tet_element` with `Complex64::new(0.0, 1.0)` for `eps_omega`
  yields purely imaginary `M_local` diagonal). `cargo doc -p yee-fem
  --no-deps` clean.
- **Verification:** `cargo test -p yee-fem --release element_complex &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on `nalgebra::SMatrix<Complex64, 6, 6>`
  not implementing some trait the v0 code uses (e.g. `RealField`) → store
  the matrix as `nalgebra::DMatrix<Complex64>` for the dispersive path and
  keep `SMatrix<f64, 6, 6>` for the real path. Performance hit at v0
  granularity is tolerable.
- **LOC:** ~180.

### Step D2 — `ComplexInverseIterEigen` peer + `SparseEigenComplex` trait

- **Brief:** Add a complex peer of the v0 `InverseIterEigen<f64>` in
  `crates/yee-fem/src/solve.rs`. New trait
  `pub trait SparseEigenComplex { fn solve(&self, k: &CsrMatrix<Complex64>,
  m: &CsrMatrix<Complex64>, num_eigs: usize, sigma: Complex64) ->
  Result<EigenpairList<Complex64>, Error>; }`. Concrete
  `ComplexInverseIterEigen` mirrors the v0 algorithm: build `(K − σM)` as a
  `faer::sparse::SparseColMat<Complex64>`, factor once, deflated complex
  inverse iteration on `(K − σM)^{-1} M`. Sort by `Re(k²)` ascending.
  Eigenvector normalisation uses the **transposed** (not Hermitian) inner
  product `e^T M e = 1` so the Hellmann–Feynman derivative in D5 lands in
  the natural form for complex symmetric pencils (spec §11). Pattern file:
  `crates/yee-fem/src/solve.rs` itself — copy the v0 algorithm and
  search-and-replace `f64 → Complex64`, `.transpose() → .transpose()`
  (already transposed, no Hermitian conjugate), `.norm() → .map(|z|
  z.norm()).sum().sqrt()` (complex 2-norm in the Frobenius sense).
- **Lane:** `crates/yee-fem/src/solve.rs`,
  `crates/yee-fem/tests/complex_inverse_iter_smoke.rs` (create).
- **Base SHA dep:** D1 (consumes the complex element-matrix path indirectly
  via assembly, but the test fixture is hand-built so D1 is not a strict
  prereq for D2 itself; the gate consumer in D6 needs both).
- **DoD:** unit tests pass — `complex_diag_pencil_recovers_eigenvalues`
  (hand-built 4×4 complex diagonal `K = diag(1+0.1j, 2+0.2j, 5+0.05j,
  10+1j)`, `M = I`, expect eigenvalues to match the diagonal to `1e-10`);
  `complex_inverse_iter_reduces_to_real_when_imag_zero` (same eigenvalues
  as v0 `InverseIterEigen` on a real pencil to `1e-12`);
  `eigenvectors_t_M_normalised` (`(e^T M e − 1).abs() < 1e-10`).
- **Verification:** `cargo test -p yee-fem --release
  complex_inverse_iter_smoke && cargo clippy -p yee-fem --all-targets --
  -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on `faer::sparse::FaerLuSolver<Complex64>`
  not in the base SHA's `faer` pin → emit `nalgebra::DMatrix<Complex64>` and
  use `nalgebra`'s dense `LU` for the inner solve. At fem-eig-002's ~2 k
  DoFs the dense complex LU runs in a few seconds; performance regression
  is acceptable for the gate. Document in `validation/README.md`.
- **LOC:** ~260.

### Step D3 — `MaterialDatabase` + `permittivity_derivative`

- **Brief:** Two sub-steps:
  - **D3a (yee-core relocation):** create
    `crates/yee-core/src/material.rs` and move the `Material` enum from
    `crates/yee-fdtd/src/material.rs`. Keep the `MaterialMap` grid-tag
    struct in `yee-fdtd` (it's FDTD-grid-shaped, not FEM-tet-shaped).
    Update `yee-fdtd/src/material.rs` to `pub use yee_core::material::Material;`
    so all existing `yee-fdtd` callers compile unchanged. Verify
    `cargo test -p yee-fdtd --release` stays green.
  - **D3b (derivative + database):** add
    ```rust
    impl Material {
        pub fn permittivity_derivative(&self, omega: f64) -> Complex64 {
            // closed-form per-pole dε/dω:
            // - Vacuum:        0
            // - Drude:   ω_p² · (2ω − jγ) / (ω² − jγω)²
            // - Lorentz: −Δε ω₀² · (−2ω + 2jδ) / (ω₀² − ω² + 2jδω)²
            // - Debye:   −Δε · jτ / (1 + jωτ)²
        }
    }
    ```
    plus `crates/yee-fem/src/dispersive.rs`'s `MaterialDatabase {
    materials: Vec<Material> }` with `new`, `permittivity_at(tet_id,
    omega) -> Complex64`, `permittivity_derivative_at(tet_id, omega) ->
    Complex64`. Index by `TetId = usize`. `From<&[f64]> for
    MaterialDatabase` lifts a v0 free-space `eps_r: Vec<f64>` to a
    `Vec<Material>` of constant `Drude { eps_inf: eps_r[i], omega_p: 0.0,
    gamma: 0.0 }` (lossless constant `ε = ε_inf`) so v0 callers keep
    compiling.
- **Lane:** `crates/yee-core/src/{lib,material}.rs`,
  `crates/yee-fdtd/src/material.rs`, `crates/yee-fem/src/dispersive.rs`
  (create), `crates/yee-fem/tests/material_database.rs` (create).
- **Base SHA dep:** none — D3a parallel-safe with D1, D2. D3b depends on
  D3a.
- **DoD:** unit tests pass — `permittivity_derivative_drude` (closed-form
  derivative against finite-difference `(eps(ω+h) − eps(ω−h))/(2h)` agrees
  to `1e-7` at `ω = 2π · 10 GHz`); same tests for `Lorentz` and `Debye`;
  `permittivity_at_zero_omega_drude_returns_finite` (sanity at the
  γ-determined limit); `vacuum_derivative_is_zero`; existing
  `cargo test -p yee-fdtd --release` stays green.
- **Verification:** `cargo test -p yee-core --release material &&
  cargo test -p yee-fdtd --release && cargo test -p yee-fem --release
  material_database` exits 0.
- **Escape hatch:** blocked > 15 min on a circular dep between yee-core and
  yee-fdtd because the `MaterialMap` grid-tag struct uses `ndarray` (a
  yee-fdtd-only dep) → keep `Material` in `yee-fdtd` and have `yee-fem`
  depend on `yee-fdtd` for the enum. Cleaner placement is in `yee-core`
  but the dep-graph result is what matters; document as a finding either
  way.
- **LOC:** ~280.

### Step D4 — `DispersiveSolver::solve_at_frequency` (single trial ω)

- **Brief:** Implement `DispersiveSolver::solve_at_frequency(&self, mesh,
  omega: Complex64) -> Result<DispersiveEigenpairs, Error>` in
  `crates/yee-fem/src/dispersive.rs`. For each tet, compute `eps_omega =
  material_db.permittivity_at(tet_id, Re(omega))` and pass to the complex
  `assemble_tet_element` (D1). Assemble complex `K(ω)`, `M(ω)` via the
  global scatter from v0 (sign-aware orientation flip, PEC row/column
  elimination — all path-unchanged from v0 because the orientation logic
  is real). Solve the *linearised* generalised eigenproblem at shift
  `σ = (omega/c)²` via `ComplexInverseIterEigen` (D2). Return the
  linearised eigenvalues `θ` and M-orthonormalised eigenvectors — this is
  **not yet a self-consistent dispersive eigenmode**; the Newton outer
  loop in D5 closes the loop.
- **Lane:** `crates/yee-fem/src/dispersive.rs`,
  `crates/yee-fem/src/assembly.rs` (complex path),
  `crates/yee-fem/tests/solve_at_frequency.rs` (create).
- **Base SHA dep:** D1 + D2 + D3 all merged.
- **DoD:** unit tests pass — `solve_at_frequency_vacuum_matches_v0`
  (with `MaterialDatabase` constructed via `From<&[f64]>` from `eps_r = 1.0
  everywhere`, the linearised eigenvalues at `omega = 2π · 9.66 GHz` match
  the v0 free-space `fem-eig-001` eigenvalues to `1e-6`);
  `solve_at_frequency_lossy_has_imaginary_eigenvalues` (with a Drude
  material the lowest eigenvalue has `Im(θ) < 0`). Wall-time < 20 s in
  `--release` on the fem-eig-002 mesh.
- **Verification:** `cargo test -p yee-fem --release solve_at_frequency &&
  cargo clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the complex assembly orientation
  signs producing a non-symmetric `K(ω)` → cross-check that the orientation
  flip is applied identically to the real and complex paths (it should be
  — the orientation is geometric, not material). Diff
  `K(ω).to_real_part()` against the v0 `K` on the same mesh.
- **LOC:** ~320.

### Step D5 — Newton-Raphson `ω`-tracker (single mode)

- **Brief:** Implement `DispersiveSolver::track_mode(&self, mesh,
  omega_warm_start: f64, e_warm_start: Option<DVector<Complex64>>) ->
  Result<DispersiveEigenpair, Error>`. Algorithm per spec §4.2 pseudocode:
  starting from `ω = Complex64::from(omega_warm_start)`, repeat
  - `(theta, e) ← solve_at_frequency(mesh, ω)` (D4) and pick the mode
    nearest to `e_warm_start` (or the lowest by `Re(θ)` on the first
    iteration);
  - `F = theta − (ω/c)²`;
  - if `|F| < tol_residual` or `|Δω| < tol_omega` or `iter > max_iter`,
    return;
  - `F'_prime = e^T (dK/dω − θ dM/dω) e / (e^T M e) − 2ω/c²`, where the
    matrix derivatives `dK/dω`, `dM/dω` are assembled per-tet via the
    complex element-matrix path with `eps_omega →
    permittivity_derivative_at(tet_id, ω)` and `mu_omega → 0` (μ_r is
    real and frequency-independent in v1 per spec §2);
  - if `|F'_prime| < tol_F'`, **fall back to bisection** on a small
    (Re, Im) box around the current ω — track the sign of `Re(F)` and
    `Im(F)` on the four corners.
  - else `ω ← ω − F / F'_prime`.
  Pattern file: any Newton-method reference; the spec carries the explicit
  pseudocode in §4.2. The `dK/dω`, `dM/dω` matrices are *not* materialised
  globally — only their action on `e` is computed, via a per-tet quadratic
  form `Σ_T e_T^T (dK_T/dω − θ dM_T/dω) e_T` where `e_T` is the local
  6-entry restriction of `e` to tet `T`.
- **Lane:** `crates/yee-fem/src/dispersive.rs`,
  `crates/yee-fem/tests/newton_tracker.rs` (create).
- **Base SHA dep:** D4 merged.
- **DoD:** unit tests pass — `track_mode_vacuum_converges_in_one_step`
  (with `MaterialDatabase` from `eps_r = 1.0`, the Newton residual is zero
  in one iteration because `dK/dω = 0`, `dM/dω = 0`); `track_mode_lossy_drude_converges_5_iter`
  (a single-tet fixture with a hand-derived analytic complex root converges
  in ≤ 5 iterations to `|F| < 1e-10`);
  `track_mode_bisection_fallback_triggers_when_F_prime_small` (forced
  fixture with `dε/dω ≈ 0` switches to bisection and still converges).
- **Verification:** `cargo test -p yee-fem --release newton_tracker && cargo
  clippy -p yee-fem --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the Hellmann–Feynman derivative
  sign / transposition (complex symmetric eigenproblems use `e^T`, not
  `e^H` — spec §11) → derive the 2×2 hand case first: `K = diag(1+jε, 2)`,
  `M = I`, with explicit `dK/dε` and confirm the derivative formula
  numerically before deploying on tet meshes. Per CLAUDE.md §5 escape-hatch
  rule, falling back to *numerical* `F'` via finite differences
  `(θ(ω+h) − θ(ω−h)) / (2h)` is acceptable for v1 ship if the analytic
  derivative remains stubbornly wrong; document as a Phase 4.fem.eig.1.0.2
  finding.
- **LOC:** ~360.

### Step D6 — fem-eig-002 validation gate (lossy-SiO₂ cavity)

- **Brief:** Implement spec §9 validation. Construct
  `TetMesh3D::cavity_uniform(0.010, 0.005, 0.020, 8, 4, 16)`. Build a
  `MaterialDatabase` of `Material::Drude { eps_inf: 3.78, omega_p:
  2π·0.4e9, gamma: 2π·2.0e9 }` on every interior tet (the cavity walls
  are PEC; the fill is uniform Drude). Compute the analytic complex
  reference per spec §9.1: solve `ω² ε_Drude(ω) / c² = (π/a)² + (π/d)²`
  numerically for ω via inner Newton on the closed-form ε(ω). Run
  `DispersiveSolver::track_mode(&mesh, omega_warm_start =
  v0_air_resonance, None)`. Assert: (1) `|Re(f_FEM) − Re(f_analytic)| /
  Re(f_analytic) ≤ 0.005`; (2) `|Im(f_FEM) − Im(f_analytic)| /
  |Im(f_analytic)| ≤ 0.05`; (3) `iterations ≤ 8`; (4) bisection fallback
  did **not** trigger. Register
  `run_fem_eig_002_lossy_sio2_cavity` in `crates/yee-validation/src/lib.rs`
  mirroring the existing per-validation drivers.
- **Lane:** `crates/yee-validation/src/lib.rs`,
  `crates/yee-validation/tests/fem_eig_002_lossy_sio2_cavity.rs` (create),
  `crates/yee-fem/validation/README.md` (modify).
- **Base SHA dep:** D5 merged (and transitively D1 + D2 + D3 + D4).
- **DoD:** test passes within the ±0.5 % Re(f) and ±5 % Im(f) bounds;
  `validation/README.md` has `fem-eig-002 (lossy-SiO₂ TE_{101})` row; CI
  workflow `ci.yml` test step picks up the new validation crate test
  automatically (`cargo test --workspace --release`). Wall-time < 90 s in
  `--release`.
- **Verification:** `cargo test -p yee-validation --release
  fem_eig_002_lossy_sio2_cavity` exits 0.
- **Escape hatch:** blocked > 15 min with Re(f) error between 0.5 % and 1 %
  → refine to (12,6,24) and retry; if still > 0.5 % the failure is not
  mesh resolution — investigate (a) Newton sign error in D5, (b) complex
  inverse-iter mode-ordering picking the wrong branch, (c) Material
  parameters not matching the analytic root. Do **not** weaken the ±0.5 %
  bound. Per CLAUDE.md §4 "no solver feature ships without a
  published-benchmark validation case"; if the gate cannot be met,
  fem-eig-002 does not ship and the failure is a spec-level decision.
- **LOC:** ~340.

### Step D7 (optional) — Python binding `yee.fem.solve_cavity_dispersive(...)`

- **Brief:** Extend `crates/yee-py/src/fem.rs` with
  `solve_cavity_dispersive(a, b, d, nx, ny, nz, materials, omega_warm_start,
  num_modes = 1) -> list[tuple[complex, np.ndarray]]`. The `materials` arg
  is a Python list of dicts with `kind: "Drude" | "Lorentz" | "Debye" |
  "Vacuum"` and the per-kind parameters; the binding constructs a
  `MaterialDatabase` Rust-side. Mirror the existing `yee.fem.solve_cavity`
  binding pattern from Phase 4 T8. Pytest case
  `crates/yee-py/tests/test_fem_dispersive.py` re-runs fem-eig-002 from
  Python with the same ±0.5 % / ±5 % tolerances.
- **Lane:** `crates/yee-py/src/fem.rs`, `crates/yee-py/tests/test_fem_dispersive.py`.
- **Base SHA dep:** D6 merged.
- **DoD:** `maturin develop -p yee-py --release` succeeds;
  `pytest crates/yee-py/tests/test_fem_dispersive.py` exits 0; Python
  returns a `complex` whose real part is `8.62e9 ± 0.5 %` and imaginary
  part is `−9.5e6 ± 5 %`.
- **Verification:** `cd crates/yee-py && maturin develop --release && pytest
  tests/test_fem_dispersive.py` exits 0.
- **Escape hatch:** blocked > 15 min on PyO3 0.28 ABI / abi3-py310 mismatch
  with `Complex64` returns → ship `(re: float, im: float)` tuples instead
  of `complex`, document as a Phase 4.fem.eig.1.0.3 finding, surface for
  resolution in a follow-up yee-py-lane PR.
- **LOC:** ~220.

## Track sequencing

Critical path: `D1 → D4 → D5 → D6` and `D2 → D4 → D5 → D6` and `D3 → D4 →
D5 → D6`.

```
                    D1 ──┐
                         │
D3 ──────────────────────┼── D4 ── D5 ── D6 ──┬── D7 (optional)
                         │
                    D2 ──┘
```

- **D1, D2, D3 run in parallel** at the start. All three branch off
  `5602609` and touch disjoint files inside `yee-fem` (D1: `element.rs`;
  D2: `solve.rs`; D3: `dispersive.rs` + `yee-core::material`).
- **D4 depends on D1 + D2 + D3** (consumes the complex element matrices,
  the complex inverse-iter, and the material database).
- **D5 depends on D4** (consumes `solve_at_frequency` as the Newton inner
  step).
- **D6 depends on D5** (production gate consumes the full pipeline).
- **D7 is optional and depends on D6**.

Within CLAUDE.md §5's "up to 5 parallel agents" envelope: peak parallelism
is at the start (D1 ‖ D2 ‖ D3 — three agents). Serial bottleneck is the
`D4 → D5 → D6` chain, ~3 agents-days end-to-end.

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| **fem-eig-002 Re(f)** — Newton-converged Re(f_TE101) vs analytic | D6 | `|Re(f) − Re(f_analytic)| / Re(f_analytic) ≤ 0.5 %` | `< 90 s` `--release` |
| **fem-eig-002 Im(f)** — Newton-converged Im(f_TE101) vs analytic | D6 | `|Im(f) − Im(f_analytic)| / |Im(f_analytic)| ≤ 5 %` (equivalently Q ±5 %) | covered by same test |
| **Newton convergence** — outer iteration budget | D6 (sub-assertion) | iterations ≤ 8 | covered by same test |
| **No bisection fallback on gate** — F'-degenerate path not triggered | D6 (sub-assertion) | fallback counter is zero | covered by same test |

Both Re(f) and Im(f) rows land in `crates/yee-fem/validation/README.md`.
Per CLAUDE.md §4 "no solver feature ships without a published-benchmark
validation case" — fem-eig-002 is the published benchmark (Pozar §3.1
lossy-dielectric extension + Bucur et al. SiO₂ permittivity reference).
Higher-application gates are scoped to later phases:

- **fem-eig-003 (lossy DRA)** — dispersive Petosa DRA validation.
  Phase 4.fem.eig.3. Out of scope here.
- **Production-scale complex eigensolve** — ≥ 100 k DoFs. Phase 4.fem.eig.2.

## Lane / file inventory

| Step | Files |
|------|-------|
| D1 | `crates/yee-fem/src/element.rs`, `crates/yee-fem/Cargo.toml`, `crates/yee-fem/tests/element_complex.rs` (create) |
| D2 | `crates/yee-fem/src/solve.rs`, `crates/yee-fem/tests/complex_inverse_iter_smoke.rs` (create) |
| D3 | `crates/yee-core/src/{lib,material}.rs` (create `material.rs`), `crates/yee-fdtd/src/material.rs`, `crates/yee-fem/src/dispersive.rs` (create), `crates/yee-fem/tests/material_database.rs` (create) |
| D4 | `crates/yee-fem/src/{dispersive,assembly,lib}.rs`, `crates/yee-fem/tests/solve_at_frequency.rs` (create) |
| D5 | `crates/yee-fem/src/dispersive.rs`, `crates/yee-fem/tests/newton_tracker.rs` (create) |
| D6 | `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/fem_eig_002_lossy_sio2_cavity.rs` (create), `crates/yee-fem/validation/README.md` |
| D7 (opt) | `crates/yee-py/src/fem.rs`, `crates/yee-py/tests/test_fem_dispersive.py` (create) |

Cross-lane consumers (`yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`,
`yee-cuda`) are not touched in 4.fem.eig.1.

## Risk register

Spec §11 risks mapped to steps:

1. **Newton convergence near branch points** (spec §11). Mitigated by
   bisection fallback when `|F'| < tol_F'`. **Materialises in Step D5** as
   the conditional branch in the inner loop. The fem-eig-002 gate is
   constructed to be safely inside the quadratic-convergence basin (Drude
   pole far below ω of interest), so the fallback is exercised only by the
   D5 unit test, not the D6 gate.
2. **Complex eigenvalue ordering ambiguity** (spec §11). For lightly-lossy
   cavities (Q > 100) ordering by `Re(k²)` is unambiguous. Mitigated by an
   eigenvector continuity check in D5: project the new eigenvector onto
   the previous; if `|<e_new, e_old>|² < 0.5`, flag and bisection-fallback.
3. **ε_∞ vs. ε(ω) confusion** (spec §11). Mitigated by D1's complex
   signature on `assemble_tet_element`: the type system rules out the
   scalar-ε_∞ shortcut on the dispersive path. v0 callers must explicitly
   `Complex64::from`.
4. **`faer` complex sparse LU surface area** (spec §11). **Materialises in
   Step D2** — pre-flight (1) confirms the API; escape hatch falls back to
   dense `nalgebra::ComplexEigen` at fem-eig-002's 2 k DoF scale.
5. **Hellmann–Feynman transposed vs. Hermitian** (spec §11). **Materialises
   in Step D5**. Two-line cross-check on a 2×2 hand fixture in the
   `newton_tracker.rs` test set is the canary.
6. **Warm-start sensitivity** (spec §11). **Materialises in Step D6**.
   fem-eig-002 warm-starts from the v0 free-space resonance (16.77 GHz)
   tracking to a converged Re(f) = 8.62 GHz — a 2× ratio is comfortably
   inside the Drude pole's basin. Other geometries get an explicit
   frequency-sweep warm-start chain via the `omega_warm_start` parameter
   on `track_mode`.

## Out of scope

Explicit non-goals for this plan, per spec §2 and §13:

- **No Beyn 2012 / contour-integral nonlinear eigensolve.** Phase
  4.fem.eig.1.5 if and when Newton-with-bisection-fallback proves
  insufficient on a published case.
- **No multi-pole dispersion expansions.** Single-pole Drude / Lorentz /
  Debye only. Phase 4.fem.eig.1.1.
- **No magnetic dispersion μ(ω).** Phase 4.fem.eig.1.2.
- **No driven 3-D FEM.** Eigenmode only.
- **No higher-order Nedelec.** First-order Whitney-1 only — same as v0.
- **No production-scale complex sparse eigensolve.** Phase 4.fem.eig.2.
- **No periodic / Floquet BCs.** Phase 4.fem.eig.4+.
- **No GPU.** CPU-only scalar complex FP64.
- **No DRA validation.** Phase 4.fem.eig.3.
- **No anisotropic / tensor ε(ω).** Scalar isotropic complex only.
- **No CLI / GUI exposure** beyond the optional Python binding (D7).

## Final verification

```bash
cargo build  -p yee-core -p yee-fdtd -p yee-fem -p yee-validation
cargo clippy -p yee-core -p yee-fdtd -p yee-fem -p yee-validation \
  --all-targets -- -D warnings
cargo test   -p yee-core --release
cargo test   -p yee-fdtd --release
cargo test   -p yee-fem  --release
cargo test   -p yee-validation --release fem_eig_002
cargo fmt    --check --all
cargo doc    --no-deps -p yee-core -p yee-fdtd -p yee-fem
```

All eight must exit 0. Every existing `crates/yee-mom/`,
`crates/yee-fdtd/`, `crates/yee-mesh/`, `crates/yee-fem/` test (including
the shipped `fem-eig-001` v0 gate) stays green — Phase 4.fem.eig.1 is a
strict extension, not a refactor. The mom-001 dipole gate and the
Phase 2.fdtd.3 dispersive FDTD gates are untouched.

## Estimated total

- LOC: ~1 960 core (D1 ~180, D2 ~260, D3 ~280, D4 ~320, D5 ~360, D6 ~340,
  D7 ~220).
- Wall-time per agent: 4–6 days end-to-end at one-engineer pace. Critical
  path `D1 → D4 → D5 → D6` is ~3 days; `D2 → D4 → D5 → D6` and `D3 → D4 →
  D5 → D6` run in parallel through D4. D7 adds ~1 day.
- Risk concentration: Step D5 (Newton tracker on complex symmetric pencils)
  is the load-bearing engineering risk per spec §11; the 2×2 hand-fixture
  unit test is the canary. Step D2 (complex inverse iteration over `faer`
  sparse Complex64 LU) is the load-bearing infrastructure risk; the
  pre-flight + dense-fallback escape hatch isolates it.
