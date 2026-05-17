# Phase 1.3.1.1 — Numerical 2-D Cross-Section Eigensolver — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan task-by-task.

**Goal:** Land a numerical 2-D FEM eigensolver inside `yee-mom` for arbitrary wave-port cross-sections, expose it via `WavePort::with_numerical_cross_section`, and validate it against the analytic TE10 mode (Phase 1.3.1.0) within 0.1% on β and 1% on Z_w.

**Companion spec:** `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`

**Tech stack additions:** `nalgebra-sparse` (sparse matrix assembly), optional `arpack-rs` (sparse shift-and-invert Arnoldi) behind a feature flag, dense fallback via `nalgebra::SymmetricEigen` (already pulled in via `nalgebra = 0.34`).

**Architecture:** Three new `pub(crate)` modules under `crates/yee-mom/src/eigensolver/`: `mesh.rs` (TriMesh2D adaptor + edge enumeration), `assembly.rs` (Nedelec edge basis + sparse A/B matrix fill), `solve.rs` (dense fallback + sparse shift-and-invert dispatch). One public type added to `crates/yee-mom/src/ports.rs`: `NumericalCrossSection`, with the matching `ModalDistribution::Numerical2D` variant and `WavePort::with_numerical_cross_section` builder.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-mom/Cargo.toml` | Modify | Add `nalgebra-sparse`; add optional `arpack-rs` behind `arpack` feature |
| `crates/yee-mom/src/lib.rs` | Modify | Declare `pub(crate) mod eigensolver;` |
| `crates/yee-mom/src/eigensolver/mod.rs` | Create | Module root + re-exports |
| `crates/yee-mom/src/eigensolver/mesh.rs` | Create | `TriMesh2D` adaptor, edge enumeration |
| `crates/yee-mom/src/eigensolver/assembly.rs` | Create | Nedelec edge basis, A & B sparse assembly |
| `crates/yee-mom/src/eigensolver/solve.rs` | Create | Dense fallback + sparse Arnoldi dispatch |
| `crates/yee-mom/src/ports.rs` | Modify | `NumericalCrossSection` + `ModalDistribution::Numerical2D` + `WavePort::with_numerical_cross_section` + `WavePort::rhs` match arm |
| `crates/yee-mom/tests/eigensolver_wr90.rs` | Create | Validation gate A (rectangular waveguide TE10 cross-check) |
| `crates/yee-mom/tests/eigensolver_loaded_wr90.rs` | Create | Validation gate B (slow-wave WR-90 + PEC septum) |
| `crates/yee-mom/validation/README.md` | Modify | New `eigensolver-001` / `eigensolver-002` rows |

---

## Step 0 — Dependencies and feature gate

**Files:** `crates/yee-mom/Cargo.toml`, root `Cargo.toml` (`[workspace.dependencies]`).

- [ ] Add to `[workspace.dependencies]`:
  - `nalgebra-sparse = "0.10"` (matches `nalgebra = 0.34` ABI; tracks the latest 0.10 minor at time of writing).
  - `arpack-rs = "0.3"` — **optional**, lives behind feature `arpack`.
- [ ] In `crates/yee-mom/Cargo.toml`:
  - `nalgebra-sparse = { workspace = true }`
  - `arpack-rs = { workspace = true, optional = true }`
  - `[features] arpack = ["dep:arpack-rs"]`. Default features stay empty.
- [ ] Run `cargo check -p yee-mom` (no-feature) and `cargo check -p yee-mom --features arpack` separately; the latter is allowed to fail in CI with a clear error if upstream `arpack-rs` doesn't build, in which case the spec's escape hatch (dense fallback only) takes over and the `arpack` feature is dropped from this milestone.

Estimated LOC: ~10. Verification: `cargo check -p yee-mom` exits 0.

---

## Step 1 — `NumericalCrossSection` stub + `ModalDistribution::Numerical2D`

**Files:** `crates/yee-mom/src/ports.rs`, `crates/yee-mom/src/eigensolver/mod.rs`.

- [ ] Add module declaration `pub(crate) mod eigensolver;` to `lib.rs`. `eigensolver/mod.rs` re-exports `TriMesh2D` and the eventual `solve_eigenproblem` entry point but is otherwise empty.
- [ ] In `ports.rs`:
  - Add `Numerical2D(NumericalCrossSection)` to `ModalDistribution`.
  - Add `NumericalCrossSection` struct with caches initialised to `NaN` / sentinel zero values; `solve(freq_hz) -> Result<()>` initially returns `yee_core::Error::Unimplemented("Phase 1.3.1.1 step 1 stub")`.
  - Add `WavePort::with_numerical_cross_section(mode)` builder.
  - Add the `ModalDistribution::Numerical2D` arm to `WavePort::rhs`; for the stub, it returns the same uniform-weighted RHS as `Uniform` so the build still passes existing tests.
- [ ] Unit test: `wave_port_numerical_stub_matches_uniform_before_solve` — verifies that without calling `solve`, the RHS is bit-for-bit identical to a `WavePort` with `Uniform` distribution.

Estimated LOC: ~80. Verification: `cargo test -p yee-mom --release ports::tests` exits 0.

---

## Step 2 — Nedelec edge-element basis and sparse A/B assembly

**Files:** `crates/yee-mom/src/eigensolver/mesh.rs`, `crates/yee-mom/src/eigensolver/assembly.rs`.

- [ ] `TriMesh2D`: a 2-D mesh struct with `Vec<[f64; 2]>` vertices, `Vec<[u32; 3]>` triangles, edge enumeration via `RwgBasis`-style shared-edge walk, and a boundary-edge predicate (Dirichlet `E_t = 0` on PEC walls).
- [ ] `LocalNedelecBasis`: per-triangle 3-edge basis (Whitney-1 forms). Provides element-stiffness and element-mass matrices for the curl-curl and ε-weighted terms; closed-form integrals on linear triangles (Pelosi/Coccioli/Selleri 2009 §3.4).
- [ ] `assemble_a_b(mesh: &TriMesh2D, eps_r: &HashMap, mu_r: &HashMap, freq_hz: f64) -> (CsrMatrix<Complex64>, CsrMatrix<Complex64>)`:
  - A holds curl-curl − k₀² ε_r mass terms.
  - B holds the ε_r mass term that multiplies β² on the RHS.
  - Both are `nalgebra-sparse::csr::CsrMatrix<Complex64>`.
- [ ] Unit tests:
  - `edge_enumeration_two_tri_mesh` — exactly 1 interior edge for a 2-triangle quad.
  - `element_stiffness_symmetric_on_unit_triangle` — local A is symmetric, eigenvalues real to within 1e-12.
  - `assembled_a_b_dimensions_match_dofs` — global matrix shape matches the interior edge count.

Estimated LOC: ~350. Verification: `cargo test -p yee-mom eigensolver::assembly::tests --release` exits 0.

---

## Step 3 — Dense fallback eigensolve on a tiny mesh

**Files:** `crates/yee-mom/src/eigensolver/solve.rs`.

- [ ] Implement `solve_dense(a: CsrMatrix<Complex64>, b: CsrMatrix<Complex64>) -> (f64, DVector<Complex64>)`:
  - Convert A and B to dense `nalgebra::DMatrix<f64>` (lossless case for the validation gates A/B).
  - Symmetric-pencil reduction: `L = chol(B)`; `M = L^{-1} A L^{-T}`; solve via `nalgebra::SymmetricEigen` on `M`.
  - Return the smallest positive eigenvalue (that's the dominant mode's `β²`) and the back-transformed eigenvector.
- [ ] Wire `NumericalCrossSection::solve` to call `assemble_a_b` and `solve_dense`; cache `β`, `z_w`, and `mode_profile`.
- [ ] Unit test on a hand-rolled 5-triangle WR-90 fixture (very coarse): the smallest-mode `β` is positive, finite, within an order of magnitude of `π / a`. Tolerance is loose — this is an assembly sanity check, not a validation gate.

Estimated LOC: ~180. Verification: `cargo test -p yee-mom eigensolver::solve::tests --release` exits 0 in `< 5 s`.

---

## Step 4 — Sparse shift-and-invert (feature-gated)

**Files:** `crates/yee-mom/src/eigensolver/solve.rs`.

- [ ] Behind `#[cfg(feature = "arpack")]`, implement `solve_sparse_shift_invert(a, b, sigma)`. Use `arpack-rs`'s generalized-eigenproblem entry point with shift `σ = (ω √ε_r,max / c)²`.
- [ ] Dispatch in `NumericalCrossSection::solve`: prefer sparse path when feature is on **and** the mesh has > 200 interior edges; else fall through to `solve_dense`.
- [ ] Bench (non-CI, informational): a 200-element mesh dense vs sparse on the same fixture. Expected: dense ~50–200 ms, sparse ~10–30 ms; both well under the 10 s budget. Capture numbers in a `// bench:` comment block.
- [ ] If `arpack-rs` fails to build at this step, drop the feature, document the escape hatch in `validation/README.md`, and proceed to Step 5 on the dense path. The DoD gates do not require sparse.

Estimated LOC: ~120 (skip if escape-hatched). Verification: `cargo check -p yee-mom --features arpack` exits 0 OR the feature is dropped with a documented rationale.

---

## Step 5 — Validation gate A (rectangular waveguide TE10 cross-check)

**Files:** `crates/yee-mom/tests/eigensolver_wr90.rs`, `crates/yee-mom/tests/fixtures/wr90_cross_section.rs`.

- [ ] Fixture: a structured `TriMesh2D` of the WR-90 cross-section (22.86 × 10.16 mm, air, eps_r = 1.0). Mesh density: 20 × 10 quads → 400 triangles → ~580 interior edges. Boundary edges tagged for PEC (Dirichlet `E_t = 0`).
- [ ] Integration test `eigensolver_wr90_te10`:
  - At `f = 10 GHz`, build `NumericalCrossSection`, call `solve(10e9)`.
  - Build a reference `RectangularWaveguideTe10 { a: 0.02286, b: 0.01016, eps_r: 1.0 }`.
  - Assert `|β_num − β_analytic| / β_analytic ≤ 0.001` (0.1%).
  - Assert `|Z_w_num.magnitude() − Z_w_analytic| / Z_w_analytic ≤ 0.01` (1%).
  - Assert L2 norm of `(E_t_num − E_t_analytic)` ≤ 1% of `||E_t_analytic||_L2`, where `E_t_analytic = sin(π x / a)` sampled at edge midpoints.
- [ ] Add row to `crates/yee-mom/validation/README.md` under a new "Wave-port mode solver" section.

Estimated LOC: ~150. Verification: `cargo test -p yee-mom --release eigensolver_wr90_te10` exits 0 in `< 30 s`.

---

## Step 6 — Validation gate B (slow-wave loaded WR-90)

**Files:** `crates/yee-mom/tests/eigensolver_loaded_wr90.rs`, `crates/yee-mom/tests/fixtures/wr90_with_septum.rs`.

- [ ] Fixture: WR-90 cross-section with a 1 mm-thick vertical PEC septum at `x = a/2`. The septum is realised by tagging the two adjacent triangle strips as Dirichlet-on-all-edges. Mesh density similar to Step 5.
- [ ] Integration test `eigensolver_wr90_septum_slow_wave`:
  - At `f = 10 GHz`, build `NumericalCrossSection`, `solve(10e9)`.
  - Compute `k_0 = 2π · 10e9 / c₀`.
  - Assert `β_num > k_0` (slow-wave inequality — the septum confines the mode to a smaller effective cross-section).
  - Regression-track `β_num` to a hard-coded value within ±1% (re-baseline once the assembly path settles).
- [ ] Add row to `crates/yee-mom/validation/README.md`.

Estimated LOC: ~120. Verification: `cargo test -p yee-mom --release eigensolver_wr90_septum_slow_wave` exits 0.

---

## Step 7 — Wire `WavePort::rhs` `Numerical2D` arm

**Files:** `crates/yee-mom/src/ports.rs`.

- [ ] Replace the Step-1 stub in `WavePort::rhs`'s `Numerical2D` arm: for each port-edge in `basis.port_basis_indices(self.tag)`, evaluate `mode.e_tangential_at(midpoint)` (a Nedelec interpolation), weight by edge length, scale by `self.voltage`.
- [ ] Doc note that the modal-profile sign convention is fixed by the dominant-eigenvector convention from `solve_dense` (positive-going wave, `β > 0`).
- [ ] Unit test on the WR-90 fixture: a `WavePort::with_numerical_cross_section` with the solved mode at 10 GHz produces an RHS whose port-edge weighting is within 1% (L2) of the equivalent `WavePort::with_rectangular_te10` RHS at the same frequency.

Estimated LOC: ~90. Verification: `cargo test -p yee-mom --release ports::tests::numerical_matches_te10_within_1pc` exits 0.

---

## Final verification

```bash
cargo build  -p yee-mom
cargo clippy -p yee-mom --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-mom
```

All five must exit 0. mom-001 dipole gate (`dipole_z_at_resonance`) must remain green — the eigensolver is opt-in and the `Uniform` / `Te10` paths are unchanged.

---

## Estimated total

- LOC: ~1100 (assembly is the bulk; tests are ~370).
- Wall-time per agent: 2–3 days at one-engineer-pace; less under TDD if Steps 2 and 3 are co-developed.
- Risk: Step 2 (Nedelec assembly) is the hardest; the escape hatch at Step 4 protects against the second-hardest (arpack-rs build).
