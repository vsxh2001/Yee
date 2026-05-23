# Phase 1.3.1.1 step 4 — in-tree block LOBPCG sparse eigensolver

**Status:** Draft
**Owner:** TBD
**Phase:** 1.3.1.1 step 4 (sparse eigensolver upgrade; ROADMAP
"Pending (high priority)" line 153).
**Depends on:** Track NNNNNN (`fb6be04`) — `SparseEigen` /
`SparseEigenComplex` traits + `InverseIterEigen` /
`ComplexInverseIterEigen` shift-invert escape-hatch implementations in
`crates/yee-fem/src/solve.rs`.
**Blocks:** robust multi-mode cross-section eigensolves for clustered /
degenerate waveguide spectra (Phase 1.3.1.1 step 5 quasi-TEM
microstrip wave-ports consume the lowest-k block).

## 1. Goal

Add a **block LOBPCG** (Locally Optimal Block Preconditioned Conjugate
Gradient, Knyazev 2001) implementation of the existing `SparseEigen`
trait, computing the `num_eigs` smallest `k²` eigenpairs of
`K e = k² M e` **simultaneously** rather than one-at-a-time. The
shift-invert operator `(K − σM)⁻¹M`, already factored once via faer
sparse LU in `build_shifted`, serves as the LOBPCG preconditioner.

The new solver is **pure-Rust, no new external dependency** — it reuses
faer (already a workspace dep) for the dense Rayleigh-Ritz subproblem
and the existing sparse LU. The `arpack-rs` path is explicitly **not**
taken (see §5).

Existing `InverseIterEigen` is retained; LOBPCG is an additive
`SparseEigen` impl selected by the caller. No validation tolerance is
weakened.

## 2. Background

### 2.1 Why the escape-hatch needs an upgrade

`InverseIterEigen` (Track NNNNNN, `solve.rs:171`) is **deflated
inverse-power iteration, one mode at a time**: mode `i` iterates
`x ← (K−σM)⁻¹Mx`, M-orthogonalising against the `i−1` already-converged
vectors each step. This is correct but has two structural weaknesses
exactly where waveguide eigenproblems live:

1. **Clustered / degenerate spectra.** Rectangular and symmetric
   cross-sections have genuine modal degeneracies (e.g. TE/TM pairs,
   `TE_{mn}`/`TE_{nm}` on a square). Sequential deflation converges
   slowly and accumulates M-orthogonality error across a cluster: each
   mode is only as orthogonal as the running Gram-Schmidt against
   floating-point-imperfect predecessors.
2. **Convergence rate.** Single-vector inverse iteration converges
   linearly at rate `|θ₂/θ₁|`; for tightly clustered shift-invert
   eigenvalues that ratio → 1 and the per-mode iteration budget
   (`max_iter = 1000`) is consumed.

Block LOBPCG addresses both: it carries an `n × b` block (`b ≥
num_eigs`, typ. `b = num_eigs + guard`), enforces M-orthonormality of
the whole block via a single Rayleigh-Ritz step per iteration, and
resolves clusters because the block subspace spans the degenerate
eigenspace directly.

### 2.2 The shift-invert preconditioner

LOBPCG minimises the Rayleigh quotient `ρ(x) = (xᵀKx)/(xᵀMx)` over a
3b-dimensional search space spanned by the current block `X`, the
preconditioned residual `W = T·R`, and the previous block `P`. With
`T = (K − σM)⁻¹M` (the same operator inverse iteration applies), the
preconditioned residual aims the search at the smallest-`k²` modes near
`σ`. faer's sparse LU of `(K − σM)` — already built by
`build_shifted` + `sp_lu` — supplies `T` as one triangular
solve per column. No second factorisation; the cost delta vs inverse
iteration is the dense `3b × 3b` Rayleigh-Ritz eigensolve per outer
iteration, negligible for `b ≤ 20`.

## 3. Approach

### 3.1 New type `LobpcgEigen`

`crates/yee-fem/src/solve.rs`:

```rust
/// Block LOBPCG (Knyazev 2001) shift-invert eigensolver.
pub struct LobpcgEigen {
    pub max_iter: usize,
    pub tol: f64,
    /// Guard columns added to the block beyond `num_eigs` for
    /// cluster robustness. Block size b = num_eigs + guard.
    pub guard: usize,
}
impl SparseEigen for LobpcgEigen { /* ... */ }
```

`solve(k, m, num_eigs, sigma)`:

1. Shape validation — identical guards to `InverseIterEigen`
   (square, matching dims, `1 ≤ num_eigs ≤ n`).
2. `let shifted = build_shifted(k, m, sigma)?;` then
   `shifted.sp_lu()` → preconditioner `T`. **Reuse the existing
   free functions** — do not duplicate.
3. Block size `b = (num_eigs + guard).min(n)`. Seed `X₀` as an
   `n × b` matrix of deterministic seeds (extend `seed_vector` to a
   block, or QR a deterministic pseudo-random block for reproducible
   CI) then M-orthonormalise via the existing `m_orthogonalize` /
   `m_normalize` building blocks (block variants).
4. LOBPCG outer loop (Knyazev Algorithm 4.1):
   - `R = K·X − M·X·Λ` (block residual; Λ = current Ritz values).
   - `W = T·R` (preconditioned residual, one LU solve per column).
   - M-orthonormalise `W`, `P` against `X`.
   - Rayleigh-Ritz on `S = [X | W | P]`: form `Sᵀ K S` and `Sᵀ M S`
     (small `3b × 3b` dense), solve the dense generalized symmetric
     eigenproblem via faer, take the `b` smallest, update
     `X ← S·C_b`, `P ← (W,P part of S·C_b)`.
   - Convergence: per-column relative residual
     `‖K xᵢ − k²ᵢ M xᵢ‖₂ / (k²ᵢ ‖M xᵢ‖₂) < tol` for the leading
     `num_eigs` columns.
5. Return the leading `num_eigs` Ritz pairs in an `EigenpairList`,
   sorted ascending by `k²`, eigenvectors M-orthonormal — **same
   postcondition contract as `InverseIterEigen`** (`solve.rs:96-102`).

### 3.2 Numerical guards

- **Rank deficiency in `[X|W|P]`.** Near convergence `P` collapses
  into `span(X)`; the `3b × 3b` Gram matrix `SᵀMS` goes
  near-singular. Standard mitigation: drop `P` columns whose
  M-norm-after-orthogonalisation falls below `√ε`, shrinking the
  search block that iteration. Document inline (Knyazev §4 "soft
  locking").
- **Shift `σ` near the spectrum.** Same failure mode as inverse
  iteration; `(K − σM)` sparse LU fails → `Error::Numerical`, identical
  message style to `solve.rs:246`.
- **Determinism.** CI must be bit-reproducible: seed the initial block
  from a fixed-seed deterministic generator (no `rand` thread RNG),
  documented in the type doc.

### 3.3 Complex arm (optional, this phase)

`SparseEigenComplex` / `ComplexInverseIterEigen` (`solve.rs:552,589`)
have the analogous one-at-a-time structure. A `ComplexLobpcgEigen` is
**out of scope for step 4** — flag as a step-4.1 follow-on. Lossy
dispersive cavities (`fem-eig-002`) keep the complex inverse-iteration
path. Spec §6 records this boundary.

## 4. Validation

LOBPCG must reproduce every existing eigensolver gate to the **same or
tighter** error, and additionally resolve a degenerate cluster the
sequential solver handles poorly.

- DoD-V1. `cargo test -p yee-fem --test lobpcg_smoke` and
  `complex_lobpcg_smoke` green (complex unchanged; real smoke run
  against `LobpcgEigen` too).
- DoD-V2. `solve.rs` unit tests
  (`recovers_smallest_eigenvalue_on_known_dense_pencil`,
  `scaled_identity_pencil`, `eigenvectors_m_orthogonal`) **duplicated
  for `LobpcgEigen`** — known 4×4 pencil eigenvalues {0.5, 1.2, 3.4,
  7.8} recovered to `1e-8`; block M-orthogonality `eᵀMe ≈ I`.
- DoD-V3. `cargo test -p yee-mom --test eigensolver_wr90` — WR-90 TE10
  cutoff still within the existing band; record LOBPCG iteration count
  vs InverseIterEigen.
- DoD-V4. `cargo test -p yee-validation --test
  fem_eig_001_rectangular_cavity` — TE_{101} 0.09% / mode-10 RMS 0.37%
  band held when the gate's solver is swapped to `LobpcgEigen`
  (or a new parametrised variant; do not regress the default).
- DoD-V5. **New degenerate-cluster unit test:** a pencil with a known
  double eigenvalue (e.g. a square-cross-section discrete Laplacian
  with `TE_{12}`/`TE_{21}` degeneracy, or a synthetic 6×6 with a
  repeated root) — `LobpcgEigen` returns both members M-orthonormal to
  `1e-6`; assert the residual on each is below tol. This is the
  capability `InverseIterEigen` is weak at and the reason step 4
  exists.

## 5. Why not `arpack-rs`

ROADMAP names "step 4 sparse arpack-rs / LOBPCG". The cross-section
eigensolver design spec (2026-05-17) §fallback already anticipated
this: `arpack-rs` binds system ARPACK (Fortran/LAPACK), which (a)
violates the workspace "feature flags default OFF for anything
requiring an external toolchain" rule (CLAUDE.md §3) if made default,
and (b) the published `lobpcg` crate was found unavailable/unusable at
Track NNNNNN time (`solve.rs:147`). An **in-tree pure-Rust block
LOBPCG** sidesteps both: no system library, no feature gate, no CI
lint-clean risk, and it is ~200 lines on top of faer primitives the
crate already links. This matches the project's pure-Rust LA ethos
(faer, no `nalgebra-lapack`). If a future need for >10⁵-DoF
cross-sections arises, an `arpack` *optional feature* can be added then
behind the same `SparseEigen` trait — the trait is the swap point and
remains so.

## 6. Definition of done

DoD-1. `LobpcgEigen: SparseEigen` lands in
`crates/yee-fem/src/solve.rs`, fully documented (`#![warn(missing_docs)]`
clean), reusing `build_shifted` / `lu_solve` / `csr_matvec` /
`m_orthogonalize` / `m_normalize` (no duplication).
DoD-2. No new entry in `Cargo.toml` dependencies (pure-Rust, faer
only).
DoD-3. DoD-V1…V5 all green.
DoD-4. `crates/yee-fem/src/lib.rs` module doc + `solve.rs` header note
the LOBPCG availability and the `arpack` deferral; the complex arm
boundary (§3.3) is documented.
DoD-5. Tutorial / theory touch: `docs/src/theory/` eigensolver note (if
present) gains a LOBPCG paragraph; `ROADMAP.md` step-4 line marked
shipped with the merge SHA. ADR-0050 records the in-tree-LOBPCG vs
arpack-rs decision.
DoD-6. Lint floor clean: `cargo fmt --check --all` +
`cargo clippy --workspace --all-targets -- -D warnings`.

## 7. Risks

(a) **LOBPCG basis ill-conditioning** near convergence (§3.2). Mitigated
by soft-locking / `P`-column dropping; the new degenerate-cluster test
(DoD-V5) is the canary.
(b) **Rayleigh-Ritz dense solve** uses faer's symmetric generalized
eigensolver on `SᵀMS` — verify faer exposes a stable generalized
symmetric path; if not, reduce via Cholesky of `SᵀMS` then standard
symmetric eigensolve of the transformed `SᵀKS`. Either is dense and
small.
(c) **Determinism in CI** (§3.2). Fixed-seed block init; assert
reproducibility by running the smoke test twice in the same process if
cheap.

## 8. References

* Knyazev, "Toward the Optimal Preconditioned Eigensolver: LOBPCG",
  SIAM J. Sci. Comput. 23(2), 2001.
* `crates/yee-fem/src/solve.rs` — `SparseEigen`, `InverseIterEigen`,
  `build_shifted`, helper functions.
* `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`
  §fallback — the arpack-rs / dense-fallback contingency this spec
  resolves.
* ADR-0022 — Phase 1.3.1.1 eigensolver spec deferred-impl.
