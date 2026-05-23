# Phase 1.3.1.1 step 5.3 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-3-direct-sparse-beta-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/{solve,assembly,mod}.rs`,
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs`, `ROADMAP.md`.
**Out of lane** (findings, not fixes): `reference.rs` (the oracle —
never edit to force a match), `crates/yee-fem/**` (read its
`build_shifted`/`sp_lu` pattern only), `crates/yee-py/**`. No
`Cargo.toml` dependency (faer is already present).

## Step ladder

### G1 — Sparse shift-and-invert on the β-direct pencil
1. Build `K = k₀²B − A` and `B₁` as faer `SparseColMat` (from the
   assembled COO/dense — reuse the assembly output; faer triplet
   build like the yee-fem `build_shifted`).
2. Physics-informed shift `σ₀ = (k₀² − k_c²)·⟨ε_r⟩` from the existing
   cutoff-pencil Stage-1 selection (the hybrid's β² estimate).
3. `sp_lu(K − σ₀B₁)`, inverse-iterate `z ← (K−σ₀B₁)⁻¹ B₁ z` to the
   eigenpair nearest σ₀; transverse-energy filter rejects a spurious
   capture (re-shift if needed). β² = Rayleigh quotient on the
   converged **β-direct** eigenvector (now exact for that mode).

**Verification:** `cargo test -p yee-mom --lib eigensolver` green; the
direct-solve β on the horizontal slab is reported.

### G2 — Mesh-refinement convergence study
Run the inhomogeneous slab at ≥3 mesh densities (8×8, 16×16, 24×24…)
with the sparse direct solve; record β → reference convergence. This is
the (a)-discretization vs (b)-eigenvector-mismatch discriminator.

**Verification:** convergence table emitted; β trend documented.

### G3 — Gate closure
1. **FR-4 (ε_r=4.4):** add/repoint the inhomogeneous gate to assert
   numerical β within ≤5% of `reference.rs` — a **failing** gate
   (the §4 published-benchmark closure at a representative contrast).
2. **ε_r=10.2:** tighten if ≤5% achieved; else document the
   discretization-limited residual + queue step-5.4.
3. No regression: uniform-fill ≤1%, ε_r=1 canary, wr90, coupling
   guards, Z_w.

**Verification:** full block below, all exit 0.

### G4 — ROADMAP + ADR-0054
ROADMAP step-5.3 line (§4 inhomogeneous gap closed at FR-4, or the
residual finding + step-5.4 queue). ADR-0054 already written —
update only if the disposition differs from the spec.

## Full verification (all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --lib eigensolver
cargo test -p yee-mom --test eigensolver_inhomogeneous
cargo test -p yee-mom --test eigensolver_wr90
cargo test -p yee-mom --test wave_port_numerical_te10 --test te10_waveport
git diff --stat -- '**/Cargo.toml'        # expect EMPTY
```

## Escape-hatch
If the sparse direct solve cannot reliably land the physical mode past
the spurious cluster within 25 min → fall back to a DENSE direct
β-direct solve with the physics-informed shift (still recovers the true
eigenvector at validation `n`, resolving (b); defer the sparse/finer-
mesh to a follow-on), document, queue step-5.4 for sparse. If even the
dense direct solve cannot beat the hybrid's β at FR-4, STOP and surface
— it would mean (a) discretization dominates even at FR-4, a different
finding. NEVER edit `reference.rs`; NEVER relax the uniform anchor or
the ε_r=1 canary.

## Out-of-scope (findings, not fixes)
* Higher-order / curl-conforming p-refinement (step-5.4 if (a) dominates).
* Lossy complex-ε_r sparse solve (Phase 1.3.1.2).
* yee-fem; yee-py.
