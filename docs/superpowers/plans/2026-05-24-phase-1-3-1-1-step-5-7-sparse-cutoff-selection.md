# Phase 1.3.1.1 step 5.7 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-24-phase-1-3-1-1-step-5-7-sparse-cutoff-selection-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/{solve,mod}.rs`,
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs`,
`crates/yee-mom/examples/`, `ROADMAP.md`, `docs/src/decisions/0058-*.md`.
**Out of lane** (findings, not fixes): `reference.rs` (NEVER edit), the
p=2 element matrices in `assembly.rs` (consume), `ports.rs` selection
contract (the public solve must behave the same), `crates/yee-fem/**`
(read the `LobpcgEigen`/`build_shifted` PATTERN only — re-implement the
small block-LOBPCG in yee-mom rather than cross-crate-coupling),
`crates/yee-py/**`. No `Cargo.toml` dependency.

## Step ladder

### M1 — Sparse cutoff shift-invert
In `cutoff_candidates`: build the cutoff pencil `A x = k_c² B x` as faer
sparse; `sp_lu(A − σB)` at a small positive σ (a fraction of the analytic
air-cutoff (π/a)², above the gradient floor); inverse-iterate / block-
LOBPCG to the few eigenpairs nearest σ → the low-cutoff physical
candidates. Reuse the step-5.3 sparse LU + the step-4 LOBPCG pattern.

**Verification:** `cargo test -p yee-mom --lib eigensolver` green.

### M2 — Dense-path agreement (DoD-1)
At the existing validation meshes, the sparse candidates must yield the
SAME dominant mode as the dense path. Keep dense `cutoff_candidates` as a
small-`n` fallback/reference (DoF-threshold or `_rq` variant); add a unit
test pinning sparse≈dense dominant β at FR-4 / homogeneous / vertical-slab.

**Verification:** FR-4 1.39%, homogeneous 1.5e-14, uniform 305.117,
vertical-slab 243.51 bit-identical (or within a tight tol if the sparse
iteration introduces benign rounding).

### M3 — Finer-mesh convergence + ε_r=10.2 closure (DoD-2)
Now-affordable finer ε_r=10.2 meshes (12×12, 16×16, 24×24…) → reference
582.95. If ≤5%, promote `loaded_beta_hi_contrast_*` to a failing gate
(§4 high-contrast closed). Else document the improved p=2 trend + queue
step-5.8. Finer points as a `--release` example if inline is too slow.

**Verification:** the high-contrast gate passes ≤5%, or the documented
finer-mesh finding.

### M4 — Mesh-scaling demo (DoD-3) + ROADMAP + ADR-0058
A timed test or release example at n well past the old ~457 cap.
ROADMAP step-5.7 line. ADR-0058 already written.

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
- If the positive-shift sparse solve cannot cleanly exclude the gradient
  cluster within ~30 min, add gradient-deflation (project against the
  discrete-gradient range — the standard Nedelec remedy). If THAT is also
  blocked, commit the sparse-solve attempt behind a feature/threshold
  with the dense path as default, document the blocker, queue step-5.8 —
  do NOT ship a sparse selection that disagrees with the dense path on
  the validation gates.
- If sparse agrees with dense (M2) + finer meshes run (M3) but ε_r=10.2
  is still > 5%, that is the DoD-2 document-and-queue branch (the
  residual is then deeper than discretization) — not a failure; do NOT
  weaken gates.
- NEVER edit `reference.rs`; NEVER weaken FR-4/uniform/ε_r=1; NEVER change
  the public `NumericalCrossSection::solve` selected mode for the
  existing cases (DoD-1 is the guard).

## Out-of-scope (findings, not fixes)
* p≥3. Lossy complex-ε_r (Phase 1.3.1.2). yee-fem; yee-py.
* The breadth-rotation tracks (mom-002 numerical-port, FDTD Q6,
  real-port) remain documented-as-grind-risky.
