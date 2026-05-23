# Phase 1.3.1.1 step 5.4 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-4-interface-graded-mesh-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/mesh.rs`,
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs`,
`crates/yee-mom/examples/` (optional release convergence example),
`ROADMAP.md`, `docs/src/decisions/0055-*.md`.
**Out of lane** (findings, not fixes): `reference.rs` (oracle — never
edit), `solve.rs` (consume the step-5.3 sparse solve read-only; if a
solver change is needed, surface as a finding), `assembly.rs` (the
first-order element matrices are reused as-is), `crates/yee-fem/**`,
`crates/yee-py/**`. No `Cargo.toml` dependency.

## Step ladder

### H1 — Graded horizontal-slab mesh
Add a `TriMesh2D` builder in `mesh.rs` that clusters `y`-grid lines
geometrically toward the interface `y = d₁` (symmetric grading, finest
at the interface), `x` uniform, a grid line placed EXACTLY at `d₁` so
the dielectric/air material partition stays sharp. Unit-test the mesh
(node count, the interface line is present, tags correct).

**Verification:** `cargo test -p yee-mom --lib eigensolver::mesh` green.

### H2 — Convergence study (ε_r=10.2)
Drive the graded mesh through the step-5.3 sparse `solve_dense_mixed` at
≥3 grading/DoF points; report β → reference 582.95. Run as a `--release`
example (`examples/`) if the DoF exceeds the inline dense-selection cap
(~12×12); a small inline point for the gate.

**Verification:** convergence table emitted; best β recorded.

### H3 — Gate disposition
- **If ≤5% reached:** add `loaded_beta_hi_contrast_matches_reference`
  (ε_r=10.2 graded) as a failing gate; §4 closed at high contrast.
- **Else (plateau):** record the best achieved β + improvement,
  keep the ε_r=10.2 reconciliation a non-failing diagnostic, queue
  step-5.5 (p-refinement) with the h-plateau evidence.
- No regression: FR-4 gate, uniform anchor, ε_r=1 canary, wr90,
  coupling guards, Z_w stay green (graded mesh is additive).

**Verification:** full block below, all exit 0.

### H4 — ROADMAP + ADR-0055
ROADMAP step-5.4 line (§4 closed at high contrast, or the h-plateau +
step-5.5 queue). ADR-0055 already written; update only if disposition
differs.

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
If graded `h`-refinement does NOT reach ≤5% within 25 min of tuning,
that is the DoD-2 "plateau" branch — do NOT keep grinding grading
parameters; record the best β + the improvement curve, keep the
diagnostic non-failing, queue step-5.5 (p-refinement) with the evidence
that first-order `h` is insufficient. A partial improvement (e.g.
16%→8%) is a valid, reportable outcome. NEVER edit `reference.rs`, NEVER
weaken FR-4 / uniform / ε_r=1 gates, NEVER change `solve.rs`/`assembly.rs`
to force a match (surface as a finding).

## Out-of-scope (findings, not fixes)
* Second-order / curl-conforming p-refinement (step-5.5 if h plateaus).
* Sparse cutoff-pencil selection (only if the inline gate demands the
  high-DoF point; otherwise the release example suffices).
* Lossy complex-ε_r (Phase 1.3.1.2). yee-fem; yee-py.
