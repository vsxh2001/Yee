# Phase 1.3.1.1 step 5.5 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-5-second-order-nedelec-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/{assembly,mesh,mod}.rs`,
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs`,
`crates/yee-mom/examples/` (optional release convergence example),
`ROADMAP.md`, `docs/src/decisions/0056-*.md`.
**Out of lane** (findings, not fixes): `reference.rs` (oracle — never
edit), `solve.rs` (the step-5.3 sparse solve is reused order-agnostically;
a solver change is a finding), `crates/yee-fem/**`, `crates/yee-py/**`.
No `Cargo.toml` dependency.

## Step ladder

### J1 — Triangle Gauss quadrature + p=2 element matrices
Add a 2-D triangle Gauss rule (≥ degree-4, e.g. 6-point) and the p=2
element matrices: second-order Nedelec for `E_t` (2 DoF/edge + 2
interior; Jin §9.4 / Webb hierarchal) and quadratic nodal for `E_z` (6
nodes); curl is non-constant so all integrals go through the quadrature.
Keep the first-order matrices intact (default path).

**Verification:** `cargo test -p yee-mom --lib eigensolver::assembly`
green; a UNIT test comparing a p=2 element-matrix entry to an
independent quadrature value (mirror the step-5 `local_b_ze`
independent-quadrature pin) — proves the element matrices before any
eigensolve.

### J2 — Higher-order DoF map + assembly + PEC elimination
Extend the DoF map (`mesh.rs`) for the p=2 edge/interior/midpoint DoFs +
orientation/sign handling; assemble the global `(A, B, B₁)` blocks at
p=2; extend PEC Dirichlet elimination. An `ElementOrder` selector;
first-order is the default.

**Verification:** `cargo test -p yee-mom --lib eigensolver` green.

### J3 — Correctness anchor: p=2 on homogeneous WR-90 (DoD-4)
Before the high-contrast case, confirm p=2 reproduces the analytic TE10
β on the homogeneous WR-90 at least as accurately as p=1 — a clean,
no-singularity correctness check that the p=2 element matrices + DoF map
+ solve are right.

**Verification:** a p=2 WR-90 unit/integration test passes at ≤ the p=1
error.

### J4 — High-contrast convergence + §4 closure (DoD-1/2)
p=2 ε_r=10.2 convergence study (≥3 meshes) → reference 582.95, showing
faster-than-p1 convergence. If ≤5%: add `loaded_beta_hi_contrast_p2_matches_reference`
failing gate (§4 closed at high contrast). High-DoF points as a
`--release` example if the dense selection caps inline meshes.

**Verification:** the new gate passes ≤5%; FR-4/uniform/ε_r=1/wr90
unchanged (no regression).

### J5 — ROADMAP + ADR-0056
ROADMAP step-5.5 line (§4 closed at high contrast, or — if p=2 still
short — the finding + step-5.6 queue). ADR-0056 already written; update
only if disposition differs.

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
This is the biggest implementation in the chain — heed the budget.
- If the p=2 element matrices do NOT pass the independent-quadrature unit
  test (J1) or the homogeneous-TE10 anchor (J3) within ~40 min, the
  element formulation is wrong — STOP, commit the partial + the failing
  anchor `#[ignore]`'d, surface the specific discrepancy (which matrix /
  which DoF), do NOT proceed to the high-contrast case on a broken
  element.
- If J1/J3 pass but J4 ε_r=10.2 still exceeds 5% after p=2 + refinement,
  that is a documented finding (the interface may have a weak
  singularity p doesn't fully resolve) — keep the diagnostic non-failing,
  record the p=2 convergence curve vs the p=1 plateau (still valuable
  evidence), queue step-5.6. Do NOT weaken any gate.
- NEVER edit `reference.rs`; NEVER weaken FR-4/uniform/ε_r=1; NEVER change
  `solve.rs` to force a match (surface as a finding).

## Out-of-scope (findings, not fixes)
* p≥3 / adaptive hp. Lossy complex-ε_r at p=2 (Phase 1.3.1.2).
* Sparse cutoff-pencil selection. yee-fem; yee-py.
