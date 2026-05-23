# Phase 1.3.1.1 step 5.6 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-6-p2-robust-mode-selection-design.md`
**Base SHA:** `<post-scoping-commit>` (set at dispatch)
**Lane:** `crates/yee-mom/src/eigensolver/{solve,mod}.rs`,
`crates/yee-mom/src/ports.rs`, `crates/yee-mom/src/eigensolver/assembly.rs`
(the p=2 uniform-fill anchor test only), `crates/yee-mom/tests/eigensolver_inhomogeneous.rs`,
`crates/yee-mom/examples/`, `ROADMAP.md`, `docs/src/decisions/0057-*.md`.
**Out of lane** (findings, not fixes): `reference.rs` (oracle — NEVER
edit), the p=2 element matrices in `assembly.rs` (validated in 5.5 —
consume; a change is a finding), `crates/yee-fem/**`, `crates/yee-py/**`.
No `Cargo.toml` dependency.

## Step ladder

### K1 — ε_eff-screened / physics-seeded mode selection
In `solve_dense_mixed`: replace "smallest transverse-dominated k_c²" with
selection of the transverse-dominated candidate maximising the β-direct
Rayleigh quotient β²=R(x) (= highest ε_eff = physical dominant quasi-TEM).
Combine with the existing transverse-energy + propagating (real, β²>0)
filters. Optionally seed σ₀ from a physics ε_eff estimate. **Must reduce
to current behaviour where it is already correct** (p1/homogeneous/FR-4).

**Verification:** `cargo test -p yee-mom --lib eigensolver` + the FR-4 /
homogeneous / vertical-slab cases stay green (bit-identical or improved).

### K2 — p=2 uniform-fill anchor (review P1-1)
Lib test (assembly.rs tests or a helper): p=2 uniformly-filled WR-90
(ε_r=2.55) β ≤1% vs analytic `√(ε_r k₀²−(π/a)²)`. Closes the ε_r≠1
`assemble_mixed_p2` coverage gap.

**Verification:** the p=2 uniform anchor passes.

### K3 — ElementOrder wiring through ports.rs
Expose `ElementOrder::Second` via `NumericalCrossSection` (constructor
opt or `solve` arg); first-order default. p=2 reachable end-to-end.

**Verification:** an end-to-end p=2 solve through `NumericalCrossSection`
works (a smoke test).

### K4 — ε_r=10.2 closure gate (DoD-3)
With K1's fixed selection, assert p=2 ε_r=10.2 β ≤5% vs reference 582.95
→ failing gate (promote the step-5.5 documented-finding study). If the
selection picks the physical mode but p=2 still > 5% at the dense-cap
mesh, run a finer `--release` example + document; close if achievable,
else the narrower finding (selection-fix validated by ε_eff recovery).

**Verification:** the new gate passes ≤5% (or the documented finer-mesh
branch).

### K5 — ROADMAP + ADR-0057
ROADMAP step-5.6 line (§4 closed across contrasts, or the finding). ADR-0057
already written; update only if disposition differs.

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
- If the ε_eff-screened selection cannot preserve the p1/FR-4/homogeneous
  gates (regresses them) within ~25 min, scope the new selection to p=2
  ONLY (order-conditional), keeping p1 bit-identical. If even that is
  blocked, commit K2/K3 (the anchor + wiring) + the selection attempt
  `#[ignore]`'d, surface the specific mis-selection, queue step-5.7.
- If selection is fixed (ε_eff recovers from 4.8 toward 8.17) but ≤5% not
  reached at the dense-cap mesh, that is the DoD-3 finer-mesh branch — a
  documented finding (selection-fix validated), not a failure; do NOT
  weaken gates.
- NEVER edit `reference.rs` or the validated p=2 element matrices; NEVER
  weaken FR-4/uniform/ε_r=1.

## Out-of-scope (findings, not fixes)
* Sparse cutoff selection beyond what K1 needs (a perf follow-on).
* p≥3. Lossy complex-ε_r (Phase 1.3.1.2). yee-fem; yee-py kwargs.
