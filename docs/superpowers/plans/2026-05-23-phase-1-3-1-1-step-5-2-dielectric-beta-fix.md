# Phase 1.3.1.1 step 5.2 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-2-dielectric-beta-fix-design.md`
**Base SHA:** `981926b` (post step-5.1 merge)
**Lane:** `crates/yee-mom/src/eigensolver/{solve,assembly,mod}.rs`,
`crates/yee-mom/src/ports.rs` (only if β/Z_w plumbing needs it),
`crates/yee-mom/tests/eigensolver_inhomogeneous.rs`, `ROADMAP.md`,
`docs/src/decisions/0053-*.md`.
**Out of lane** (findings, not fixes): `reference.rs` (the verified
oracle — do NOT change it to make the solver "match"; if the fix matches
the uniform-fill analytic but not the reference, surface as a finding),
`crates/yee-fem/**`, `crates/yee-py/**`.

## Step ladder

### F1 — Confirm the bug with a uniform-fill test (cheapest isolation)

Add a uniformly-filled-guide test to `eigensolver_inhomogeneous.rs`:
WR-90 fully filled with ε_r=2.55 (a tag-uniform mesh), analytic
`β = √(ε_r k₀² − (π/a)²)` (≈305.16 at 10 GHz per the step-5.1
fully-filled anchor). Run the CURRENT solver against it and record the
(expected) failure — this isolates the β-extraction bug from
inhomogeneity + coupling. Commit the test `#[ignore]`'d-or-failing with
the measured wrong value as the bug witness.

**Verification:** test runs; record numerical vs analytic β.

### F2 — Fix the β-extraction (spec §3 option A)

Reformulate `solve_dense` + `solve_dense_mixed` to solve
`(k₀² T_ε − S) x = β² T_1 x` (eigenvalue = β², RHS = **unweighted**
mass `T_1 = ∫N·N`), or the minimal equivalent that makes β correct for
ε_r≠1. Requires `assemble_*` to expose / build the unweighted `T_1`
alongside the ε_r-weighted `T_ε`. Update mode selection to "largest
valid β²" and re-check the spurious-mode floor under the new operator.
Keep the ε_r=1 path bit-identical.

**Verification:** `cargo test -p yee-mom --lib eigensolver` green;
`cargo test -p yee-mom --test eigensolver_wr90` green (ε_r=1 canary
unchanged).

### F3 — Re-validate against analytic + reference; tighten the gate

1. F1 uniform-fill test now PASSES (≤1% vs analytic) — un-ignore it.
2. `eigensolver_inhomogeneous.rs`: the reconciliation diagnostic becomes
   a **failing gate** — numerical β matches `reference.rs` ≤5% for the
   horizontal (ε_r=10.2) and vertical (ε_r=2.2) slabs. Replace the stale
   regression values (180.23 / 201.52 — they were wrong) with the
   corrected ones.
3. Coupling guards (‖E_z‖/‖E_t‖ horizontal slab, zero-`B_tz` delta) +
   Z_w reduction must stay green.

**Verification:** full block below, all exit 0.

### F4 — Docs + ROADMAP + ADR-0053

ROADMAP: step-5.2 line — bug fixed, §4 published-benchmark gap CLOSED.
ADR-0053: record the root cause (β=k₀²−k_c² only valid for ε_r=1) + the
reformulation + the gate tightening. Note ADR-0051/0052 superseded on
the β-extraction point.

**Verification:** `mdbook build docs/` exit 0.

## Full verification (before done; all exit 0)
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

If the reformulation does NOT make β match the uniform-fill analytic
(DoD-1) within 25 min, the root cause is not (only) the β-extraction —
STOP, commit the uniform-fill bug-witness test `#[ignore]`'d, surface
the measurements + what you ruled out, and queue step-5.3. Do NOT
weaken the reference, do NOT touch `reference.rs`, do NOT relax the
homogeneous canary. If β matches the uniform analytic but still not the
inhomogeneous reference, that is a narrower finding (inhomogeneous-only
residual) — ship the uniform fix + analytic gate, document the
inhomogeneous gap, queue step-5.3.

## Out-of-scope (findings, not fixes)
* `reference.rs` (the oracle — never edit to force a match).
* Lossy/heterogeneous complex ε_r β-extraction (Phase 1.3.1.2).
* Sparse mixed solve; yee-fem; yee-py.
