# Phase 1.3.1.1 step 5.1 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-1-transcendental-reference-design.md`
**Base SHA:** `6f716df` (post step-5 merge + ADR reconciliation)
**Lane:** `crates/yee-mom/src/eigensolver/**` (a new `reference.rs` or a
test helper), `crates/yee-mom/tests/eigensolver_inhomogeneous.rs`,
`ROADMAP.md`, `docs/src/decisions/0052-*.md`.
**Out of lane** (surface as findings, do NOT fix): the step-5 mixed
solver internals (`assemble_mixed`/`solve_dense_mixed` — only consume
them; if the reconciliation exposes a solver bug, surface it as a
finding with the measurement, do NOT silently patch the solver here),
`crates/yee-fem/**`, `crates/yee-py/**`.

## Step ladder

### R1 — Reference dispersion + independent verification

1. Implement `slab_loaded_beta(a, b, d1, eps_r, freq_hz) -> f64`: the
   LSE transverse-resonance transcendental for a horizontally-stratified
   rectangular guide, bracketed root-find (bisection or Brent) on the
   verified dispersion (spec §3). **Confirm the LSE-vs-LSM family
   against Pozar §3.6 first** — match it to the numerical dominant
   mode's field orientation.
2. **Independent unit test (DoD-1):** assert `slab_loaded_beta`
   reproduces a textbook-tabulated slab-loaded-guide β (cite source +
   page) — this isolates a reference bug from a solver bug.

**Verification:** `cargo test -p yee-mom --lib eigensolver::reference`
(or wherever the helper lands) exit 0; the textbook value matches.

### R2 — Reconcile numerical vs reference

1. In `eigensolver_inhomogeneous.rs`, compute `slab_loaded_beta` for the
   horizontal-slab geometry (ε_r=10.2) and compare to the numerical β
   (201.52). Report both + the relative error at 2-3 mesh densities.
2. Branch:
   - **Agree ≤5%:** make the transcendental comparison the primary gate
     (keep the bracket as secondary). Tighten if agreement is better.
   - **Disagree** (reference independently verified in R1): keep the
     V2′ bracket, add the reference comparison as a **reported but
     non-failing** diagnostic with the measured gap, root-cause notes
     (mesh? mode family? β-extraction?), and queue step-5.2. Surface the
     solver-vs-reference gap as a finding to the orchestrator.

**Verification:** `cargo test -p yee-mom --test eigensolver_inhomogeneous`
exit 0 (gate either tightened-and-passing or bracket-retained with the
diagnostic).

### R3 — Docs + ROADMAP + ADR-0052

1. Test docstrings: state whether the §4 published-benchmark gap is
   closed (R2 agreed) or still open with the measured gap (R2 disagreed).
2. `ROADMAP.md`: step-5.1 line — gap closed, or the reconciliation
   finding + step-5.2 queued.
3. ADR-0052 `docs/src/decisions/0052-phase-1-3-1-1-step-5-1-transcendental-reference.md`:
   record the LSE/LSM family decision, the reference source, and the
   reconciliation outcome.

**Verification:** `mdbook build docs/` exit 0;
`grep -n "step 5.1\|step-5.1" ROADMAP.md`.

## Lint floor (every commit)
```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
```

## Full verification (before declaring done; all exit 0)
```
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p yee-mom --lib eigensolver
cargo test -p yee-mom --test eigensolver_inhomogeneous
cargo test -p yee-mom --test eigensolver_wr90
git diff --stat -- '**/Cargo.toml'        # expect EMPTY
```

## Escape-hatch

If the reference cannot be made to match a textbook-tabulated value
within R1 (the dispersion derivation itself is the blocker) > 25 min →
commit the reference attempt behind `#[ignore]` with the derivation
notes, keep the V2′ bracket gate intact, and surface the blocker —
do NOT ship an unverified reference (an unverified reference is worse
than the honest bracket gate). If the reference IS verified but
disagrees with the solver, that is the R2 "disagree" branch — a
documented finding, not a blocker; ship it as a non-failing diagnostic.

## Out-of-scope (surface as findings, do not fix)

* Patching the step-5 mixed solver if the reconciliation exposes a bug
  (surface it; a solver fix is its own reviewed change).
* The vertical-slab LSM dual if it needs a separate dispersion (note as
  step-5.2).
* Sparse mixed solve; CPW multi-conductor.
