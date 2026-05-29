# ADR-0085: Filter Phase F0.1 — register synthesis gates in the aggregator

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0084 (F0 synthesis core), ADR-0082/0083 (`yee validate --list`)

---

## Context

ADR-0084 shipped `yee-synth` + `yee-filter` with gates `synth-001` (g-values),
`synth-002` (coupling/Qe), `filt-001` (ideal response meets mask) as crate
tests. ADR-0084 deferred registering them in the `yee-validation` aggregator
(so they surface in `yee validate --list`) to this follow-on, Phase F0.1, to
keep F0's lane off the aggregator file.

## Decision

Register the three synthesis gates in the single-source `case_registry()`:

- `yee-validation` gains driver wrappers `run_synth_001()`, `run_synth_002()`,
  `run_filt_001()` returning `CaseResult` (each reproduces the F0 gate check by
  calling `yee_synth`/`yee_filter` and comparing to the published reference;
  `Passed` iff within tolerance). All three are `ExecutionPolicy::Run` (fast,
  pure-math). New dep edges: `yee-validation → yee-synth, yee-filter`.
- A new `Solver::Synth` variant categorizes them (ids `synth-001`, `synth-002`,
  `filt-001`).
- `yee-cli` gains `ValidateTarget::Synth` matching `synth-*` and `filt-*`, the
  `run_validate_list` `Solver::Synth => "Synth"` arm, and the synopsis/doc
  update. `yee validate synth` runs/lists only the synthesis cases.

## Consequences

**Ships:** 3 new `Run` cases in `run_all()`/`list_cases()`; `Solver::Synth`;
`ValidateTarget::Synth`; they appear in `yee validate --list[ --json]` and
`yee validate synth`. Gates: a `yee-validation` test asserting the three drivers
return `Passed`; a `yee-cli` test that `yee validate synth --list` shows
`synth-001` + `filt-001`. `run_all()` ordering: append the synth cases after the
fem-eig block (registry is the single source — no drift).

**Not in scope:** new synthesis math (F0 shipped it); EM; layout.

**No new external dependency** (only internal crate edges). Lane:
`crates/yee-validation/**`, `crates/yee-cli/**`, root `Cargo.toml` (yee-validation
dep entries).

---

## References
ADR-0084; `docs/superpowers/specs/2026-05-29-filter-f0-1-aggregator-reg-design.md`;
`docs/superpowers/plans/2026-05-29-filter-f0-1-aggregator-reg.md`.
