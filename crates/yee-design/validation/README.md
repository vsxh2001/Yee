# yee-design — Validation

Phase 3.nl.0 ships the natural-language design surface as a walking
skeleton (spec
`docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`).
The `nl-001` production gate enforces spec §9 — every canonical prompt
parses to a typed `DesignIntent`, satisfies the spec §7 schema, emits a
deterministic `yee.toml`, and round-trips byte-identically through
`<out>.intent.json`.

## Canonical references

- Spec §9 — production-gate composition (schema / round-trip / solver /
  offline).
- Plan
  `docs/superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md`
  step R6 — DoD, escape hatch, lane.
- CLAUDE.md §4 — validation-gate policy and the `mom-001` precedent for
  separating fast lint-floor checks from multi-minute solver gates.
- CLAUDE.md §10 — `MultilayerGreens` placeholder; the existing patch
  driver (`mom-003`) is itself `Skipped` pending Phase 1.1.1.

## Cases — Phase 3.nl.0

The `nl-001` gate is composed of four sub-gates per spec §9. Every
sub-gate is run for each of the ten canonical prompts shipped in
[`crates/yee-design/validation/prompts.toml`][prompts]; the row below
records the wiring location, the always-on / `#[ignore]`'d disposition,
and the relevant tolerance.

| ID                            | Tolerance / Assertion                          | Default CI | `#[ignore]` |
|-------------------------------|------------------------------------------------|------------|-------------|
| `nl-001 (offline)`            | `yee_design::parse_offline(prompt)` succeeds   | yes        | no          |
| `nl-001 (schema)`             | spec §7 schema (frequency / substrate / enums) | yes        | no          |
| `nl-001 (round-trip)`         | `emit → intent.json → emit` byte-identical     | yes        | no          |
| `nl-001 (solver, ±5% f)`      | `|f_min − f_target| / f_target ≤ 0.05`        | no         | yes         |

The first three sub-gates run by default in every `cargo test -p
yee-validation --release` invocation; their wall-time is sub-second per
prompt (pure-Rust regex parse, serde round-trip, deterministic emit).
The solver sub-gate is `#[ignore]`'d for the reasons in the row notes
below.

## Solver sub-gate disposition (CLAUDE.md §10 inheritance)

Per CLAUDE.md §10, `MultilayerGreens` is a Phase 1.1.0 placeholder
(one-image DCIM only). The closest existing patch-resonance driver in
this repo is `mom-003`, which is itself `Skipped` per
`crates/yee-validation/src/lib.rs::run_mom_003` until Phase 1.1.1 lands
the real Sommerfeld-integral / multi-image DCIM extraction.

R6 inherits that posture: the `nl-001 (solver, ±5% f)` sub-gate is
plumbed end-to-end (the emitted `yee.toml` is a valid input to
whichever patch driver lands first, and
`yee_validation::nl_001::run_nl_001_solver_gate` returns a structured
`CaseResult` per prompt) but the test entry point at
`crates/yee-validation/tests/nl_001_canonical_prompts.rs` is
`#[ignore]`'d. The slow-gate body currently returns
`CaseStatus::Skipped` with a notes string spelling out the upstream
dependency; when Phase 1.1.1 ships, the body switches to the real
solve without touching the test layout.

## Running

```bash
# Fast sub-gates (A + B + C) — always-on in default CI.
cargo test -p yee-validation --release --test nl_001_canonical_prompts

# Slow solver sub-gate (D) — opt-in via `--include-ignored`. Will
# return `Skipped` per prompt until Phase 1.1.1 lands the real
# MultilayerGreens driver; then re-runs the actual ±5 % gate.
cargo test -p yee-validation --release --test nl_001_canonical_prompts \
    -- --include-ignored
```

The structured `CaseResult` payloads — one per prompt for sub-gate D —
are intended to be folded into the same JSON / Markdown report shape as
the `mom-001` / `mom-002` / `fem-eig-001` cases; the
`yee_validation::nl_001` module is `pub` so a downstream
[`Report`][report]-rollup landing in 3.nl.4 (production validation
matrix) can include all ten rows without re-implementing the sub-gate
composition.

[prompts]: ./prompts.toml
[report]: ../../yee-validation/src/lib.rs
