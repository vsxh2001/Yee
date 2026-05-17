# ADR-0031: Phase 3.nl.0 NL design surface implementation plan

## Status

Accepted — 2026-05-18 (plan only; track execution deferred to
follow-up agents).

## Context

ADR-0028 locked the Phase 3.nl.0 scope. The next step is ~1 800
LOC of plan across a new `crates/yee-design/` crate, a Python
sidecar, a CLI subcommand, and a production gate. Track QQQQQ
(merge `6c55644`) lands the plan. Pre-flight finding: CLAUDE.md
§2's "yee-surrogate has not landed" is stale at base SHA
`f217df0`; the plan routes Stage 4 through the existing crate's
no-op refiner and surfaces the CLAUDE.md update as a docs-lane
finding.

## Decision

The plan splits into **seven tracks**:

- **R1 — `yee-design` scaffold + `DesignIntent` types.**
  Workspace cross-lane on root `Cargo.toml` members; substrate
  library loaded at build time; `serde_json` byte-identical
  round-trip as a property test.
- **R2 ‖ R3 ‖ R4 parallel post-R1.** R2: Balanis Ch. 14
  synthesis, gated by Example 14.1 (±0.5% on `W`, `L`, `y_0`).
  R3: deterministic `toml_edit` emitter (sorted keys, fixed
  `{:.6e}`) + byte-identical regen test (spec §8 determinism
  gate). R4: Python sidecar wrapping Anthropic Messages API
  with tool-use + `input_schema`; the schema JSON at
  `crates/yee-design/src/intent_schema.json` is a deliberate
  cross-lane shared artefact so Rust (`include_str!`) and
  Python (`jsonschema`) read the same source.
- **R5 — `yee design` CLI subcommand.** Depends on R1+R2+R3; R4
  optional since `--offline` is the default-CI path.
- **R6 — nl-001 production gate.** Four sub-gates per prompt:
  schema validity, byte-identical regen, solver `±5%`, offline
  success. Hardware-gated `#[ignore]` if 30 min overruns.
- **R7 (optional) — mdBook tutorial.** Parallel with R6.

Four load-bearing decisions locked here:

- **Offline parser is the CI default, not a fallback.** Default
  CI has no `ANTHROPIC_API_KEY`; live-API tests are `pytest`-
  marked (`-m anthropic`).
- **JSON schema is the Rust ↔ Python seam.** Single source of
  truth; no duplication.
- **Secrets policy is explicit.** Key read from env, never
  logged, never in `Provenance`, never on disk.
- **No new Rust LLM client.** Anthropic via Python sidecar;
  direct Rust client would be a `TECH_STACK.md` change, deferred.

Critical path: `R1 → R2 → R3 → R5 → R6`; R4 parallel with R2/R3;
R7 parallel with R6. Peak active set is three — within
CLAUDE.md §5's envelope.

## Consequences

- **CI stays green without credentials.** `anthropic` in
  optional `[llm]` extra; default `maturin develop` is
  unaffected.
- **Schema duplication is closed off.** One file edit, both
  sides re-read.
- **Loose-tolerance posture is inherited cleanly.** ±5% checks
  the surface; the underlying solver inherits CLAUDE.md §10's
  placeholder caveat until Phase 1.1.1.2 (ADR-0025).
- **CLAUDE.md §2 staleness about `yee-surrogate` is a docs-lane
  finding,** not fixed in-lane.

## References

- `docs/superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md`
- `docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`
- Track QQQQQ merge commit `6c55644`.
- ADR-0028 — Phase 3.nl.0 scope lock (this plan's parent).
- Anthropic Messages API — `tools` + `input_schema` structured
  output.
- CLAUDE.md §3, §4, §5, §6, §10.
