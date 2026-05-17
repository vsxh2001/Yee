# ADR-0028: Phase 3.nl.0 natural-language design surface scope

## Status

Accepted — 2026-05-18 (spec only; implementation deferred to
follow-up tracks — see ADR-0031).

## Context

`ROADMAP.md` Phase 3 carves out a natural-language design surface:
an engineer asks for "a 2.4 GHz patch on RO4003C with ≥ 100 MHz
bandwidth" and gets a runnable Yee project file back. Today that
needs hours of Balanis Ch. 14 hand-arithmetic for a calculation
that is fundamentally textbook. The spec is unusual for the
project in that it touches an LLM — introducing determinism
(`ROADMAP.md` line 150: *"all interactions are reproducible
script — the natural-language layer is convenience, not magic"*)
and CI-without-credentials questions that scope must resolve up
front. Track NNNNN (merge `1a6ea9b`) lands the spec.

## Decision

Phase 3.nl.0 ships a one-shot prompt → project-file pipeline
constrained on four axes:

1. **One geometry family: rectangular inset-fed microstrip patch**
   (Balanis Ch. 14). Synthesis equations are unambiguous and
   `mom-003` already anchors against a published benchmark. Other
   families are 3.nl.2.
2. **Five-stage pipeline with exactly one non-deterministic
   stage.** Stage 1 (LLM parser) is the only stage that touches
   a model; Stages 2–5 (geometry resolution, Balanis synthesis,
   no-op surrogate refinement, emission) are pure functions of
   the typed `DesignIntent`. Stage 4 BO lands in 3.nl.1.
3. **The TOML is the truth, not the prompt.** Every invocation
   writes a sibling `yee.intent.json` (structured intent + LLM
   provenance + substrate-library version). Re-running with the
   saved intent regenerates the TOML **byte-identically** —
   sorted keys, fixed `{:.6e}`, pinned `toml_edit`.
4. **Two parser modes: LLM via Anthropic Messages API tool-use,
   plus an offline deterministic regex grammar.** The offline
   mode handles all 10 canonical prompts so the gate runs in
   default CI without credentials. The LLM uses structured-
   output (`tools` + `input_schema`); never proposes `W` / `L`,
   only a schema-validated `DesignIntent`. Local CLI / Python
   only; web-facing deferred (3.nl.3 carries prompt-injection).

Validation gate **nl-001** — 10 canonical prompts, each producing
a TOML whose `|S_11|` minimum lands within ±5% of the requested
frequency under `yee run`. Tolerance is deliberately loose: the
gate checks the **surface**, not the Green's accuracy (inherits
CLAUDE.md §10's `mom-003` posture). Lane: new `crates/yee-design/`
plus touches on `yee-py`, `yee-cli`, `yee-validation`.

## Consequences

- **A one-sentence intent becomes a runnable project file.**
  `yee run` accepts the emitted TOML unchanged.
- **Determinism is enforced by construction.** Without sorted
  keys, pinned serializer, and byte-identical regen, the surface
  fails the `ROADMAP.md` invariant.
- **CI default has no LLM dependency.** Offline parser carries
  the gate; `anthropic` in `[llm]` extra; live-API tests are
  `pytest`-marked.
- **Loose-tolerance posture inherited, not weakened.** Tightening
  waits on Phase 1.1.1.2 (ADR-0025) and 3.nl.1's surrogate loop.
- **Defers four sub-projects:** surrogate refinement (3.nl.1),
  more families (3.nl.2), multi-turn agent (3.nl.3), held-out
  matrix (3.nl.4).

## References

- `docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`
- Track NNNNN merge commit `1a6ea9b`.
- ADR-0031 — Phase 3.nl.0 implementation plan (companion).
- C. A. Balanis, *Antenna Theory*, 4th ed., Wiley 2016, Ch. 14.
- Anthropic Messages API — `tools` + `input_schema`.
- `ROADMAP.md` Phase 3 line 150 — reproducibility invariant.
- CLAUDE.md §3, §4, §10.
