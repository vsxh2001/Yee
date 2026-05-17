# Phase 3.nl.0 — Natural-Language Design Surface (Walking Skeleton)

**Status:** Draft
**Owner:** TBD
**Phase:** 3.nl.0
**Depends on:** Phase 1.cli.1 (`yee validate`, `yee bench`, project-file plumbing), Phase 1.frontend.0 (yee-py PyO3 surface), Phase 3.gp.0 (surrogate framework, shipped — optional refinement hook only).
**Blocks:** Phase 3.nl.1 (surrogate-refinement loop), Phase 3.nl.2 (additional geometry families), Phase 3.nl.3 (interactive agent loop), Phase 3.nl.4 (production validation suite).

## 1. Motivation

The target user is a working RF / antenna engineer who already knows what they want: a 2.4 GHz patch on RO4003C with ≥ 100 MHz bandwidth and ≥ 6 dBi gain. Today, getting from that one-sentence intent to a Yee project file that the solver will accept requires (a) opening Balanis Ch. 14 to recover the inset-fed patch synthesis equations, (b) computing `W`, `L`, `y_0`, and the feed-line dimensions by hand, (c) writing the corresponding `yee.toml` by hand, and (d) discovering at solve time that one of the substrate parameters was mistyped. We have repeatedly observed this loop consume **two engineer-days for what is fundamentally a thirty-second textbook calculation**.

Phase 3.nl.0 closes that gap. The natural-language design surface accepts a free-form design intent, parses it into a structured `DesignIntent`, applies textbook synthesis equations to produce initial dimensions, and emits a `yee.toml` that the existing `yee run` pipeline accepts unchanged. The surface is a **convenience layer over an existing pipeline**, not a new solver.

The failure mode the surface prevents is hand-coding a starting geometry that a closed-form formula already gives. The failure modes the surface does **not** prevent — and explicitly leaves to other phases — are listed in §2.

## 2. Non-goals

The NL surface is explicitly **not**:

- A solver, optimizer, or simulator. It emits inputs; it does not run physics.
- A free-form CAD generator. It selects from a closed enum of supported geometry families and parameterizes them; it does not synthesize novel topologies.
- A substitute for engineering judgement. The emitted project file is a **starting point**; the user is expected to review it before running.
- A closed-loop optimizer. Surrogate-driven refinement to actually hit the bandwidth / gain targets is Phase 3.nl.1's responsibility; Phase 3.nl.0 emits the textbook initial estimate and stops.
- A free-form chat interface. Phase 3.nl.0 is one-shot: prompt in, project file out, no follow-up turn. The interactive agent loop is Phase 3.nl.3.
- A long-lived stateful service. Each invocation is stateless; reproducibility comes from logging the parsed `DesignIntent`, not from a session.
- Web-facing. Local CLI / local Python only; the prompt-injection threat model in §10 is the reason.

## 3. Scope decision

Phase 3.nl.0 is a walking skeleton in the sense of CLAUDE.md §3: the minimum end-to-end pipe through every layer, accuracy floor where the textbook puts it.

In scope:

- **One geometry family:** rectangular inset-fed microstrip patch antenna (Balanis Ch. 14). Chosen because the synthesis equations are unambiguous, the validation case `mom-003` already targets a 2.4 GHz patch on FR-4, and the Phase 1 closed-form references give us a published-benchmark anchor for the validation gate in §9.
- **One design-intent grammar** (§7): target frequency, substrate (named-preset or explicit `{eps_r, h_mm, loss_tangent}`), optional gain target in dBi, optional fractional bandwidth target.
- **One emitter:** `yee.toml` project file matching the format consumed by `yee run` today (TOML, per `TECH_STACK.md` §"Config / project files").
- **CLI subcommand:** `yee design "<prompt>" --output project.toml`.
- **Python entrypoint:** `yee.design.from_prompt(prompt: str) -> ProjectFile`.
- **Offline-mode fallback:** when no LLM credentials are available, parse a restricted English subset via a deterministic regex / template grammar so the validation gate (§9) can run in CI without network access.

Deferred to 3.nl.1+ (see §12): surrogate-refinement loop, additional geometry families, interactive Claude-as-tool agent, multi-prompt sessions, the production validation matrix.

## 4. Interaction model

Three modes, in increasing complexity. All three preserve the **reproducibility invariant** from `ROADMAP.md` line 150: "Underneath this surface, all interactions are reproducible script — the natural-language layer is convenience, not magic."

- **(a) One-shot NL → project file.** `yee design "..." -o p.toml` — runs the pipeline once, writes `p.toml`, exits. The prompt is preserved as a comment header in the emitted file (`# nl-prompt: ...`) so re-running `yee run p.toml` is fully deterministic and the original intent is recoverable from the artifact alone.
- **(b) Interactive refinement.** Phase 3.nl.3; out of scope here. Stub the API surface so a later agent loop can call the same pipeline as a tool.
- **(c) Batch / scripted.** The emitted `yee.toml` plus its sibling `yee.intent.json` (the structured `DesignIntent`, see §8) **are** the artifacts of record. Re-running the pipeline with the saved `DesignIntent` and an empty prompt regenerates the same project file byte-identically (§8). The LLM stage is short-circuited; the textbook-synthesis stage is pure.

The invariant: **the YAML/TOML is the truth, not the prompt**. The prompt is convenience input.

## 5. Pipeline architecture

The pipeline is five separable stages. Each stage is a function with a known input / output type, testable independently of the LLM:

```
NL prompt
   │
   ▼  ┌─────────────────────────────────────────────────────────┐
   │  │ Stage 1: Intent parser                                  │
   │  │   - LLM call with structured-output schema (§7)         │
   │  │   - or offline regex/template fallback                  │
   │  └─────────────────────────────────────────────────────────┘
   ▼
DesignIntent  (typed, validated, reproducible)
   │
   ▼  ┌─────────────────────────────────────────────────────────┐
   │  │ Stage 2: Geometry-family resolver                       │
   │  │   - DesignIntent.family enum → GeometryFamily impl      │
   │  │   - Phase 3.nl.0: only RectangularPatch supported       │
   │  └─────────────────────────────────────────────────────────┘
   ▼
GeometryFamily + DesignIntent
   │
   ▼  ┌─────────────────────────────────────────────────────────┐
   │  │ Stage 3: Initial-dimension calculator                   │
   │  │   - Balanis Ch. 14 closed-form synthesis                │
   │  │   - Pure function, deterministic, unit-tested standalone│
   │  └─────────────────────────────────────────────────────────┘
   ▼
InitialEstimate  (W, L, y_0, feed-line W, substrate stack-up)
   │
   ▼  ┌─────────────────────────────────────────────────────────┐
   │  │ Stage 4: (optional) Surrogate-refinement hook           │
   │  │   - Phase 3.nl.0: NO-OP. Hook exists; passes through.   │
   │  │   - Phase 3.nl.1 plugs in BO / NSGA-II loop.            │
   │  └─────────────────────────────────────────────────────────┘
   ▼
RefinedEstimate
   │
   ▼  ┌─────────────────────────────────────────────────────────┐
   │  │ Stage 5: Project-file emitter                           │
   │  │   - Renders yee.toml + yee.intent.json                  │
   │  │   - Deterministic; key ordering stable                  │
   │  └─────────────────────────────────────────────────────────┘
   ▼
yee.toml  (+ yee.intent.json)
```

Stages 2–5 contain no LLM calls and no I/O outside the emit step; they are unit-testable without network or credentials. Stage 1 is the only stage that touches a model.

## 6. API sketch

New crate: `yee-design` (Rust library, lives under `crates/yee-design/`). Avoiding overload of `yee-cli` keeps the LLM dependency feature-gated.

```rust
/// Structured representation of a parsed natural-language design intent.
///
/// Phase 3.nl.0: one family (rectangular patch). Round-trips to JSON for
/// reproducibility (see §8).
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
pub struct DesignIntent {
    pub family: GeometryFamily,
    pub target_frequency_hz: f64,
    pub substrate: Substrate,
    pub gain_target_dbi: Option<f64>,
    pub bandwidth_target_mhz: Option<f64>,
    /// Verbatim original prompt; preserved for the project-file header.
    pub source_prompt: String,
    /// LLM model id + temperature, or "offline" for the fallback parser.
    pub provenance: Provenance,
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq, Eq)]
pub enum GeometryFamily {
    RectangularPatch,
    // 3.nl.2: Wilkinson, HairpinFilter, ...
}

#[derive(serde::Serialize, serde::Deserialize, Debug, Clone)]
pub enum Substrate {
    Named(NamedSubstrate),  // FR4, RO4003C, RO5880, AluminaTC, ...
    Explicit { eps_r: f64, h_mm: f64, loss_tangent: f64 },
}

pub struct InitialEstimate { /* family-specific dimension struct */ }

pub struct ProjectFile {
    pub toml: String,
    pub intent_json: String,
}

/// One-shot pipeline. Returns the rendered project file.
pub fn from_prompt(prompt: &str, opts: &DesignOptions) -> Result<ProjectFile, Error>;

/// Pure path: skip the LLM, use a pre-parsed intent.
pub fn from_intent(intent: &DesignIntent) -> Result<ProjectFile, Error>;
```

Python bindings live in `yee-py`:

```python
import yee
proj = yee.design.from_prompt(
    "I need a 2.4 GHz inset-fed patch on RO4003C with at least 100 MHz "
    "bandwidth and gain over 6 dBi"
)
proj.write("patch_2g4.toml")
```

CLI: `yee design "<prompt>" --output project.toml [--model <id>] [--offline]`.

## 7. Structured-output schema

The LLM stage emits one JSON object conforming to the schema below. Anything else is a hard reject and triggers one re-prompt; after the second failure the pipeline aborts with a non-zero exit code. Schema-as-tool is the standard Anthropic Messages API pattern (tool use with `input_schema`).

```json
{
  "type": "object",
  "required": ["family", "target_frequency_hz", "substrate"],
  "properties": {
    "family": { "enum": ["rectangular_patch"] },
    "target_frequency_hz": { "type": "number", "minimum": 1.0e6, "maximum": 1.0e12 },
    "substrate": {
      "oneOf": [
        { "type": "object", "required": ["named"],
          "properties": { "named": { "enum": ["FR4", "RO4003C", "RO5880", "AluminaTC"] } } },
        { "type": "object", "required": ["eps_r", "h_mm", "loss_tangent"],
          "properties": {
            "eps_r": { "type": "number", "minimum": 1.0, "maximum": 100.0 },
            "h_mm": { "type": "number", "minimum": 0.05, "maximum": 10.0 },
            "loss_tangent": { "type": "number", "minimum": 0.0, "maximum": 0.5 }
          } }
      ]
    },
    "gain_target_dbi": { "type": "number", "minimum": -10.0, "maximum": 40.0 },
    "bandwidth_target_mhz": { "type": "number", "minimum": 0.0 }
  },
  "additionalProperties": false
}
```

The substrate library (Phase 3.nl.0: four entries) lives in `crates/yee-design/substrates.toml` and is read at build time so adding a substrate is a no-code change. The library is **versioned**: the version string ends up in `DesignIntent.provenance` so a re-run against a moved substrate is detectable.

## 8. Determinism + reproducibility

The pipeline is split deliberately so the **only** non-deterministic stage is Stage 1 (LLM intent parsing). All other stages are pure functions of `DesignIntent`.

Every `yee design` invocation writes two files:

- `<out>.toml` — the project file. First line is `# nl-prompt: <verbatim prompt>`. Second line is `# yee-design: <intent-hash> <provenance>`.
- `<out>.intent.json` — the serialized `DesignIntent` including the LLM model id, temperature, schema version, and substrate-library version.

Re-running with the saved `DesignIntent` (`yee design --from-intent <out>.intent.json -o <out>.toml`) regenerates `<out>.toml` **byte-identically**. This is testable as a property in CI: round-trip every canonical prompt in §9, hash the emitted TOML, re-run from `intent.json`, assert hash equality. This is the gate that distinguishes "convenience over a script" from "magic" — without it the NL surface fails the `ROADMAP.md` invariant.

Stage-5 emitter must therefore:

- Sort TOML keys lexicographically within each table.
- Use a fixed numeric format (`{:.6e}` for floats; one canonical form).
- Pin the TOML serializer (`toml_edit`) version in the workspace `Cargo.toml`.

## 9. Validation gate

Phase 3.nl.0 closes when the following 10 canonical prompts round-trip through the pipeline and the emitted project file, when run through `yee run`, produces an `|S11|` minimum within **±5%** of the requested frequency. Tolerance is loose deliberately: this validation is checking the **surface**, not the surrogate (Phase 3.nl.1) and not the multilayer Green's accuracy (Phase 1.1.1, still placeholder per CLAUDE.md §10).

Canonical prompts (`crates/yee-design/validation/prompts.toml`):

1. `"2.4 GHz patch on FR4"`
2. `"5.8 GHz patch on RO4003C with 200 MHz bandwidth"`
3. `"915 MHz inset-fed patch on FR4, gain over 5 dBi"`
4. `"3.5 GHz patch antenna on Rogers RO5880"`
5. `"a 1.575 GHz GPS patch, FR4 substrate, 1.6 mm thick"`
6. `"design a 2.45 GHz ISM-band patch with 100 MHz bandwidth on RO4003C"`
7. `"24 GHz patch on alumina, 0.5 mm substrate, εr 9.8"`
8. `"5 GHz WiFi patch on FR4, gain at least 6 dBi"`
9. `"a patch for 868 MHz LoRa on FR4"`
10. `"explicit substrate eps_r 3.0 h 0.508 mm tan_delta 0.0027, design a 10 GHz patch"`

Gate composition:

- **Schema gate.** All 10 prompts parse to a valid `DesignIntent` against §7 (offline-mode tested in CI; LLM-mode tested in the nightly job that has credentials).
- **Round-trip gate.** For each prompt, `from_prompt → from_intent` emits the byte-identical TOML.
- **Solver gate.** For each prompt, `yee run <emitted.toml>` (loose-tolerance microstrip / patch path; see CLAUDE.md §10 on placeholder Green's) finds `|S11|` minimum at `f_min` with `|f_min − f_target| / f_target ≤ 0.05`.
- **Offline gate.** All 10 prompts must succeed with `--offline` (the deterministic fallback parser) so the gate runs in default CI without API credentials.

Gate non-goals: the gate does **not** assert the gain or bandwidth target is met (that requires the surrogate-refinement loop in 3.nl.1). The gate asserts the surface produces a runnable, on-frequency starting point.

## 10. Risks / open questions

- **LLM hallucinated dimensions.** Mitigated: the LLM is constrained to emit `DesignIntent` only (schema-validated, §7). Dimensions are computed downstream by deterministic Stage 3 code. The LLM never proposes `W` or `L`.
- **Substrate-library drift.** A substrate-library version moves under a saved `intent.json`. Mitigation: substrate-library version is in `provenance`, and Stage 2 errors loudly on a version mismatch.
- **Prompt injection.** Not a concern in Phase 3.nl.0 (local CLI only). Becomes a concern in Phase 3.nl.3 (agent loop) and any web-facing variant; non-goal here, surfaced now so 3.nl.3's spec inherits it.
- **Eval-set drift.** The 10 canonical prompts are the entire eval set. They will overfit. Mitigation: 3.nl.4 extends to a held-out matrix; 3.nl.0 explicitly accepts overfit-to-eval as a known limit.
- **Hidden coupling to a specific model version.** The model id is logged in `provenance`. A future re-validation against a new model is a separate gate.
- **Schema-output failures cascading.** The Messages-API tool-use path can refuse to call the tool. Mitigation: one re-prompt with the schema embedded in the user turn; second failure → error out, suggest `--offline`.
- **TBD: licensing of named-substrate datasheet parameters.** Rogers / Taconic / Isola publish `eps_r` curves; one-frequency literals are fine, full curves may not be. Phase 3.nl.0 uses literals only; TBD: confirm with `THIRD_PARTY_LICENSES.md` author before shipping any substrate beyond FR4.

## 11. Dependencies

- **Phase 3.gp.0 surrogate framework** (shipped) — referenced by Stage 4 as the optional refinement hook. Phase 3.nl.0 itself uses a no-op pass-through; the integration point exists so 3.nl.1 is a wiring change, not an architectural one.
- **Phase 1.cli.1 / `yee run`** — emitted `yee.toml` must conform to the project-file format `yee run` accepts. Phase 3.nl.0 lane is read-only on the project-file format; if a missing field is needed, surface as a finding to the `yee-cli` lane.
- **LLM SDK choice.** Recommend **Anthropic Messages API** with tool use, invoked from a Python sidecar (the user already has `anthropic` installed per the Phase 3.bo.0 notebook tooling). Rust-side dependency is `reqwest` only; the LLM call is shelled out via `yee-py` so the Rust crate is `cuda`-style feature-gated (`--features llm`) and builds green with no network / no credentials in default CI. **Offline mode** parses a restricted English subset via a hand-written regex / template grammar (`regex` crate, no new dep) and is the gate path in default CI.
- **TOML serializer:** `toml_edit` (already in the workspace lock).

## 12. Phase numbering ladder

- **3.nl.0** — Walking skeleton: one family (rectangular patch), textbook init, no surrogate, offline-mode fallback, 10-prompt gate. **This spec.**
- **3.nl.1** — Surrogate-refinement loop. Stage 4 becomes a real BO call against the Phase 3.gp.0 / 3.bo.0 stack; the surface starts hitting bandwidth / gain targets, not just frequency.
- **3.nl.2** — Additional geometry families: Wilkinson divider, hairpin filter, microstrip line (the existing `mom-002` / `mom-004` / `mom-006` cases). Each new family is a `GeometryFamily` variant plus a textbook synthesis routine.
- **3.nl.3** — Interactive Claude-as-tool agent loop. Multi-turn refinement, the solver becomes a tool the agent calls, the project file evolves across turns. Prompt-injection threat model lives here.
- **3.nl.4** — Production validation matrix: held-out prompts, multi-family stress, the `ROADMAP.md` Phase-3 gate ("End-to-end: text prompt → working design that meets stated specs to within 10% on at least 5 canonical antenna / filter classes").

## References

- Balanis, *Antenna Theory: Analysis and Design*, 4th ed., Wiley 2016, Ch. 14 (rectangular microstrip patch synthesis: `W`, `L`, `y_0`, edge / inset feed).
- Pozar, *Microwave Engineering*, 4th ed., Wiley 2012, §3.8 (microstrip line characteristic impedance — feed-line sizing).
- Anthropic, *Messages API — Tool use and structured outputs* (constrained JSON output via `tools` + `input_schema`).
- `ROADMAP.md` Phase 3 §"Natural-language design surface" (line 150).
- `TECH_STACK.md` §"Config / project files" (TOML for `yee.toml`).
- CLAUDE.md §3 (walking-skeleton first), §4 (validation-gate requirement), §10 (placeholder Green's tolerance posture inherited here).
