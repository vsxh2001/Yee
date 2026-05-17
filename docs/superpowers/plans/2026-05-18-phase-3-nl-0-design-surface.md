# Phase 3.nl.0 — Natural-Language Design Surface (Walking Skeleton) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan track-by-track.

**Companion spec:** `docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`
**Base SHA:** `f217df0` (Track PPPPP merge).
**Target phase:** 3.nl.0 only. 3.nl.1–3.nl.4 are explicitly deferred — see §"Out of scope".
**Tech-stack additions:** `toml_edit` (already in lock for `yee-cli`); `regex` (already a transitive dep, promote to direct in `yee-design`); `serde_json` (already in lock); Python sidecar adds `anthropic>=0.40` to `yee-py`'s optional dev extras (gated behind a `pytest` mark).

---

## Goal

Phase 3.nl.0 ships a one-shot natural-language → project-file surface: a free-form prompt like *"2.4 GHz inset-fed patch on RO4003C with ≥ 100 MHz bandwidth"* parses into a typed `DesignIntent`, applies the Balanis Ch. 14 inset-fed-patch synthesis equations to produce starting dimensions, and emits a `yee.toml` that the existing `yee run` pipeline accepts unchanged. One geometry family (rectangular inset-fed patch), one design-intent grammar, two execution modes (LLM-mediated via Anthropic Messages API; deterministic offline regex fallback), and a 10-prompt validation gate enforcing `|S11|` minimum within ±5% of the requested frequency. No surrogate refinement, no optimizer loop, no multi-turn agent — those are 3.nl.1+ per the spec.

## Pre-flight — `yee-surrogate` status reconciliation

CLAUDE.md §2 (timestamp SHA `33e28db`) states the `yee-surrogate` crate "has not landed." At the current base SHA `f217df0`, **`crates/yee-surrogate/` does exist** and `yee-py` already re-exports `yee.surrogate.GaussianProcess` (Track NNNNN finding #5 confirmed). The plan therefore routes Stage 4's optional refinement hook through the existing `yee-surrogate` crate (no-op pass-through in 3.nl.0; 3.nl.1 wires the BO call). It does **not** assume a new `yee-design` consumer of `yee-surrogate` types — Stage 4 is a `trait Refiner { fn refine(&self, est: InitialEstimate) -> InitialEstimate; }` whose default impl is the identity. CLAUDE.md should be updated after merge; that is **out of lane** here and surfaced as a finding.

## File structure

| File | Action | Responsibility |
|------|--------|----------------|
| `Cargo.toml` (workspace root) | Modify | Add `yee-design` workspace member; promote `regex` / `toml_edit` to `[workspace.dependencies]` if not already pinned. |
| `crates/yee-design/Cargo.toml` | Create | New crate, defaults: `serde`, `serde_json`, `toml_edit`, `regex`, `thiserror`. No LLM dep — Stage 1 LLM call lives in the Python sidecar. |
| `crates/yee-design/src/lib.rs` | Create | Crate root, `#![forbid(unsafe_code)]`, `#![warn(missing_docs)]`, re-exports. |
| `crates/yee-design/src/intent.rs` | Create | `DesignIntent`, `GeometryFamily`, `Substrate`, `NamedSubstrate`, `Provenance` types (spec §6). |
| `crates/yee-design/src/estimate.rs` | Create | `InitialEstimate::from_intent` — Balanis Ch. 14 synthesis equations. |
| `crates/yee-design/src/emit.rs` | Create | `ProjectFile` + deterministic `toml_edit` emitter (sorted keys, fixed `{:.6e}` floats). |
| `crates/yee-design/src/offline.rs` | Create | Deterministic regex / template-grammar parser used by `--offline` and by default in CI. |
| `crates/yee-design/substrates.toml` | Create | Named-substrate library (FR4, RO4003C, RO5880, AluminaTC) — versioned at top of file. |
| `crates/yee-design/validation/prompts.toml` | Create | 10 canonical prompts (spec §9). |
| `crates/yee-design/tests/balanis_example.rs` | Create | Unit test against Balanis published example. |
| `crates/yee-design/tests/emit_roundtrip.rs` | Create | Emit → re-parse round-trip + byte-identical regen from `intent.json`. |
| `crates/yee-py/src/lib.rs` | Modify | Add `mod design;`, register `yee.design` submodule. |
| `crates/yee-py/src/design.rs` | Create | Python sidecar: Anthropic Messages tool-use call → `DesignIntent` JSON; exposed as `yee.design.from_prompt_llm`. |
| `crates/yee-py/tests/test_design_llm.py` | Create | `pytest.mark.anthropic` test gated by `ANTHROPIC_API_KEY`. |
| `crates/yee-cli/src/main.rs` | Modify | Wire `yee design "<prompt>" --output <p.toml> [--offline] [--model <id>]` subcommand. |
| `crates/yee-cli/tests/design_offline.rs` | Create | Integration test exercising the 10 canonical prompts via `--offline`. |
| `crates/yee-validation/src/lib.rs` | Modify | `run_nl_001_canonical_prompts` driver; new `nl-001` cases. |
| `crates/yee-validation/tests/nl_001_canonical_prompts.rs` | Create | The production gate (10-prompt end-to-end). |
| `crates/yee-design/validation/README.md` | Create | `nl-001 (schema)`, `nl-001 (round-trip)`, `nl-001 (solver)`, `nl-001 (offline)` rows. |
| `docs/src/tutorials/04-nl-design-surface.md` | Create | mdBook tutorial (Step 7, optional). |

No changes to `yee-core`, `yee-fdtd`, `yee-mom`, `yee-mesh`, `yee-cuda`, `yee-gui`, `yee-plotters`, `yee-surrogate`. `yee-design` is a new crate; the `unsafe_code` floor is preserved.

## Step ladder

### Step 1 (Track R1) — `yee-design` crate scaffold + `DesignIntent` types

- **Brief:** Create `crates/yee-design/` workspace member with `Cargo.toml`, `src/lib.rs`, `src/intent.rs`. Implement the spec §6 type surface (`DesignIntent`, `GeometryFamily`, `Substrate`, `NamedSubstrate`, `Provenance`) with `serde::{Serialize, Deserialize}` derives. Load `substrates.toml` at build time via `include_str!` + `toml::from_str` into a `static OnceLock<SubstrateLibrary>`. Lane crosses into the workspace root `Cargo.toml`; that is **explicit cross-lane** — call it out in the agent's report. No other crate touched.
- **Lane:** `Cargo.toml` (root, members array only), `crates/yee-design/**`.
- **Base SHA dep:** none — branches off `f217df0` directly.
- **DoD:** crate builds clean (`cargo check -p yee-design`); `DesignIntent` round-trips through `serde_json::{to_string, from_str}` byte-identically (property test, 100 random samples); `cargo clippy -p yee-design --all-targets -- -D warnings` exits 0; `cargo doc -p yee-design --no-deps` is `missing_docs`-clean.
- **Verification:** `cargo check -p yee-design && cargo clippy -p yee-design --all-targets -- -D warnings && cargo test -p yee-design --release intent && cargo doc -p yee-design --no-deps` exits 0.
- **Escape hatch:** blocked > 15 min on workspace `Cargo.toml` resolver / lock churn → revert lock changes, re-run `cargo check --workspace`, take `--theirs` per CLAUDE.md §5. Do not hand-merge `Cargo.lock`.
- **LOC:** ~220.

### Step 2 (Track R2) — Balanis Ch. 14 initial-estimate calculator

- **Brief:** Implement `crates/yee-design/src/estimate.rs::InitialEstimate::from_intent(&DesignIntent) -> Result<InitialEstimate, Error>` as a pure function. For `GeometryFamily::RectangularPatch`: compute `W` (Balanis 14-6), effective `eps_reff` (Balanis 14-1), edge-effect `ΔL` (Balanis 14-2), physical `L` (Balanis 14-3), inset-feed offset `y_0` (Balanis 14-20a using 50 Ω feed-line `Z_0`), and feed-line width from Pozar §3.8 synthesis at the substrate stack-up. `InitialEstimate` is a flat struct of named `f64` dimensions in metres + a copy of the resolved `Substrate`. Unit-test against Balanis Example 14.1 (`f = 10 GHz`, `ε_r = 2.2`, `h = 1.588 mm` → published `W ≈ 11.86 mm`, `L ≈ 9.06 mm`). Tolerance ±0.5% on each dimension.
- **Lane:** `crates/yee-design/src/estimate.rs`, `crates/yee-design/tests/balanis_example.rs`.
- **Base SHA dep:** Step 1 merged.
- **DoD:** Balanis Example 14.1 test exits 0 within ±0.5% on `W`, `L`, `y_0`; one additional spot test at 2.4 GHz on FR-4 (`ε_r = 4.4`, `h = 1.6 mm`) verifies finite, positive, plausible-magnitude dimensions; pure-function property — `from_intent` called twice with the same input returns bit-identical output (no `f64::NAN` paths).
- **Verification:** `cargo test -p yee-design --release estimate` exits 0; `cargo clippy -p yee-design --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on the inset-offset transcendental (Balanis 14-20a uses `cos²(π y_0 / L) = R_in / R_edge`) → start with `y_0 = 0.3 · L` as a closed-form lower bound, mark `// TBD: tighten when Balanis 14-20a numeric root agreed`, file the inset-offset accuracy as Phase 3.nl.0.1. The frequency gate (spec §9) does not need a precise `y_0` — `|S11|` minimum vs frequency is dominated by `L`.
- **LOC:** ~250.

### Step 3 (Track R3) — deterministic project-TOML emitter + round-trip

- **Brief:** Implement `crates/yee-design/src/emit.rs::emit(&InitialEstimate, &DesignIntent) -> ProjectFile`. Use `toml_edit::DocumentMut` so key insertion order is controllable; sort keys lexicographically inside each table; format every `f64` via `format!("{:.6e}", x)`. The first two lines of the emitted TOML are spec-§8 metadata comments (`# nl-prompt: ...` then `# yee-design: <sha256(intent_json)> <provenance>`). `ProjectFile { toml: String, intent_json: String }`. Round-trip test: emit, re-parse the TOML via the same `toml` crate `yee run` uses, assert all geometry / substrate / frequency fields survive; second round-trip test: re-emit from the saved `intent.json` and assert **byte-identical** output (the spec §8 determinism gate).
- **Lane:** `crates/yee-design/src/emit.rs`, `crates/yee-design/tests/emit_roundtrip.rs`.
- **Base SHA dep:** Step 2 merged.
- **DoD:** round-trip test passes for at least two intents (2.4 GHz FR-4, 5.8 GHz RO4003C); byte-identical regeneration test passes (`emit` twice → `assert_eq!`); `serde_json` of `DesignIntent` matches the emitted `intent.json` field-for-field. The TOML schema field names must align with whatever `yee-cli`'s `Run` subcommand consumes today — if `yee run` is still the Phase-0 stub at this base SHA, document the field names in a `// TBD: confirm with yee-cli once Run parser lands` comment and pick a forward-compatible shape (`[geometry]`, `[substrate]`, `[frequency]`, `[ports]`).
- **Verification:** `cargo test -p yee-design --release emit` exits 0.
- **Escape hatch:** blocked > 15 min on TOML field-name conflict with the `yee run` parser → check the `examples/patch-2g4` example for the existing project-file shape; mirror it. If `yee run` is still a stub, lock in the forward shape and surface the cross-lane mismatch as a finding for `yee-cli`. Do not extend the `yee-cli` parser in this lane.
- **LOC:** ~280.

### Step 4 (Track R4) — Python sidecar: Anthropic Messages API tool-use call

- **Brief:** Implement `crates/yee-py/src/design.rs` exposing `yee.design.from_prompt_llm(prompt: str, model: str | None = None, api_key: str | None = None) -> dict`. Calls the Anthropic Messages API with a `tools` array containing one tool whose `input_schema` is the JSON schema from spec §7 (load from `crates/yee-design/src/intent.rs` via a `pub const INTENT_SCHEMA: &str = include_str!("intent_schema.json")` so both Rust and Python read the same artefact). On tool-use response, validate the returned JSON against the schema (use `jsonschema` Python package, listed under `yee-py`'s `[project.optional-dependencies] llm`). On schema-validation failure, re-prompt once with the schema appended to the user message; second failure → raise `yee.design.SchemaRejectedError`. Secrets policy: read `ANTHROPIC_API_KEY` from env unless the caller passes `api_key=`; never log it; never include it in `Provenance`. `Provenance` records model id + temperature + schema version, no secrets.
- **Lane:** `crates/yee-py/src/lib.rs` (one `mod design;` + registration), `crates/yee-py/src/design.rs`, `crates/yee-py/tests/test_design_llm.py`, `crates/yee-design/src/intent_schema.json` (the schema file; created here so the Rust crate's `include_str!` resolves — spec §7 makes this an intentional cross-lane shared artefact).
- **Base SHA dep:** Step 1 merged (uses `DesignIntent` shape + schema file).
- **DoD:** `pytest -m anthropic crates/yee-py/tests/test_design_llm.py` passes against the live API when `ANTHROPIC_API_KEY` is set (three prompts: a clean one, an under-specified one that should re-prompt, a hostile one that should `SchemaRejectedError`); `pytest -m "not anthropic"` skips the LLM tests cleanly; `maturin develop -m crates/yee-py/Cargo.toml` succeeds without `anthropic` installed. CI default does **not** run the `anthropic`-marked tests.
- **Verification:** `cd crates/yee-py && pytest -m "not anthropic"` exits 0; with `ANTHROPIC_API_KEY` set, `pytest -m anthropic` exits 0.
- **Escape hatch:** blocked > 15 min on the Anthropic SDK tool-use response shape → fall back to a plain Messages call with the schema in the system prompt and a `json.loads(message.content[0].text)` parse + `jsonschema.validate`. Mark `// TBD: migrate to native tool_use once SDK shape confirmed`. Do not pin a specific SDK minor version in `yee-py`'s required deps — keep the `anthropic` dep in the `llm` extra so default `maturin develop` is unaffected.
- **LOC:** ~310 (Python ~180, Rust binding glue ~80, schema JSON ~50).

### Step 5 (Track R5) — `yee design` CLI subcommand wiring R1–R3 (+ R4 optional)

- **Brief:** Add the `Design { prompt: String, output: PathBuf, offline: bool, model: Option<String> }` arm to `yee-cli`'s `Command` enum. Default path: if `--offline` or if `ANTHROPIC_API_KEY` is unset, invoke `yee_design::offline::parse(&prompt)?`; else shell out to the Python sidecar via `pyo3` in-process (re-using the existing `yee-py` embedding pattern — confirm one exists; if not, fall back to `std::process::Command::new("python") -c "import yee; ..."`). Both paths produce a `DesignIntent`; pass it to `yee_design::estimate::InitialEstimate::from_intent` and then `yee_design::emit::emit`. Write `<output>` and `<output>.intent.json`. Print the resolved dimensions to stdout for the engineer to eyeball. The 10 canonical prompts (`crates/yee-design/validation/prompts.toml`) must succeed end-to-end with `--offline`.
- **Lane:** `crates/yee-cli/src/main.rs`, `crates/yee-cli/tests/design_offline.rs`, `crates/yee-design/src/offline.rs`.
- **Base SHA dep:** Steps 1, 2, 3 merged. R4 is optional — `--offline` path must work without it.
- **DoD:** `yee design "2.4 GHz patch on FR4" -o /tmp/p.toml --offline` writes `/tmp/p.toml` and `/tmp/p.toml.intent.json` (file existence + non-empty + parses as TOML); all 10 canonical prompts succeed; integration test asserts the emitted TOML contains the requested centre frequency within ±0.1% (the offline parser's job; solver gate is R6).
- **Verification:** `cargo test -p yee-cli --release design_offline` exits 0; `cargo clippy -p yee-cli --all-targets -- -D warnings` exits 0.
- **Escape hatch:** blocked > 15 min on `pyo3` in-process embedding from `yee-cli` → drop in-process call, shell out via `std::process::Command::new("python3").args(["-c", ...])`. Document the choice as a `// TBD: revisit in-process embed in 3.nl.0.1`. Do not introduce a new Rust LLM client in `yee-cli` — that is a tech-stack change requiring a `TECH_STACK.md` update, out of scope here.
- **LOC:** ~280 (CLI arm ~80, offline parser ~150, integration test ~50).

### Step 6 (Track R6) — `nl-001` production validation gate (10-prompt end-to-end)

- **Brief:** Implement spec §9 in `crates/yee-validation/`. For each of the 10 canonical prompts: run `yee design --offline` to produce a project TOML, run `yee run <project.toml>` to obtain `|S11|(f)`, find the frequency `f_min` at which `|S11|` is minimum across the swept band, assert `|f_min − f_target| / f_target ≤ 0.05`. Schema gate: every prompt parses to a valid `DesignIntent` against the spec §7 schema (use the same `jsonschema` validation crate via FFI or a pure-Rust `jsonschema` crate — prefer the latter to keep CI dep-free of Python). Round-trip gate: hash-equal regen from `intent.json`. **Inherits the loose-tolerance posture from CLAUDE.md §10** — the `MultilayerGreens` placeholder means the solver gate is checking the surface, not the accuracy floor.
- **Lane:** `crates/yee-validation/src/lib.rs`, `crates/yee-validation/tests/nl_001_canonical_prompts.rs` (create), `crates/yee-design/validation/README.md` (create).
- **Base SHA dep:** Step 5 merged.
- **DoD:** all 10 prompts pass the four sub-gates (schema, round-trip, solver, offline); wall-time `< 30 min` `--release` (10 prompts × ~3 min each on the existing `mom-003` patch path); hardware-gate behind `#[ignore]` if overrun. Validation README rows: `nl-001 (schema)`, `nl-001 (round-trip)`, `nl-001 (solver, ±5% f)`, `nl-001 (offline)`.
- **Verification:** `cargo test -p yee-validation --release nl_001_canonical_prompts` exits 0.
- **Escape hatch:** blocked > 15 min with > 5% frequency error on any prompt → first verify the offline parser extracted the correct `target_frequency_hz` (it almost certainly did — the regex is deterministic); then check the Balanis `L` formula in R2. If still > 5%, mark the affected prompt `#[ignore]` with a `// TBD: investigate <prompt> frequency drift` comment and surface as a finding. Do **not** loosen the 5% gate — spec §9 is firm.
- **LOC:** ~280.

### Step 7 (Track R7, optional) — mdBook tutorial `04-nl-design-surface.md`

- **Brief:** Write `docs/src/tutorials/04-nl-design-surface.md` mirroring the structure of `01-microstrip-line.md` and `02-dipole-from-python.md`: motivation, prompt → CLI invocation → emitted TOML → `yee run` → `yee plot` for `|S11|`. Cover both `--offline` and LLM paths. Surface the spec §10 prompt-injection caveat in a callout. Link to the spec from the tutorial header and from `docs/src/SUMMARY.md`. No code beyond shell invocations; this is the engineer-facing on-ramp.
- **Lane:** `docs/src/tutorials/04-nl-design-surface.md` (create), `docs/src/SUMMARY.md` (modify).
- **Base SHA dep:** Step 5 merged.
- **DoD:** `mdbook build docs/` exits 0; tutorial renders; the worked-example prompt + emitted TOML round-trip when copy-pasted (smoke-test by hand once before commit).
- **Verification:** `mdbook build docs/` exits 0; no broken links per `mdbook-linkcheck` if installed.
- **Escape hatch:** blocked > 15 min on the LLM-path screenshot / output capture → omit the LLM section, link to `crates/yee-py/tests/test_design_llm.py` instead, ship the offline-only walkthrough.
- **LOC:** ~180 (Markdown).

## Track-letter sequencing

R1 must land first — it defines `DesignIntent` for everything downstream. R2 and R3 are sequential on R1 (R3 consumes `InitialEstimate` from R2). R4 is independent of R2/R3 — it only needs `DesignIntent` from R1, so it can run in parallel with R2 + R3. R5 depends on R1 + R2 + R3 (R4 is optional at the CLI layer because `--offline` is the default-CI path). R6 depends on R5. R7 depends on R5 (the tutorial demonstrates the working CLI).

Critical path: `R1 → R2 → R3 → R5 → R6` (5 sequential merges) + `R4` parallel to R2/R3 + `R7` parallel to R6. Within CLAUDE.md §5's "up to 5 parallel agents" envelope: at any instant the active set is ≤ 3 (R2 ‖ R3 ‖ R4 after R1; then R5 alone; then R6 ‖ R7).

```
R1 ──┬── R2 ── R3 ── R5 ──┬── R6
     ├── R4 ───────────── ┘    └── R7
     └─ (R4 also unblocks R5's LLM path; R5 ships with offline as default)
```

## Validation rollup

| Gate | Step | Tolerance | Run-time |
|------|------|-----------|----------|
| Schema — 10 prompts parse to valid `DesignIntent` | R6 | spec §7 schema, hard reject | `< 1 s` |
| Round-trip — `emit` → `intent.json` → `emit` byte-identical | R6 (also covered by R3 unit test) | byte-identical | `< 1 s` |
| **`nl-001` solver gate** — `|S11|` minimum within ±5% of requested freq | R6 | `|f_min − f_target| / f_target ≤ 0.05` | `< 30 min` `--release`, hardware-gate if overrun |
| Offline parser — all 10 prompts succeed without `ANTHROPIC_API_KEY` | R6 | exit 0 on every prompt | `< 1 s` |

Per CLAUDE.md §4 "No solver feature ships without a published-benchmark validation case" — `nl-001` is the gate. The published benchmark *for the surface* is the spec §9 prompt set; the **solver** underneath continues to inherit the loose-tolerance posture from `mom-003` until Phase 1.1.1 lands the real `MultilayerGreens`. This is documented at the README row.

## Lane / file inventory

| Step | Files |
|------|-------|
| R1 | `Cargo.toml` (root, cross-lane — workspace members only), `crates/yee-design/Cargo.toml` (create), `crates/yee-design/src/lib.rs` (create), `crates/yee-design/src/intent.rs` (create), `crates/yee-design/substrates.toml` (create) |
| R2 | `crates/yee-design/src/estimate.rs` (create), `crates/yee-design/tests/balanis_example.rs` (create) |
| R3 | `crates/yee-design/src/emit.rs` (create), `crates/yee-design/tests/emit_roundtrip.rs` (create) |
| R4 | `crates/yee-py/src/lib.rs` (modify), `crates/yee-py/src/design.rs` (create), `crates/yee-py/tests/test_design_llm.py` (create), `crates/yee-design/src/intent_schema.json` (create — shared artefact) |
| R5 | `crates/yee-cli/src/main.rs` (modify), `crates/yee-cli/tests/design_offline.rs` (create), `crates/yee-design/src/offline.rs` (create), `crates/yee-design/validation/prompts.toml` (create) |
| R6 | `crates/yee-validation/src/lib.rs` (modify), `crates/yee-validation/tests/nl_001_canonical_prompts.rs` (create), `crates/yee-design/validation/README.md` (create) |
| R7 | `docs/src/tutorials/04-nl-design-surface.md` (create), `docs/src/SUMMARY.md` (modify) |

Cross-lane consumers: R1 touches the workspace root `Cargo.toml`; R4 creates a shared schema JSON under `crates/yee-design/` from the `yee-py` lane (explicit because both Rust and Python read the same file — keep the schema as the single source of truth). Both crossings are documented in their step briefs; out-of-lane edits beyond these must be surfaced as findings, not fixed.

## Risk register

1. **LLM hallucination produces invalid `DesignIntent`** (spec §10). Mitigation: the LLM is constrained to emit JSON conforming to the spec §7 schema (`tools` + `input_schema`); a schema-rejection path re-prompts once, then aborts with `SchemaRejectedError` and a `--offline` suggestion. **Surfaces in R4** as the schema-validation step in `from_prompt_llm`. The downstream Stage 3 dimension computation is deterministic — the LLM never proposes `W` or `L`.
2. **Substrate-library drift under a saved `intent.json`** (spec §10). The substrate library has a version string at the top of `substrates.toml`; `Provenance` records it; Stage 2 errors loudly on a version mismatch. **Surfaces in R1** as the `SubstrateLibrary::current_version()` constant + the R3 round-trip test (which is a regression on version-string preservation).
3. **API-key handling and secrets policy** (spec §10, prompt-injection footnote). Mitigation: never log `ANTHROPIC_API_KEY`; never include it in `Provenance`; never persist it to disk. Tests use a dummy key (`os.environ["ANTHROPIC_API_KEY"] = "test-..."`) but only the live-API tests behind `pytest -m anthropic` actually transact. **Surfaces in R4** as the explicit env-var read + the `pytest` mark gating CI exposure.
4. **CI runs without network access** — the default CI lane has no `ANTHROPIC_API_KEY`. Mitigation: the offline parser handles all 10 canonical prompts deterministically; the solver gate (R6) runs `--offline` so it works headless. **Surfaces in R5** as the default behaviour and in R6 as the gate path.
5. **`yee run` not yet a complete parser at base SHA `f217df0`** (CLAUDE.md §10 records it as the Phase-0 stub at SHA `33e28db`; status at `f217df0` is unverified by this plan). Mitigation: R3 picks a forward-compatible TOML shape and surfaces any `yee run`-parser mismatch as a cross-lane finding; R6 may need to invoke the solver via `yee_validation::run_mom_003_patch`-style helpers rather than the CLI entry point until `yee run` lands. **Surfaces in R3 and R6** as documented `// TBD` markers if hit.
6. **Eval-set drift / overfit-to-10-prompts** (spec §10). Phase 3.nl.0 explicitly accepts overfit — that is what "walking skeleton" means in CLAUDE.md §3. The 3.nl.4 held-out matrix is the long-term mitigation. **Documented in §"Out of scope"**, not a 3.nl.0 risk.
7. **CLAUDE.md §2 contradicts the actual workspace contents** re: `yee-surrogate`. The plan's pre-flight identifies this. Mitigation: route Stage 4 (no-op refiner) through `yee_surrogate`'s public surface; surface "update CLAUDE.md §2" as a cross-lane finding after merge. **Materialises in pre-flight** and persists as a documentation finding.

## Out of scope

Explicit non-goals for this plan, per spec §2 and §12:

- **No surrogate-refinement loop** — Stage 4 in the spec is a no-op pass-through; the real BO call lands in Phase 3.nl.1.
- **No additional geometry families** — only rectangular inset-fed patch. Wilkinson, hairpin, microstrip line are Phase 3.nl.2.
- **No interactive / multi-turn agent loop** — one-shot only; Phase 3.nl.3.
- **No web-facing UI** — local CLI + local Python only. Prompt-injection threat model is deferred with the agent loop.
- **No held-out eval matrix** — the 10 canonical prompts overfit by design. Phase 3.nl.4.
- **No bandwidth / gain target hitting** — the surface emits the textbook *starting point*; meeting bandwidth / gain is the surrogate-refinement job in 3.nl.1. The R6 gate checks frequency only.
- **No CLAUDE.md edit in-lane** — the §2 "yee-surrogate has not landed" stale statement is surfaced as a finding for the docs lane to fix post-merge.
- **No new Rust LLM client** — Anthropic calls go through the Python sidecar in `yee-py`. A direct Rust `reqwest` client would be a tech-stack change requiring `TECH_STACK.md` review; deferred.

## Pre-flight installs

```bash
# Already in the standard toolchain — confirm before R4:
uv pip install anthropic jsonschema pytest

# R7 (optional):
cargo install mdbook --locked
```

`anthropic` and `jsonschema` live under `yee-py`'s `[project.optional-dependencies] llm`; default CI does not install them. Local development does (one-line `uv pip install -e crates/yee-py[llm]`).

## Final verification

```bash
cargo build  -p yee-design -p yee-cli -p yee-validation
cargo clippy -p yee-design -p yee-cli -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-design --release
cargo test   -p yee-cli --release
cargo test   -p yee-validation --release nl_001
cargo fmt    --check --all
cargo doc    --no-deps -p yee-design
cd crates/yee-py && pytest -m "not anthropic"
```

All eight must exit 0. Existing `mom-001` / `mom-003` / FDTD gates stay green — `yee-design` is opt-in and no existing crate's behaviour changes.

## Estimated total

- LOC: ~1 800 (R1 ~220, R2 ~250, R3 ~280, R4 ~310, R5 ~280, R6 ~280, R7 ~180).
- Wall-time per agent: 4–6 days end-to-end at one-engineer pace. Critical path R1→R2→R3→R5→R6 is ~4 days; R4 and R7 add ~1 day in parallel.
- Risk concentration: R4 (Anthropic SDK shape) and R6 (10-prompt frequency gate). Both have explicit escape hatches that preserve merge throughput.
