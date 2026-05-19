# Natural-language design surface

This tutorial walks through Yee's Phase 3.nl.0 natural-language design
surface end-to-end. You will turn a free-form prompt such as
`"2.4 GHz patch on FR4"` into a structured
[`DesignIntent`][yee-design], apply Balanis Ch. 14 + Pozar §3.8
closed-form synthesis to compute starting dimensions, emit a
deterministic `yee.toml` project file plus its `intent.json` sidecar,
then hand the project file to `yee run` and plot `|S11|`.
The Phase 3.nl.0 walking skeleton is the convenience layer in front
of the existing solver pipeline — see ADR-0028 and ADR-0031 for the
scope decision and implementation plan.

## Goal

Reduce the engineer's path from "I want a 2.4 GHz patch on FR4" to a
Yee project file from a hand-edited TOML (Balanis page-flipping
required) to a one-line CLI invocation. The surface is reproducible:
the prompt is preserved as a header comment, the structured intent is
written next to the TOML, and re-running from `intent.json`
regenerates the TOML byte-identically.

## What's shipping in Phase 3.nl.0

The walking-skeleton surface under [`crates/yee-design/`][yee-design]
covers:

- **One geometry family** — rectangular inset-fed microstrip patch
  (Balanis Ch. 14 synthesis equations 14-1 through 14-20a; Pozar §3.8
  for the 50 Ω feed-line width).
- **One design-intent grammar** — target frequency, substrate
  (named-preset `FR4 / RO4003C / RO5880 / AluminaTC` or explicit
  `{eps_r, h_mm, loss_tangent}`), optional gain and bandwidth targets.
- **Two parse modes** — a deterministic offline regex / template
  parser (`yee_design::offline::parse`, the default-CI path) and an
  Anthropic Messages API tool-use path served by the `yee-py` sidecar
  at `yee.design.from_prompt_llm`.
- **Deterministic emitter** — `yee_design::emit` renders TOML with
  lexicographic key ordering, fixed `{:.6e}` float format, and a
  `# nl-prompt: …` / `# yee-design: <sha256> source=<mode>` header
  pair. Round-trip through `intent.json` is byte-identical.
- **`nl-001` production gate** — 10 canonical prompts, four
  sub-gates (offline / schema / round-trip / solver). The first three
  run in default CI; the solver sub-gate is `#[ignore]`'d pending
  Phase 1.1.1 (see CLAUDE.md §10 on the `MultilayerGreens`
  placeholder).

Compared to authoring `yee.toml` by hand: no opening Balanis to
recover `W`, `L`, `y_0`, and the feed-line width; no risk of
mistyping the substrate `eps_r`; the project file is reproducible
from the saved `intent.json` even after the surface itself moves.

## Prerequisites

- Rust 1.92+ (`rust-toolchain.toml` pin).
- A `yee` binary on `PATH` — `cargo build --release -p yee-cli` and
  add `target/release/` to your `PATH`, or invoke via
  `cargo run --release -p yee-cli --`.
- *Optional, for the LLM path:* Python 3.10+ with the `yee-py`
  bindings installed with the `llm` extra. From the repo root:

  ```bash
  uv venv .venv
  source .venv/bin/activate
  uv pip install maturin pytest numpy
  uv pip install -e crates/yee-py[llm]
  ```

  The `[llm]` extra pulls in `anthropic` and `jsonschema`; default
  CI does not install them, so the offline path always works without
  this step.

## Offline path

The offline parser is the default-CI path and works without network
or credentials. Invoke `yee design` with `--offline`:

```bash
yee design "2.4 GHz patch on FR4" -o /tmp/patch-2g4.toml --offline
```

Expected stdout:

```text
Wrote /tmp/patch-2g4.toml
Wrote /tmp/patch-2g4.toml.intent.json
Resolved design:
  center_frequency = 2.4000 GHz
  substrate eps_r  = 4.4000
  width            = 38.0100 mm
  length           = 29.4216 mm
  inset_offset     = 10.8207 mm
  feed_width       = 3.0590 mm
```

The emitted `/tmp/patch-2g4.toml` is the spec §8 deterministic TOML:

```toml
# nl-prompt: 2.4 GHz patch on FR4
# yee-design: ccb54c85a05c0e688306979793ae362a2010da173e651666f4056c9b992433f3 source=offline

[frequency]
center_hz = 2.400000e9
span_hz = 4.800000e8
sweep_points = 201

[geometry]
feed_width_m = 3.058975e-3
inset_offset_m = 1.082071e-2
length_m = 2.942159e-2
type = "rectangular_inset_patch"
width_m = 3.800997e-2

[substrate]
eps_r = 4.400000e0
h_m = 1.600000e-3
loss_tangent = 2.000000e-2

[[ports]]
id = 1
inset_offset_m = 1.082071e-2
kind = "delta_gap"
location = "feed"
z0_ohm = 5.000000e1
```

The sidecar `/tmp/patch-2g4.toml.intent.json` records the typed
intent so re-running the emitter from this file regenerates the TOML
byte-for-byte:

```json
{"family":"rectangular_patch","target_frequency_hz":2400000000.0,"substrate":{"name":"FR4"},"source_prompt":"2.4 GHz patch on FR4","provenance":{"source":"offline","schema_version":"1","substrate_library_version":"1"}}
```

Swap the prompt for `"5.8 GHz patch on RO4003C"` to land on Rogers
RO4003C at 5.8 GHz; the offline parser handles the ten canonical
prompts in [`crates/yee-design/validation/prompts.toml`][prompts]
deterministically.

## LLM path

Drop `--offline` and set `ANTHROPIC_API_KEY` in the environment to
route Stage 1 through the Anthropic Messages API. The Python sidecar
in `yee-py` (`yee.design.from_prompt_llm`) handles the tool-use call
with the spec §7 JSON schema as the tool's `input_schema`:

```bash
export ANTHROPIC_API_KEY=...
yee design "2.4 GHz patch on FR4 for IoT" -o /tmp/patch-iot.toml
```

> **Note.** The CLI's LLM-path wiring is the Phase 3.nl.0.1 follow-up;
> at this base SHA `yee design` without `--offline` prints a pointer
> to the sidecar and exits non-zero per the R5 escape hatch. Call the
> sidecar from Python today:
>
> ```python
> import yee.design
> intent = yee.design.from_prompt_llm("2.4 GHz patch on FR4 for IoT")
> ```
>
> See [`crates/yee-py/tests/test_design_llm.py`][test-llm] for the
> end-to-end example, gated by the `pytest -m anthropic` mark so the
> default test suite skips it cleanly.

Stage 1 is the only non-deterministic stage. Stages 2–5 (geometry
resolve, Balanis synthesis, optional refinement no-op, emit) are pure
functions of the `DesignIntent` — see spec §5 for the pipeline
diagram.

## Prompt-injection caveat

> If a prompt contains adversarial substrings designed to fool the
> LLM into emitting a different `DesignIntent` (for example a
> frequency far from the requested band), the schema validator will
> reject it and the call will raise `yee.design.SchemaRejectedError`.
> The schema is the source of truth, not the prompt.

The full threat-model discussion lives in spec §10. Phase 3.nl.0 is
local-CLI / local-Python only; web-facing exposure is explicitly out
of scope and lands with the interactive agent loop in Phase 3.nl.3.

## Running the solve

Feed the emitted TOML to `yee run`:

```bash
yee run /tmp/patch-2g4.toml --output /tmp/patch-2g4.s1p
```

> **Status note.** Per CLAUDE.md §10, `yee run` is the Phase-0 stub
> at this point in the roadmap — it acknowledges the project file but
> does not yet drive the planar-MoM solver end-to-end through the
> emitted TOML schema. The `nl-001` solver sub-gate inherits the
> same posture: it is wired but `#[ignore]`'d until Phase 1.1.1
> lands the real `MultilayerGreens`. In the meantime, drive the
> patch solve from Python with `yee.PlanarMoM().run(...)` per
> [Tutorial 2 — Half-wave dipole from Python](02-dipole-from-python.md),
> using the dimensions echoed by `yee design` as inputs.

## Plotting

Once the Touchstone file lands, the `yee plot` subcommand emits a
PNG or SVG of `|S11|` in dB versus frequency:

```bash
yee plot /tmp/patch-2g4.s1p --format db --output /tmp/patch-2g4-s11.png
```

Other kinds: `--format smith` for a Smith chart, `--format phase` for
phase in degrees, `--format both` for the dB and Smith charts in one
invocation. The output extension (`.png` / `.svg`) picks the backend.

## Validation

The Phase 3.nl.0 production gate is `nl-001`, wired at
[`crates/yee-validation/tests/nl_001_canonical_prompts.rs`][nl-001].
The 10 canonical prompts in
[`crates/yee-design/validation/prompts.toml`][prompts] are exercised
against four sub-gates per the spec §9 composition:

| Sub-gate          | Tolerance / assertion                          | Default CI | `#[ignore]` |
|-------------------|------------------------------------------------|------------|-------------|
| offline           | `yee_design::parse_offline(prompt)` succeeds   | yes        | no          |
| schema            | spec §7 schema (frequency / substrate / enums) | yes        | no          |
| round-trip        | `emit → intent.json → emit` byte-identical     | yes        | no          |
| solver (±5 % f)   | `|f_min − f_target| / f_target ≤ 0.05`         | no         | yes         |

The first three sub-gates are sub-second per prompt and always run.
The solver sub-gate is `#[ignore]`'d for the reasons in CLAUDE.md §10
(the `MultilayerGreens` placeholder) and the Phase 1.1.1 dependency
recorded in the validation README. Run the full set with:

```bash
cargo test -p yee-validation --release --test nl_001_canonical_prompts \
    -- --include-ignored
```

See [`crates/yee-design/validation/README.md`][val-readme] for the
sub-gate disposition table and the Phase 1.1.1 retirement plan.

## What's next

Phase 3.nl.0 is deliberately a walking skeleton — one geometry
family, one prompt grammar, no surrogate refinement. The roadmap
sketches the immediate follow-ups:

- **Phase 3.nl.0.1** — wire the LLM path into `yee-cli` directly via
  PyO3 in-process embedding (deferred from R5 because it requires
  opting `yee-cli` into PyO3 and a `python3` toolchain at link time,
  a tech-stack change out of scope for the walking skeleton).
- **Phase 3.nl.1** — Stage 4 of the pipeline (currently a no-op
  refinement pass-through) becomes a real Bayesian-optimization loop
  against the Phase 3.gp.0 / 3.bo.0 stack; the surface starts hitting
  bandwidth and gain targets, not just frequency.
- **Phase 3.nl.2** — additional geometry families (Wilkinson divider,
  hairpin filter, microstrip line — the existing
  `mom-002` / `mom-004` / `mom-006` cases each get a `GeometryFamily`
  variant plus textbook synthesis).
- **Phase 3.nl.3** — interactive Claude-as-tool agent loop, where the
  solver itself becomes a tool the agent invokes. The
  prompt-injection threat model lives here.
- **R6 D-gate retirement** — when Phase 1.1.1 lands the real
  Sommerfeld-integral / multi-image DCIM `MultilayerGreens`, the
  `nl-001 (solver, ±5 % f)` sub-gate flips from `#[ignore]`'d to
  always-on without any change to the test layout.

## References

- **Spec** —
  [`docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`](../../superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md)
  (pipeline architecture, schema, determinism contract,
  validation gate, threat model).
- **Plan** —
  [`docs/superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md`](../../superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md)
  (track-by-track R1-R7 breakdown).
- **ADR-0028** — Phase 3.nl.0 NL design surface scope.
- **ADR-0031** — Phase 3.nl.0 implementation plan.
- **ADR-0035** — Berenger Huygens-surface subgridding (related FDTD
  follow-on context).
- **ADR-0036 / ADR-0037 / ADR-0038** — `mom-002` validation-strategy
  decisions that drive the loose-tolerance posture inherited by the
  `nl-001` solver sub-gate.
- **Balanis, D. M.**, *Antenna Theory: Analysis and Design*, 4th ed.,
  Wiley 2016 — Ch. 14 (rectangular microstrip patch synthesis: `W`,
  `L`, `y_0`, edge / inset feed).
- **Pozar, D. M.**, *Microwave Engineering*, 4th ed., Wiley 2011 —
  §3.8 (microstrip-line characteristic impedance, feed-line sizing).
- **CLAUDE.md §10** — `MultilayerGreens` placeholder status and the
  loose-tolerance posture the `nl-001` solver sub-gate inherits.

[yee-design]: https://github.com/your-org/yee/tree/main/crates/yee-design
[prompts]: https://github.com/your-org/yee/blob/main/crates/yee-design/validation/prompts.toml
[nl-001]: https://github.com/your-org/yee/blob/main/crates/yee-validation/tests/nl_001_canonical_prompts.rs
[val-readme]: https://github.com/your-org/yee/blob/main/crates/yee-design/validation/README.md
[test-llm]: https://github.com/your-org/yee/blob/main/crates/yee-py/tests/test_design_llm.py
