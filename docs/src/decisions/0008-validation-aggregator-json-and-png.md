# ADR-0008: Validation aggregator emits both a JSON report and per-case PNG artifacts

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

The `yee-validation` crate is the project's cross-cutting regression
gate. Phase 1.validation.0 shipped the aggregator skeleton with four
skip-placeholder cases. Phase 1.validation.1 wired `mom-001` (the
half-wave dipole) into it and added a `yee-validate` binary that
prints a human-readable summary or a machine-readable JSON document.
Phase 1.validation.2 added per-case PNG plot emission and a CI job
that uploads the resulting `validation/results/` directory as a
workflow artifact.

The aggregator has two distinct consumers, and the question is what
output formats serve them:

- **Downstream tooling.** A future regression dashboard, a release-
  notes generator, a "did this PR regress mom-001 by more than 1%?"
  bot — all want a stable, structured, diff-able summary they can
  parse without rendering plots. JSON is the obvious choice; the
  project already standardises on `serde_json` and `Touchstone v1.1`
  is the only competing external format, which is too domain-specific
  for a cross-crate aggregator.
- **Human reviewers.** A PR reviewer or a contributor running the
  validation locally wants to *see* what the solver produced. For
  S-parameter cases (`mom-001` today, `mom-002` / `mom-003` later,
  and any future cross-validation case) the canonical visualisations
  are an S₁₁-magnitude trace in dB and a Smith chart. A pure-text
  JSON summary buries the failure mode behind a single
  pass/fail/tolerance number.

Three decision points sat in front of the team:

1. **JSON only, no plots.** Simplest. Defers plotting to a
   downstream tool that reads the JSON. Rejected because the
   downstream tool does not exist yet, and the time-to-insight for a
   reviewer is gated on its existence.
2. **Plots only, no JSON.** Rejected because plots are not
   diff-able. A regression in `mom-001`'s real part by 0.4% is
   invisible in a Smith chart and obvious in JSON.
3. **Both, written side-by-side under `validation/results/`.** The
   aggregator's `Report` is `serde::Serialize` and the binary writes
   it as `validation/results/report.json`. Each case that opts into
   plotting writes one or more PNGs to
   `validation/results/<case>/<plot>.png`. CI uploads the whole
   `validation/results/` tree as a workflow artifact with 14-day
   retention.

The plotting backend is `yee-plotters`, the existing static-export
crate that wraps `plotters` with a Yee-flavoured S-parameter and
Smith-chart API. This adds `yee-plotters` (and its transitive
dependency on `plotters` + `libfontconfig1-dev` on Linux) to
`yee-validation`'s dependency graph.

The compile-time cost of `yee-plotters` was weighed against the
benefit. The `plotters` tree is substantial (font handling, bitmap
rendering, SVG), and adding it to `yee-validation` means anyone
running `cargo test -p yee-validation` pays the compile cost. The
mitigation considered was a `yee-plotters` feature gate on
`yee-validation` so a non-feature build skips the plotter
dependency. The decision was to **not** feature-gate it: the
aggregator's whole purpose is to produce reviewer-facing artefacts,
and a configuration where the aggregator silently skips plots is a
bug-shaped attractor. The compile cost is paid once per build and
amortised across every validation run.

The CI side has its own decision. The naive choice was to publish
the PNGs to GitHub Pages alongside the mdBook (`docs/`). Pages
deployment is **gated by a repo-settings action the maintainer must
take** (Source: GitHub Actions), and that action is already a known
"first-time setup" failure mode for the docs workflow. Coupling
plot publication to Pages would tie validation artefact visibility
to a separate manual step. The chosen path uploads
`validation/results/` as a workflow artifact, which:

- Works on any GitHub repository out of the box, no Pages setting
  required.
- Has a configurable retention window (14 days, matching the default
  practice for non-release artefacts).
- Is downloadable directly from the run page by a reviewer, without
  any HTML routing layer.

A Pages-published validation gallery can be added later as a separate
follow-up that does not block the artifact path.

## Decision

`yee-validation::Report::run_all()` returns a `Report` struct that is
`#[derive(serde::Serialize, serde::Deserialize)]`, with one
`CaseResult` per case. The schema is:

```rust
#[derive(Serialize, Deserialize)]
pub struct Report {
    pub yee_version: String,
    pub git_sha: String,
    pub generated_at: String,           // RFC 3339 UTC
    pub cases: Vec<CaseResult>,
}

#[derive(Serialize, Deserialize)]
pub struct CaseResult {
    pub id: String,                     // "mom-001"
    pub status: Status,                 // Pass | Fail | Skip
    pub metrics: serde_json::Value,
    pub tolerance: serde_json::Value,
    pub plot_paths: Vec<PathBuf>,       // paths under validation/results/
}
```

The `yee-validate` binary writes the serialised `Report` to
`validation/results/report.json`. The same binary delegates to each
case's `emit_plots(out_dir)` hook; for `mom-001` that produces
`validation/results/mom-001/s11_db.png` and
`validation/results/mom-001/smith.png`.

`yee-validation`'s `Cargo.toml` takes a hard dependency on
`yee-plotters` (no feature gate). The crate's integration test
asserts that the PNG files exist and are non-trivial in size
(`> 1 KiB`), which is a cheap proxy for "the plot actually
rendered" without parsing PNG internals.

The CI job uploads `validation/results/` as a workflow artifact with
14-day retention. Pages deployment of validation plots is **not**
wired up; it is deferred to a separate follow-up.

## Consequences

**What becomes easier:**

- Machine-readable summary for downstream tooling. A future
  regression bot can fetch `report.json` from the artifact and diff
  it against the previous run without any rendering pipeline.
- Human-readable plots for review. A reviewer can click through to
  the workflow artifact, unzip, and see the S₁₁ trace and Smith
  chart for any case that opts in.
- Plot publication is decoupled from Pages activation. The
  validation pipeline works on a fresh fork with zero repo-settings
  configuration.
- Adding a new validation case is uniform: implement
  `Case::run() -> CaseResult`, optionally implement
  `Case::emit_plots(out_dir)`, register in `run_all()`.

**What becomes harder:**

- `yee-validation` now compile-depends on `yee-plotters`, which
  pulls `plotters` into the dependency graph. Local builds without
  `libfontconfig1-dev` will fail on this crate. The CI already
  installs `libfontconfig1-dev` and `pkg-config`; local setups need
  the same packages (documented in the root `CLAUDE.md` §7).
- Cold `cargo build -p yee-validation` is slower than the pre-Phase
  1.validation.2 baseline by the cost of compiling `plotters` and
  its transitive font / image dependencies.
- Workflow artifacts have a 14-day retention by default. Older
  artefacts are not addressable. A long-term archive needs the
  separate Pages-deploy follow-up.

**What's now closed off:**

- Feature-gating the plot output on `yee-validation`. The aggregator
  always tries to emit plots for cases that opt in; a build without
  plotting support is not a supported configuration.
- Publishing validation plots through the docs Pages deploy in the
  same workflow as the artifact upload. The two channels are
  intentionally separate.

## References

- `crates/yee-validation/` — the aggregator crate.
- `crates/yee-validation/src/report.rs` — `Report` and `CaseResult`
  schema.
- `crates/yee-validation/tests/png_artifacts.rs` — integration test
  asserting PNGs exist and exceed the 1 KiB floor.
- `crates/yee-validate/` — the binary that runs the aggregator and
  writes JSON + PNGs.
- `crates/yee-plotters/` — the `plotters`-backed static-export crate
  used for the S₁₁ dB and Smith-chart renderings.
- `.github/workflows/ci.yml` — the `validation-artifacts` job that
  uploads `validation/results/` (retention 14 d).
- ADR-0001 — license; `plotters` is Apache-2.0 / MIT, GPL-3.0
  compatible.
