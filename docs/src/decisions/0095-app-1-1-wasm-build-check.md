# ADR-0095: App.1.1 — verify `yee-studio` builds to `wasm32` + CI wasm-readiness gate

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0089 (app architecture), ADR-0092 (App.1.0 `desktop`-feature gating)

---

## Context

App.1.0 (ADR-0092) gated `yee-studio`'s eframe shell behind a default `desktop`
feature so `StudioState` builds eframe-free. The claim "App.1 / the web app is
reachable" rests on `yee-studio --no-default-features` actually compiling to the
`wasm32-unknown-unknown` target — which has **not** been verified (the target was
not installed). App.1.1 verifies it and locks it in with a CI gate, before the
heavier full eframe-web + `trunk` deploy (App.1.2).

## Decision

- Install `wasm32-unknown-unknown` (`rustup target add`) and verify
  `cargo build -p yee-studio --no-default-features --target wasm32-unknown-unknown`
  exits 0 — i.e. the WASM-safe flow logic (`StudioState` + `yee-synth`/`yee-filter`)
  genuinely compiles to WASM.
- Add a CI job (a step in `.github/workflows/ci.yml`, or a small dedicated
  workflow) that installs the target and runs that build, gating future
  regressions of the light-flow WASM-safety on every push.

If the target cannot be installed in this environment (offline), surface that
(escape hatch) — the CI gate is still the deliverable, and the local build is
confirmed on CI.

## Consequences

**Ships:** a CI wasm-build-check (the App.1 web-readiness gate) and a confirmed
`wasm32` build of `yee-studio --no-default-features`. Gate: the workflow step is
present + well-formed; locally `cargo build … --target wasm32-unknown-unknown`
exits 0 (or, if the target won't install here, the agent surfaces it and the CI
job is the gate).

**Not in scope (App.1.2):** the full eframe *web* build (egui canvas/WebGPU
backend for `mod app`), `trunk`/`index.html`, and the static deploy. This
increment only proves + gates the light-flow WASM compile.

**No new runtime dependency.** Lane: `.github/workflows/**` (+ a short build note
in `crates/yee-studio/` docs if useful).

---

## References
- ADR-0089 (desktop+web app), ADR-0092 (`desktop` feature). `rust-toolchain.toml`
  (pinned 1.92); CLAUDE.md §8 (CI layout).
- `docs/superpowers/plans/2026-05-29-app-1-1-wasm-build-check.md`.
