# App.1.1 — wasm32 build check + CI gate — Plan

**ADR:** ADR-0095 · **Date:** 2026-05-29 · (combined spec+plan — small CI/build increment)

## Lane
`.github/workflows/**` (+ optionally a short build note under `crates/yee-studio/`
docs). Out of lane (yee-studio source, any other crate) → finding, not fix.

## Base
Worktree `worktrees/wasm`, branch `feature/app-1-1-wasm-build`, base `bde5cfe`.

## Pattern files
- `.github/workflows/ci.yml` — the existing Rust CI matrix; add a `wasm-build`
  job (or step) mirroring its `rustup`/`cargo` style + the `actions/checkout` +
  toolchain setup it uses.

## Steps
1. Locally: `rustup target add wasm32-unknown-unknown`, then
   `cargo build -p yee-studio --no-default-features --target wasm32-unknown-unknown`.
   Confirm exit 0 (StudioState + yee-synth/yee-filter compile to WASM, no eframe).
   If the target install fails (offline) → escape hatch: surface it; the CI job is
   still the deliverable.
2. Add a CI job `wasm-build` in `.github/workflows/ci.yml`: checkout, install Rust
   1.92 + `wasm32-unknown-unknown` target, run
   `cargo build -p yee-studio --no-default-features --target wasm32-unknown-unknown`.
   Match the file's existing job/style; do NOT alter the existing jobs' behaviour.
3. (Optional) a one-line note in `crates/yee-studio/README` or lib doc on the
   wasm build command — only if it stays in-lane/light.

## Verify
```
rustup target add wasm32-unknown-unknown
nice -n 19 cargo build -p yee-studio --no-default-features --target wasm32-unknown-unknown --jobs 2   # exit 0
# YAML well-formed: python3 -c "import yaml,sys; yaml.safe_load(open('.github/workflows/ci.yml'))"
```
(The CI job itself runs on push — locally only the cargo build + YAML validity are checkable.)

## Escape hatch
`rustup target add` fails / no network → STOP the local-build step, still add the
CI job (the gate), and surface that the local wasm build could not be confirmed
here. Do NOT touch yee-studio source to force it.

## Done when
`cargo build … --target wasm32-unknown-unknown` exits 0 (or the install-failure is
surfaced); the CI `wasm-build` job is added + YAML valid; existing jobs unchanged;
`git diff --stat bde5cfe..HEAD` shows only `.github/workflows/**` (+ optional doc)
+ the 2 committed docs.
