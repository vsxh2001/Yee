# Filter — `yee filter synth --gerber` export wiring — Design Spec

**ADR:** ADR-0102 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Write a single-copper-layer Gerber from `yee filter synth` via the shipped
`yee_export::layout_to_gerber`. CLI wiring; mirrors the existing `--layout-svg`
path. No FDTD, no new physics.

## Changes (`crates/yee-cli/**` ONLY)
- `Cargo.toml`: add `yee-export = { workspace = true }`.
- The `filter synth` clap subcommand (`src/main.rs`): add `--gerber <PATH>`
  (optional); thread into `run_synth`.
- `src/filter.rs::run_synth`: it already builds `substrate` + (for `--layout-svg`)
  the `Layout` via `dimension_edge_coupled_layout`. When `--gerber` is given,
  write `yee_export::layout_to_gerber(&layout, &yee_export::GerberOptions::default())`
  to the path. Reuse the same `layout` value the `--layout-svg` branch computes
  (compute once; don't duplicate the call if both flags are set). The dims `Err`
  path already returns non-zero before reaching here — unchanged.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-cli` exit 0 (no FDTD; fast).
4. **`cli_gerber` test (`crates/yee-cli/tests/`):** invoke the synth path with
   `--gerber <tmpfile>` for the committed Chebyshev 0.5 dB N=5 fixture (FR-4
   default substrate); assert exit success and that the written file's contents
   contain `%FSLAX46Y46*%` and `M02*`. Use the crate's existing CLI test idiom
   (the `cli_dims` / `--layout-svg` test is the precedent).

## Out of scope
Board outline / drill / multi-layer / KiCad / STEP (F1.4.1+); studio export
button; any change to `layout_to_gerber` (consume it as-is).
