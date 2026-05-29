# Filter — `yee filter synth` physical-dimensions + layout-SVG output — Design Spec

**ADR:** ADR-0098 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Surface F1.2.0's physical dimensions (and the layout SVG) from `yee filter
synth`. CLI wiring over the shipped `yee_filter::dimension_edge_coupled` /
`dimension_edge_coupled_layout`. No FDTD, no new physics.

## Changes (`crates/yee-cli/**` ONLY)
- `Cargo.toml`: add `yee-layout = { workspace = true }`.
- The `filter synth` clap subcommand (in `src/main.rs` or wherever the subcommand
  is defined): add `--eps-r <f64>` (default `4.4`), `--h-mm <f64>` (default
  `1.6`), `--layout-svg <PATH>` (optional). Thread them into `run_synth`.
- `src/filter.rs::run_synth`: extend the signature to accept the substrate
  params + optional layout-svg path. After the existing synthesis:
  1. Build `yee_layout::Substrate` from `eps_r` + `h = h_mm·1e-3` (read the
     `Substrate` field names/units from `yee-layout`).
  2. `match dimension_edge_coupled(&project, &substrate)`:
     - `Ok(dims)` → print a block: `line_width`, `resonator_length`, and each
       `gaps_m[i]` with `target_k[i]` (show SI metres + mm).
     - `Err(e)` → print a clear diagnostic (e.g. "cannot realize coupling k=… on
       εr=…, h=… (gap out of [..]) ") and return a non-zero `ExitCode`.
  3. If `--layout-svg` given, write
     `dimension_edge_coupled_layout(&project, &substrate)?.to_svg()` to the path.
- Keep the existing Touchstone + `--plot` behaviour unchanged.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-cli` exit 0 (no FDTD; fast).
4. **`cli_dims` test (`crates/yee-cli/tests/`):** invoke the dims path for the
   committed Chebyshev 0.5 dB N=5 fixture spec (FR-4 default substrate) — assert
   the produced `EdgeCoupledDimensions` (or the formatted output) has
   `line_width_m > 0`, `resonator_length_m > 0`, every `gaps_m[i] > 0`, and that
   writing the layout SVG yields a string containing `<svg` and `</svg>`. Use the
   existing yee-cli test harness/pattern (assert_cmd, or a direct call into a
   `run_synth`-style entry — match what the crate already does).

## Out of scope
Gerber/KiCad (F1.4); `qe`→feed gap (F1.2.1); graphical preview; any change to the
`dimension_*` algorithm (consume it as-is).
