# F2.2-cli — lumped-board export in `yee filter synth` — implementation plan

**Spec:** [2026-06-04-f2-2-cli-lumped-board-export-design.md](../specs/2026-06-04-f2-2-cli-lumped-board-export-design.md)
**ADR:** [ADR-0158](../../src/decisions/0158-f2-2-cli-lumped-board-export.md)
**Fork from:** `main` (4e444e7 or later).

## Steps

1. Read `crates/yee-cli/src/filter.rs` (`run_synth`, the export block ~line 189, how `Substrate` is
   built from `--eps-r`/`--h-mm`) and `main.rs` (how `filter synth` flags are declared/threaded).
   Read `crates/yee-filter/src/board.rs` (`lumped_board`, `Footprint`, `LumpedBoard.layout`) and
   `crates/yee-filter/src/lumped.rs` (`synthesize_lumped`, `LumpedError`).
2. `main.rs`: add `--lumped` (flag) + `--footprint <0402|0603|0805>` (default `0603`, parse to
   `yee_filter::Footprint`) to the `filter synth` subcommand; pass into `run_synth`.
3. `filter.rs::run_synth`: new params `lumped: bool, footprint: Footprint`. In the export block:
   if `lumped`, `let ladder = synthesize_lumped(&proj)?;` (map `LumpedError` into the CLI error
   context), `let board = lumped_board(&ladder, &substrate, footprint);`, set `layout = board.layout`
   — then the SAME existing `--gerber`/`--kicad-pcb`/`--layout-svg` writers run on it. The planar
   branch is unchanged (the `if lumped { … } else { <planar dims→layout> }` selects the source
   `Layout`; the writers below are shared). Keep the single-`Layout` invariant.
4. Gate `crates/yee-cli/tests/cli_lumped_export.rs` (mirror the planar `cli_gerber` gate): synth a
   small bandpass spec with `--lumped --gerber <tmp>`; assert non-empty + valid RS-274X markers +
   `G36` region count ≥ 2·N + differs from the planar `--gerber` output for the same spec. Fast,
   non-`#[ignore]`.
5. Verify (boxed, but light — no FEM): `cargo fmt --check`; `cargo clippy -p yee-cli --all-targets
   -- -D warnings`; `cargo test -p yee-cli` (incl. the new gate). All exit 0.

## Dispatch

- ONE agent, worktree off main, lane `crates/yee-cli/src/{filter.rs,main.rs}` +
  `crates/yee-cli/tests/cli_lumped_export.rs`. The increment is pure-compute (no FEM) → the agent
  can run its own gate (no code/gate-run split needed). Reviewer (never self-review) after; I verify
  boxed + merge `--no-ff`.

## CI

- The new gate runs in the default `cargo test --workspace` (fast, non-ignored). Confirm yee-cli's
  test deps are present; no new heavy/release gate job needed.
