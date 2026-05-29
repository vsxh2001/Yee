# Filter — `yee filter synth` dims + layout-SVG output — Plan

**Spec:** `2026-05-30-filter-f1-2-cli-dims-output-design.md` · **ADR:** ADR-0098

## Lane
`crates/yee-cli/**` ONLY (`Cargo.toml`, `src/filter.rs`, the clap subcommand in
`src/main.rs`, `tests/`). Do NOT edit `yee-filter`/`yee-layout`/any other crate —
consume their public API. Out of lane → finding.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-f1-2-cli-dims`.

## Pattern files
- `crates/yee-cli/src/filter.rs` — `run_synth` (extend it; mirror its style for
  the `--plot` optional-output handling).
- `crates/yee-cli/src/main.rs` — the `filter synth` clap subcommand (how
  `--output`/`--plot` flags are declared + passed) — add `--eps-r`/`--h-mm`/
  `--layout-svg` the same way.
- `crates/yee-filter/src/dimension.rs` — `dimension_edge_coupled` /
  `dimension_edge_coupled_layout` / `EdgeCoupledDimensions` / `DimError` (the API
  to call; read the field names + error variants).
- `crates/yee-layout/src/lib.rs` — `Substrate` field names/units + `Layout::to_svg`.
- existing `crates/yee-cli/tests/` — the CLI test harness/pattern to imitate.

## Steps
1. `Cargo.toml`: add `yee-layout = { workspace = true }`.
2. main.rs: add the three flags to `filter synth`; pass through.
3. filter.rs: extend `run_synth`; build `Substrate`; call `dimension_edge_coupled`;
   print the dims block (SI + mm); handle `Err` with a clear message + non-zero
   exit; write the SVG when `--layout-svg` is given.
4. tests: `cli_dims` per DoD 4.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-cli --jobs 2
```
Do NOT run `cargo test --workspace`, FDTD, or mom-001.

## Escape hatch
Blocked > 15 min — the clap subcommand wiring fights the existing arg structure,
or the committed fixture's couplings can't be realized on FR-4 (so the dims path
errors) → STOP and surface: the exact blocker + (if a fixture issue) the
computed target_k vs the achievable range. Do NOT weaken the gate or fabricate
dims. Do NOT edit `yee-filter`/`yee-layout`.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` shows only `crates/yee-cli/**`
(+ `Cargo.lock`) + the 3 committed docs.
