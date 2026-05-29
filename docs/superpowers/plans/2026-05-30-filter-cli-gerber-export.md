# Filter — `yee filter synth --gerber` export wiring — Plan

**Spec:** `2026-05-30-filter-cli-gerber-export-design.md` · **ADR:** ADR-0102

## Lane
`crates/yee-cli/**` ONLY (`Cargo.toml`, `src/main.rs`, `src/filter.rs`, `tests/`).
Do NOT edit `yee-export`/`yee-filter`/`yee-layout`/other crates — consume their
public API. Out of lane → finding.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-cli-gerber`.

## Pattern files
- `crates/yee-cli/src/filter.rs` — `run_synth`; specifically the `--layout-svg`
  branch (where `dimension_edge_coupled_layout` is called + a file written) — add
  the `--gerber` write the same way, reusing the computed `layout`.
- `crates/yee-cli/src/main.rs` — the `filter synth` subcommand flags
  (`--layout-svg` declaration) to mirror for `--gerber`.
- `crates/yee-export/src/lib.rs` — `layout_to_gerber` + `GerberOptions::default`.
- `crates/yee-cli/tests/cli_filter.rs` — the `cli_dims` test (which already
  exercises `--layout-svg`) is the precedent for `cli_gerber`.

## Steps
1. `Cargo.toml`: add `yee-export = { workspace = true }`.
2. main.rs: add `--gerber <PATH>` to `filter synth`; pass through.
3. filter.rs: compute the `Layout` once (it may already be computed for
   `--layout-svg`); when `--gerber` is set, write `layout_to_gerber(&layout,
   &GerberOptions::default())` to the path.
4. test: `cli_gerber` per DoD 4.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-cli --jobs 2
```
Do NOT run `cargo test --workspace`, FDTD, mom-001.

## Escape hatch
Blocked > 15 min — the clap wiring or the shared-`layout` refactor in `run_synth`
fights the existing structure → STOP and surface the exact blocker. Do NOT
duplicate the dimensioning call in a way that diverges the SVG vs Gerber layout.
Do NOT edit yee-export/yee-filter.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` = only `crates/yee-cli/**`
(+ `Cargo.lock`) + the 3 committed docs.
