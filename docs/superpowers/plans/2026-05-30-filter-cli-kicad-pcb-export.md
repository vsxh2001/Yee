# Filter — `yee filter synth --kicad-pcb` export wiring — Plan

**Spec:** `2026-05-30-filter-cli-kicad-pcb-export-design.md` · **ADR:** ADR-0106

## Lane
`crates/yee-cli/**` ONLY (`src/filter.rs`, `src/main.rs`, `tests/`). Do NOT edit
`yee-export`/`yee-layout`/other crates — consume the public API
(`yee_export::{layout_to_kicad_pcb, KicadPcbOptions}`). Out of lane → finding.

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-cli-kicad-pcb`.

## Pattern files (MIRROR exactly)
- `crates/yee-cli/src/filter.rs` — the `--gerber` path: the `gerber:
  Option<&Path>` param on `run_synth`, the `layout_svg.is_some() ||
  gerber.is_some()` export block, the `layout_to_gerber` + `std::fs::write` +
  `println!` branch (~lines 183–201). Add the `--kicad-pcb` branch identically.
- `crates/yee-cli/src/main.rs` — the clap `--gerber` arg definition + its
  thread-through into `run_synth(..)`. Add `--kicad-pcb` the same way.
- The existing `cli_gerber` integration test (find it: `grep -rn cli_gerber
  crates/yee-cli/tests`) — clone it to `cli_kicad_pcb`, asserting the written
  file starts with `(kicad_pcb` and contains `Edge.Cuts`.

## Steps
1. `src/filter.rs`: add `kicad_pcb: Option<&Path>` to `run_synth`; extend the
   export-guard; add the KiCad write branch; update the doc-usage line.
2. `src/main.rs`: add the `--kicad-pcb <PATH>` clap arg; pass it into `run_synth`.
3. Add the `cli_kicad_pcb` test (mirror `cli_gerber`).

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-cli --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-cli --jobs 2
```
This box is MEMORY-CONSTRAINED. `yee-cli` is already built from prior work, so
this is an INCREMENTAL rebuild (one crate + relink) — keep `--jobs 2`. Do NOT run
`cargo test --workspace`, FDTD, mom-001, or build the whole workspace. If
`cargo test -p yee-cli` itself OOMs (exit 137 / SIGKILL) — note that yee-cli's
test target may pull heavy bins — fall back to `cargo build -p yee-cli` +
`cargo test -p yee-cli --test <the_cli_synth_test_file>` (the single test file),
and SURFACE the OOM in your report rather than running the whole suite.

## Escape hatch
Blocked > 15 min (the `synth` clap args live somewhere unexpected, the `cli_gerber`
test idiom doesn't clone cleanly, or the build OOMs) → STOP and surface the exact
blocker + what you tried. Do NOT edit yee-export/yee-layout; do NOT change
`layout_to_kicad_pcb`; do NOT widen scope to studio buttons.

## Done when
DoD 1–4 pass; `git diff --stat <base>..HEAD` = only `crates/yee-cli/**` (+ the 3
committed docs); `--kicad-pcb` writes a `(kicad_pcb …)` file sharing the single
`Layout` with `--gerber`/`--layout-svg`.
