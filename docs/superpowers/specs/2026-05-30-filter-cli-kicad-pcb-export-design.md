# Filter — `yee filter synth --kicad-pcb` export wiring — Design Spec

**ADR:** ADR-0106 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Make the shipped `yee_export::layout_to_kicad_pcb` (F1.4.1b, ADR-0105) reachable
from the CLI, exactly as `--gerber` (ADR-0102) exposes `layout_to_gerber`. After
this, `yee filter synth … --kicad-pcb out.kicad_pcb` writes a KiCad 7 board file
the user opens directly in the KiCad PCB editor — the goal's "kicad export"
endpoint, reachable end-to-end (spec → synth → dims → layout → `.kicad_pcb`).

## Changes (`crates/yee-cli/**` ONLY)
- `src/filter.rs`:
  - `run_synth` gains a `kicad_pcb: Option<&Path>` parameter (beside the existing
    `gerber: Option<&Path>`).
  - In the existing optional-layout-export block (the one gated on
    `layout_svg.is_some() || gerber.is_some()`), extend the guard to also fire on
    `kicad_pcb.is_some()`, and add a branch mirroring the `--gerber` one:
    `layout_to_kicad_pcb(&layout, &KicadPcbOptions::default())` → `std::fs::write`
    → a `println!("  wrote KiCad PCB: …")`. The SAME single `layout` is reused
    (so SVG / Gerber / KiCad can never diverge).
  - Update the function's doc-usage line to include `[--kicad-pcb <out.kicad_pcb>]`.
- `src/main.rs` (or wherever the `filter synth` clap args live): add a
  `--kicad-pcb <PATH>` arg (optional, `PathBuf`), mirroring `--gerber`, and thread
  it into the `run_synth(..)` call.

## DoD (machine-checkable)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-cli --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-cli` exit 0 — including a new `cli_kicad_pcb` test mirroring
   the existing `cli_gerber` test: run `synth` with `--kicad-pcb <tmp>` (a small
   spec), assert the file exists and its contents start with `(kicad_pcb` and
   contain `Edge.Cuts`.
4. `yee filter synth --help` lists `--kicad-pcb` (verifiable via the clap arg
   being present; the test in DoD 3 covers the functional path).

## Out of scope
Studio export buttons (yee-studio lane — separate); drill/multi-layer/footprints
(F1.4.1c); any change to `layout_to_kicad_pcb` itself (it shipped in F1.4.1b).
Keep the change a thin mirror of the `--gerber` path — one `Layout`, three
optional writers.
