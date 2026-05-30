# ADR-0106: Filter — `yee filter synth --kicad-pcb` export wiring

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0105 (F1.4.1b `layout_to_kicad_pcb`), ADR-0102 (`--gerber`
wiring — the pattern this mirrors), ADR-0100/0103 (Gerber emitters),
`FILTER-DESIGN-ROADMAP.md`

---

## Context

F1.4.1b (ADR-0105) shipped `yee_export::layout_to_kicad_pcb`, but no user can
invoke it: only the library API exists. The Gerber emitter got its CLI handle in
ADR-0102 (`yee filter synth --gerber <out.gbr>`). The KiCad `.kicad_pcb` emitter
should get the symmetric handle so the "kicad export" goal endpoint is reachable
end-to-end from a single command.

## Decision

Add `--kicad-pcb <PATH>` to `yee filter synth`, mirroring `--gerber` exactly:
`run_synth` gains a `kicad_pcb: Option<&Path>` parameter; the existing
optional-layout-export block (already shared by `--layout-svg` and `--gerber`,
all computing one `Layout` so they cannot diverge) gains a branch that calls
`layout_to_kicad_pcb(&layout, &KicadPcbOptions::default())` and `std::fs::write`s
the result. The change is confined to `crates/yee-cli/**`; `layout_to_kicad_pcb`
itself is unchanged.

## Consequences

**Ships:** `yee filter synth … --kicad-pcb out.kicad_pcb` writes a KiCad 7 board
file. Combined with `--gerber` and `--layout-svg`, one synth invocation can emit
the SVG preview, the fab Gerber, and the editable KiCad board — all from the same
`Layout`.

**Gate:** `cargo test -p yee-cli` passes, including a new `cli_kicad_pcb` test
(mirror of `cli_gerber`) that runs `synth --kicad-pcb <tmp>` and asserts the file
starts with `(kicad_pcb` and contains `Edge.Cuts`.

**Not in scope:** the studio export button (yee-studio lane — a follow-on);
drill/multi-layer/footprints (F1.4.1c); any change to the emitter.

---

## References
- ADR-0102 (the `--gerber` wiring this mirrors); ADR-0105 (the emitter).
- `docs/superpowers/specs/2026-05-30-filter-cli-kicad-pcb-export-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-cli-kicad-pcb-export.md`.
