# ADR-0098: Filter — `yee filter synth` emits physical dimensions + layout SVG

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0088 (`yee filter synth --plot`), ADR-0097 (F1.2.0 dimensional
synthesis), ADR-0086 (`yee-layout` geometry + SVG), `FILTER-DESIGN-ROADMAP.md`

---

## Context

F1.2.0 (ADR-0097) added `yee_filter::dimension_edge_coupled` — a `CouplingMatrix`
→ physical edge-coupled microstrip dimensions mapping — but nothing user-facing
calls it. `yee filter synth` today stops at the abstract coupling matrix + the
spec-mask verdict (+ optional `|S21|` plot). The natural next step is to surface
the **physical dimensions** (line width, resonator length, inter-resonator gaps)
and the **layout SVG** the studio/F1.0 path already knows how to render.

## Decision

Extend `yee filter synth` (in `yee-cli`) to compute and print the F1.2.0
dimensions and optionally write the layout SVG:

- New optional flags: `--eps-r <f64>` and `--h-mm <f64>` (substrate permittivity
  and dielectric height) with **FR-4 defaults** (`εr = 4.4`, `h = 1.6 mm`), and
  `--layout-svg <path>` to write the `yee_layout::Layout` SVG.
- After synthesis, call `dimension_edge_coupled(&project, &Substrate{eps_r, h})`
  and print a dimensions block: line width, resonator length, and each gap with
  its target `k` (SI + mm). If `dimension_edge_coupled` returns `Err` (e.g.
  `GapNotBracketed` for an unrealizable coupling on the chosen substrate), print
  a clear diagnostic and exit non-zero — do NOT silently skip.
- With `--layout-svg`, write `dimension_edge_coupled_layout(...).to_svg()`.

`yee-cli` gains a `yee-layout` dependency (already WASM-irrelevant — the CLI is a
native binary).

## Consequences

**Ships:** dimensions printout + optional SVG from `yee filter synth`. Closes the
"synthesis → user-visible geometry" loop on the CLI surface.

**Gate (`yee-cli` tests):** a test that for the committed Chebyshev 0.5 dB N=5
fixture (FR-4 default substrate) the dimension block is produced with sane values
(width/length/gaps > 0, gaps in mm range) and `--layout-svg` writes a non-empty,
well-formed SVG (`<svg` … `</svg>`). No FDTD; sub-second.

**Not in scope:** Gerber/KiCad export (F1.4); the `qe`→feed gap (F1.2.1); any
graphical preview. Pure CLI wiring over the shipped `dimension_*` API.

---

## References
- ADR-0097 (`dimension_edge_coupled` / `dimension_edge_coupled_layout`).
- `docs/superpowers/specs/2026-05-30-filter-f1-2-cli-dims-output-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-f1-2-cli-dims-output.md`.
