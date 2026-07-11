# FS.3 — layout import: Gerber first

**Date:** 2026-07-08
**Track:** FULL-SUITE-ROADMAP FS.3. Import is the "bring your own board"
door every commercial tool has; our own writer (`yee_export`, ADR-0198
byte-checked in the studio) defines the first dialect to support.

## Decomposition

- **FS.3.0 (walking skeleton, this spec): RS-274X region-subset import.**
  `yee_export::import::gerber_to_polygons(&str) -> Result<Vec<Polygon>>`
  parsing exactly the dialect our writer emits — `%FSLAX46Y46*%` (any
  `AXmnYmn` accepted), `%MOMM*%` (imperial rejected with a clear error),
  `G04` comments, `%AD…%`/`D<code>*` aperture bookkeeping (regions ignore
  them), `G36*`/`G37*` region contours from modal `X…Y…D02/D01` words,
  `M02*`. Contours drop the explicit closing vertex (Polygon convention);
  coordinates are exact (`fixed46 → metres` is `n·1e-9`).
  `gerber_to_layout(gerber, substrate, ports)` wraps the polygons (Gerber
  carries no stackup/ports — the caller provides them; the studio's
  import flow will ask).
- **Gate `gerber-rt-001`** (unit, instant): for real generator layouts
  (hairpin BPF, inset patch, quasi-Yagi, 2×1 array):
  1. `import(export(L))` vertex counts match and every vertex equals the
     original to ≤ 1 nm (the 4.6 quantum);
  2. **byte-stability**: `export(import(export(L))) == export(L)`
     byte-identical — the house artifact philosophy, and the strongest
     cheap proof the parse is lossless over the writer dialect.
  Plus error-path units: imperial units, unclosed region, draw before
  move, junk words.
- **FS.3.1 (queued):** stroked-path import (the outline layer), arbitrary
  aperture flashes (D03), DXF; then the studio import → verify → export
  flow (FULL-SUITE gate: an imported reference board measures within
  tolerance of its native twin).

## Non-goals (FS.3.0)

Arc segments (G02/G03), polarity (%LP), step-repeat, macro apertures,
inches. All rejected with explicit errors, never silently mis-parsed.
