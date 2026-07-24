# FS.3.3 — DXF import (closes FS.3)

**Date:** 2026-07-24 · **Track:** FS.3 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-export/**` (+ docs)
**Predecessors:** FS.3.0/3.1/3.2 Gerber import chain (ADR-0209/0217/0220/0229) —
subset parser with named rejections, `gerber_to_outline`/`gerber_to_layout`,
byte-identical round-trip, bit-identical measurement twin. FS.3's last item.

## Approach — mirror the Gerber discipline

DXF (ASCII R12+ group-code format) import for trace outlines, with the same
subset-plus-named-rejections philosophy that made the Gerber importer reviewable:

- **Supported subset**: `LWPOLYLINE` (closed, straight segments — the trace-outline
  workhorse) and closed `POLYLINE`/`VERTEX` chains (R12 fallback); `$INSUNITS`
  header respected for mm/inch (reject unitless files or document the assumed
  default with a named lenient mode — pick ONE, document it); entities on any
  layer, with an optional layer filter argument.
- **Bulge (arc) segments**: tessellate at the same pinned 1 µm chord tolerance as
  FS.3.2b Gerber arcs (reuse the tessellation helper if it is reusable; else mirror
  its contract + gate idiom).
- **Named rejections** (typed errors, each tested): open polylines, CIRCLE/ARC/
  ELLIPSE/SPLINE entities, TEXT/DIMENSION, blocks/INSERT, 3-D (nonzero Z /
  elevation), unsupported `$INSUNITS`.
- **Entry points**: `dxf_to_outline(bytes, opts) -> Outline...` mirroring
  `gerber_to_outline`'s output type exactly, so `gerber_to_layout`'s downstream
  (outline→Layout) is shared, not duplicated — the twin-gate chain then covers DXF
  for free at the geometry level.

## Deliverables

1. `yee_export::import_dxf` (or module `dxf`) per above, with unit tests: exact
   vertex parsing (hand-authored minimal DXF fixtures in the test file or
   `tests/data/` — follow where gerber-rt keeps fixtures), bulge tessellation
   sagitta ≤ 1 µm (CW+CCW), units handling, and the full rejection matrix.
2. **Gate `dxf-rt-001`**: the S.6 stub-board trace geometry hand-authored as a
   DXF fixture → `dxf_to_outline` → structural identity vs the native generator's
   polygons at the gerber-rt tolerance (0.5 nm for straight segments; tessellated
   bulges at their own pinned tolerance). Geometry-only — no FDTD (the FS.3.2c
   twin gate already proved outline→measurement; do not re-run a 5-min solve).
3. **ADR-0230** + FS.3 roadmap row → **FS.3 COMPLETE**.

## Constraints

- No new dependency for parsing (group-code pairs are `lines().tuples()`-simple;
  a hand parser matches the repo's Gerber precedent) — if the implementer believes
  a crate is genuinely warranted, that is a STOP-and-surface, not a unilateral add
  (TECH_STACK.md discipline).
- Existing gerber gates + studio import gate unmodified/green.
- Rejections are typed and each carries a test — the subset boundary stays explicit.

## Non-goals

DXF *export*; INSERT/block expansion; 3-D; studio panel wiring (panel accepts
Gerber today; DXF wiring is a small follow-on once the library seam exists);
DXF-sourced FDTD twin (covered transitively).
