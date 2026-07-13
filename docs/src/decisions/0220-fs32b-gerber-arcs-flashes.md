# ADR-0220: FS.3.2b — Gerber import: circular arcs + aperture flashes

**Date:** 2026-07-13 · **Status:** accepted · **Track:** FS.3
**Spec:** `docs/superpowers/specs/2026-07-13-fs32-gerber-arcs-flashes-design.md`

## Decision

Extend `yee_export::import` (FS.3.0, ADR-0209) — the IMPORT side only —
with the two constructs every real CAD export uses that the FS.3.0
subset rejected by name:

1. **Circular-arc region segments.** `G75*` multi-quadrant mode is
   modal state; `G02`/`G03` set CW/CCW circular interpolation
   (standalone words or inline `G02X…I…J…D01` prefixes; `G01` returns
   to linear). An in-region circular `D01` reads modal `X`/`Y` plus
   non-modal `I`/`J` centre offsets (0 when omitted, per Ucamco G75)
   and is tessellated into contour vertices. `start == end` is a full
   360° circle (Ucamco §4.7.2). Arc **endpoints are exact** (the 4.6
   fixed-point word, `n·1e-9` m); only interior vertices are
   synthesized, at uniform angular steps with the radius linearly
   interpolated between the (quantized) start/end radii.
2. **Flashed apertures (`D03`)** outside regions. `%AD…%` definitions
   are parsed per D-code (`C,<dia>` / `R,<x>X<y>`, no holes); `D<code>*`
   (code ≥ 10) selects. Rect flashes convert **exactly** (4 corners,
   CCW from lower-left); circle flashes tessellate (vertex 0 on the +x
   axis, CCW). Flash polygons interleave with regions in file order.

**Chord tolerance pinned:** `pub const ARC_CHORD_TOL_M = 1.0e-6` (1 µm).
Max angular step `φ_max = 2·acos(1 − tol/r)` (chord sagitta
`r(1 − cos(φ/2)) ≤ tol`), `n = ceil(sweep/φ_max)` segments (min 1;
circles min 4).

**New named rejections** (the FS.3.0 philosophy — never mis-parse):
`G74` single-quadrant mode and arcs before `G75` →
`UnsupportedCommand`; `D03` inside a region → `FlashInRegion` (spec
forbids it); flash of an unselected/undefined D-code →
`UnknownAperture`; flash of obround/polygon/macro/holed/zero-size
apertures → `UnsupportedAperture` (bookkept as `Other`, rejected only
if flashed); zero-radius or start/end radii differing > 10·tol →
`BadArc`. Inches, polarity, and stroked draws (now linear *or*
circular) stay rejected; `gerber_to_outline` is unchanged (arcs there
queued with DXF).

**Writer untouched.** `layout_to_gerber` emits no `G02`/`G03`/`D03`
(verified), so there is nothing to round-trip and `gerber-rt-001`
byte-stability is unaffected. Emitting arcs would lose the polygon
representation — a different feature, out of FS.3.2b scope.

## Verification

Gate `gerber-rt-003` (`crates/yee-export/tests/gerber_arcs_flashes.rs`,
instant, non-ignored, GREEN): pinned numbers —

- quarter arc r = 1 mm: **n = 18** segments
  (`φ_max(1 mm) = 0.0894502 rad`), endpoints exact ≤ 0.5 nm, every
  vertex on-circle ≤ 1 nm, measured sagitta ≤ 1 µm, CW/CCW orientation
  both checked;
- full-circle region r = 1 mm (start == end): **n = 71** vertices;
- circle flash ⌀1 mm: **n = 50** vertices, vertex 0 at (r, 0);
- rect flash `R,2X1` at (5, 5) mm: 4 exact corners;
- region/flash file-order interleaving; the full rejection matrix
  above plus the FS.3.0 rejections re-asserted.

`gerber-rt-001`'s stale "arcs are rejected" assertion was updated to
pin `G74` instead (arcs are now in-subset).

## Lesson (cheap but real)

Interpolation mode is **modal**: a linear closing edge after a `G03`
arc must be preceded by `G01`, or its `D01` is an arc with a
defaulted (0, 0) centre offset — which this importer correctly flags
as `BadArc` (zero radius) rather than silently drawing a line. The
first cut of the gate hit exactly this; real CAD files always emit the
`G01`, and the failure mode is a named error, not corruption.
