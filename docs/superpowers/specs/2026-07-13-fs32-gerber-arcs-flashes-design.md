# FS.3.2 — Gerber import: circular arcs + aperture flashes

**Date:** 2026-07-13
**Track:** FULL-SUITE-ROADMAP FS.3.2 (import side). FS.3.0 (ADR-0209)
parses exactly the region dialect our writer emits and rejects, by name,
everything else — including arcs (`G02`/`G03`) and flashes (`D03`). Real
CAD exports use both constantly: KiCad rounds corners with region arcs
and emits every pad as a flash. This spec extends the IMPORT side only;
the writer stays region-only (extending it is out of scope, and there is
nothing to round-trip — see Non-goals).

## Decomposition

- **FS.3.2b (this spec): arcs + C/R flashes in `yee_export::import`.**

  1. **Circular-arc region segments.** `G75*` (multi-quadrant mode)
     accepted as modal state; `G02*`/`G03*` set clockwise/counter-
     clockwise circular interpolation (standalone or as a `G02X…`/`G03X…`
     prefix on the coordinate word); `G01*` returns to linear. Inside a
     region, a `D01` with active circular interpolation reads modal
     `X`/`Y` plus **non-modal** `I`/`J` centre offsets (each defaulting
     to 0 when omitted, per the Ucamco G75 rule) and is **tessellated**
     into polygon vertices at a documented chord tolerance (below).
     `start == end` in multi-quadrant mode is a full 360° circle (Ucamco
     §4.7.2). Arc endpoints are placed **exactly** (the 4.6 fixed-point
     word, `n·1e-9` m); only the interior vertices are synthesized.
  2. **Flashed apertures (D03)** outside regions, for the two standard
     templates our target boards need:
     - `C,<dia>` — circle, tessellated at the same chord tolerance,
       vertex 0 on the +x axis, CCW;
     - `R,<x>X<y>` — rectangle, converted **exactly** to its 4 corners
       CCW from lower-left.
     `%ADD…%` definitions are now parsed and bookkept per D-code;
     `D<code>*` (code ≥ 10) selects the current aperture. Flash polygons
     land in the output `Vec<Polygon>` in file order, interleaved with
     regions.

- **Chord tolerance (pinned):** `ARC_CHORD_TOL_M = 1.0e-6` (1 µm — two
  orders below any λ/20 cell this suite meshes; segment counts stay
  double-digit for mm-scale features). The maximum angular step for
  radius `r` is `φ_max = 2·acos(1 − tol/r)` (sagitta of a chord
  subtending `φ` is `r(1 − cos(φ/2))`); an arc sweeping `θ` uses
  `n = ceil(θ/φ_max)` uniform segments (min 1; circles min 4). This is a
  `pub const` so gates and downstream code pin the same number.

- **Named rejections (new).** Everything below fails with an explicit
  error, never a silent mis-parse:
  - `G74*` (single-quadrant arc mode — legacy, ambiguous): rejected as
    `UnsupportedCommand`;
  - an arc `D01` before `G75*`: `UnsupportedCommand` (single-quadrant
    semantics would apply);
  - `D03` inside a `G36` region (forbidden by the Gerber spec):
    `FlashInRegion`;
  - `D03` with no aperture selected / an undefined D-code:
    `UnknownAperture`;
  - `D03` of a non-C/R template (obround, polygon, macro), a holed
    aperture (`C,1X0.5`), or a degenerate (zero-size) one:
    `UnsupportedAperture`;
  - geometrically broken arcs (zero radius, start/end radii differing
    by > 10·tol): `BadArc`.

- **Existing rejections stay named:** inches (`%MOIN%` →
  `ImperialUnits`), polarity (`%LP…%` → `UnsupportedCommand`), stroked
  draws outside a region in the copper importer (`UnsupportedCommand` —
  linear *or* circular). `gerber_to_outline` keeps rejecting arcs,
  flashes and regions unchanged.

## Gate `gerber-rt-003` (unit, instant, non-ignored)

Hand-written Gerber snippets with analytically predictable geometry:

1. **Quarter arc** (r = 1 mm, 90° CCW `G03`, I/J offsets): endpoint
   vertices exact (≤ 0.5 nm, the existing gate idiom); segment count
   equals the pinned `n = 18`; every synthesized vertex on the circle to
   ≤ 1 nm; measured max sagitta ≤ `ARC_CHORD_TOL_M`.
2. **Full-circle region** (start == end, r = 1 mm): closes, `n = 71`
   segments, sagitta bound holds. Clockwise `G02` quarter arc mirrors
   case 1 (orientation check).
3. **Rect flash** `R,2X1` at (5, 5) mm: exactly 4 vertices equal to the
   half-extent corners (≤ 0.5 nm).
4. **Circle flash** `C,1` at origin: `n = 50` vertices, vertex 0 at
   (r, 0), all on the circle, sagitta bound; flash + region in one file
   import in file order.
5. **Rejection paths:** every bullet in "Named rejections" above, plus
   the FS.3.0 rejections re-asserted (inches, polarity, stroked draws)
   so the subset boundary cannot silently widen.

## Non-goals (FS.3.2b)

- **Writer arcs/flashes.** `layout_to_gerber` emits regions with linear
  segments only (verified — no `G02`/`G03`/`D03` anywhere in the
  emitter), so there is nothing to round-trip; `gerber-rt-001`
  byte-stability is untouched. Emitting arcs would *lose* the polygon
  representation and is a different feature.
- Obround/polygon/macro apertures, aperture holes, `%LP` polarity,
  step-repeat, inches, stroked-draw copper — still rejected by name.
- Arc support in the outline importer (`gerber_to_outline`) — the
  Edge.Cuts dialect our writer emits is rectangular; queued with DXF.
