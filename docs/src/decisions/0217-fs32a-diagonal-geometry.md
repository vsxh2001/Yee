# ADR-0217: FS.3.2a — diagonal-edge geometry under full-wave test

**Date:** 2026-07-12 · **Status:** accepted · **Track:** FS.3 (geometry generality)
**Spec:** `docs/superpowers/specs/2026-07-12-fs32a-diagonal-geometry-design.md`

## Context

`yee_voxel::point_in_polygon` is a general even-odd raycast, so arbitrary
polygons have always rasterized — but every EM gate in the repo drove
axis-aligned rectangles, leaving the geometry-generality claim behind
Gerber/DXF import untested where it matters: in a measurement. The
classic physics-referenced testcase for a non-axis-aligned edge is the
mitered 90° microstrip bend (Douville & James): chopping the outer corner
removes the bend's excess capacitance and its reflection.

## Decision

1. **`yee_layout::double_jog(substrate, w, run_x, gap_x, jog_dy, MiterStyle)`**
   — a four-bend through line with both ports x-facing at equal y, so the
   uniform *and* graded two-port fixtures apply unchanged, and
   `reference_through_line` degenerates correctly (equal port y).
   `MiterStyle::Square` vs `MiterStyle::Mitered { f }` (45° cut, legs
   `f·w`; default 0.7 ≈ the published optimum for w/h ≈ 1.9).
   Construction: straights overlap corners by `0.2·w` (point-sampling
   seam robustness — the quasi-Yagi lesson) with the overlap capped clear
   of the cut (`f ≤ 0.8` enforced); each corner is its **own polygon** so
   the automesh rulebook's per-polygon AABBs refine every bend (a single
   outline polygon would present one blob AABB).
2. **Gate `voxel-poly-001`** (instant, GREEN): a 45°-cut square
   rasterizes to the exact predicted staircase (cut line chosen off cell
   centres — a boundary-exact centre is even-odd-implementation-defined,
   measured on the first run); the mitered jog masks strictly fewer PEC
   edges than the square jog on the identical grid, by the right area
   (~196 mask edges ≈ 4 cut triangles × 2 mask arrays).
3. **Gate `engine-miter-001`** (release, dedicated step in the graded CI
   job; `graded_` test-name prefix keeps it out of the blanket step):
   square vs mitered double-jog through the graded fixture, double-ratio
   |S21| on the ADR-0216 criterion band (3.5–5.75 GHz). Asserts:
   (a) mitered band-mean |S21| ≥ square; (b) mitered worst in-band |S21|
   ≥ −6 dB; (c) the miter advantage does not shrink with frequency
   (excess-C reflection ∝ f), 0.02 linear slack.

## Measured

(pinned from the first green run — see the gate's `--nocapture` output)

## Non-goals

Optimal-miter sweep, y-directed ports (fixture lane), Gerber arcs/flashes
(FS.3.2 proper), DXF, MoM cross-check.
