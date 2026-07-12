# FS.3.2a — diagonal-edge geometry under full-wave test (mitered bends)

**Date:** 2026-07-12 · **Track:** FS.3 (geometry generality)
**Plan:** `docs/superpowers/plans/2026-07-12-fs32a-diagonal-geometry.md`

## Problem

`yee_voxel::point_in_polygon` is a general even-odd raycast, so arbitrary
polygons already rasterize — but every EM gate in the repo drives
axis-aligned rectangles. The geometry-generality claim behind Gerber
import (FS.3) and future DXF/GDSII needs a full-wave gate on a
**non-axis-aligned edge**. The classic, physics-referenced testcase is
the mitered 90° bend: chopping the outer corner (a 45° cut — Douville &
James) reduces the bend's excess capacitance and its reflection.

## Design

1. **`yee_layout::double_jog`** generator: port 1 → x-run → 90° up-jog →
   y-run (Δy) → 90° down-jog → x-run → port 2, both ports at the same y
   (fixture-compatible: both feeds x-directed). `MiterStyle::Square`
   builds each corner from overlapping rects; `MiterStyle::Mitered { f }`
   replaces each corner square with a 5-vertex polygon whose outer corner
   is cut by a 45° edge (cut legs `f·w`; default f = 0.7, near the
   Douville–James optimum for w/h ≈ 1.9 on FR-4).
2. **Gate `voxel-poly-001`** (instant): a 45°-cut polygon rasterizes to
   the expected staircase on a hand-computed grid; the mitered corner
   masks strictly fewer cells than the square corner; single-rect layouts
   stay bit-identical.
3. **Gate `engine-miter-001`** (release): the double-jog board (FR-4,
   w = 3 mm, four corners) through the certified graded fixture
   (`two_port_board_jobs_graded`) — three runs on per-DUT grids
   (square DUT, mitered DUT, each with its thru reference): launch-
   normalized double-ratio |S21| over 3.5–5.75 GHz (the ADR-0216
   criterion band). Asserts: (a) mitered mean |S21| ≥ square mean |S21|
   (the miter physics, in linear magnitude); (b) mitered worst-case
   in-band |S21| above a floor pinned from the first green run; (c) the
   square-vs-mitered gap grows with frequency (reflection ∝ f for the
   corner's excess C — assert gap at the top decade of the band ≥ gap at
   the bottom).

## Non-goals

Optimal-miter sweep (f is fixed), y-directed ports (fixture lane),
arcs/flashes in Gerber (FS.3.2 proper), MoM cross-check.
