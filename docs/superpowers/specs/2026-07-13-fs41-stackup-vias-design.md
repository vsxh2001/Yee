# FS.4.1 — vias through multilayer stackups (design)

**Date:** 2026-07-13 · **Track:** FS.4 (`FULL-SUITE-ROADMAP.md`) · **ADR:** 0221

## Problem

FS.4.0 (ADR-0215) gave the FDTD flow N-layer stackups with buried traces
and lids, but the only via primitive is the single-layer
`yee_voxel::with_via_at_cell(model, i, j, k_top)` — a ground-to-trace
`E_z` PEC column starting at `k = 0` (R.1, ADR-0194). A multilayer board
needs two more shapes:

- **through-via** — a column through the *entire* stack (ground → lid on
  a lidded board);
- **blind via** — a column between two arbitrary z node-planes (e.g.
  layer-1 top to layer-2 top), touching neither outer plate.

The engine protocol already carries the `pec_mask_ez` field end-to-end
(R.1), so — exactly like FS.4.0 — the gap is purely the voxel-side
helper, not the solver.

## Design

Keep the shipped cell-index post-processing idiom (`with_via_at_cell` /
`truncate_ground_at_cell`): callers map layout-frame coordinates to grid
cells through the same `x₀ = bbox.min − margin` origin the voxelizer
used, and the helper mutates `pec_mask_ez` in place.

1. **`with_via_between(model, i, j, k_lo, k_hi)`** — the generalization:
   `E_z` edges `k = k_lo .. k_hi` at grid column `(i, j)` become PEC,
   i.e. a metal post spanning node-plane `k_lo` to node-plane `k_hi`.
   `k_lo`/`k_hi` are **grid cell indices**, not stackup layer indices —
   the caller quantizes layer heights exactly the way `voxelize_stackup`
   did (`round(h/dx).max(1)` cumulative bands). Panics on a malformed
   range (caller bug, not data).
2. **`with_through_via_at_cell(model, i, j)`** — convenience for the
   full stack: `with_via_between(model, i, j, 0, nz)`. On a lidded
   stackup this lands on the lid plane (node `nz` *is* the lid); on an
   open stack it reaches the domain top (top boundary PEC), which is the
   degenerate open-domain reading — the gates only pin the lidded case.
3. **`with_via_at_cell(model, i, j, k_top)`** is re-expressed as
   `with_via_between(model, i, j, 0, k_top)` — bit-identical mask, same
   panic condition; the R.1 gate `engine-via-001` and the existing unit
   test keep passing untouched.

## Gates

1. **`voxel-stackup-002`** (instant, `crates/yee-voxel/tests/`): on the
   FS.4.0 3-layer lidded stack (2+3+2 cells, trace at k = 5, lid at
   k = 7): a through-via masks exactly `E_z` `k = 0..7` at its column; a
   blind via `with_via_between(…, 2, 5)` masks exactly `k = 2..5`;
   **total set-cell count over the whole mask equals the sum of the two
   columns** (nothing else was touched — the exact-count idiom is
   stronger than spot-checking neighbours); the four neighbour columns
   are spot-checked anyway; and `with_via_at_cell` ≡
   `with_via_between(0, k_top)` bit-identical on a microstrip model.
2. **`engine-stackup-via-001`** (release, ignored,
   `crates/yee-engine/tests/`): differential three-run fixture on ONE
   grid (the proven `engine-via-001` structure) moved onto a symmetric
   **stripline** built by `voxelize_stackup`:
   - control: mid-line λ/4 **open** stub → deep |S21| notch at its
     design f₀;
   - DUT: same stub **shorted to the bottom ground by a through-via**
     at its far end → the λ/4-open resonance is gone at that frequency
     (a shorted λ/4 stub is an input *open*; its own short-stub notch
     sits near 2 f₀, outside the measured band);
   - reference: bare line, shared by both ratios
     (`sparams::transmission_db` — notch-location/depth-shaped asserts,
     which is the ADR-0204-sanctioned use of the single ratio).

   Hygiene carried over from `engine-stripline-eeff-001` (ADR-0215) and
   re-derived for this fixture in the test comments: **b ≥ 16 cells**
   (confined lidded mode resolution); lateral PEC box modes are
   **obviated by CPML side walls** (the box-mode cutoff rule applied to
   the absorbing-wall fixture); no time gate (absorbing terminations +
   ring-down-length run instead of the PEC-box pulse-tail rule).
   Measured first with `--nocapture`, asserts pinned from the measured
   numbers.

## Non-goals (queued)

Finite via *barrel* diameter (> one cell column) and antipads; via pads;
buried-via layer-pair convenience taking `Stackup` layer indices
directly; via inductance extraction vs closed form; multi-trace-layer
layouts (still one `trace_layer` per voxelization); automesh awareness
of vias (FS.4.2).

## Lane

`crates/yee-voxel/**`, `crates/yee-engine/tests/**` (new files),
`docs/superpowers/specs+plans`, ADR-0221, one `SUMMARY.md` line.
