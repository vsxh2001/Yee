# ADR-0205: FS.1a.0 — truncated ground: the quasi-Yagi unblock

**Status:** Accepted
**Date:** 2026-07-08
**Related:** FULL-SUITE-ROADMAP FS.1 (antenna catalog), ADR-0194 (R.1
`with_via_at_cell`, the post-processing idiom this follows), ADR-0108 (the
ground/dielectric z-stack this leaves untouched).
**Spec:** `docs/superpowers/specs/2026-07-08-fs1a-quasi-yagi-design.md`

## Decision

`yee_voxel::truncate_ground_at_cell(model, i_ground_end)` clears the
`k = 0` tangential-E PEC beyond the first `i_ground_end` grid columns, with
an exact, documented edge rule (`Ex` nodes kept for `i < i_ground_end`,
`Ey` nodes for `i ≤ i_ground_end` — the edge plane is
`x = x₀ + i_ground_end·dx`); `i_ground_end ≥ nx` is a no-op. The 2026-07-07
antenna assessment found the quasi-Yagi blocked *solely* by the voxelizer's
full-ground assumption; this is the smallest change that removes the block.
Post-processing on `MicrostripModel` (the R.1 `with_via_at_cell` idiom) was
chosen over a `VoxelOptions`/`Layout` field because ~20 construction sites
build those types by struct literal — additive API, zero churn, and the
default path stays bit-identical by construction.

## Gate `voxel_002` (unit, milliseconds)

- Full ground pre-truncation: every `k = 0` node PEC (hand-counted
  populations `nx·(ny+1)` / `(nx+1)·ny`).
- Truncation at `g = nx/2`: populations exactly `g·(ny+1)` / `(g+1)·ny`,
  the edge checked node-by-node on both sides for every `j`, trace-layer
  masks untouched.
- The documented no-op leaves both masks **bit-identical** to an
  untruncated build.

## Queued (same spec)

- **FS.1a.1**: `yee_layout::quasi_yagi` generator (Kaneda/Deal topology on
  the FR-4-class stack; the published X-band reference's 0.09 mm CPS gaps
  are uniform-dx infeasible until FS.0b) + full-wave S11 gate
  `engine-antenna-005`.
- **FS.1a.2**: end-fire pattern + front-to-back gate via the A.2 NTFF
  machinery — the wrong-phase-balun detector.
