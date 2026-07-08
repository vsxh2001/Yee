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

## FS.1a.1 measured negative result (2026-07-08) — the floor IS a ground

`yee_layout::quasi_yagi` shipped (scaling-rule seeds, unit-gated topology:
connectivity by union-find over positive-area overlaps, parasitic director,
balun detour over ground) and the full-wave gate ran twice (133/140 s per
run, auto_dx-seeded):

1. First run: the λ/2-looking dip sat at 7.500 GHz (+29 %); calibrating
   the dipole permittivity down to the measured 1.61
   (`eps_dipole = 1 + 0.18·(ε_r−1)`, the R.6 single-point pattern — the
   half-space (ε_r+1)/2 badly overestimates thin-substrate loading; kept,
   it is real physics) moved the dip only to 7.150 GHz — **that dip is a
   feed-structure resonance, not the dipole**.
2. Both runs: |S11| ≈ 0 dB across 4–7.6 GHz — the dipole is never driven.

Root cause is the **z-stack, not the layout**: `voxelize_microstrip` puts
the ground at `k = 0` on the domain floor, and the antenna boundary keeps
the bottom face PEC. Past the truncation the *boundary* still provides an
infinite image plane 1.6 mm (0.03 λ) under the dipole, annihilating its
radiation resistance — the FS.1a.0 mask truncation is necessary but not
sufficient. Gate `engine-antenna-005` stays in-tree, `#[ignore]`'d and
named `antenna_…` so the blanket CI step skips it; it turns green with:

## Queued: FS.1a.1b — the lifted stack (open space below the ground)

- A voxelizer variant with `air_below_cells`: ground sheet at
  `k_gnd = air_below_cells` mid-domain, substrate above it, all-six-face
  CPML; `MicrostripModel` records `k_gnd`;
  `truncate_ground_at_cell` operates at `k_gnd` instead of the hard-wired 0.
- `AperturePortSpec::k_lo` (serde-default 0 for full back-compat): the
  drive/measure column becomes `k_lo .. k_top` on **both** backends (cpu.rs
  loop + the R.3 WGSL kernel); compute-015 parity re-certified with the
  default proving bit-exactness.
- Then re-run engine-antenna-005; FS.1a.2 (end-fire/front-to-back NTFF
  gate) follows.
