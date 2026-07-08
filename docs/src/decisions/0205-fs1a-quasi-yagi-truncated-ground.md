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

## FS.1a.1b — the lifted stack: SHIPPED, gate GREEN (same day)

- `voxelize_microstrip_open(layout, opts, air_below_cells)`: ground sheet
  at `k_gnd = air_below_cells` mid-domain, free air (and a bottom
  absorber's room) beneath, dielectric `k_gnd..k_top`;
  `MicrostripModel.k_gnd`; `truncate_ground_at_cell` operates at `k_gnd`.
  `air_below_cells = 0` ≡ the classic stack (which delegates).
- `AperturePortSpec::k_lo` (serde default 0, validated `< k_top`): the
  drive/measure column spans `k_lo..k_top`. The compute-side
  `AperturePort` was already cell-list based, so **no kernel change on
  either backend** — the assumption lived only in the engine translation
  and the voxelizer. All 15 construction sites pass `k_lo: 0`; every
  existing path is bit-identical.
- Gate re-run on the lifted stack + all-six-face CPML (the only PEC left
  is the masked ground sheet + traces): **dip 5.950 GHz / −20.9 dB vs the
  designed 5.80 GHz → 2.6 % error** (gate ≤ 10 %, depth pinned ≤ −10 dB),
  broadband |S11| baseline −3…−4 dB (real radiated power), matched band
  ≈ 5.8–6.4 GHz. The ε = 1.61 dipole calibration verified blind — the
  mode landed where it predicted once the antenna could radiate.

## FS.1a.2 — end-fire pattern gate: SHIPPED, GREEN first run

`engine-antenna-006` (azimuth NTFF cut at the measured 5.95 GHz, lifted
stack — the box's bottom face sits in free air below the ground sheet):
**F/B = 12.3 dB** (pinned ≥ 6), main lobe toward the director (−1.9 dB at
±30°, −6.4 dB at 60°), minimum exactly over the reflector at φ = 180°.
The balun's 180° split verified by radiation physics, not just match.
FS.1a is complete: truncated ground → generator → S11 → pattern, all
gated, both antenna gates in the antenna CI job.
