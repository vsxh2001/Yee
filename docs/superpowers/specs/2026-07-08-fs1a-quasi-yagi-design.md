# FS.1a — quasi-Yagi antenna: truncated ground + generator + gates

**Date:** 2026-07-08
**Track:** FULL-SUITE-ROADMAP FS.1 (antenna catalog). The 2026-07-07
antenna-support assessment found the quasi-Yagi *blocked* solely by
`voxelize_microstrip`'s full-ground assumption — the topology is otherwise
planar rectangles, exactly what the certified pipeline handles.

## Why quasi-Yagi first

It is the highest-value unlock per unit work: a directive, end-fire,
moderate-gain printed antenna (the workhorse for point-to-point links and
arrays), all-planar (no vias in the classic Deal/Kaneda/Qian/Itoh form),
and its defining structural feature — the **truncated ground plane acting
as the reflector element** — is a one-option change to the voxelizer.

## Decomposition (walking-skeleton first)

- **FS.1a.0 — truncated ground (this increment):**
  `VoxelOptions::ground_x_max_m: Option<f64>` — when set, the k = 0
  tangential-E PEC ground covers only cells whose x ≤ the given coordinate
  (grid frame = layout frame, same `x0 − margin` origin as the traces);
  `None` (default) is the full ground plane, **mask-identical** to today.
  Unit gates: default masks byte-equal to the pre-option output on a
  reference layout; truncated masks PEC exactly where specified; the
  engine/board fixture passes the option through untouched.
- **FS.1a.1 — generator + S11 gate:** `yee_layout::quasi_yagi(f0,
  substrate, z0)` emitting the Kaneda/Deal topology (50 Ω feed →
  microstrip T-junction balun with half-guide-wavelength delay arm → CPS
  section → driven dipole arms + one director; ground truncated ~λg/4
  behind the driven element), dimensions from the published scaling rules
  (driven ≈ 0.45–0.5 λ_diel, director ≈ 0.3 λ0, element spacings 0.15–0.3
  λ0). Gate `engine-antenna-005`: full-wave directional |S11| dip within
  the A-track tolerance of the design f0 at ≤ −10 dB.
- **FS.1a.2 — end-fire pattern gate:** NTFF cut (the A.2 machinery):
  main lobe toward the director (end-fire), front-to-back ratio ≥ a
  measured-then-pinned floor. This is the quasi-Yagi's *purpose* made
  machine-checkable; also the first pattern gate on a non-broadside
  radiator.

## Constraints and risks

- The published X-band references (Deal et al., εr = 10.2, h = 0.635 mm,
  0.09 mm CPS gaps) are **uniform-dx infeasible** (auto_dx feature rule →
  45 µm cells). FS.1a.1 therefore designs on the FR-4-class stack the rest
  of the repo certifies (mm-scale features) and gates against the design
  frequency + pattern physics, not against a paper's exact curve. A
  paper-exact reproduction becomes feasible with FS.0b's graded grid.
- The balun's phase balance is the classic failure point; the S11 gate
  will catch gross failure, the FS.1a.2 pattern gate catches polarity
  errors (a wrong-phase balun still matches but radiates broadside/split).

## Gates

- `voxel_002` (FS.1a.0, unit, instant): default = mask-identical;
  truncation boundary exact.
- `engine-antenna-005` (FS.1a.1, release): S11 dip location + depth.
- `engine-antenna-006` (FS.1a.2, release): end-fire lobe + F/B floor.
