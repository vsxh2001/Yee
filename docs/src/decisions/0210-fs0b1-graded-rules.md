# ADR-0210: FS.0b.1 — graded mesh rules + graded voxelization: the payoff step

**Date:** 2026-07-11. **Status:** accepted. **Phase:** FS.0b.1
(`FULL-SUITE-ROADMAP.md`). **Spec:**
`docs/superpowers/specs/2026-07-11-fs0b1-graded-rules-design.md`. **Plan:**
`docs/superpowers/plans/2026-07-11-fs0b1-graded-rules.md`.

## Context

FS.0b.0 (ADR-0208) certified the graded CPU kernel (compute-018 bit-exact
on uniform, compute-019 taper reflection −52.7 dB) but nothing generated
spacings from a `Layout` and `yee-voxel` was uniform-only. The measured
motivation is ADR-0204's: the uniform convergence residual
(max linear Δ|S| = 0.198) sits entirely in staircase-limited feature
regions, and the next uniform pass costs ~2.4 h (19 M cells) — refining
everywhere because you cannot refine somewhere.

## Decision

- **Rule generator** `yee_engine::automesh::auto_spacings(layout, f_max,
  &GradedMeshOptions) -> AutoSpacings` (per-axis primal widths + the
  origin/z-stack metadata the voxelizer needs):
  - **Coarse ceiling** = `auto_dx_bulk` — λ/20-in-dielectric and h/3,
    clamped to [1 µm, 1 mm], deliberately **without** the uniform
    rulebook's feature/2 term: features refine **locally**, so a narrow
    gap no longer caps the whole domain (the graded payoff, unit-gated).
  - **Fine bands** at `fine = min(min_feature/2, coarse/2)` inside
    `edge ± guard` around every trace-AABB edge and across every
    inter-trace axis gap; `guard` defaults to the substrate height (the
    fringing length scale).
  - **Geometric taper** between fine and coarse with cell-to-cell ratio
    ≤ 1.3 including junction steps (the compute-019-certified regime;
    growth > 1.3 is rejected).
  - **z:** substrate = `ceil(h/(coarse/2))` cells of exactly `h/n_sub`
    (the ADR-0108 no-air-gap stack), then geometric growth into the air.
  - **Absorbers:** the `npml` outermost x/y cells stay bit-equal coarse
    (the FS.0b.0 `validate_cpml_layers` scope rule); a fine band reaching
    an absorber is an **error**, not a silent clamp. `JobSpec::dx_m`
    stays `coarse`, the nominal spacing the CPML σ_max recipe assumes.
- **Graded voxelizer** `yee_voxel::voxelize_microstrip_graded(layout,
  &GradedVoxelGrid) -> GradedMicrostripModel`: the `voxelize_inner`
  rasterization (ground sheet, exact substrate fill, cell-centre
  point-in-polygon trace masks) against true per-axis coordinates;
  returns raw eps/PEC `Array3`s + node coordinates, not a `YeeGrid`
  (whose scalar dx/dt are meaningless graded). Port cells by coordinate
  lookup (`partition_point` over nodes).
- **Coordinate generation is run-wise, not a running sum.** Coordinates
  are computed per maximal run of identical widths as `origin + (base +
  m·d)` / `origin + (base + (m+0.5)·d)`, so constant arrays reproduce the
  uniform voxelizer's `x0 + i·dx` / `x0 + (i+0.5)·dx` **bit-exactly** by
  construction. This was forced by measurement, not taste: the first
  implementation used a naive cumulative sum, and gate `voxel-graded-001`
  caught the Ex PEC masks diverging — trace edges sit exactly on node
  planes (margins are whole cells off the bbox), and the ~1 ulp
  cumulative-sum offset flips the even-odd point-in-polygon test there.
  The FS.0b.0 lesson (choose expressions that degenerate bit-exactly)
  applies to rasterization too.

## Rule iterations (measured)

- **Iteration 0 — the feature-rule-only fine spacing (rejected without a
  solve).** With `fine = min_feature/2` alone, the FS.0a stub board gets
  **no refinement at all** (min_feature/2 = 1.5 mm > coarse 0.533 mm):
  the generated mesh is identical to the uniform pass-0 grid, whose notch
  ADR-0204 already measured at **5.100 GHz — 5.2 % from the converged
  4.850 GHz**, outside the ±2 % gate. Rejected by that existing
  measurement; the `coarse/2` term in the fine rule is the fix (one
  halving at the staircase-limited edges is exactly what the FS.0a
  uniform trajectory 0.533 → 0.267 mm showed sufficient).
- **Iteration 1 — the shipped rules (first solved iteration): PASS, no
  retuning needed.** `fine = min(min_feature/2, coarse/2)`, guard = h =
  1.6 mm, substrate dz = coarse/2-class. Measured 2026-07-11 on the S.6
  stub board (f_max 6 GHz): coarse 0.533 mm, fine 0.267 mm, grid
  282×110×41 = **1,271,820 cells**, k_top = 6, dt = 4.622e-13 s,
  10 125 steps. Notch **4.900 GHz @ −37.2 dB** → **1.03 %** from the
  uniform-converged 4.850 GHz (gate ±2 %). For calibration: 4.900 GHz
  is exactly the uniform pass-1 (dx = 0.377 mm) notch, one 50 MHz bin
  from the pass-2 answer — refined-mesh physics, not the coarse pass-0
  5.100 GHz. Neither the guard width nor the fine spacing needed a
  second solved iteration.

## Measured results (gates)

- **voxel-graded-001** (`voxel_graded_001_uniform_bitexact.rs`, fast):
  constant arrays vs `voxelize_microstrip` — eps, Ex/Ey PEC masks, dims,
  and port cells exactly equal (plus graded-z-stack and nonuniform
  coordinate-lookup sanity). PASS.
- **automesh unit gates** (fast): growth ratio ≤ 1.3 everywhere on all
  three axes (junctions included); every trace edge ± guard and every
  inter-trace gap covered by fine cells; absorber layers bit-equal
  coarse; total length covers the domain with < one coarse cell of
  overshoot and `x0` exact; single-rect layout near-uniform (bit-equal
  coarse away from the four edge bands, coarse cells dominate);
  growth > 1.3 and absorber-colliding guards are errors. PASS.
- **engine-graded-001** (`engine_graded_notch.rs`, `#[ignore]`, release):
  the S.6 stub-notch board on the `auto_spacings` grid — JobSpec built
  directly on the graded voxelizer; DUT + through-line reference on the
  same DUT-derived grid; double-ratio |S21| over 3.5–6 GHz / 50 MHz;
  probe triples on scanned uniform-coarse stretches (12 coarse cells =
  the FS.0a 6.4 mm — `fit_standing_wave` needs equal spacing, so a
  triple must never straddle a taper). **Measured 2026-07-11: notch
  4.900 GHz @ −37.2 dB → err 1.03 % vs the uniform-converged 4.850 GHz
  (assert ≤ 2 %); depth ≤ −20 dB with 17 dB to spare;
  cells_graded / cells_uniform = 1,271,820 / 6,679,200 = 0.190** — the
  uniform comparison grid is the dx = 0.267 mm pass-2 fixture
  (506×176×75) built with the FS.0a loop's exact rescaled options (not
  solved). **Ratio pinned at 0.25** (~24 % margin). Runtime: reference
  solve 152.4 s, both solves 306.1 s release on 4 cores at the same
  10 125-step window and dt as the uniform pass — the graded pair
  replaces a ~5×-costlier uniform pass-2 pair at matching physics. PASS
  on the first solved rule iteration.

## Consequences

- The graded pipeline is end-to-end real: rules → voxelization → the
  FS.0b.0 kernel — but only the gate drives it. FS.0b.2 should rewire
  `board.rs` / `converge_two_port` so every consumer (studio, Python,
  WS) gets graded passes, and take the GPU kernel.
- `voxelize_finite_board` / lifted-ground variants stay uniform-only;
  the graded voxelizer covers the classic floor-ground stack.
- Fine spacing is a single global value per mesh; per-band fine values
  (e.g. a 0.1 mm gap next to a 0.5 mm gap) are a later increment.
