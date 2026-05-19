# Phase 2.fdtd.7 — Subgridding (Walking Skeleton)

**Status:** Draft
**Owner:** TBD
**Phase:** 2.fdtd.7.0
**Depends on:** Phase 2.fdtd.0/1/2 (Yee grid + CPML + NTFF, shipped). No strict dependency on dispersive ADE (Phase 2.fdtd.3) or TF/SF (Phase 2.fdtd.5), but composes with them in later sub-phases.
**Blocks:** End-to-end PCB / antenna validation that requires sub-cell features — thin laminate substrates, slot antennas, fillets, fine conductor edges — without globally refining the uniform Yee grid.

## Motivation

The current `YeeGrid` is uniform: `dx = dy = dz`, set once at construction and applied across the whole computational volume. That works for the resonant-cavity and free-space-dipole gates we ship today, where the smallest feature (a quarter wavelength) is comfortably ≥ 10 cells. It does **not** work for the production targets called out in `ROADMAP.md` Phase 2:

- **Thin laminate substrates.** A 2.4 GHz patch on 0.508 mm RO4003C sits over a dielectric `≈ λ_eff / 80` thick. Resolving it as ≥ 10 cells across forces `dx ≤ 50 µm` globally; at that resolution a 100 × 100 mm test domain is `2000³ ≈ 8 · 10⁹` cells — out of reach on a single workstation and not a productive use of multi-GPU either.
- **Slot antennas and fine apertures.** Slot widths of `λ / 200` are routine. The CPML reflection gate (Phase 2.fdtd.1) is unaffected; what's affected is whether the slot field is resolved at all.
- **Conductor edges and fillets.** Field singularities at conductor corners benefit super-linearly from local refinement; staircasing them on a coarse uniform grid corrupts the radiated pattern at the few-percent level that NTFF validation cares about.

Without subgridding we accept one of three losses: (a) coarse uniform grid + lost accuracy on the thin feature, (b) fine uniform grid + intractable memory / runtime, (c) conformal techniques alone (Dey-Mittra, Phase 2 deliverable) which help with conductor-edge staircasing but do not help with the substrate-thickness problem. Conformal and subgridding are complementary, not substitutes.

This spec defines the **walking-skeleton** scope (Phase 2.fdtd.7.0) — the smallest end-to-end nest that exercises the coarse / fine interface contract — and the deferred phases that flesh it out.

## Scope decision

Phase 2.fdtd.7.0 (this spec) ships:

- **Single nested region.** One axis-aligned, cuboidal fine sub-region inside one coarse parent `YeeGrid`. Multiple nests and nest-inside-nest are deferred to 7.1+.
- **2× refinement ratio, isotropic.** `dx_fine = dx_coarse / 2`, ditto `dy`, `dz`. Other integer ratios (3×, 5×) and anisotropic ratios are deferred to 7.4.
- **Time-subcycling at 2×.** The fine grid steps twice per coarse step, with `dt_fine = dt_coarse / 2`, so the fine CFL is satisfied without throttling the coarse step. The Chevalier / Okoniewski temporal-interpolation pattern (Okoniewski 1997, Chevalier 1997) supplies fine-grid boundary `E` values at the half-step.
- **Linear spatial interpolation** of coarse → fine boundary `E` values (Chevalier 1997 §III). Fine → coarse update uses area-averaging of fine `H` to the coarse face. Higher-order (cubic) interpolation is deferred to 7.4.
- **Non-dispersive, isotropic, lossless materials inside the fine region.** The fine grid carries the same `eps_r`, `mu_r` scalars as today's `YeeGrid`. Dispersive ADE inside the fine region is deferred to 7.2 explicitly because the ADE auxiliary fields must be re-staggered and re-clocked at the interface and that is its own piece of work.
- **No CPML or TF/SF box inside the fine region.** CPML lives only on the outer faces of the coarse parent; TF/SF, if present, lives entirely in the coarse grid. The fine region is assumed to sit interior to both. Co-location of a fine region with a CPML or TF/SF face is a documented runtime error in 7.0 and a feature in 7.3.
- **CPU-only, single-threaded, scalar FP64.** Same execution model as today's `WalkingSkeletonSolver`. GPU lands later in the standard Phase-2 GPU pass.

Out of scope for 7.0 (numbered phases below):

- 7.1 — multiple nests, nest-inside-nest.
- 7.2 — dispersive ADE inside the fine region.
- 7.3 — co-location with CPML / TF/SF faces.
- 7.4 — higher refinement ratios and higher-order spatial interpolation.
- 7.5 — production validation case (see §6).

## Interface design

A 2× nest places the fine grid origin on a coarse `E`-node and extends `(2 · N_x_fine, 2 · N_y_fine, 2 · N_z_fine)` fine cells. The coarse / fine interface is the six faces of the cuboid. Per face, two update paths cross the boundary:

- **Coarse → fine (`E`-field driver).** At the start of each fine step, the fine grid's outer `E_t` (tangential `E`) needs a Dirichlet-style value sampled from the coarse grid. The coarse-grid `E_t` lives on coarse edges that lie inside the fine-grid boundary plane. Where a fine boundary edge coincides with a coarse edge, we copy. Where it falls between two coarse edges (every other fine edge does, by the 2× ratio), we **linearly interpolate** between the two flanking coarse `E_t` samples. Temporal interpolation: a fine step lands at `t = n · dt_coarse + k · dt_fine` for `k ∈ {1, 2}`; for `k = 1` the coarse `E_t` is interpolated linearly in time between its values at the bracketing coarse steps.
- **Fine → coarse (`H`-field consistency).** After the fine grid has completed both sub-steps for the coarse interval, the coarse grid's `H` components that lie on the interface need a consistent value. The Chevalier scheme replaces the coarse-grid `H_t` on the interface face with the **area-average** of the four fine-grid `H_t` cells covering that coarse face. This is the step that closes the energy balance: without it, the coarse stencil sees a stale `H` on its interior face and the late-time energy drifts.

Field-component ownership at the interface (2× nest, axis-aligned face normal to `+x`):

- The coarse `E_y`, `E_z` on the interface face are **owned by the fine grid**: the coarse `E_t` arrays are overwritten with the area-/edge-average of the fine `E_t` after each fine sub-step. (Symmetric to Berenger 2003 §III.)
- The coarse `H_x` on the interface face is owned by the coarse grid (it lives on a coarse face, no fine counterpart on the same plane).
- The fine `E_y`, `E_z` on the interface face are **owned by the coarse grid** via the interpolation above.
- The fine `H_x` on the interface face is owned by the fine grid as usual.

This ownership pattern (E owned by the finer grid, tangential coupling is one-way per field) is the Chevalier 1997 prescription and is what keeps the interface reciprocal at second order in `dx_coarse`.

Time-stepping pattern per coarse step `n → n + 1`, factoring CPML and source updates as no-ops on the fine grid in 7.0:

```text
1.  coarse: update_h          (H^{n+1/2}_c)
2.  coarse: cpml_h, source_h
3.  fine sub-step k = 1:
     a. interpolate coarse E_t at t = n + 1/4 onto fine boundary E_t
     b. fine: update_h        (H^{n+1/4}_f)
     c. fine: update_e        (E^{n+1/2}_f)
4.  fine sub-step k = 2:
     a. interpolate coarse E_t at t = n + 3/4 onto fine boundary E_t
     b. fine: update_h        (H^{n+3/4}_f)
     c. fine: update_e        (E^{n+1}_f)
5.  coarse: update_e          (E^{n+1}_c), using H^{n+1/2}_c
6.  coarse: cpml_e, source_e
7.  overwrite coarse E_t on the interface face with the area-average of
    fine E_t (close the reciprocity loop)
```

References:

- Okoniewski, Okoniewska, Stuchly, "Three-dimensional subgridding algorithm for FDTD", IEEE T-AP 1997, 45(3) — the foundational temporal-subcycling pattern.
- Chevalier, Luebbers, Cable, "FDTD local grid with material traverse", IEEE T-AP 1997, 45(3) — the spatial-interpolation prescription and the energy-balance closure we use.
- Berenger, "A Huygens subgridding for the FDTD method", IEEE T-AP 2006, 54(12) — the late-time-stability analysis that motivates the area-average fine → coarse coupling.
- Taflove & Hagness, *Computational Electrodynamics* 3rd ed., §13 — survey and stability discussion.

## API sketch

The walking skeleton adds one new type and one new solver wrapper. The existing `WalkingSkeletonSolver` is unchanged; the subgridded variant composes it.

```rust
/// Axis-aligned, cuboidal sub-region nested at 2× resolution inside a
/// parent [`YeeGrid`].
///
/// Phase 2.fdtd.7.0: single nest, isotropic 2× refinement, time-subcycling
/// at 2×, non-dispersive isotropic materials, no co-location with CPML or
/// TF/SF faces.
pub struct SubgridRegion {
    /// Coarse-cell indices of the nest corner (inclusive lo, exclusive hi).
    pub lo: (usize, usize, usize),
    pub hi: (usize, usize, usize),
    /// The fine grid. `dx_fine = dx_coarse / 2`; `dt_fine = dt_coarse / 2`.
    fine: YeeGrid,
    /// Cached coarse `E_t` snapshots at the bracketing coarse steps,
    /// used for temporal interpolation during fine sub-steps.
    e_t_snapshots: InterfaceSnapshots,
}

impl SubgridRegion {
    /// Build a 2× nest covering coarse cells `lo..hi`. The fine grid
    /// inherits `eps_r`, `mu_r` from the parent unless overridden after
    /// construction.
    pub fn new(parent: &YeeGrid, lo: (usize, usize, usize), hi: (usize, usize, usize)) -> Self;

    /// Borrow the fine grid (e.g. to set per-cell materials before the
    /// first step).
    pub fn fine_grid(&self) -> &YeeGrid;
    pub fn fine_grid_mut(&mut self) -> &mut YeeGrid;
}

/// Subgridded driver that wraps a [`WalkingSkeletonSolver`] and one
/// [`SubgridRegion`].
///
/// Composes with CPML on the parent's outer faces (configure on the
/// inner [`WalkingSkeletonSolver`] as today). Does NOT compose with a
/// CPML or TF/SF face that intersects the nest — that's Phase 2.fdtd.7.3.
pub struct SubgriddedSolver {
    inner: WalkingSkeletonSolver,
    region: SubgridRegion,
}

impl SubgriddedSolver {
    pub fn new(inner: WalkingSkeletonSolver, region: SubgridRegion) -> Self;

    /// Advance one coarse step, performing two fine sub-steps in between.
    pub fn step(&mut self);
}
```

`SubgriddedSolver::step` implements the seven-stage sequence in the interface-design section. The parent `WalkingSkeletonSolver::grid_and_cpml_mut` split-borrow (already exists in `lib.rs` at line 158) is the primitive used to drive coarse `update_h` / `update_e` without re-implementing the CPML wiring. The fine-grid kernels reuse `update::update_h` / `update::update_e` unchanged — the only new code is the interface module that ships interpolation and area-averaging.

## Stability and reciprocity considerations

Late-time instability is the classical subgridding failure mode (Berenger 2003, §IV). The mechanism is asymmetric coupling: if the coarse → fine path adds energy at a different rate than the fine → coarse path extracts it, an unphysical eigenmode at the interface grows exponentially over `O(10⁴)` time steps even when the per-step coupling error is `< 10⁻⁶`.

Two regression tests gate the interface:

- **Round-trip energy test.** Initialise a Gaussian-modulated sinusoidal pulse inside the fine region, propagate forward through the interface into the coarse region, reflect off PEC walls, propagate back through the interface, and integrate `∫ ε_0 |E|² + µ_0 |H|² dV` over the whole domain at `t = 0` and `t = N · T_round_trip` with `N ≥ 50`. Energy must not drift by more than **0.5%** of the initial energy. This catches both the asymmetric-coupling failure mode and gross reciprocity breakage.
- **Plane-wave traversal contrast.** Drive a TF/SF plane wave (Phase 2.fdtd.5) through the coarse grid such that the fine region sits inside the TF box; measure the scattered field outside the TF box. With the fine region empty (vacuum), the scattered field must be **≥ 60 dB below** the incident field. Anything worse means the interface is radiating spuriously.

Both tests are run for `N ≥ 10000` coarse steps to surface late-time growth.

The reciprocity guarantee comes from the area-average fine → coarse coupling closing the discrete energy balance to second order in `dx_coarse`, cf. Chevalier 1997 §IV. We are not attempting Berenger's Huygens-surface scheme (Berenger 2006), which has stronger stability guarantees but a much heavier implementation cost; if the tests above fail at high cell counts we revisit that choice.

## Validation gate

**fdtd-007 — dielectric-loaded thin slot antenna ([TBD verify which Maloney-Smith paper] — see ADR-0041).**

Geometry: a thin slot of width `w = 0.5 mm` and length `L = 30 mm` cut into an infinite PEC ground plane, with a dielectric slab (`ε_r = 2.2`, thickness `h = 1.524 mm`) backing the slot. Driven by a delta-gap voltage at the slot centre. Coarse grid `dx_coarse = 1 mm`, fine grid `dx_fine = 0.5 mm` covering a `(40 × 6 × 4) mm` box centred on the slot and substrate.

- **Reference: [TBD — citation under review.](../../src/decisions/0041-fdtd-007-reference-correction.md)** The original draft cited Maloney & Smith, "A study of transient radiation from the Wu-King resistive monopole — FDTD analysis and experimental measurements", IEEE T-AP 1993, 41(5), pp. 668–676, Fig. 9. **That paper is a cylindrical-monopole paper and does not contain the cited slot geometry** (Track UUUUUUUU finding, commit `d56c460`, verified by Track XXXXXXXX against IEEE Xplore document 222286). The `f_res = 8.9 GHz` and `|S_11| = -22 dB` figures encoded in `yee_validation::FDTD_007_FRES_REF_HZ` / `FDTD_007_S11_DB_REF` are therefore unverified; the physics gates remain `#[ignore]`'d pending `fdtd-007.1` resolution. See ADR-0041 for the chain of evidence and the candidate follow-ups.
- Tolerance: resonant frequency within **±2%**, `|S_11|` at resonance within **±1 dB** — applies once the reference is resolved.
- Comparator: the **same problem run on a globally uniform `dx = 0.5 mm` grid** must produce a result within 0.3% / 0.3 dB of the subgridded result. This is the internal sanity check; the (resolved) external reference is the external gate.

Run time budget: `< 30 min` wall-time in `--release` on the CI Linux runner. If it overruns, the test is hardware-gated like the cuSOLVER tests (Phase 1.5 precedent).

## Risks / open questions

- **Late-time instability at the interface.** Standard subgridding failure mode; mitigated by the Chevalier area-average and gated by the 10k-step round-trip energy test. If 2× linear interpolation proves unstable we fall back to the Huygens-surface variant (Berenger 2006), which is a strictly bigger change.
- **Fine CFL throttling the global step.** Time-subcycling (2× sub-steps) precisely avoids this in 7.0. At higher refinement ratios (7.4) the sub-step count grows linearly and the fine-grid update cost can dominate; that's a budgeting concern, not a correctness one.
- **Co-location with CPML.** Out of scope for 7.0; documented as a runtime error. The CPML auxiliary `Ψ` arrays would need their own interface-coupling story (they don't satisfy the same energy balance as the bare Yee stencil), which is real work and lands in 7.3.
- **Interaction with NTFF.** NTFF samples on a closed Huygens surface that today is grid-aligned on the coarse grid. The walking-skeleton requirement is that the NTFF surface not intersect a fine region; if it does, the NTFF DFT bins are sampled inconsistently and the radiated pattern is corrupted. Surface this as a runtime error in 7.0.
- **Conformal techniques (Dey-Mittra) inside the fine region.** Not a 7.0 concern (we are still on staircase geometry workspace-wide) but called out as a future composition point: Dey-Mittra cell-cut coefficients are local to each cell and should compose cleanly with a finer Yee grid.

## Dependencies

No strict new external dependency. The implementation uses `ndarray` (already in the dep graph) for the fine-grid arrays and reuses `update::update_h` / `update::update_e` unchanged.

The one internal refactor that makes this cleaner — but is not strictly required — is splitting `WalkingSkeletonSolver::step` into per-stage helpers (`update_h_only`, `apply_cpml_h`, `update_e_only`, `apply_cpml_e`) so the subgridded driver can interleave fine sub-steps without re-implementing CPML wiring. The `grid_and_cpml_mut` split-borrow at `crates/yee-fdtd/src/lib.rs:158` already provides the lower-level primitive; the per-stage helpers are a quality-of-life improvement and not blocking.

## Phase numbering

- **Phase 2.fdtd.7.0** — walking skeleton: single nest, 2×, axis-aligned, time-subcycling, non-dispersive, interior to CPML/TF-SF. **This spec.**
- **Phase 2.fdtd.7.1** — multiple disjoint nests; nest-inside-nest.
- **Phase 2.fdtd.7.2** — dispersive ADE materials inside the fine region (Drude / Lorentz / Debye, composing with Phase 2.fdtd.3).
- **Phase 2.fdtd.7.3** — co-location of fine region with a CPML or TF/SF face.
- **Phase 2.fdtd.7.4** — higher refinement ratios (3×, 5×) and higher-order spatial interpolation (cubic).
- **Phase 2.fdtd.7.5** — production validation: openEMS cross-check on a 2.4 GHz inset-fed patch on RO4003C with the substrate inside a fine region, to within ±0.1 dB on the in-band `S_11` minimum.

## Lane

`crates/yee-fdtd/**` only. New module `crates/yee-fdtd/src/subgrid.rs` plus integration tests under `crates/yee-fdtd/tests/`. No changes to `yee-core`, `yee-mesh`, `yee-cli`, or `yee-py` in 7.0.

## References

- Okoniewski, M., Okoniewska, E., Stuchly, M. A., "Three-dimensional subgridding algorithm for FDTD", *IEEE Trans. Antennas Propag.* 45(3), 1997, pp. 422–429.
- Chevalier, M. W., Luebbers, R. J., Cable, V. P., "FDTD local grid with material traverse", *IEEE Trans. Antennas Propag.* 45(3), 1997, pp. 411–421.
- Berenger, J.-P., "A Huygens subgridding for the FDTD method", *IEEE Trans. Antennas Propag.* 54(12), 2006, pp. 3797–3804.
- Maloney, J. G., Smith, G. S., "A study of transient radiation from the Wu-King resistive monopole — FDTD analysis and experimental measurements", *IEEE Trans. Antennas Propag.* 41(5), 1993, pp. 668–676. **[TBD verify — this paper is a cylindrical-monopole paper, not a slot-antenna paper; the original `fdtd-007` attribution is disputed. See [ADR-0041](../../src/decisions/0041-fdtd-007-reference-correction.md).]**
- Taflove, A., Hagness, S. C., *Computational Electrodynamics: The Finite-Difference Time-Domain Method*, 3rd ed., Artech House 2005, §13.
