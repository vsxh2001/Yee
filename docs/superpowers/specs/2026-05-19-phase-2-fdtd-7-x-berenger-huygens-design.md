# Phase 2.fdtd.7.x — Berenger 2006 Huygens-surface subgridding (closure rewrite)

**Status:** Draft
**Owner:** TBD
**Phase:** 2.fdtd.7.x
**Depends on:** Phase 2.fdtd.7 Q1–Q5 + Q4.1 helper (merged through `618b49b`). Reuses the `SubgridRegion`, `SubgriddedSolver`, coarse→fine interpolation surface (Q3), and the `snapshot_fine_h_mid_step` time-centering helper (Q4.1).
**Blocks:** Q5 strict 0.5%-of-peak plane-wave-traversal gate, Q6 10 000-step round-trip energy gate, Q7 fdtd-007 Maloney-Smith production gate. All three are currently `#[ignore]`'d or unimplemented because the spec §3 Okoniewski-style direct-copy closure is discrete-energy-balance unstable.

## Motivation

Phase 2.fdtd.7.0 (spec `2026-05-18-phase-2-fdtd-7-subgridding-design.md`) chose a Chevalier 1997 / Okoniewski 1997 style **direct-field-copy** closure for the fine↔coarse interface:

- Coarse → fine: linear spatial + temporal interpolation of `E_t` (Q3).
- Fine → coarse: area-averaged `H_t` overwrite + edge-averaged `E_t` overwrite (Q4).

The spec §6 risk register flagged late-time interface instability as the load-bearing risk and explicitly named **Berenger 2006's Huygens-surface scheme** as the documented fallback if 2× linear interpolation proved unstable.

Track VVVVVV (merge `a2abb4c`, commit `72c825c`) attempted to clear the strict 0.5%-of-peak Q5 traversal gate by time-centering the fine `H_t` average (adding the `snapshot_fine_h_mid_step` helper) and diagnosed the failure as **fundamental, not a phase bug**:

> The `overwrite_coarse_e_from_fine` reads fine boundary `E_t` cells that are Dirichlet-set by interpolation (not updated by `update_fine_e`), round-tripping an interpolation result rather than a fresh fine evolution; the H closure simultaneously writes the time-averaged fine `H` back to coarse, but coarse `update_e` already used the stale pre-closure coarse `H` to update coarse `E_t` at the interface plane.

That is: the two-direction direct-copy closure cannot be made discrete-energy-balanced while the interface plane is simultaneously a Dirichlet boundary for the fine grid and a sample source for the coarse grid. The interface plane is over-determined.

This spec specifies the Berenger 2006 fallback. The fine↔coarse closure is replaced by a one-directional **equivalent-current re-radiation** scheme that decouples the two grids at the field level: the fine grid evolves freely, and its boundary fields are translated into equivalent `J` / `M` surface currents that drive the coarse grid. The coarse grid never overwrites fine-grid storage and never reads from a Dirichlet plane.

## Goal

Replace the spec `2026-05-18-…-subgridding-design.md` §3 stage-6 (`average_fine_h_to_coarse`) and stage-7 (`overwrite_coarse_e_from_fine`) direct-copy closure with a Berenger 2006 Huygens-surface scheme that:

1. **Retires the Q5 strict 0.5%-of-peak agreement gate** over the first 500 steps of the canonical `(64, 32, 32)` coarse + `(16, 16, 16)`-coarse-cell fine plane-wave traversal test (`crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs` — un-`#[ignore]` it).
2. **Restores discrete energy balance on the Q6 10 000-step round-trip case** (`|W(10000) − W(0)| / W(0) ≤ 0.5%`).

This is a **closure-only rewrite.** The Q1 step-refactor, Q2 `SubgridRegion` scaffold, Q3 coarse→fine `E_t` interpolation, and Q4.1 `snapshot_fine_h_mid_step` helper are all kept; only the spec §3 stage-6/7 fine→coarse direction is replaced.

## Background

Berenger, J.-P., "A Huygens subgridding for 3-D FDTD", *IEEE Trans. Antennas Propag.* 54(12), 2006, pp. 3797–3804 — DOI `10.1109/TAP.2006.886504`.

The key idea: the equivalence principle says any field exterior to a closed surface `S` can be reproduced by suitable equivalent surface currents on `S`. Concretely, given the field `(E, H)` on the fine side of the interface, the surface currents

```text
J = +n̂ × H     (electric current on S)
M = −n̂ × E     (magnetic current on S)
```

placed on `S` re-radiate exactly the fine-side fields into the exterior region. The coarse grid, surrounding the fine grid, sees only these surface currents — never the fine-grid storage itself. No coarse-side field is ever overwritten by a fine-side sample; the coupling is one-directional per face per step:

- **Fine → coarse:** equivalent currents `J`, `M` on the Huygens surface drive the coarse `update_e` / `update_h` as right-hand-side source terms (analogous to TF/SF, Phase 2.fdtd.5).
- **Coarse → fine:** unchanged from spec `2026-05-18` Q3 — linear spatial + temporal interpolation supplies the fine grid's outer `E_t` Dirichlet boundary. This direction has never been an instability source because the fine grid is naturally a closed system once its boundary is specified.

The closure of the energy balance comes from the equivalence principle itself, not from a discrete area-average. There is no round-trip overwrite, so there is no over-determined interface plane.

## Mathematical formulation

### Total-field / scattered-field convention

The fine subdomain carries the **total field** `(E_tot, H_tot)`. The coarse grid carries the **scattered field** `(E_sc, H_sc)` *inside* the fine box's coarse-cell footprint and the **total field** *outside* it. (This is symmetric to the TF/SF slab convention from Phase 2.fdtd.5: the Huygens surface separates a TF region from an SF region.) The choice is documented as part of the API; the alternative (TF outside, SF inside) is equivalent up to sign and is not adopted here.

### Surface currents

On each of the six Huygens faces of the fine box, with outward unit normal `n̂` pointing from the fine (TF) region into the coarse (SF) region:

```text
J_S(r, t) = +n̂ × H_tot(r, t)        on S, sourced by coarse update_e
M_S(r, t) = −n̂ × E_tot(r, t)        on S, sourced by coarse update_h
```

`E_tot` is sampled on the fine grid's outer-layer `E_t` edges (the same edges Q3's interpolation writes — but at the post-`update_fine_e` time level, not the Dirichlet-set time level — see §4). `H_tot` is sampled on the fine grid's outer-layer `H_t` faces at the **time-centered** mid-step level supplied by the Q4.1 `snapshot_fine_h_mid_step` helper.

### Time-stepping pattern

Per coarse step `n → n + 1`:

```text
1.  region.snapshot_coarse_e_t(parent)             (Q3, start-of-step)
2.  coarse: inner.update_h_only                    (H^{n+1/2}_c, no source yet)
3.  fine sub-step k = 1:
     a. region.interpolate_coarse_e_to_fine(0.25)  (Q3)
     b. region.snapshot_fine_h_mid_step            (Q4.1, captures H^{n+1/4}_f)
     c. region.update_fine_h                       (H^{n+3/4}_f after this — see Q5 note)
     d. region.update_fine_e                       (E^{n+1/2}_f)
4.  region.snapshot_coarse_e_t_end(parent)         (Q3, post-coarse-H)
5.  fine sub-step k = 2:
     a. region.interpolate_coarse_e_to_fine(0.75)
     b. region.update_fine_h                       (H^{n+3/4}_f confirmed)
     c. region.update_fine_e                       (E^{n+1}_f)
6.  inner.update_e_only                            (E^{n+1}_c, NO Huygens source yet)
7.  region.inject_equivalent_currents_to_coarse(parent)
     — adds (J_S × dt / ε_0) to coarse E_t on the Huygens surface
     — adds (M_S × dt / μ_0) to coarse H_t on the Huygens surface
8.  inner.apply_cpml_h, inner.apply_cpml_e         (outer-face PML)
9.  inner.advance_clock
```

Stage 7 is the only new closure call. Stages 1–6 reuse Q1/Q3/Q4.1 helpers verbatim; stages 8–9 reuse the existing CPML / clock pipeline.

### Discretisation of `n̂ × E` and `n̂ × H` at corners

Each of the six fine-box faces enumerates a set of Yee edges (for `E_t`) and faces (for `H_t`) on its interior outer layer. Edges shared by two adjacent fine-box faces (the 12 edges of the cuboid) and corners shared by three (the 8 corners) are handled by **enumerating contribution per face independently and summing into the coarse grid**. The equivalence-principle integrand has no corner singularity at second order, so the per-face sum is the discrete analogue of the closed-surface integral. (Berenger 2006 §III; see also Taflove & Hagness §13.7 for the equivalent TF/SF corner treatment.)

## Why this retires the instability

The spec `2026-05-18` Q4 closure was **bidirectional**: fine → coarse `H_t` overwrite *and* fine → coarse `E_t` overwrite *and* coarse → fine `E_t` Dirichlet — three field-level coupling channels sharing the same interface plane. The plane was over-determined, and VVVVVV's diagnosis showed the consequence: `overwrite_coarse_e_from_fine` reads a Dirichlet value (Q3 interpolation result, not a Maxwell-evolved fine `E`) and round-trips it back to coarse, while coarse `update_e` has already advanced using the stale pre-closure coarse `H`. No discrete-energy-balanced fix exists at the closure layer alone.

Berenger 2006's scheme is **one-directional per face per step**:

- Coarse → fine is still Dirichlet (Q3 interpolation, unchanged).
- Fine → coarse is **not a copy** — it is an equivalent-current source term added to the coarse grid's `update_e` / `update_h` right-hand side. The fine grid's storage is read but never written by the coarse grid; the coarse grid's storage is updated by a known source term, not overwritten.

The interface plane is no longer over-determined. The fine grid evolves Maxwell's equations with its Q3-Dirichlet boundary; the coarse grid evolves Maxwell's equations with a Huygens-surface source. The two are coupled by the equivalence principle, which guarantees consistency to the order of the per-grid stencil — second-order in `dx_coarse` (Berenger 2006 §IV).

## Public API delta

### Kept (no change)

- `SubgridRegion::snapshot_coarse_e_t`, `snapshot_coarse_e_t_end` (Q3).
- `SubgridRegion::interpolate_coarse_e_to_fine` (Q3).
- `SubgridRegion::update_fine_h`, `update_fine_e` (Q5).
- `SubgridRegion::snapshot_fine_h_mid_step` (Q4.1) — now repurposed as the time-centered `H_tot` source for the `J_S` current.
- `WalkingSkeletonSolver::update_h_only` / `update_e_only` / `apply_cpml_h` / `apply_cpml_e` / `advance_clock` (Q1).

### Replaced

- `SubgridRegion::average_fine_h_to_coarse(&mut parent)` — **removed from the step pipeline** but retained as `pub fn` for posterity / future diagnostic use (see ADR-0035). Marked `#[doc(hidden)]` and documented as "Phase 2.fdtd.7 Q4 closure, replaced by Berenger 2006 in 7.x — see ADR-0035."
- `SubgridRegion::overwrite_coarse_e_from_fine(&mut parent)` — same treatment.

### Added

```rust
impl SubgridRegion {
    /// Berenger 2006 Huygens-surface closure. Injects equivalent surface
    /// currents `J = +n̂ × H_tot` (sampled from `snapshot_fine_h_mid_step`)
    /// and `M = -n̂ × E_tot` (sampled from the outer-layer fine `E_t` after
    /// `update_fine_e`) onto the six interface faces of the parent coarse
    /// grid.
    ///
    /// One-directional, post-coarse-update: must be called after
    /// `WalkingSkeletonSolver::update_e_only` for the coarse step and before
    /// `apply_cpml_h` / `apply_cpml_e`. The fine grid's storage is read
    /// only; the coarse grid's `E` and `H` arrays are mutated in-place at
    /// the Huygens surface only.
    ///
    /// References: Berenger 2006, *IEEE T-AP* 54(12), §III.
    pub fn inject_equivalent_currents_to_coarse(&self, parent: &mut YeeGrid);
}
```

`SubgriddedSolver::step` is rewritten to the stage list in §3 above. `step_with_gaussian_source_ez` follows the same refactor.

## Validation

| Gate | Test | Tolerance | Status before 7.x | Status after 7.x |
|------|------|-----------|-------------------|------------------|
| **Q5 strict plane-wave traversal** | `subgrid_plane_wave_traversal::strict_05pct_peak_over_500_steps` | `max\|E_z_sub − E_z_ref\| / peak ≤ 0.5%` over first 500 steps, 5 probes downstream of nest | `#[ignore]`'d (closure unstable) | **un-`#[ignore]`'d, passes** |
| **Q6 round-trip energy drift** | `subgrid_energy_balance::round_trip_10000_steps` (new) | `\|W(10000) − W(0)\| / W(0) ≤ 0.5%` | Not implemented (Q6 blocked on Q5) | **passes** |
| **Q7 fdtd-007 Maloney-Smith** | `yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (new) | `f_res ±2%`, `\|S_11(f_res)\| ±1 dB`, 0.3% / 0.3 dB sanity vs uniform-fine | Not implemented (Q7 blocked on Q5) | **passes (hardware-gate if > 30 min)** |

Q5 and Q6 are the spec-level retiring criteria for Phase 2.fdtd.7.x. Q7 is the published-benchmark validation gate per CLAUDE.md §4 ("no solver feature ships without a published-benchmark validation case") and is forward-ported from the original 7.0 plan unchanged.

## Risks / open questions

1. **TF/SF accounting bookkeeping at the Huygens surface.** The total-field / scattered-field convention chosen in §3 puts TF inside the fine box and SF outside; sign errors in the `J = +n̂ × H` / `M = −n̂ × E` formulas, or in the outward-normal direction, are the most likely first-cut failure mode. Mitigation: re-use the sign-convention discipline from Phase 2.fdtd.5 TF/SF (ADR-0026) and write the per-face signs in a single helper table; verify on the Q5 strict gate before Q6.
2. **Equivalent-current discretization at the 12 cuboid edges and 8 corners.** Berenger 2006 §III handles this by **per-face enumeration with no special-casing** (the per-face integrand is regular at second order). Risk: corner-cell contributions double-count if the per-face loops are written naively to include the cuboid's edge cells in two faces simultaneously. Mitigation: define each face's edge / face index range as **half-open in the tangential directions** (`lo ≤ i < hi`, with the cuboid edges assigned to exactly one of the two adjacent faces by a deterministic axis-ordering rule). Verified by a unit test that sums a constant `J_S` over the closed surface and recovers the correct discrete divergence.
3. **Time-centering mismatch between `J_S` and `M_S`.** `J_S = n̂ × H_tot` is sampled at `t = n + 1/2` via Q4.1's `snapshot_fine_h_mid_step`; `M_S = −n̂ × E_tot` is sampled at `t = n + 1` via the post-`update_fine_e` outer fine `E_t`. The coarse `update_e` consumes the `J_S` source at `t = n + 1/2` (correct) and the coarse `update_h` (next step's stage 2) consumes the `M_S` source at `t = n + 1` (correct in Yee staggering). Risk: an off-by-half-step error here is silent on the Q5 short-time gate but surfaces on the Q6 10 000-step energy-balance gate. Mitigation: spec the time levels explicitly in the API doc; the Q6 gate is the canary.
4. **Q7 hardware-gating.** Maloney-Smith was already budgeted at `< 30 min` `--release` in the original Q7 brief; the Berenger closure adds two extra per-step traversals of the six Huygens faces (the J and M injections) which is `O(N²)` work — a few percent of total cost. If Q7 still overruns, hardware-gate behind `#[ignore]` per Phase 1.5 / mom-001 precedent (CLAUDE.md §4).

## Dependencies

No new external dependencies. The implementation reuses `ndarray`, the Q3 snapshot / interpolation surface, the Q4.1 mid-step helper, and the Q1 `WalkingSkeletonSolver` per-stage helpers. The `#![forbid(unsafe_code)]` floor on `yee-fdtd` is preserved.

## Phase numbering

- **Phase 2.fdtd.7.x** — Berenger 2006 Huygens-surface closure replaces spec `2026-05-18` §3 stage-6/7. **This spec.** Retires Q5 strict gate, Q6 energy gate, Q7 Maloney-Smith gate.
- Phase 2.fdtd.7.1 onward — unchanged from spec `2026-05-18` §9 (multi-nest, ADE-in-nest, CPML co-location, higher ratios, RO4003C patch openEMS cross-check).

## Lane

`crates/yee-fdtd/**` only. Modifies `crates/yee-fdtd/src/subgrid.rs` and `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs`; adds `crates/yee-fdtd/tests/subgrid_energy_balance.rs`. Q7 touches `crates/yee-validation/src/lib.rs` and adds `crates/yee-validation/tests/fdtd_007_maloney_smith_slot.rs` (forward-ported unchanged from the original Q7 brief).

## References

- Berenger, J.-P., "A Huygens subgridding for 3-D FDTD", *IEEE Trans. Antennas Propag.* 54(12), 2006, pp. 3797–3804, DOI `10.1109/TAP.2006.886504`. **Primary.**
- Spec being amended: `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md` §6 (risks register, Berenger fallback flag).
- VVVVVV diagnosis: commits `a2abb4c` (merge), `72c825c` (Q4.1 helper), and the `#[ignore]` block at `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs:111`.
- ADR-0035 — decision record for the fallback (companion to this spec).
- Taflove, A., Hagness, S. C., *Computational Electrodynamics: The Finite-Difference Time-Domain Method*, 3rd ed., Artech House 2005, §13.7 (TF/SF corner treatment, analogous to Huygens corner treatment).
- Chevalier, M. W., Luebbers, R. J., Cable, V. P., "FDTD local grid with material traverse", *IEEE Trans. Antennas Propag.* 45(3), 1997, pp. 411–421 — the closure being replaced.
- Maloney, J. G., Smith, G. S., *IEEE Trans. Antennas Propag.* 41(5), 1993, Fig. 9 — fdtd-007 reference (forward-ported).
