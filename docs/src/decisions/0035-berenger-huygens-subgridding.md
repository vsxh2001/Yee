# ADR-0035: Berenger 2006 Huygens-surface subgridding replaces Okoniewski direct-copy closure

## Status

Accepted — 2026-05-19 (spec + plan only; implementation deferred
to track ZZZZZZ — see companion files).

## Context

Phase 2.fdtd.7.0 (ADR-0027, ADR-0030) shipped a Chevalier 1997 /
Okoniewski 1997 style direct-copy closure for the fine↔coarse
subgridding interface: linear coarse → fine `E_t` interpolation
(Q3) plus area-averaged `H_t` overwrite + edge-averaged `E_t`
overwrite (Q4) for the fine → coarse direction. The Q5
plane-wave-traversal integration test failed the strict
0.5%-of-peak agreement gate over 500 steps and was `#[ignore]`'d
pending diagnosis (commit `fab286a`, merge `426a36c`).

Track VVVVVV (commit `72c825c`, merge `a2abb4c`) added the
`snapshot_fine_h_mid_step` time-centering helper (Q4.1) on the
hypothesis that the failure was a half-step phase bug in the
fine → coarse `H` average. With Q4.1 the rel-err improved but
did not clear 0.5%. VVVVVV's report diagnosed the residual as
fundamental:

> The `overwrite_coarse_e_from_fine` reads fine boundary `E_t`
> cells that are Dirichlet-set by interpolation (not updated by
> `update_fine_e`), round-tripping an interpolation result rather
> than a fresh fine evolution; the H closure simultaneously
> writes the time-averaged fine `H` back to coarse, but coarse
> `update_e` already used the stale pre-closure coarse `H` to
> update coarse `E_t` at the interface plane.

The interface plane is **over-determined** by the bidirectional
closure: three field-level coupling channels (coarse → fine `E_t`
Dirichlet, fine → coarse `H_t` overwrite, fine → coarse `E_t`
overwrite) share the same Yee edges/faces, and no
discrete-energy-balanced fix exists at the closure layer alone.

Spec `2026-05-18-phase-2-fdtd-7-subgridding-design.md` §6 risks
register explicitly flagged Berenger 2006's Huygens-surface
scheme as the documented fallback for this exact failure mode.

## Decision

Replace the spec `2026-05-18` §3 stage-6 + stage-7 bidirectional
direct-copy closure with the **Berenger 2006 Huygens-surface
scheme** (Berenger, *IEEE T-AP* 54(12), 2006, pp. 3797–3804,
DOI `10.1109/TAP.2006.886504`):

1. **Coarse → fine direction unchanged.** Q3 linear spatial +
   temporal interpolation of `E_t` continues to supply the fine
   grid's outer Dirichlet boundary. This direction has never
   been an instability source.
2. **Fine → coarse direction becomes equivalent-current
   re-radiation.** Compute `J = +n̂ × H_tot` and
   `M = −n̂ × E_tot` on the six Huygens faces (TF inside the
   fine box, SF outside — convention documented in the spec).
   Inject as RHS source terms to the coarse `update_e` and
   `update_h` after the coarse field-update stages, before
   CPML. The fine grid's storage is read but never written
   by the coarse grid; the coarse grid's storage is updated
   by a source term, not overwritten.
3. **No closure round-trip; interface plane is no longer
   over-determined.** The fine grid evolves Maxwell with a
   Q3-Dirichlet boundary; the coarse grid evolves Maxwell
   with a Huygens-surface source; the equivalence principle
   guarantees second-order consistency.
4. **Q4.1 `snapshot_fine_h_mid_step` helper is repurposed**
   as the time-centered `H_tot` source for `J_S`. Q1 step
   helpers and Q3 snapshot/interpolation surface unchanged.
5. **Old Q4 helpers retained `#[doc(hidden)]`** —
   `average_fine_h_to_coarse` and `overwrite_coarse_e_from_fine`
   stay in `subgrid.rs` for posterity and potential future
   diagnostic use, but are no longer wired into the step
   pipeline. Their removal is a separate spec amendment.

Validation gates retired by this fallback: Q5 strict 0.5%-of-peak
plane-wave traversal (un-`#[ignore]`'d), new Q6 10 000-step
round-trip energy-drift gate (`≤ 0.5%`), and the original Q7
fdtd-007 Maloney-Smith production gate (forward-ported unchanged).

## Consequences

- **Phase 2.fdtd.7.x replaces Q4–Q7 of the original 7.0 plan
  in the step ladder.** Q4 (`average_fine_h_to_coarse` etc.)
  and Q4.1 (`snapshot_fine_h_mid_step`) **remain in code**
  because the Q4.1 helper is reused as the `J_S` source; only
  the call site changes. Q5, Q6, Q7 become B3, B4, B5 in the
  new plan with identical tolerances.
- **One additional public function** on `SubgridRegion`:
  `inject_equivalent_currents_to_coarse`. No new dependencies;
  `#![forbid(unsafe_code)]` floor preserved.
- **TF/SF sign-convention discipline carries over from
  Phase 2.fdtd.5** (ADR-0026). The Huygens-surface accounting
  is structurally identical to the TF/SF slab; mistakes
  surface in B3's short-time Q5 gate before reaching B4's
  long-time energy gate.
- **Maloney-Smith Fig. 9 reference is still hand-digitised**
  per ADR-0030 — that escape hatch is unchanged.
- **CLI / Python / GUI exposure stays deferred** to
  Phase 2.fdtd.7.0.1.
- **Old Q4 closure stays compileable but `#[doc(hidden)]`**
  to keep the diff against `main` reviewable and to leave a
  rollback option open if the Berenger pipeline surfaces a
  new failure mode.

## References

- `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
- `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-x-berenger-huygens.md`
- Spec being amended: `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md` §6
- ADR-0027 — Phase 2.fdtd.7.0 scope (this ADR's parent)
- ADR-0030 — Phase 2.fdtd.7 implementation plan (this ADR's
  parent)
- VVVVVV diagnosis: commit `72c825c`, merge `a2abb4c`,
  `#[ignore]` block at
  `crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs:111`
- J.-P. Berenger, "A Huygens subgridding for 3-D FDTD",
  *IEEE Trans. Antennas Propag.* 54(12), 2006, pp. 3797–3804,
  DOI `10.1109/TAP.2006.886504`
- CLAUDE.md §3, §4, §5
