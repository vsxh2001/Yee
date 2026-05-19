# ADR-0038: Berenger M-side closure — compensating-source amendment to Phase 2.fdtd.7.x

## Status

Accepted — 2026-05-19 (spec + plan only; implementation deferred to a
follow-on track — see companion files).

## Context

Phase 2.fdtd.7.x (ADR-0035) shipped the Berenger 2006 Huygens-surface
subgridding closure to replace the Phase 2.fdtd.7.0 Okoniewski direct-copy
closure. The 7.x implementation landed in four steps:

- **B1** (Track EEEEEEE) — `inject_equivalent_currents_to_coarse` skeleton.
- **B2** (Track FFFFFFF) — `SubgriddedSolver::step` Berenger pipeline rewrite.
- **B2.1** (Track LLLLLLL) — split injection: `inject_j_to_coarse_e` +
  `inject_m_to_coarse_h` with the time-centring fix
  (`snapshot_fine_e_end_of_step` defers M by one coarse step).
- **B2.2** (Track OOOOOOO, commit `464c7ba`) — coarse-ghost subtraction
  on the **J side only**: `J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)`,
  per Berenger 2006 §III canonical equivalent-current form.

B2.2's 100-step canary peak `|E_z|_fine` dropped from 31 V/m (B2 baseline,
un-ghosted) to 2.75 V/m (J-ghosted). The 500-step divergence onset moved
from step ~98 to step ~137 — strict improvement, but the strict Q5 gate
(0.5% peak agreement, 500 steps), the 500-step Berenger canary, and the
Q6 long-time energy gate (10 000 steps) remain `#[ignore]`'d.

The expected follow-up — applying coarse-ghost subtraction
**symmetrically** to the M source per Berenger 2006 §III pure form
(`M = -n̂ × (E_TF_fine − E_SF_coarse_ghost)`) — was tried by OOOOOOO and
**rejected**:

> The "symmetric" M ghost candidate (coarse E at the surface plane) is
> Dirichlet-tied to the fine grid by Q3 `interpolate_coarse_e_to_fine`,
> so `fine_E_surface − coarse_E_surface ≈ 0` by construction. Enabling
> M ghost subtraction empirically REGRESSED the 100-step canary from
> 2.75 V/m to ≈ 1 kV/m peak `|E_z|_fine`.

Root cause: Phase 2.fdtd.7's Q3 coarse→fine `E_t` interpolation is a
**Dirichlet boundary condition** on the fine grid. To first order in
`dt_fine`, the post-Q3 / pre-`update_fine_e` fine surface E *equals* the
time-interpolated coarse surface E — the M canonical-form difference
vanishes by construction. The J side does not suffer this because Q3
only writes E_t; the fine `H_t` is updated freely by `update_fine_h`
from fine E curls, so `H_TF_fine − H_SF_coarse_ghost` is genuinely
non-zero.

## Decision

Amend the Phase 2.fdtd.7.x spec §3 fine→coarse closure mathematics for
the M side. The J-side coarse-ghost subtraction (B2.2) is **kept
unchanged** — it is the load-bearing improvement of 7.x and is
preserved.

Recommended option (per the companion spec
`docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md`
§3): **Option β — compensating M source**.

1. **Capture two fine-`E_t` snapshots per coarse step**:
   `E_pre` = fine outer-layer E_t after Q3 Dirichlet writes, before
   `update_fine_e`; `E_post` = same E_t after sub-step 2's
   `update_fine_e` completes.
2. **Compute the compensating M source** as
   `M = -n̂ × (E_post − E_pre)`. The difference is the Maxwell-evolved
   part of fine E that escapes the Q3 Dirichlet tie; the
   Dirichlet-tied part nullifies in the subtraction.
3. **Inject `M` to coarse `H_t` via the existing
   `inject_m_to_coarse_h` pipeline.** No coarse-ghost subtraction is
   applied to the compensating source (the differencing is already
   done in the `E_post − E_pre` sample; re-subtracting a coarse ghost
   would double-count and regenerate the OOOOOOO failure mode). The
   per-face helper carries a debug-assert guard against the
   combination.

Held as **escape hatch** if Option β fails the strict gates: **Option α
— drop the Q3 Dirichlet** and replace it with a second-order Mur
absorbing BC on the fine grid's outer `E_t` layer. Restores Berenger's
canonical M source at the cost of a Mur-floor accuracy ceiling
(~−60 dB on plane waves at normal incidence; the strict Q5 0.5%
tolerance may need re-spec to 1% if Option α is taken).

## Consequences

- **Phase 2.fdtd.7.x J-side ghost subtraction (B2.2) is preserved
  unchanged.** Only the M side changes.
- **Q5 strict, 500-step canary, and Q6 energy-drift gates retire
  *conditionally*** — pending Steps C1–C4 of the companion plan
  passing. If Option β degenerates to the current B2.2 state, Step
  C5 pivots to Option α with a spec amendment.
- **One additional public function on `SubgridRegion`:**
  `snapshot_fine_e_pre_update`. One existing public function
  (`snapshot_fine_e_post_update`, currently named
  `snapshot_fine_e_end_of_step`) gains broadened semantics. No new
  dependencies; `#![forbid(unsafe_code)]` floor preserved.
- **B2.1 `snapshot_fine_e_end_of_step` retained `#[doc(hidden)]`** for
  one release cycle as a rollback option. Its removal is a separate
  spec amendment.
- **Existing 100-step `berenger_step_propagates_without_divergence`
  canary stays green** at peak `|E_z|_fine ≤ 2.75` V/m — Step C2's
  acceptance criterion is "no short-time regression vs B2.2".
- **TF/SF sign-convention discipline carries over** from Phase 2.fdtd.5
  (ADR-0026) and Phase 2.fdtd.7.x (ADR-0035).

## References

- `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md`
- `docs/superpowers/plans/2026-05-19-phase-2-fdtd-7-y-m-coupling.md`
- ADR-0035 — Berenger Huygens-surface subgridding (Phase 2.fdtd.7.x
  parent decision).
- ADR-0027, ADR-0030 — Phase 2.fdtd.7.0 scope and implementation plan
  (grandparent decisions).
- OOOOOOO M-side ghost regression: commit `464c7ba` and the pinned
  diagnosis at `crates/yee-fdtd/tests/berenger_traversal.rs:103-107`
  (B2.2 500-step canary `#[ignore]` block).
- J.-P. Berenger, "A Huygens subgridding for 3-D FDTD",
  *IEEE Trans. Antennas Propag.* 54(12), 2006, pp. 3797–3804,
  DOI `10.1109/TAP.2006.886504` — §III canonical equivalent-current
  form.
- CLAUDE.md §3, §4, §5.
