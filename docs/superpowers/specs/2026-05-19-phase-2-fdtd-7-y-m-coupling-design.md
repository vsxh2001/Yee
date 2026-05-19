# Phase 2.fdtd.7.y — Berenger M-side closure amendment (compensating-source)

**Status:** Draft
**Owner:** TBD
**Phase:** 2.fdtd.7.y
**Depends on:** Phase 2.fdtd.7.x B1 / B2 / B2.1 / B2.2 (merged through `464c7ba`).
Reuses `SubgridRegion::inject_j_to_coarse_e`, `inject_m_to_coarse_h`, the Q3 coarse →
fine `E_t` interpolation surface, the Q4.1 `snapshot_fine_h_mid_step` helper, and the
B2.1 `snapshot_fine_e_end_of_step` helper.
**Blocks:** Phase 2.fdtd.7.x Q5 strict 0.5%-of-peak plane-wave-traversal gate, the
500-step Berenger traversal canary (`berenger_traversal::berenger_step_propagates_without_divergence_500_steps`),
and the Q6 10 000-step round-trip energy-drift gate — all three currently `#[ignore]`'d.

## Motivation

Phase 2.fdtd.7.x B2.2 (Track OOOOOOO, commit `464c7ba`) shipped coarse-ghost
subtraction on the **J side** of the Berenger equivalent-current closure, taking
`J = +n̂ × (H_TF_fine − H_SF_coarse_ghost)` per Berenger 2006 §III canonical
form. The 100-step canary peak `|E_z|_fine` dropped from 31 V/m (B2 baseline,
un-ghosted) to 2.75 V/m (B2.2, J-ghosted), and the 500-step divergence onset
moved from step ~98 to step ~137 — a strict improvement, but the strict Q5 gate
and the Q6 long-time energy gate remain `#[ignore]`'d.

The natural follow-up — applying the same coarse-ghost subtraction
**symmetrically** to the M source per Berenger 2006 §III pure form
(`M = -n̂ × (E_TF_fine − E_SF_coarse_ghost)`) — was tried and **rejected**.
Track OOOOOOO's report:

> The "symmetric" M ghost candidate (coarse E at the surface plane) is
> Dirichlet-tied to the fine grid by Q3 `interpolate_coarse_e_to_fine`, so
> `fine_E_surface − coarse_E_surface ≈ 0` by construction. Enabling M ghost
> subtraction empirically REGRESSED the 100-step canary from 2.75 V/m to
> ≈ 1 kV/m peak `|E_z|_fine`.

This spec amends the Phase 2.fdtd.7.x §3 fine→coarse closure mathematics for
the M side. The J-side ghost subtraction is **kept unchanged** — it is the
load-bearing improvement of 7.x B2.2 and the underlying asymmetry between the
J and M closures is precisely what motivates the amendment.

## Goal

Retire the two strict gates that the Phase 2.fdtd.7.x M-equivalent-current
approach cannot retire on its own:

1. **Q5 strict 0.5%-of-peak plane-wave-traversal gate** (`crates/yee-fdtd/tests/subgrid_plane_wave_traversal.rs:111`,
   `#[ignore]`'d) plus the 500-step Berenger traversal canary
   (`crates/yee-fdtd/tests/berenger_traversal.rs:103`, `#[ignore]`'d).
2. **Q6 round-trip energy-drift gate** (`|W(10 000) − W(0)| / W(0) ≤ 0.5%`).

This is a **closure-amendment-only** scope. The J-side `inject_j_to_coarse_e`
is unchanged. The Q3 coarse→fine Dirichlet interpolation is unchanged. Only
the M sampling time level and the M ghost-subtraction term are amended.

## Background

Berenger 2006 §III gives the canonical Huygens-surface equivalent-current form:

```text
J_S(r, t) = +n̂ × (H_TF(r, t) − H_SF(r, t))     on S
M_S(r, t) = -n̂ × (E_TF(r, t) − E_SF(r, t))     on S
```

In a textbook implementation with two independent grids that both evolve
Maxwell freely, `H_TF − H_SF` and `E_TF − E_SF` are non-zero on `S` because the
two grids diverge after `t = 0` (their `t = 0` initial conditions agree, but
their per-grid CFL clocks and stencils accumulate distinct numerical solutions
of the same continuum field). Berenger's claim is that the **difference** is
exactly the equivalent surface current that re-radiates the TF region's field
into the SF region.

The Yee subgridding pipeline in Phase 2.fdtd.7 violates the "two free grids"
assumption on the **E side only**:

- **H side (J source):** the coarse `H_t` slot on the surface plane is
  updated by the coarse `update_h_only` from coarse `E` curls **without any
  reach into the fine grid**. The fine `H_t` is similarly updated by
  `update_fine_h` from fine `E` curls **without overwriting from coarse**.
  The two are genuinely independent — `H_TF_fine − H_SF_coarse_ghost` is
  non-zero and Berenger's canonical formula works (B2.2 J-side ghost
  subtraction).
- **E side (M source):** the fine `E_t` on the surface plane is **set by
  Q3 Dirichlet interpolation from coarse** before each fine sub-step's
  `update_fine_e`. To first order in `dt_fine`, the fine surface E
  *equals* the time-interpolated coarse surface E. Any difference
  `fine_E_surface − coarse_E_surface` is `O(h_fine²)` from the
  truncation of the Q3 spatial interpolation and `O(dt_fine²)` from the
  fine sub-step's curl correction — both quantities much smaller than
  the leading-order field. Berenger's canonical M formula effectively
  vanishes, and enabling it actively destabilises the closure (the
  small residual is dominated by per-grid round-off noise and amplifies
  through the leapfrog feedback loop).

This asymmetry is the root cause of the Q5 / Q6 gate ignorance.

## Mathematical formulation

### Root cause restated

The Q3 coarse→fine E_t interpolation is a **Dirichlet boundary condition** on
the fine grid: at the start of each fine sub-step, the outer-layer fine `E_t`
is overwritten with a coarse-derived value. After Q3 writes, but *before*
`update_fine_e` runs, the fine surface E literally equals the coarse surface E
to interpolation order. The M source

```text
M_S = -n̂ × (E_TF_fine − E_SF_coarse_ghost)
```

samples both fields at the same physical surface plane. If both are taken at
the **post-Q3, pre-update_fine_e** time level, the difference is zero modulo
interpolation truncation — i.e. the canonical M source is *eliminated by
construction* before `update_fine_e` has any chance to add Maxwell evolution
to the fine E.

The J side does not suffer this because the fine `H_t` is **never
Dirichlet-set** by Q3 — Q3 only writes E_t. The fine `H_t` is updated freely
by `update_fine_h` from fine E curls.

### Candidate fix Option α — drop Q3 Dirichlet

Replace the Q3 coarse→fine `E_t` interpolation with an **absorbing /
Sommerfeld radiation boundary** on the fine grid's outer layer. The fine E
then evolves Maxwell freely, and `fine_E_surface − coarse_E_surface` is
non-zero in general. Berenger's canonical M formula recovers its physical
content.

**Trade-off:** the absorbing / Sommerfeld BC on a Cartesian Yee fine grid
introduces its own reflection at the well-known `~ −40 dB` floor for
first-order Mur or `~ −60 dB` for second-order Mur on plane waves at normal
incidence. That floor becomes the new accuracy ceiling for the coupled
sub-grid system — `0.5%` of peak ≈ `−46 dB`, so first-order Mur is
insufficient; second-order Mur is borderline; a CPML on the fine grid is
overkill and would re-couple to the J side via the fine `H_t` PML auxiliaries.

### Candidate fix Option β — compensating M source (RECOMMENDED)

Keep the Q3 Dirichlet interpolation. Sample the M source at a time level
where the fine E has had a chance to depart from the Q3-Dirichlet value.

Concretely: capture two snapshots of the outer-layer fine `E_t`:

1. **`E_pre`** — the fine `E_t` immediately after Q3 writes the Dirichlet
   value and **before** `update_fine_e` runs. By construction
   `E_pre ≈ coarse_E_interpolated`.
2. **`E_post`** — the fine `E_t` immediately after `update_fine_e`
   completes its sub-step. The fine sub-step's curl-of-H update has added
   a Maxwell-evolved correction whose magnitude is
   `(dt_fine / ε_0) · (∇ × H_fine)_t` evaluated on the surface plane —
   physically non-zero whenever fine H carries the propagating wave.

Define the **compensating M source** as

```text
M_S^compensating = -n̂ × (E_post − E_pre)
```

This is the discrete analogue of the temporal derivative of the fine E on
the surface plane, weighted by `dt_fine`. It captures exactly the
Maxwell-evolved part of the fine E that escapes the Q3 Dirichlet
constraint, and discards the Dirichlet-tied part that would nullify against
the coarse ghost.

`M_S^compensating` is then injected to the coarse `H_t` slot per the
existing `inject_m_to_coarse_h` pipeline. No coarse-ghost subtraction is
applied to this term — the differencing is already done in the
`E_post − E_pre` sample, so re-subtracting a coarse ghost would
double-count (and reintroduces the failure mode that OOOOOOO documented).

**Trade-off:** the compensating-source magnitude is `O(dt_fine)` smaller
than the J source magnitude `O(1)`. The M-side injection then has a much
lower signal-to-round-off ratio than the J side — if the round-off
amplification through the leapfrog feedback loop dominates, Option β
degenerates to the current B2.2 state (M source vanishes, only J ghost
works). This is the dominant residual risk.

### Recommended option

**Option β.** Lower implementation risk: retains the spec §3 stage structure
(Q3 stays, only the M sampling time level changes), preserves the Phase
2.fdtd.7.x B2.2 J-side ghost-subtraction code path, and is reversible if
empirically it degenerates to the current state (the new snapshot becomes
unused).

Option α is held as the **escape hatch** for Step C5 if Option β fails.

## Public API delta

### Kept (no change)

- `SubgridRegion::inject_j_to_coarse_e` — unchanged Phase 2.fdtd.7.x B2.2
  J-side coarse-ghost subtraction.
- `SubgridRegion::interpolate_coarse_e_to_fine`, `snapshot_coarse_e_t`,
  `snapshot_coarse_e_t_end` (Q3) — unchanged.
- `SubgridRegion::snapshot_fine_h_mid_step` (Q4.1) — unchanged J-side feed.
- `SubgridRegion::snapshot_fine_e_end_of_step` (B2.1) — **retained** as the
  `E_post` source, but its semantics broaden (see below).

### Added

```rust
impl SubgridRegion {
    /// Snapshot the fine `E_t` on the outer Huygens layer *before*
    /// `update_fine_e` runs (i.e. immediately after the Q3 Dirichlet
    /// interpolation has been applied). Pair with
    /// [`Self::snapshot_fine_e_post_update`] taken *after* `update_fine_e`
    /// completes; [`Self::inject_m_to_coarse_h`] then computes the
    /// compensating M source `-n̂ × (E_post − E_pre)`.
    ///
    /// Phase 2.fdtd.7.y Option β; see spec
    /// `2026-05-19-phase-2-fdtd-7-y-m-coupling-design.md` §3.
    pub fn snapshot_fine_e_pre_update(&mut self);

    /// Snapshot the fine `E_t` on the outer Huygens layer *after*
    /// sub-step 2's `update_fine_e` completes. Replaces / supersedes
    /// [`Self::snapshot_fine_e_end_of_step`] as the `E_post` source.
    pub fn snapshot_fine_e_post_update(&mut self);
}
```

### Modified

- `SubgridRegion::inject_m_to_coarse_h` — switches the M source from
  `M = -n̂ × E_fine_end_of_step` (Phase 2.fdtd.7.x B2.1 form) to
  `M = -n̂ × (E_post − E_pre)` (Phase 2.fdtd.7.y Option β form). No
  coarse-ghost subtraction is applied to the compensating source. The
  per-face inner helper `inject_*_face` gains a single boolean
  `use_compensating_source` flag that selects which of the two cached
  fine-E arrays to read.
- `SubgriddedSolver::step` and `step_with_gaussian_source_ez` —
  rewired to call `snapshot_fine_e_pre_update` immediately after the
  fine sub-step 2's `interpolate_coarse_e_to_fine` and before its
  `update_fine_e`, and `snapshot_fine_e_post_update` immediately
  after that `update_fine_e` completes.

## Validation

| Gate | Test | Tolerance | Status before 7.y | Status after 7.y |
|------|------|-----------|-------------------|------------------|
| **Berenger 500-step canary** | `berenger_traversal::berenger_step_propagates_without_divergence_500_steps` | `peak \|E_z\|_fine < 1e3` V/m over 500 steps | `#[ignore]`'d (peak ≈ 1.035e3 at step ~137 after B2.2) | **un-`#[ignore]`'d, passes** |
| **Q5 strict plane-wave traversal** | `subgrid_plane_wave_traversal::strict_05pct_peak_over_500_steps` | `max\|E_z_sub − E_z_ref\| / peak ≤ 0.5%` over 500 steps, 5 probes downstream | `#[ignore]`'d | **un-`#[ignore]`'d, passes** |
| **Q6 round-trip energy-drift** | `subgrid_energy_balance::round_trip_10000_steps` | `\|W(10 000) − W(0)\| / W(0) ≤ 0.5%` | `#[ignore]`'d (blocked on Q5) | **un-`#[ignore]`'d, passes** |

Phase 2.fdtd.7.x's published-benchmark gate (fdtd-007 Maloney-Smith Fig. 9)
remains under Phase 2.fdtd.7.x scope; Phase 2.fdtd.7.y is an amendment to the
closure mathematics only and does not need its own published-benchmark gate.

## Risks / open questions

1. **Round-off noise floor.** Option β's compensating source is
   `O(dt_fine)` in magnitude relative to a hypothetical canonical-form M
   source. On a 64-bit float Yee grid the round-off floor is `~ 1e-15`
   relative; `dt_fine ~ 1e-12` s on the 7.x canonical geometry, so the
   ratio is `~ 1000:1` — comfortable, but not infinite. **Surfaces in
   Step C3** as the strict Q5 gate either passing or remaining at the
   `~ 1%` level (compensating source too small to bridge). Mitigation: if
   C3 sees `0.5% < rel_err < 5%`, the compensating source is the right
   sign and order of magnitude but undersized; pivot to Option α.
2. **Compensating-source degeneration.** If `E_post − E_pre` is
   numerically indistinguishable from zero (Q3 + `update_fine_e` lock
   the surface E to the coarse value bit-exact within a sub-step),
   Option β literally reduces to the current B2.2 state. **Surfaces in
   Step C3** as the strict Q5 gate failing at exactly the same level as
   the current B2.2 (peak `|E_z| ≈ 1e3` V/m at step ~137). Mitigation:
   Step C5 — switch to Option α.
3. **Inadvertent double-counting.** If a future maintainer adds
   coarse-ghost subtraction to `inject_m_to_coarse_h` while Option β is
   live, the closure regresses to the failure mode OOOOOOO documented
   (the compensating source already does the differencing). **Mitigation:**
   `inject_m_to_coarse_h`'s docstring explicitly forbids coarse-ghost
   subtraction in the Option β regime; the per-face helper carries a
   `debug_assert!(!do_ghost || !use_compensating_source)` guard.

## Dependencies

No new external dependencies. Reuses `ndarray` array clones for the
pre/post snapshots. The `#![forbid(unsafe_code)]` floor on `yee-fdtd` is
preserved.

## Phase numbering

- **Phase 2.fdtd.7.y** — Option β compensating-source M closure
  amendment to Phase 2.fdtd.7.x. **This spec.** Retires the 500-step
  Berenger canary, the Q5 strict 0.5% gate, and the Q6 energy-drift
  gate.
- Phase 2.fdtd.7.y.α (escape hatch) — Option α absorbing fine BC.
  Conditional on Option β failing Step C3.
- Phase 2.fdtd.7.1 onward — unchanged from Phase 2.fdtd.7.x §9.

## Lane

`crates/yee-fdtd/**` only. Modifies `crates/yee-fdtd/src/subgrid.rs` and
the two `#[ignore]`'d traversal test files. Adds the Q6 energy-balance
test only if not already created under Phase 2.fdtd.7.x B4.

## References

- Berenger, J.-P., "A Huygens subgridding for 3-D FDTD", *IEEE Trans.
  Antennas Propag.* 54(12), 2006, pp. 3797–3804, DOI
  `10.1109/TAP.2006.886504`. **Primary**, §III canonical equivalent-current
  form.
- Spec being amended:
  `docs/superpowers/specs/2026-05-19-phase-2-fdtd-7-x-berenger-huygens-design.md`
  §3 fine→coarse closure mathematics.
- ADR-0035 — Berenger Huygens-surface subgridding decision record for
  Phase 2.fdtd.7.x.
- ADR-0038 — Berenger M-side closure compensating-source amendment
  (companion to this spec).
- Track OOOOOOO empirical regression: commit `464c7ba`, J-side ghost
  subtraction merged, M-side ghost tried and rejected. See the
  `#[ignore]` block at
  `crates/yee-fdtd/tests/berenger_traversal.rs:103-107` for the
  pinned diagnosis.
- Taflove, A., Hagness, S. C., *Computational Electrodynamics: The
  Finite-Difference Time-Domain Method*, 3rd ed., Artech House 2005,
  §7.6 (Mur absorbing BC on Cartesian grids), §13.7 (TF/SF analogue of
  Huygens corner treatment).
