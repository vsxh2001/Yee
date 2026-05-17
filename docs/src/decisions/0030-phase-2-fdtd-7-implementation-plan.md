# ADR-0030: Phase 2.fdtd.7 subgridding implementation plan

## Status

Accepted — 2026-05-18 (plan only; track execution deferred to
follow-up agents).

## Context

ADR-0027 locked the Phase 2.fdtd.7.0 scope: single 2× axis-aligned
nest, time-subcycled, EM-only, no co-location with CPML / TF-SF,
fdtd-007 Maloney-Smith gate plus a 10 000-step round-trip energy-
drift gate. The next step is the task ladder for ~1 580 LOC of new
code (`SubgridRegion`, `SubgriddedSolver`, six-face interpolation
/ area-averaging, two gates) without churning the existing
`crates/yee-fdtd/tests/` regressions. Track PPPPP (merge `f217df0`)
lands the plan.

## Decision

The plan splits into **one pre-flight refactor plus seven
sequential / parallel tracks**:

- **Step 1 (Q1) — `WalkingSkeletonSolver::step` refactor into
  composable helpers.** The spec called this "quality-of-life,
  not blocking." The plan disagrees: every fine sub-step needs
  `update_h_only` / `update_e_only` *without* the CPML / PEC /
  clock-advance side effects bundled into the existing `step`
  family. Bundling refactor + new feature would make the diff
  hostile to review. Q1 lands the refactor alone, byte-identical
  against every existing `yee-fdtd/tests/` regression.
- **Q1 → Q2 → Q3 → Q4 → Q5 sequential critical path.** Q2:
  `SubgridRegion` + `SubgriddedSolver` scaffold. Q3: coarse →
  fine `E_t` interpolation. Q4: fine → coarse area-averaging
  closure (Chevalier 1997 §IV). Q5: seven-stage time-subcycling
  + plane-wave-traversal integration test. Q3 / Q4 are siblings
  serialised to keep `subgrid.rs` diffs reviewable (CLAUDE.md §5).
- **Q6 ‖ Q7 parallel post-Q5.** Q6: 10 000-step round-trip
  energy-drift stability gate (≤ 0.5%). Q7: fdtd-007 Maloney-
  Smith production gate (resonance ±2%, `|S_11|` ±1 dB, plus
  0.3% / 0.3 dB internal sanity check).
- **Validation rollup with explicit run-time budgets.** Q6 `< 5
  min` `--release`; Q7 `< 30 min`, hardware-gated `#[ignore]`
  if it overruns (Phase 1.5 / mom-001 precedent). Maloney-Smith
  Fig. 9 values are hand-digitised with an escape hatch (±5% by
  eye, `// TBD`, do not invent reference numbers).

Critical path: five serial merges + one parallel pair, within
CLAUDE.md §5's "up to 5 parallel agents." Lane: `crates/yee-fdtd/
**` plus targeted touches on `crates/yee-validation/`. The
`#![forbid(unsafe_code)]` floor is preserved throughout.

## Consequences

- **Pre-flight refactor is the right call.** Q1 alone is ~80 LOC
  of mechanical extraction against byte-identical regression;
  bundled into Q2 it would have been ~260 LOC mixing behaviour-
  preserving refactor with behaviour-changing new feature.
- **Late-time instability surfaces in Q6, not Q7.** The
  Chevalier area-average is the closure; asymmetric coupling
  can still grow over `O(10⁴)` steps. Q6 records drift between
  0.5% and 5% as `// regression-tracked` and surfaces as a
  finding — the 0.5% gate is not weakened without a spec
  amendment.
- **fdtd-007 reference is hand-digitised.** Q7's escape hatch
  documents the ±5% precision with a `// TBD`.
- **CLI / Python / GUI exposure deferred to 7.0.1.**

## References

- `docs/superpowers/plans/2026-05-18-phase-2-fdtd-7-subgridding.md`
- `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`
- Track PPPPP merge commit `f217df0`.
- ADR-0027 — Phase 2.fdtd.7.0 scope lock (this plan's parent).
- J.-P. Berenger, *IEEE Trans. Antennas Propag.* 54(12), 2006,
  §IV — asymmetric-coupling stability analysis.
- J. G. Maloney, G. S. Smith, *IEEE Trans. Antennas Propag.*
  41(5), 1993, Fig. 9 — fdtd-007 reference.
- CLAUDE.md §3, §4, §5.
