# ADR-0113: Filter Phase F2.4 — tolerance / yield (Monte-Carlo)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0111 (F2.0 LC ladder), ADR-0112 (F2.1 BOM + E-series tolerances),
the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

The lumped-LC goal names **tolerance consideration**. F2.1 records each part's
series tolerance (E24 ±5 %, E96 ±1 %); spec compliance of the as-built filter is a
distribution, not a point. Nothing yet quantifies it.

## Decision

Add `crates/yee-filter/src/tolerance.rs`: `monte_carlo_yield(&LumpedLadder,
ESeries, &SpecMask, n_samples, seed) -> YieldResult`. It snaps each L/C to the
chosen E-series value (F2.1), then for M seeded trials perturbs every value within
`±tolerance` (uniform), rebuilds the ladder, evaluates `ladder_s21` (F2.0) against
the spec mask, and reports `yield_fraction` + worst-case in-band RL / stopband
rejection. A tiny in-module seeded PRNG (SplitMix64) keeps it dep-free, WASM-safe,
and reproducible — no `rand` crate.

Gate `yield_001` (cheb N=5): determinism (same seed → same yield); nominal-realized
sanity; and the robust invariant **`yield(E96) ≥ yield(E24)`** (tighter parts never
reduce yield) — avoiding brittle magic-number assertions. `yield_fraction ∈ [0,1]`.

## Consequences

**Ships:** yield/tolerance analysis for a lumped filter — the goal's named
"tolerance consideration." Feeds the UI (a yield readout + part-tolerance
sensitivity) and informs E24-vs-E96 part choice.

**Gate:** `cargo test -p yee-filter` green incl. `yield_001`. Pure-math,
sub-second.

**Not in scope:** sensitivity ranking, sigma/Cpk, gaussian part distributions,
parasitics (F2.1b), FDTD-based yield (F2.3), UI, a `rand` dependency.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-4-tolerance-yield-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-4-tolerance-yield.md`.
