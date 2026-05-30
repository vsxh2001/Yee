# Filter Phase F2.4 — tolerance / yield (Monte-Carlo) — Design Spec

**ADR:** ADR-0113 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

The lumped-LC goal names **tolerance consideration**. Real parts have tolerances
(E24 ±5 %, E96 ±1 %); the as-built filter's spec compliance is a *distribution*.
F2.4 = **Monte-Carlo yield**: sample each component within its tolerance, rebuild
the LC ladder, evaluate the realized response against the spec mask, and report
the **fraction that passes** (yield) + worst-case margins. Pure-math, WASM-safe,
no FDTD. Builds only on F2.0 (`LumpedLadder` + `ladder_s21`) and F2.1
(`ESeries::nearest` + per-series tolerance).

## Method

1. **Realize:** for each `LcResonator`, snap `L`,`C` to the chosen `ESeries`
   value (`F2.1 nearest`).
2. **Sample (M trials, seeded):** perturb each realized value by a random factor
   in `[1−tol, 1+tol]` (uniform; `tol` = series tolerance). Use a tiny in-module
   **seeded PRNG** (SplitMix64 / xorshift — NO `rand` dependency, keeps WASM-safe +
   dep-free + reproducible).
3. **Evaluate:** rebuild the perturbed `LumpedLadder`, compute `ladder_s21` over
   the band, apply the **same mask verdict** as `filt`/`lumped` (passband ripple,
   in-band RL, stopband rejection).
4. **Aggregate:** `yield_fraction = passes/M`; also worst in-band RL and worst
   stopband rejection across trials.

## Changes (`crates/yee-filter/**` ONLY)

- New `crates/yee-filter/src/tolerance.rs`:
  - `pub struct YieldResult { yield_fraction: f64, n_samples: usize,
    worst_inband_rl_db: f64, worst_stopband_rej_db: f64 }` (+ serde).
  - `pub fn monte_carlo_yield(ladder: &LumpedLadder, series: ESeries,
    mask: &SpecMask, n_samples: usize, seed: u64) -> YieldResult`.
  - private seeded PRNG + the per-sample perturb/evaluate.
- Re-export from `lib.rs`. (Reuse `ladder_s21` (doc-hidden pub) + the mask-verdict
  helper; if the verdict logic lives only in yee-cli, lift a small shared
  `mask_pass(s21_db_at(f)…)` into yee-filter — minimal, in-lane.)

## DoD (machine-checkable; pure-math, NO FDTD)

1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 — incl. gate `yield_001` (cheb N=5 fixture):
   - **Determinism:** same `seed` → identical `yield_fraction` (reproducible).
   - **Zero-tolerance sanity:** a hypothetical 0 % tolerance (or n_samples with a
     degenerate ±0 factor) yields **1.0** (the nominal realized ladder passes — or
     document if the E24-quantized nominal already fails, in which case assert the
     nominal-pass separately).
   - **Monotonicity invariant (robust, the key check):** `yield(E96) ≥ yield(E24)`
     for the same seed + M — tighter parts never reduce yield. (Avoids brittle
     exact-number assertions.)
   - `yield_fraction ∈ [0,1]`; `n_samples` honored; M ≥ 200 runs sub-second.

## Out of scope

Correlated/√-N statistics, sensitivity ranking per component, sigma/Cpk,
non-uniform (gaussian) part distributions (uniform only for the skeleton),
parasitics (F2.1b), FDTD-based yield (would use F2.3), UI. A `rand`-crate
dependency (use the in-module PRNG).

## Why now

Pure-math, light, dispatchable immediately (depends only on shipped F2.0+F2.1, not
on the heavier F2.2 PCB / F2.3 FDTD), and directly delivers the goal's named
"tolerance consideration." The monotonicity gate is robust (no brittle magic
numbers) yet physically meaningful.
