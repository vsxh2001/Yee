# Filter Phase F2.4 — tolerance / yield (Monte-Carlo) — Plan

**Spec:** `2026-05-30-f2-4-tolerance-yield-design.md` · **ADR:** ADR-0113

## Lane
`crates/yee-filter/**` ONLY (new `src/tolerance.rs`, `lib.rs` re-export, `tests/`).
Consume F2.0 `LumpedLadder`/`ladder_s21` + F2.1 `ESeries`. Out of lane → finding.
WASM-safe: pure `f64` + an in-module seeded PRNG (NO `rand` dep) + serde.

## Base
New worktree off current `main` (re-fetch first). Branch
`feature/filter-f2-4-tolerance`.

## Pattern files
- `crates/yee-filter/src/lumped.rs` — `LumpedLadder`/`LcResonator`/`ladder_s21`
  (`#[doc(hidden)] pub`); reuse `ladder_s21` for the per-sample response. Mirror
  module-doc/serde/re-export style.
- `crates/yee-filter/src/parts.rs` — `ESeries::nearest`/`tolerance_pct` for the
  realize step.
- `crates/yee-cli/src/filter.rs` `check_mask` (+ `crates/yee-filter/tests/
  filt_001_mask.rs` / `lumped_001.rs`) — the mask-verdict logic to reuse; if it
  only lives in yee-cli, lift a tiny `mask_verdict` into yee-filter (in-lane).
- `tests/lumped_001.rs` — the cheb_bpf fixture setup; clone for `yield_001`.

## Steps
1. `src/tolerance.rs`: seeded PRNG (SplitMix64: `state = state.wrapping_add(0x9E3779B97F4A7C15); z=...; → f64 in [0,1)`), `YieldResult` (+serde), `monte_carlo_yield(ladder, series, mask, n_samples, seed)`: realize (nearest each L/C) → per sample perturb ×(1±tol·(2u−1)) → rebuild ladder → `ladder_s21` sweep → mask verdict → tally yield + worst RL/rejection.
2. `lib.rs`: re-export `YieldResult, monte_carlo_yield` (+ the lifted `mask_verdict` if added).
3. `tests/yield_001.rs`: determinism (same seed → same yield), zero-tol/nominal sanity, **`yield(E96) ≥ yield(E24)`** (same seed+M), yield∈[0,1], M=200+.

## Verify (exit 0; pure-math, fast — NO FDTD)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
```

## Escape hatch
Blocked > 25 min (the mask-verdict logic isn't reachable/liftable cleanly; the
E24-quantized nominal already fails the mask so the "zero-tol→1.0" sanity needs
rethinking; PRNG determinism issues) → STOP + surface the nominal realized verdict
+ the E24/E96 yields. Do NOT add `rand`; do NOT touch F2.0/F2.1 logic; do NOT
weaken the monotonicity gate.

## Done when
DoD 1–3 pass; `yield_001` green (incl. E96 ≥ E24 monotonicity); `git diff --stat
<base>..HEAD` = only `crates/yee-filter/**`; F2.0/F2.1 logic untouched; WASM-safe.
