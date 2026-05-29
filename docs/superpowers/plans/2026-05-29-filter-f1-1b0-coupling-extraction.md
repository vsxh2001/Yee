# Filter Phase F1.1b.0 — coupling/Qe extraction — Implementation Plan

**Spec:** `2026-05-29-filter-f1-1b0-coupling-extraction-design.md` · **ADR:** ADR-0093

## Lane
`crates/yee-filter/**` ONLY. Out of lane (any other crate) → finding, not fix.

## Base
Worktree `worktrees/extract`, branch `feature/filter-f1-1b0-extract`, base `9aa79d7`.

## Pattern files
- `crates/yee-filter/src/lib.rs` — house style, `check_mask`/`ideal_response`
  signatures, the `MaskReport` doc idiom; add the `extract` items in the same file
  or a `mod extract` re-exported at the crate root.
- `crates/yee-fdtd` `cavity_q.rs` decay-fit (the log-linear τ fit to mirror) —
  read for the method only; do NOT depend on yee-fdtd.

## Steps
1. `crates/yee-filter/src/lib.rs` (or `src/extract.rs` + `pub mod extract;` +
   `pub use extract::*;`): add `CouplingExtraction`, `extract_coupling`,
   `extract_q_ringdown` per spec. Pure f64; no new dep. Doc every public item.
2. Tests (`tests/extract.rs` or `#[cfg(test)]`): `extract_001_coupling` +
   `extract_002_q_ringdown` per spec §DoD 4–5, including the two negative controls.

## Math notes
- Split freqs: `f_lo = f0/√(1+k)`, `f_hi = f0/√(1−k)` (Hong-Lancaster §8) — use
  these to BUILD the synthetic signal, then `extract_coupling` inverts to
  `k = (f_hi²−f_lo²)/(f_hi²+f_lo²)`; the round-trip recovers `k_true`.
- Lorentzian: `1 / (1 + (2·Q·(f−fc)/fc)²)` per peak; pick Q so the two peaks are
  clearly separated at the chosen `k` (e.g. Q≈100 with k=0.04).
- Decay fit: ordinary least squares on `(t_k, ln|env_k|)`; `τ = −1/slope`.
  `Q = π f0 τ`.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
```
Pure math — sub-second. Do NOT run `cargo test --workspace`.

## Escape hatch
Blocked >15 min — peak-finder doesn't resolve the two Lorentzians at the chosen
k/Q (widen separation or sweep resolution in the TEST, not by loosening the gate
below 1e-2), or the decay fit misses τ → STOP, surface computed-vs-expected.
Do NOT add a dependency or touch yee-fdtd.

## Done when
DoD 1–5 pass; `git diff --stat 9aa79d7..HEAD` shows only `crates/yee-filter/**`
+ the 3 committed docs.
