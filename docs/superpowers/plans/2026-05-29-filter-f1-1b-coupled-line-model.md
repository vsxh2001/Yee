# Filter Phase F1.1b.gate — coupled-microstrip even/odd model — Plan

**Spec:** `2026-05-29-filter-f1-1b-coupled-line-model-design.md` · **ADR:** ADR-0094

## Lane
`crates/yee-layout/**` ONLY. Out of lane → finding. Keep `yee-layout` WASM-safe
(pure f64, no native dep) per ADR-0089. No new external dependency.

## Base
Worktree `worktrees/cline`, branch `feature/filter-f1-1b-gate-coupled-line`,
base `bde5cfe`.

## Pattern files
- `crates/yee-layout/src/lib.rs` — house style; the existing HJ single-line
  `microstrip_width` (~line 284) + `eps_eff` (~line 318) to reuse + extend.
- `crates/yee-synth/tests/` — published-table gate style (cite source + numbers).

## Steps
1. `src/lib.rs` (or `src/coupled.rs` + `pub mod coupled; pub use coupled::*;`):
   add `CoupledMicrostrip`, `coupled_microstrip(w,s,h,eps_r)`,
   `coupling_coefficient`. Implement a CITED static even/odd model; document the
   model + accuracy inline. Reuse the single-line helpers where applicable.
2. Tests: `coupled_001_vs_published` (cite source + published Z0e/Z0o; within the
   stated tol) and `coupled_002_monotonic` (z0e>z0o>0; k>0 strictly decreasing
   with gap). Use WebSearch to find/confirm a published data point; cite it.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-layout --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-layout --jobs 2
```
Pure math — sub-second. Do NOT run `cargo test --workspace`.

## Escape hatch
Blocked >15 min — can't find/verify a published Z0e/Z0o data point to gate
against (do NOT fabricate one), or the chosen model diverges >10% from every
reference you find → STOP, surface the model + what references you checked + the
computed-vs-published numbers, and propose either a different cited model or a
looser-but-justified tolerance. Do NOT loosen the gate silently or invent a ref.

## Done when
DoD 1–5 pass; `git diff --stat bde5cfe..HEAD` shows only `crates/yee-layout/**`
+ the 3 committed docs; `yee-layout` still has no FDTD/native dep.
