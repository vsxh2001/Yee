# Phase 1.plotting.4 — spec-mask overlay — Implementation Plan

**Spec:** `2026-05-29-spec-mask-overlay-design.md` · **ADR:** ADR-0087

## Lane
`crates/yee-plotters/**` only. Out of lane (yee-gui, yee-filter, any other) →
finding, not fix. Do NOT touch yee-gui (avoids the wgpu build).

## Base
Worktree `worktrees/plotting4`, branch `feature/plotting-4-spec-mask-overlay`,
base `53df105`.

## Pattern files
- `crates/yee-plotters/src/**` — read the existing multi-trace S-parameter
  magnitude draw fn (ADR-0063) and the VSWR work (ADR-0081): reuse its error
  type, colour cycle, axis/range logic, and the render-smoke-test style.

## Steps
1. Add `MaskKind`, `MaskRegion` (doc'd; `#[derive(Debug, Clone, Copy)]`, and
   `PartialEq` for tests).
2. `mask_violations(freqs_hz, trace_db, regions) -> Vec<usize>` — pure, per spec.
3. `draw_sparam_with_mask(path, freqs_hz, traces, regions, title)` — reuse the
   existing magnitude-plot scaffolding; shade each region's forbidden side
   (translucent red) before drawing traces. Match the crate's existing draw-fn
   signature/error type exactly.
4. Tests: a `mask_violations` unit test (Ceiling + Floor + compliant cases) and
   a render smoke test writing to `CARGO_TARGET_TMPDIR`/tempfile, asserting the
   output file is non-empty (mirror the ADR-0081 VSWR render test).

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-plotters --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-plotters --jobs 2
```
`yee-plotters` needs `libfontconfig1-dev` + `pkg-config` to link (CLAUDE.md §10);
they are present in this environment.

## Escape hatch
Blocked >15 min — the existing draw-fn API doesn't compose cleanly with a mask
overlay, or fontconfig link failure → STOP, surface the exact error. Do NOT pull
in yee-gui or a new dependency.

## Done when
DoD 1–5 pass; `git diff --stat 53df105..HEAD` shows only `crates/yee-plotters/**`
+ the 3 committed docs.
