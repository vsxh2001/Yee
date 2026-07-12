# FS.0b.2b — graded convergence loop (`converge_two_port_graded`)

**Date:** 2026-07-12 · **Track:** FS.0 (top priority, `FULL-SUITE-ROADMAP.md`)
**Builds on:** FS.0a `converge_two_port` (ADR-0204), FS.0b.1 `auto_spacings`
(ADR-0210), FS.0b.2a `two_port_board_jobs_graded`, FS.0b.2-GPU (ADR-0214).
**Plan:** `docs/superpowers/plans/2026-07-12-fs0b2b-graded-convergence.md`

## Problem

The FS.0a adaptive-pass loop refines a **uniform** grid: cells ×2^1.5 per
pass, and the finest pass dominates (6.68 M cells on the stub board).
FS.0b.1 measured the graded grid reproducing the uniform-converged notch at
0.190× the cells. FS.0b.2b closes FS.0b: the convergence loop itself runs
on graded grids.

## Design

1. **`GradedMeshOptions.scale: f64`** — resolution multiplier in `(0, 1]`
   applied by `auto_spacings` to both the coarse ceiling and the fine band
   spacing (1.0 = the rulebook values). This is the loop's one refinement
   knob; everything else in the options is already metre-denominated, so
   the ADR-0204 constant-physics hygiene holds by construction.
2. **`ConvergencePass.cells: usize`** — both loops report per-pass cell
   counts (the graded payoff must be measurable, not asserted).
3. **`converge_two_port_graded(dut, f_max_hz, opts, freqs, tol, max_passes)`**
   in `yee_engine::automesh`, returning the existing `Converged`:
   - per pass: `coarse = auto_dx_bulk · scale`; rescale in cells what the
     fixture denominates in coarse cells — `mesh.npml ← round(npml_m/coarse)`
     (absorber keeps physical thickness) and
     `spacing_cells ← round(spacing_m₀/coarse)` (probe βd stays put);
     build both jobs via `two_port_board_jobs_graded` (one DUT-derived
     grid), measure the launch-normalized double ratio, judge **linear**
     Δ|S21| (identical criterion to FS.0a); `scale ← scale/√2`.
   - `n_steps: None` keeps the engine-graded-001 physical-window rule
     (`9000·0.3 mm / fine`), which scales with the pass automatically.

## Gate

**`engine-automesh-002`** (release, `#[ignore]`, blanket yee-engine CI
step): the S.6 stub board, no hand-set dx anywhere. Asserts: loop verdict
converged; final-pass notch within 5 % of the uniform-converged 4.85 GHz
(ADR-0204 trajectory) at ≥ 20 dB depth; every graded pass ≤ 0.35× the
cells of the equivalent-resolution uniform grid (computed arithmetically
from the pass's own domain extents at its fine spacing).

Expected economics (why this should converge a pass early): graded pass-0
already resolves trace edges at coarse/2 = 0.267 mm — the resolution where
the uniform trajectory converged — so pass 1 (fine 0.189 mm) should move
|S21| very little.

## Non-goals

GPU-backend selection in the gate (CI has no adapter; the graded GPU path
is ADR-0214's lane and the nightly's job); multi-DUT convergence caching;
studio exposure (FS.5c lane).
