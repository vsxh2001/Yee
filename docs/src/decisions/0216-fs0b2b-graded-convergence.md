# ADR-0216: FS.0b.2b — the convergence loop on graded grids

**Date:** 2026-07-12 · **Status:** accepted · **Track:** FS.0 (close-out)
**Spec:** `docs/superpowers/specs/2026-07-12-fs0b2b-graded-convergence-design.md`
**Builds on:** ADR-0204 (`converge_two_port`), ADR-0210 (`auto_spacings`),
FS.0b.2a (`two_port_board_jobs_graded`), ADR-0214 (graded on the GPU).

## Decision

1. **`GradedMeshOptions.scale`** — one resolution multiplier in `(0, 1]`,
   applied by `auto_spacings` to the *unscaled* rule outputs (coarse and
   fine together, so a pass refines feature bands and bulk alike; deriving
   fine from the scaled coarse would let a static `min_feature/2` term pin
   the fine bands while the bulk refines). Validated alongside the other
   knobs; `for_board` defaults it to 1.0, so every existing caller is
   unchanged.
2. **`ConvergencePass.cells`** — both loops report per-pass cell counts;
   the graded payoff is measured, not asserted from folklore.
3. **`converge_two_port_graded(dut, f_max, opts, freqs, tol, max_passes)`**
   — the FS.0a loop with `scale ← scale/√2` as the refinement step. The
   metre-denominated mesh options are simply never touched (the ADR-0204
   hygiene by construction); the two coarse-cell-denominated fixture knobs
   are rescaled per pass — `mesh.npml` (absorber keeps physical thickness)
   and `spacing_cells` (probe βd stays put). `n_steps: None` keeps the
   fixture's physical-window rule, which follows the pass's fine spacing.
   Identical linear-ΔS criterion, honest unconverged verdict.

## Gate — `engine-automesh-002` (release, dedicated graded CI job)

The S.6 stub board, no hand-set dx anywhere. Asserts: loop verdict
converged; notch within 5 % of the designed 5 GHz at ≥ 20 dB depth; and
**every pass ≤ 0.35× the cells of the equivalent-resolution uniform grid**
(computed from the pass's own reconstructed spacings — the gate also
cross-checks the reconstruction against the loop's reported `cells`).

Measured (this machine, release): see the gate's `--nocapture` output
pinned in the test comments; the headline numbers live in the roadmap row.

## Why this should converge a pass earlier than uniform

Graded pass 0 already resolves trace edges at `coarse/2` = 0.267 mm — the
resolution at which the uniform trajectory converged (ADR-0204). Pass 1
(fine 0.189 mm) is expected to move |S21| little; the loop's criterion
then stops at 2 passes where the uniform loop needed 3.

## Non-goals

GPU-backend selection inside the gate (CI has no adapter; the nightly
covers ADR-0214), multi-fidelity pass reuse, studio exposure (FS.5c).
