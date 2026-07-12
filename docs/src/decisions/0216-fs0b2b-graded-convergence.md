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

## Gate — `engine-automesh-002` (release, `heavy-weekly.yml`)

The S.6 stub board, no hand-set dx anywhere. Asserts: loop verdict
converged; notch within 5 % of the designed 5 GHz at ≥ 20 dB depth; and
**every pass ≤ 0.35× the cells of the equivalent-resolution uniform grid**
(computed from the pass's own reconstructed spacings — the gate also
cross-checks the reconstruction against the loop's reported `cells`).

Measured (this machine, release, 2026-07-12; full per-bin dump in the
gate's `--nocapture` output):

| pass | coarse | fine | cells | ×uniform-eq | notch |
|------|--------|------|-------|-------------|-------|
| 0 | 0.533 mm | 0.267 mm | 1.27 M | 0.189 | 4.900 GHz @ −37.2 dB |
| 1 | 0.377 mm | 0.189 mm | 3.52 M | 0.185 | 5.100 GHz @ −37.5 dB |
| 2 | 0.267 mm | 0.133 mm | 9.87 M | 0.185 | 5.050 GHz @ −37.7 dB |

Converged at pass 2 (criterion band): final linear Δ|S| = **0.1351**
(tol 0.20); converged notch err **1.0 %** at −37.7 dB; every pass at
**≤ 0.19×** the equivalent-resolution uniform cells. Wall time ≈ 90 min
release on a fast box — hence the gate lives in `heavy-weekly.yml`, not
per-PR CI (per-PR coverage stays with engine-graded-001 + the loop's
fast unit tests).

## The band-edge lesson (measured, one red run + one forensics run)

The first run judged convergence on the full 3.5–6.0 GHz request and
failed with Δ|S| = 0.8423. The per-bin forensics localized the delta
entirely to **5.85–6.0 GHz**: a spurious dip wandering 5.95 → 5.90 GHz
between the fine passes, and a **non-physical +1.05 dB** bin at exactly
6.0 GHz (proof of measurement artifact, not physics). 6.0 GHz is
`f_max` — the frequency the λ/20 mesh rule is *designed to*, i.e. the
edge of mesh validity; the uniform ADR-0204 dump flagged the same region
as the staircase-limited skirt. **Rule: convergence-criterion bins must
stop at ~0.96·f_max** — bins at the mesh-design edge are the least
trustworthy by construction and must not drive the refinement verdict.
The gate's criterion band is 3.5–5.75 GHz. Root-causing the moving
band-edge dip is queued (FS.0b.2c candidate; low priority — it sits
outside the band any design flow should read).

## Non-goals

GPU-backend selection inside the gate (CI has no adapter; the nightly
covers ADR-0214), multi-fidelity pass reuse, studio exposure (FS.5c),
the band-edge dip root-cause (queued above).
