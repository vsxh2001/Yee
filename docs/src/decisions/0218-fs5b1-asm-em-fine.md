# ADR-0218: FS.5b.1 — space mapping with the engine as the fine model

**Date:** 2026-07-12 · **Status:** accepted · **Track:** FS.5
**Spec:** `docs/superpowers/specs/2026-07-12-fs5b1-asm-em-fine-design.md`
**Builds on:** ADR-0213 (`yee_surrogate::spacemap`, closed-form warp),
FS.0b.2a graded fixture, ADR-0216 (criterion-band rule).

## Decision

Close the coarse-model ↔ EM loop: `space_map` drives the S.6 stub board's
**measured** notch (graded engine two-port, double-ratio |S21|, 25 MHz
bins, parabolic sub-bin refinement, band capped per ADR-0216) onto an
off-design 5.3 GHz target, with the S.6 TL formula as the coarse model
and the stub length as the single knob. The fine closure logs every
`(length, frequency)` pair, so the gate reports the whole trajectory —
seed error included — at zero extra solve cost.

One library change WAS needed, and it is the increment's real finding:
**`GradedMeshOptions.snap_edges`** (default off). The first run measured
the fine model as a **staircase** in the design variable — three stub
lengths spanning 34 µm all read the identical 5.3530 GHz notch (the
point-sampled rasterization quantizes each edge to the local fine cell,
~2.5 % frequency steps here), and Broyden oscillated inside one step,
never converging (5 evals, final err 0.999 %). Snapping shifts the
nearest grid node exactly onto every trace-AABB edge (nudge spread over
4 nodes/side, max width change fine/8, junction ratio ≤ 9/7 < 1.3 — the
compute-019-certified regime; only the 8 affected cells are rewritten so
untouched coarse cells stay bit-identical for the fixture's
bit-equal-coarse probe placement), making the rasterized geometry track
the requested geometry continuously. **Any gradient/Broyden design loop
on a rastered solver needs this** — the quantization defeats it
otherwise.

## Gate — `sm-em-001` (`yee-filter`, release, dedicated CI step)

Asserts: converged within 5 fine evals at tol 0.005 (≈ 0.5 % in
frequency via df/dl ≈ −f/l); `n_fine_evals` ≤ 4; final measured notch
within 0.75 % of 5.3 GHz; final error strictly beats the coarse-only
seed error.

## Measured (GREEN, 531 s release)

| eval | stub | measured notch | err vs 5.3 GHz |
|------|------|----------------|-----------------|
| 0 (= coarse optimum) | 7.1322 mm | 5.3524 GHz | 0.988 % |
| 1 | 7.2080 mm | 5.3064 GHz | **0.121 %** |

Converged in **2 fine evals** (misalignment 0.00132, tol 0.005). Bonus
measurement: snapping alone improved the *seed* — the identical
7.1322 mm stub measured 5.1822 GHz (2.22 % err) un-snapped vs
5.3524 GHz (0.99 %) snapped; half the "coarse-model bias" was actually
rasterization error, and the remaining ~1 % is the true TL-formula
bias the mapping then absorbs in one Broyden step.

## Non-goals

Multi-knob R.4 BPF scenario (FS.5b.2), BO-vs-ASM at EM cost (the
ADR-0213 closed-form comparison stands), studio exposure (FS.5c).
