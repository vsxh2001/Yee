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

No library change was needed: the ADR-0213 API takes plain
`&dyn Fn(&[f64]) -> Vec<f64>` models, and the engine side is the
certified FS.0b.2a fixture. This ADR is the validation record.

## Gate — `sm-em-001` (`yee-filter`, release, dedicated CI step)

Asserts: converged within 5 fine evals at tol 0.005 (≈ 0.5 % in
frequency via df/dl ≈ −f/l); `n_fine_evals` ≤ 4; final measured notch
within 0.75 % of 5.3 GHz; final error strictly beats the coarse-only
seed error.

## Measured

(pinned from the first green run)

## Non-goals

Multi-knob R.4 BPF scenario (FS.5b.2), BO-vs-ASM at EM cost (the
ADR-0213 closed-form comparison stands), studio exposure (FS.5c).
