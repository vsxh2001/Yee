# FS.5b.1 — aggressive space mapping with the engine as the fine model

**Date:** 2026-07-12 · **Track:** FS.5 (optimization maturity)
**Builds on:** FS.5b.0 `yee_surrogate::spacemap` (ADR-0213, closed-form
warp testcase), the graded fixture (FS.0b.2a) + ADR-0216 band rule.
**Plan:** `docs/superpowers/plans/2026-07-12-fs5b1-asm-em-fine.md`

## Problem

FS.5b.0 proved the ASM machinery on a closed-form warp (0.00143 % in 4
fine evals vs BO's 44.8 % at the same budget). The FS.5 roadmap gate
needs the fine model to be a **real EM solve**: close the loop
coarse-model ↔ engine.

## Design

Scenario: the S.6 open-stub notch board; one design knob (stub length
`l_s`), one response (measured notch frequency).

- **Coarse model** (instant): the TL formula
  `f = c / (4 (l_s + ΔL) √ε_eff)` — the same model the S.6 synthesis
  uses, known biased vs the engine by ~1 % (ADR-0216 trajectory).
- **Fine model** (~2 release solves ≈ 5 min): `two_port_board_jobs_graded`
  double-ratio |S21|, notch located on a 25 MHz grid inside the
  ADR-0216-safe band (mesh `f_max` = 6.5 GHz, band ≤ 6.2 GHz) and refined
  sub-bin by parabolic interpolation through the notch's three dB bins.
- **Target**: notch at **5.3 GHz** (off the 5.0 design point).
  `z_star` = the coarse inverse; `x0 = z_star` (classic ASM start).
- The fine closure records every `(l_s, f)` pair, so the gate reports the
  full trajectory (seed error = eval 0) at zero extra solve cost.

## Gate — `sm-em-001` (`yee-filter`, release, dedicated CI step)

Asserts: ASM `converged` within `max_fine_evals = 5` (tol 0.005 scaled by
the nominal length ≈ 0.5 % in frequency via `df/dl ≈ −f/l`); final
measured notch within **0.75 %** of 5.3 GHz; final error strictly smaller
than the coarse-seed error (the mapping earned its keep); `n_fine_evals`
≤ 4. Numbers pinned from the first green run.

## Non-goals

Multi-knob (the R.4 BPF scenario — FS.5b.2), BO-vs-ASM at EM cost (the
5b.0 closed-form comparison stands), studio exposure (FS.5c).
