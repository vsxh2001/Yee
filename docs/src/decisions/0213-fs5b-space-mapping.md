# ADR-0213: FS.5b.0 — aggressive space mapping (Broyden ASM, deterministic)

**Date:** 2026-07-11 · **Status:** accepted · **Track:** FS.5 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-11-fs5b-space-mapping-design.md`

## Context

Yee's model pair is the textbook space-mapping setup: cheap closed forms
(coarse) + full-wave EM (fine). Direct BO on the fine model (the R.4 loop)
treats every EM solve as a black-box sample; ASM (Bandler) aligns each fine
response against the coarse model and typically converges in 3–6 fine
evaluations.

## Decision

`yee_surrogate::spacemap`:

1. **Parameter extraction** = Gauss–Newton on `‖coarse(z) − y_fine‖²`,
   central-difference Jacobian (the coarse model is cheap by definition),
   step-halving line search. Started from `z_star` each iteration — the
   best-informed initial point and the one that keeps extraction
   single-basin for well-posed warps.
2. **ASM loop**: `e_k = extract(coarse, fine(x_k)) − z_star`; Broyden
   update `B += ((Δe − B h)hᵀ)/(hᵀh)`; step `x ← x − B⁻¹e`. `B₀ = I`
   (assume aligned spaces — the classic start).
3. **Fully deterministic** — no RNG anywhere, so results reproduce
   bit-for-bit (test-pinned), matching the FS.5a yield-module policy.

## Measured gate — `surrogate-sm-001` (instant, GREEN first run)

Patch two-mode testcase: coarse `εe = (εr+1)/2`; fine adds
Hammerstad–Jensen-style `εe(W/h)` **and** Hammerstad fringing ΔL — a
physically shaped warp, not a toy offset. Spec (2.45, 3.10) GHz, budget
5 fine evaluations:

- **ASM: 0.00143 % spec error in 4 fine evals** (converged, tol 1e-4 on
  the scaled misalignment).
- **Direct BO (same budget, generous ±50 % box): 44.8 %.**
- Asserts pinned at ≤ 0.1 % ASM / ≥ 5× BO ratio — far inside the measured
  ~31 000× margin. Determinism gate: two runs bit-identical.

## Consequences

- **FS.5b.1** (queued): fine = the R.4 engine scenario — wire
  `space_map` into the hairpin-BPF refine loop and gate "fewer EM solves
  than direct BO" with real solves. The engine closure just wraps
  `verify` the way the BO loop already does.
- Extraction currently starts from `z_star`; strongly non-monotone warps
  may need multi-start extraction — deferred until a real case shows it.
