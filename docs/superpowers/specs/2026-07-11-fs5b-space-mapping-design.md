# FS.5b.0 — Aggressive space mapping walking skeleton (design)

**Date:** 2026-07-11 · **Track:** FS.5 (`FULL-SUITE-ROADMAP.md`) · **ADR:** 0213

## Problem

The textbook space-mapping setup is exactly Yee's pair: cheap closed forms
(coarse) + full-wave EM (fine). Direct BO on the fine model (R.4) spends
10–30 EM solves; aggressive space mapping (Bandler et al.) typically needs
3–6 because each fine evaluation is *aligned* against the coarse model
rather than treated as a black-box sample. FS.5b.0 ships the ASM engine
with a closed-form fine stand-in; the EM-fine gate on the R.4 scenario is
FS.5b.1.

## Design

`yee_surrogate::spacemap`:

- **Parameter extraction** `extract(coarse, target, z0, cfg)`: align the
  coarse model to an observed fine response — Gauss–Newton on
  `‖c(z) − y‖²` with central-difference Jacobian (coarse is cheap; FD is
  fine), plain normal equations via nalgebra, step-halving line search.
- **ASM loop** `space_map(fine, coarse, z_star, x0, cfg)`:
  `z_star` = coarse-optimal design (caller solves the cheap problem);
  iterate `y_k = fine(x_k)` → `z_k = extract(coarse, y_k)` →
  `e_k = z_k − z_star`; stop when `‖e_k‖` (relative to a caller scale)
  < tol; else Broyden-update `B` and step `x_{k+1} = x_k − B⁻¹ e_k`.
  `B₀ = I` (the classic ASM start: assume the spaces are aligned).
- Deterministic: no randomness anywhere (extraction and Broyden are both
  deterministic); every run reproduces bit-for-bit.

## Gates (`crates/yee-surrogate/tests/spacemap.rs`, instant)

1. **Unit-ish:** extraction recovers a known coarse design from its own
   response; fine = coarse converges in 1 fine evaluation.
2. **`surrogate-sm-001`** (roadmap FS.5 gate, closed-form edition): patch
   two-mode testcase — responses `[f₁₀, f₀₁] = c/(2L√εe), c/(2W√εe)`;
   coarse εe = (εr+1)/2, fine adds Hammerstad–Jensen-style W/h dependence
   plus a fringing length extension ΔL (a physically shaped warp, not a
   toy offset). Spec (2.45, 3.10) GHz. Assert ASM meets the spec to
   ≤ 0.1 % in ≤ 5 fine evaluations **and** direct BO (`bo::minimize` on
   the fine mismatch, same fine-evaluation budget) lands ≥ 5× worse —
   measured numbers pinned in the test.

## Lane

`crates/yee-surrogate/**`, `docs/**`, `FULL-SUITE-ROADMAP.md`.
