# ADR-0159: F1.2.1.1 — EM coupling design-curve gap correction (full-wave filter S21 toward the mask)

**Status:** Accepted
**Date:** 2026-06-04
**Related:** ADR-0147 (#1 goal: a validated full-wave filter S21), ADR-0154 N3 (the filter floored
~−27 dB), ADR-0155/K2 (k_imp≠k_eps — the dimensioning error), ADR-0157 (F1.2.1.0 single-gap BO
mechanism), ADR-0156 (Qe deferred), the full-wave-filter-S21 scope + the gap-mesh-REFUTED probe,
[[reference-em-in-loop-space-mapping]] (the research: ASM / Hong-Lancaster / OSS),
[[fem-driven-sweep-s21-viable]], [[project-filter-design-final-goal]].

---

## Context

The maintainer chose the full-wave 3-pole filter S21 (ADR-0147 #1). Two probes settled the approach:

1. **Gap-mesh lever REFUTED** (cheap probe): |S21|(F0) stays flat ~−27 to −41 dB across 10k→105k
   tets (dx 1.6→0.3 mm) — no convergence toward the mask. The ~30 dB floor is **NOT** mesh
   resolution.
2. **It is DIMENSIONING-bound:** the dimensioner sized the gaps with the analytic impedance-k
   (`coupling_coefficient`), but K2 showed it diverges ~37 % from the realized resonant-k at this
   filter's tight gaps (S/W≈0.85). The synthesized geometry does not realize the target coupling, so
   no mesh makes it pass the mask. The fix is to **correct the geometry** using the EM-measured
   coupling — the EM-in-loop refinement.

**Research (per [[reference-em-in-loop-space-mapping]], the maintainer's "implement known methods,
don't reinvent"):** the established methods are (a) the **full-wave coupling design-curve** approach
(Hong-Lancaster: EM-extract K(gap) from isolated resonator pairs, read off the gap realizing each
target K) and (b) **Aggressive Space Mapping** (Bandler: coarse synthesis + ~5 fine-EM evals,
Broyden mapping) for the coupled multi-parameter case. The k-formula `k=(f_o²−f_e²)/(f_o²+f_e²)` is
the universal one and is exactly `yee_fem::coupled_resonator_k`.

**Smoothness probe (12 boxed FEM solves) — decisive + revises F1.2.1.0:** the FEM K(gap) is a
**smooth, monotone-decreasing curve** (var-box `probe_with_gap`: 0.0611/0.0519/0.0481/0.0433/0.0394/
0.0306 over gap 1.5→3.0 mm). The F1.2.1.0 "box_w-staircase" was **over-attributed** — K(gap) is not
a wild staircase; F1.2.1.0's 13.4 % was the coarse BO budget / an extraction outlier (its
k(1.587 mm)=0.0346 does not reproduce — it interpolates to ~0.056). A fixed-box_w "fix" is NOT needed
(it actually adds a tiny non-monotone wiggle). **A smooth monotone K(gap) ⇒ a simple, robust
per-gap root-find converges to a target K in a few evals** — no fancy optimizer required for the
(near-)separable edge-coupled case.

## Decision

F1.2.1.1 = **EM coupling design-curve gap correction.** Replace the analytic-k bisection that
`dimension_edge_coupled` uses (ADR-0097, bisects `coupling_coefficient`/impedance-k) with a
**FEM-k root-find**: for each synthesis target coupling `K_target[i] = FBW·m[i][i+1]`, root-find the
gap `g_i` such that the FEM-measured `coupled_resonator_k(g_i) = K_target[i]` (secant/bisection on
the smooth monotone K(gap)). The corrected gaps realize the resonant coupling the analytic model
mis-sized. Then re-grade the corrected 3-pole filter's full-wave S21 vs the Cheb mask.

Ordered bricks, each machine-checkable:

| # | Brick | Gate | Risk | Cost |
|---|-------|------|------|------|
| **B1** | FEM-k per-gap corrector (`correct_gap_fem_k`): secant/bisection on `coupled_resonator_k(gap)` → the gap realizing a target K | NON-circular: a deliberately-mis-dimensioned seed gap converges so `\|K_fem(g*)−K_target\|/K_target < 8 %` (the F1.2.1.0 bar, now reachable on the smooth curve) in ≤ ~5 FEM evals; the corrected gap differs from the analytic-k seed | eng (de-risked: K(gap) smooth+monotone) | ~5 FEM evals (~25 min) |
| B2 | Apply B1 to all the 3-pole filter's gaps → re-build → re-grade the corrected geometry's full-wave S21 vs the Cheb mask | HONEST: record the corrected filter's in-band peak / mask margin vs the −27 dB analytic-geometry baseline; assert the measured IMPROVEMENT + the asymmetry; if it clears the mask, assert it | research-open (gap interactions: the isolated-pair K vs the in-filter coupling) | a full filter sweep |
| B3 (if B2 short) | ASM (Bandler) multi-D Broyden over all gaps using the full-filter response — handles gap interactions the per-gap root-find ignores | the corrected filter clears (or approaches) the mask | research-open | ~5 ASM iters |

**Start per-gap (B1/B2)** — the smooth separable K(gap) makes it the simplest correct method
(Hong-Lancaster design-curve); escalate to ASM (B3) only if gap interactions keep B2 short of the
mask. Qe-feed dimensioning stays deferred (ADR-0156; FEM Qe is numerical-Q-floor-limited).

## Consequences

**B1 is a ~1-2 day eng brick** (a robust 1-D root-find on the proven-smooth FEM K(gap)), high
confidence. **B2 is the real test** of whether correcting the gaps lifts the full-wave S21 toward the
mask — research-open (gap interactions / the isolated-pair-vs-in-filter-coupling approximation). If
B2 lifts the filter substantially but stops short, that is an **honest documented result** (the
corrected-geometry S21 + its mask margin + the residual cause), and B3 (ASM) is the escalation — not
a fake/weaken. This pursues the ADR-0147 #1 goal via the literature-standard full-wave coupling-design
flow.

**Honesty / non-circular:** `K_target` is the synthesis spec (`FBW·m`); `K_fem` is the EM
measurement; the root-find moves the gap. Nothing derives the target from the EM. The reviewer
enforces `gate_is_real`. Heavy FEM runs boxed (`scripts/yee-box.sh` ≤14 g/3 cpu); the agent writes
the code, the orchestrator runs the heavy gates (misfire-split). No gate merges until green.

**Not in scope / do NOT reopen:** finer-mesh as the lever (REFUTED), the fixed-box_w "fix" (not
needed — K(gap) is smooth), Qe-feed dimensioning (ADR-0156), the FDTD-resonant route (ADR-0108),
mom-002/003, fem-eig-006.

---

## Outcome (SHIPPED — B1 merge `b09db3a`, B2 merge `71c0262`)

**B1 — `correct_gap_fem_k` (the corrector) — VALIDATED.** A 20.3 %-off seed gap (2.000 mm,
k_fem=0.0481) converged to 2.0 % off (2.125 mm, k_fem=0.0392) in 3 FEM evals on the smooth K(gap)
(gate `fem-coupling-correct-001`). The EM-in-loop geometry-correction mechanism works. Reviewer no P0.

**B2 — apply to the filter + re-grade — the thesis is REFUTED (honest NO-GO, gate GREEN).** Two
orchestrator fixes were needed first: (i) `coupled_resonator_k`'s sweep band was **hardcoded
2.10–2.70 GHz** (the K1/K2 probe's f0=2.4 GHz) → the filter's f0=2.0 GHz resonances (~1.9–2.05 GHz)
fell below it → every eval `!peaks_resolvable`; fixed to **track f0** via the analytic even/odd split
(`f_lo=f_e_lo−0.20·fc`, `f_hi=f_e_hi+0.12·fc`) — **K1/K2 regression-confirmed GREEN** (k_fem
0.0612/0.0464/0.0295 monotone). (ii) the design-curve base box was oversized (~160k tets, OOM) →
matched the shipped `probe_with_gap` config (2.5·h clearance, ~80k tets).

Measured (boxed, `--release`, self-verified + reviewer APPROVE / gate-honesty sound):

- **`target_k=0.0756` is UNREACHABLE as a resonant-split k** on the 1.912 mm / FR-4 line — the
  resonant-split k **saturates ≈0.064** at the tight-gap floor (evals 1.25→0.0476, 0.875→0.0625,
  0.594→0.0641, 0.523→0.0570). The correction drives the gaps **1.622 → 0.594 mm** (best 15.2 % off,
  converged=false). The **impedance-k synthesis target over-specifies the coupling** vs what a
  resonant-split can deliver on this geometry.
- The corrected filter **LIFTS the in-band peak +5.83 dB** over N3 (−27.38 → −21.55 dB) — stronger
  coupling genuinely helps — but **still MISSES the strict Cheb mask by ~26 dB** and remains a
  ~−22 dB flat shelf, with the asymmetry margin degraded to +0.70 dB (N3 was +2.10 dB) by the
  over-tightened gaps.

**DECISIVE CONCLUSION: dimensioning-correction is a real but MINOR lever (~+6 dB); the dominant floor
is the aperture-coupling PORT fidelity (the ADR-0154 N3 finding), NOT gap dimensioning.** Neither this
per-gap correction nor a multi-D ASM over gaps (**B3 — NOT pursued**, same port wall) clears the mask
without a higher-fidelity aperture-coupling port. **The F1.2.1.x EM-in-loop dimensioning track is
port-floored and CLOSED.** A mask-clearing full-wave filter S21 needs the **deferred higher-fidelity
port** (aperture/frill reciprocity or a numerical-eigenmode aperture port) — a multi-week,
maintainer-funded track; do NOT auto-start it, and do NOT reopen B3.

**What durably landed:** the `correct_gap_fem_k` corrector (B1, reusable) + the **f0-tracking band fix
to `coupled_resonator_k`** (a genuine bug fix — the function silently only worked at f0=2.4 GHz; now
correct for any filter f0, K1/K2-regression-green) + the honest B2 gate (`#[ignore]`'d, run boxed
on-demand; NOT in routine CI — the mechanism is covered by B1 + K1/K2).

---

## References
- The EM objective (fine model, smooth K(gap) confirmed): `yee_fem::coupled_resonator_k`
  (`crates/yee-fem/src/coupled_resonator_k.rs`); the smoothness sweep data in
  the F1.2.1.1.0 probe run (`/tmp/fixedbox_probe.log`, spike `6a22824`, now reaped).
- The coarse model corrected: `yee_filter::dimension_edge_coupled` (`crates/yee-filter/src/dimension.rs`),
  `target_k = FBW·m`.
- Method + formulas: [[reference-em-in-loop-space-mapping]] (Hong-Lancaster design-curve, Bandler ASM,
  `k=(f_o²−f_e²)/(f_o²+f_e²)`, scikit-rf/openEMS).
- Spec: `docs/superpowers/specs/2026-06-04-f1-2-1-1-em-coupling-design-curve-design.md`;
  plan: `docs/superpowers/plans/2026-06-04-f1-2-1-1-em-coupling-design-curve.md`.
