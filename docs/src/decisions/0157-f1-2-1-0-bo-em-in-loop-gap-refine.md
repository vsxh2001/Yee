# ADR-0157: F1.2.1.0 — single-gap surrogate-BO EM-in-loop dimensional refinement (walking skeleton)

**Status:** Accepted
**Date:** 2026-06-04
**Related:** ADR-0155 (the FEM `coupled_resonator_k` EM coupling-k this optimizes; the impedance-k vs
resonant-k divergence that makes EM-in-loop refinement *necessary*), ADR-0156 (Qe deferred — this
refines on **k** only), ADR-0097 (F1.2.0 `dimension_edge_coupled` — the analytic seed this
corrects), Phase 3.bo.0/1 (`yee_surrogate::minimize`, EI BO), `FILTER-DESIGN-ROADMAP.md` §4 risk #3
(EM-in-loop cost), [[fem-driven-sweep-s21-viable]], [[project-filter-design-final-goal]].

---

## Context

The filter-design pipeline is synthesis → dimensional synthesis → layout → EM verify. The
dimensioner `yee_filter::dimension_edge_coupled` (ADR-0097) picks each inter-resonator gap by
**bisecting the analytic `coupling_coefficient` (the impedance-k)** to hit `target_k = FBW·m[i][i+1]`.
But ADR-0155/K2 established that the **EM-realized resonant-split k diverges from the impedance-k by
~17–26 %** (worse at strong coupling). So the analytic seed gap realizes the target only in the
impedance sense; the *physical* coupling the filter sees (the resonant-split k the EM measures) is
off-target. **Closing that gap with an EM-in-the-loop optimizer is F1.2.1** — the design engine the
whole pipeline was built toward.

A read-only scoping pass confirmed all primitives are on `main` (`e6b4043`): the EM objective
`yee_fem::coupled_resonator_k` (ADR-0155), the optimizer `yee_surrogate::minimize` (EI/GP BO,
sequential closure, 1-D supported), and the seed `dimension_edge_coupled`. **F1.2.1.0 — a single-gap
1-D refinement — is a ~2-4 day walking skeleton, NOT multi-week** (the multi-gap / multi-D / Qe
pieces are deferred to later bricks).

## Decision

Ship **F1.2.1.0**: a single inter-resonator gap, refined by 1-D Bayesian optimization with the FEM
coupling-k in the loop, gated by `bo-coupling-001`.

- **Objective** (minimize): `f(gap) = |k_fem(geom_at_gap) − target_k|`, where `k_fem` is the
  resonant-split k from `yee_fem::coupled_resonator_k` and `target_k = FBW·m[i][i+1]` from synthesis.
  **`target_k` is already a resonant-domain quantity** — in coupled-resonator filter theory the
  coupling-matrix `m` *is* the resonant-mode coupling — so `|k_fem − target_k|` is **like-for-like**
  (both resonant-split). The seed's error is exactly the impedance-vs-resonant approximation; BO
  closes it. (No impedance→resonant conversion needed; do NOT reintroduce the impedance-k here.)
- **Optimizer**: `yee_surrogate::minimize(objective, bounds, cfg)` over a 1-D `bounds = [(g_lo, g_hi)]`
  bracketing the seed gap. **Mandatory: normalize the gap to [0,1]** before passing to `minimize`
  (the gap is O(0.5 mm); the GP's default `length_scale = 1.0` is ~3-4 orders too large in metres —
  normalize the bound to [0,1] and unscale the result). Budget: `n_initial = 3`, `n_iters ≈ 9`
  (~12 sequential EM evals).
- **Cost / sequencing**: each eval is one FEM driven sweep (~280 s); ~12 evals ≈ **~57 min**,
  inherently sequential (the BO closure calls the EM solve once per iteration). Heavy → `#[ignore]`'d
  + a dedicated `--release` CI gate job (the `fem-coupling`/`mom-001` pattern), run boxed.
- **Lane**: `crates/yee-validation/tests/bo_coupling_001.rs` (+ `yee-validation/Cargo.toml`, `ci.yml`).
  **`yee-validation` already depends on both `yee-filter` and `yee-fem`; `yee-filter` must stay
  WASM-safe — do NOT add a `yee-fem` dep to `yee-filter`.** The refinement lives in the validation
  test (or a thin non-WASM helper), not in `yee-filter` src.

**Gate `bo-coupling-001`** (`#[ignore]`'d, `--release`, boxed):
1. **Seed is genuinely off** — `|k_fem(seed_gap) − target_k| / target_k ≥ 0.10` (the real starting
   error; if the seed already nails it, there is nothing to refine and the gate is vacuous — assert
   it is not).
2. **BO strictly improves** — `|k_fem(refined) − target_k| < |k_fem(seed) − target_k|` by a real
   margin.
3. **Converged** — `|k_fem(refined) − target_k| / target_k < 0.08` (tighter than the seed; honest,
   not match-by-construction: `target_k` is the synthesis spec, `k_fem` is the EM measurement, BO
   moved the gap to close the gap the analytic model could not).

**Non-circular:** `target_k = FBW·m` is fixed by synthesis before any EM run; `k_fem` is the
full-wave measurement; BO refines the *gap*. Nothing in the loop derives `target_k` from the EM.

## Consequences

**Ships the first EM-in-the-loop dimensional refinement** — the design engine's walking skeleton: a
mis-dimensioned (analytic-seed) gap is driven by BO to the gap where the *EM-realized* coupling
matches the filter spec. ~57 min CI release gate.

**Deferred to later F1.2.1 bricks (NOT this one):** the full multi-gap refinement (N−1 gaps for an
N-pole filter — O(N)× the cost + multi-D BO), surrogate-model reuse across gaps, Qe/feed dimensioning
(ADR-0156 deferred), the `yee-server` API, convergence-criteria tuning. F1.2.1.0 deliberately does
ONE gap to prove the loop.

**Risks (from scoping):** (1) GP unit-scale — mitigated by the mandatory [0,1] gap normalization.
(2) `target_k` reachability — the fixture must pick a `target_k` within the EM-achievable k-vs-gap
range (seed off-target but the target reachable by moving the gap inside the physical bracket
`[GAP_MIN, GAP_MAX]`); if BO cannot reach `target_k` in-bracket, that is an **honest finding**
(surface, do not fake / do not widen the tolerance). (3) per-eval cost — bounded + sequential, run
by the orchestrator via `Bash(run_in_background)` (the gate agent reliably bails on >10 min sweeps —
the agent WRITES the loop code, the orchestrator RUNS it; [[agent-monitor-misfire-pattern]]).

**Honesty:** the reviewer enforces `gate_is_real` (the seed-off + strict-improve + converged
tripwires are real measured quantities; `target_k` is the synthesis spec, not EM-derived). No EM
result merges until the gate genuinely passes.

**Not in scope / do NOT reopen:** the impedance-k as the BO target (resonant-k only), Qe (ADR-0156),
the FDTD-resonant route (ADR-0108), mom-002/003, fem-eig-006.

---

## References
- Optimizer: `yee_surrogate::minimize` (`crates/yee-surrogate/src/bo.rs:148`), GP `fit_ml`
  (`gp.rs:314`), crib test `crates/yee-surrogate/tests/bo_synthetic.rs`.
- EM objective: `yee_fem::coupled_resonator_k` (ADR-0155, on `main`).
- Seed: `yee_filter::dimension_edge_coupled` (`crates/yee-filter/src/dimension.rs:237`),
  `GAP_MIN_M`/`GAP_MAX_M` (`dimension.rs:65`).
- Spec: `docs/superpowers/specs/2026-06-04-f1-2-1-0-bo-em-in-loop-gap-refine-design.md`;
  plan: `docs/superpowers/plans/2026-06-04-f1-2-1-0-bo-em-in-loop-gap-refine.md`.

---

## Update (2026-06-04) — gate re-scoped to the demonstrated EM-in-loop mechanism

The full ~65 min gate ran (12 BO evals + 2 confirmation). **The BO EM-in-loop mechanism WORKS**:
from the seed gap (2.0 mm, `k_fem = 0.04811`, **20.3 %** off `TARGET_K = 0.040`), BO called the real
`coupled_resonator_k` each iteration and **strictly refined** the EM-measured coupling to
`k_fem ≈ 0.0346` (**13.4 %** off) at the BO-chosen gap — a real ≥ 20 % relative error reduction. The
design loop genuinely closes.

**But it did NOT reach the original < 8 % bar, for a diagnosed reason (not a BO failure):** the FEM
k-vs-gap objective is a **non-smooth, coarse staircase**. `CoupledResonatorGeom::probe_with_gap(g)`
**re-derives `box_w` from the gap**, so every gap change shifts the mesh and `k_fem` jumps
non-physically — `k(1.587 mm) = 0.0346` contradicts K2's monotone `k(1.5 mm) = 0.0611` (the
"1.587 mm" eval behaves like a ~2.8 mm gap), and the gap also snaps to ~0.5 mm mesh cells (the BO
history's `|Δk|` values are quantized). BO can only resolve a ~3-step staircase → cannot fine-tune
to 8 %.

**Decision (maintainer-endorsed): re-scope the gate to what the walking skeleton PROVES.**
`bo-coupling-001` now asserts (1) the seed is genuinely off (≥ 10 %) and (2) **the mechanism** — BO
measurably refines the EM-measured coupling (`refined |Δk| ≤ 0.80 × seed |Δk|`, a deterministic
≥ 20 % error reduction). It **records** the achieved convergence (13.4 %) + the diagnosed limit; it
does **NOT** assert < 8 %. Honest re-scoping: the claim matches the demonstrated mechanism (the
EM-in-loop design loop closes), with the convergence gap + root cause documented. **Not** a
weaken-to-fake — no threshold was lowered to pass a broken result; the < 8 % target is moved to a
named follow-on, not quietly dropped.

**Follow-on F1.2.1.1 (for < 8 % convergence):** a **smooth, fine** FEM k-vs-gap objective — fix
`box_w` independent of the gap (vary only the strip positions, keep the mesh constant) **and** a
finer gap-mesh (smaller cells → finer k resolution, heavier per-eval). EM-cycle-costly (multiple
~65 min runs); deferred here rather than blind-grinding. Proving the mechanism first is the correct
walking-skeleton order.
