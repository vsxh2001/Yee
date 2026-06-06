# ADR-0170: Top-C finite-Q response + route the studio lumped response to the chosen topology (T6/T7)

**Status:** Accepted (T6/T7 — close the ADR-0169 display/export coherence gap)
**Date:** 2026-06-06
**Related:** ADR-0169 (T5 — studio auto-routes the JLCPCB export to top-C, but the design view + response stay
ladder-based), ADR-0160 (`ladder_s21_lossy` finite-Q + Cohn gate `lumped-q-001`), ADR-0165 (T1 — `top_c_s21`
lossless), ADR-0163 (studio finite-Q overlay), [[project-filter-design-final-goal]].

---

## Context

After T5, when the studio's orderable JLCPCB export auto-routes to **top-C** (the narrow-band rescue), the
studio still shows the **alternating-ladder** response + board — honestly badged ("may differ from the
displayed ladder"), but the displayed response does NOT match the manufactured top-C board. The lumped flow
is tied to `LumpedLadder` because the response uses `ladder_s21` / `ladder_s21_lossy` (ADR-0160/0163) and
top-C has only the **lossless** `top_c_s21` (T1) — no finite-Q variant. Closing the coherence gap needs
(T6) a top-C finite-Q response, then (T7) routing the studio's lumped response to the auto-chosen topology.

## Decision

**T6 (this brick — engine, `yee-filter`):** `top_c_s21_lossy(net: &TopCNetwork, f_hz, z0_ohm, q_unloaded) ->
Complex64`. Identical ABCD cascade to `top_c_s21`, but each **shunt resonator** carries its unloaded-Q loss
as a parallel conductance — `Y = G + jωC + 1/(jωL)`, `G = ω₀·C/Q_u` (`ω₀ = 2π·net.f0_hz`) — exactly mirroring
`ladder_s_params_lossy`'s `Shunt` branch. The `N+1` **series coupling caps stay lossless** (resonator-only Q,
like the ladder; a separate coupling-cap `Q_c` is a documented follow-on). Guard `q_unloaded ≤ 0` or `+∞` →
`inv_q = 0` → `G = 0` → **bit-identical to `top_c_s21`** (the lossless cascade). Re-export from `lib.rs`.

**Gate `top-c-q-001`** (`crates/yee-filter/tests/`, pure-compute, non-`#[ignore]`'d, non-circular), mirroring
`lumped-q-001`:
1. **Lossless limit:** `top_c_s21_lossy(net, f, z0, f64::INFINITY)` equals `top_c_s21(net, f, z0)` bit-identical
   over a frequency sweep (and for `Q_u ≤ 0`).
2. **Cohn dissipation:** for a representative top-C BPF at `Q_u = 100`, the realized midband insertion loss
   `−20·log₁₀|S21(f₀)|` matches Cohn's `IL₀ ≈ 4.343·Σgₖ/(Q_u·FBW)` (Σ over the reactive prototype elements
   `g₁..g_N`) within the narrowband tolerance (≤ 15 %). Cohn's formula is prototype-based / topology-independent
   (it depends only on the g-values, `Q_u`, FBW), so it is the correct independent reference for top-C too.
   Non-circular (Cohn from Σg/Q_u/FBW; the IL from `top_c_s21_lossy`).

**T7 (the follow-on — studio, `yee-studio-web`):** route the lumped flow's response + (ideally) board view to
the auto-chosen `BoardTopology`. When the orderable realization is top-C, compute the ideal response from
`top_c_s21` and the finite-Q overlay from `top_c_s21_lossy` (reuse the ADR-0163 overlay infra); when it is the
ladder, keep `ladder_s21[_lossy]`. This makes the displayed response match the manufactured board — closing
the coherence gap. T7 is its own brick (separate lane/ADR-update once T6 lands).

## Consequences

- T6 gives the missing top-C finite-Q response — a clean, gated engine brick that unblocks T7.
- After T7 the studio's lumped response is topology-correct (matches the orderable export), closing the
  ADR-0169 coherence gap.
- Scope T6: `crates/yee-filter/src/{top_c.rs, lib.rs}` + `crates/yee-filter/tests/`. Pure-math, WASM-safe.
- **Not in scope (T6):** the studio routing (T7); a separate coupling-cap `Q_c` (resonator-only Q mirrors the
  ladder); a lossy `(S11,S21)` pair (T6 provides `top_c_s21_lossy → S21`, matching `top_c_s21`'s shape — a
  full pair is a trivial later add if needed).

## References
- Pattern: `yee_filter::lumped::{ladder_s_params_lossy, ladder_s21_lossy}` (ADR-0160), gate `lumped-q-001`
  (Cohn `4.343·Σg/(Q_u·FBW)`).
- Extends: `yee_filter::top_c::{top_c_s21, TopCNetwork (f0_hz, z0_ohm, shunt, coupling_caps_farad)}` (ADR-0165).
- Consumed by T7: `yee-studio-web` lumped response (ADR-0163 overlay), the ADR-0169 `orderable_upload` topology.
