# ADR-0222: FS.5c — studio exposure of yield analysis

**Date:** 2026-07-13 · **Status:** accepted · **Track:** FS.5 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-13-fs5c-studio-yield-design.md`

## Context

FS.5a (ADR-0211) shipped the yield primitive —
`yee_surrogate::yield_estimate`, a model-agnostic pass-closure Monte-Carlo
over Gaussian tolerances, deterministic in the seed (splitmix64 +
Box-Muller), Wilson 95 % CI — certified by gate `surrogate-yield-001` on a
closed-form patch-resonance testcase. It was reachable from Rust and
Python but not from the studio, where the commercial-parity question
("what fraction of manufactured boards meets spec?") is actually asked.

## Decision

Walking skeleton, exactly the existing studio command/panel idiom
(ADR-0198/0203):

1. **Command `yield_estimate`** (`studio/src-tauri/src/yield_mc.rs`,
   thin `#[tauri::command]` in `lib.rs`; `yee-surrogate` added as a path
   dependency next to `yee-engine`). The request is **design-centric**:
   `f0_hz`, `eps_r` (default 4.4), `sigma_l_m`, `sigma_eps_r`,
   `spec_halfwidth_hz`, `n_samples`, `seed`. The nominal patch length is
   derived from the ADR-0211 closed form, `L₀ = c/(2 f₀ √ε_eff)`,
   `ε_eff = (ε_r+1)/2`; a sample passes iff its resonance lands within
   ±half-width of f₀. Response carries `yield_frac`, **explicit clamped
   Wilson bounds** `ci95_lo`/`ci95_hi` (the UI shows an interval, not a
   half-width), `n_pass`, `n_samples`, and the derived `length_nominal_m`
   so the user sees what dimension σ_L perturbs. Synchronous (no
   progress events): the closed form is instant even at the 10⁷-sample
   validation cap. Non-physical draws (L ≤ 0, ε_eff ≤ 0) fail spec rather
   than panic — the pass closure is total.
2. **`YieldPanel`** in `studio/src/App.tsx` with the ADR-0211 defaults in
   UI units (2.45 GHz, ε_r 4.4, σ_L 0.1 mm, σ_εr 0.05, spec ±40 MHz,
   n = 10⁴, seed 20260711 — the gate seed, so the panel shows the same
   numbers every run), an Estimate-yield button, and a result line
   (yield %, CI [lo, hi] %, n_pass/n, L₀ mm).
3. **Gate `studio-yield-dom-001`** (`studio/src/yield.test.tsx`): default
   form values, exact invoke arguments (GHz→Hz, mm→m, MHz→Hz
   conversions), rendered yield/CI, and the error path. This is the
   studio suite's **first command-mocking test** — `@tauri-apps/api/core`
   is replaced via `vi.hoisted` + `vi.mock`, pinning the request shape
   the Rust command deserializes (earlier panel tests only exercised
   forms). Rust-side unit tests pin determinism, the ADR-0211-regime
   yield at the gate seed (range-asserted, not bit-pinned: the derived
   L₀ round-trips through f₀ and may differ from the gate's 29 mm by an
   ULP), and validation rejections.

## Deferred

- **FS.5c.1 — space-mapping exposure in the studio.** Needs a
  progress-streaming command (fine evals are engine solves, minutes) and
  a target-spec form; explicitly out of this skeleton.
- Yield over engine-verified responses (GP trained on engine samples
  behind the same command shape) — the ADR-0211 composition.

## Consequences

- The studio now answers the Optimetrics-class question end-to-end with
  numbers a gate can pin (deterministic seed ⇒ reproducible screenshots
  and support conversations).
- The command shape (tolerances in, yield + Wilson interval out) is the
  contract any future GP/engine-backed pass closure slots behind without
  touching the panel.
