# FS.5c — Studio exposure of yield analysis (design)

**Date:** 2026-07-13 · **Track:** FS.5 (optimization maturity, `FULL-SUITE-ROADMAP.md`) · **ADR:** 0222

## Problem

FS.5a shipped the yield primitive (`yee_surrogate::yield_mc`, ADR-0211):
model-agnostic Monte-Carlo pass-counting over Gaussian dimension tolerances,
deterministic in the seed, Wilson 95 % CI. It is reachable from Rust and (via
`yee.surrogate`) Python, but a studio user has no way to ask "what fraction
of manufactured boards meets spec?" — the question every commercial suite
(HFSS Optimetrics, ADS yield) answers with a form and a button. FS.5c is the
walking skeleton of that form.

## Non-goals

- Yield over *engine-verified* responses (GP trained on engine samples — the
  ADR-0211 composition). The skeleton runs the same closed-form
  patch-resonance testcase the `surrogate-yield-001` gate certified; wiring a
  GP/engine-backed pass closure behind the identical command shape is the
  follow-on.
- Space-mapping exposure in the studio. Deferred as **FS.5c.1** (it needs a
  progress-streaming command like `verify_filter` — fine evals are engine
  solves — and a target-spec form; the yield skeleton does not).
- Correlated / non-Gaussian tolerances, importance sampling (FS.5a
  follow-ons; same API).

## Design

Follow the existing studio command/panel idiom exactly (ADR-0198/0203
pattern: `*_impl` in a module, thin `#[tauri::command]` in `lib.rs`, one
React panel in `App.tsx`, DOM gates in vitest).

### Command: `yield_estimate` (studio/src-tauri/src/yield_mc.rs)

Wraps `yee_surrogate::yield_estimate` around the ADR-0211 closed-form
patch-resonance testcase `f = c / (2 L √ε_eff)`, `ε_eff = (ε_r + 1)/2`.
The request is design-centric (the user thinks in f₀, not L): the nominal
patch length is derived as `L₀ = c / (2 f₀ √ε_eff)` and perturbed by
`σ_L`; ε_r is perturbed by `σ_εr`; a sample passes iff its resonance lands
within ±`spec_halfwidth_hz` of f₀.

- **Request** (`YieldRequest`): `f0_hz`, `eps_r` (serde default 4.4),
  `sigma_l_m`, `sigma_eps_r`, `spec_halfwidth_hz`, `n_samples`, `seed`.
- **Response** (`YieldResponse`): `yield_frac`, explicit Wilson bounds
  `ci95_lo` / `ci95_hi` (clamped to [0, 1] — the UI shows an interval, not a
  half-width), `n_pass`, `n_samples`, `length_nominal_m` (so the user sees
  the derived dimension the tolerance applies to).
- **Validation** is `Result<_, String>` like every studio command: f₀,
  spec half-width > 0; σ ≥ 0; ε_r > 1; 0 < n_samples ≤ 10⁷ (the command is
  synchronous — the closed form at 10⁷ samples is still sub-second, and the
  cap keeps a typo from freezing the shell).
- **Non-physical samples fail spec** rather than panic: a draw with
  `L ≤ 0` or `ε_eff ≤ 0` counts as a fail (unreachable at sane σ; the
  closure must still be total).
- `yee-surrogate` joins `studio/src-tauri/Cargo.toml` as a path dependency,
  mirroring `yee-engine`.

### Panel: `YieldPanel` (studio/src/App.tsx)

Numeric inputs with the ADR-0211 defaults translated to the studio's usual
UI units: f₀ 2.45 GHz, ε_r 4.4, σ_L 0.1 mm, σ_εr 0.05, spec ±40 MHz,
n = 10 000, seed 20260711 (the gate seed — deterministic ⇒ the same numbers
every run). Run button, error line, result line with yield %, the Wilson
95 % CI as `[lo, hi] %`, `n_pass / n_samples`, and the derived L₀ in mm.

## Gates

- **`studio-yield-dom-001`** (`studio/src/yield.test.tsx`): the panel
  renders its 7 numeric fields with the ADR-0211 defaults; clicking Run
  fires `invoke("yield_estimate", { req })` with the correctly
  unit-converted arguments (GHz→Hz, mm→m, MHz→Hz); the resolved response's
  yield, CI bounds, and pass count appear in the DOM. The tauri `invoke` is
  mocked via `vi.mock("@tauri-apps/api/core")` — first use of a command
  mock in the studio suite (earlier panel tests only exercise forms).
- **Rust unit tests in-module**: determinism (same request ⇒ identical
  response), the ADR-0211-regime yield at the gate seed (assert within
  (0.95, 0.99) — the derived-L₀ round-trip may differ from the gate's
  pinned 29 mm by an ULP, so the exact 0.9721 is not asserted), and input
  validation rejections.
- Existing studio gates stay green: `npm run build`, `npx vitest run`,
  `cargo check --manifest-path studio/src-tauri/Cargo.toml`.

## CI

The existing studio CI job already runs the three commands above; no
workflow change needed (out of this track's lane regardless — reported,
not edited).
