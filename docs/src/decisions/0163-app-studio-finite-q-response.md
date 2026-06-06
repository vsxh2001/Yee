# ADR-0163: App — surface the finite-Q (realistic) response in the Dioxus studio

**Status:** Accepted
**Date:** 2026-06-06
**Related:** ADR-0160 (finite-Q lumped response `ladder_s21_lossy` + Cohn gate), ADR-0130 (Dioxus
`yee-studio-web` is THE studio), ADR-0143 (App.2.6 response-overlay), [[project-filter-design-final-goal]]
(the app is the final deliverable; "see the parameters at the end").

---

## Context

The live Dioxus filter studio (`yee-studio-web`, deployed at `/Yee/studio/`) computes its response via
`yee_filter::ideal_response` / `ladder_s21` — the **lossless** closed-form curve (flat-top, 0 dB midband
insertion loss). ADR-0160 shipped `yee_filter::ladder_s21_lossy` — the **realistic finite-Q** lumped
response (per-resonator unloaded-Q loss; validated against Cohn's dissipation formula to 2.2 %) — but it
is NOT surfaced in the studio. A user designing on the lumped track sees only the ideal curve, not the
insertion loss / rounded corners a built filter actually exhibits. Surfacing the realistic response is
directly the app's "see the parameters at the end" value, and `ladder_s21_lossy` is pure-math / WASM-safe.

## Decision

Add a **finite-Q realistic-response overlay** to the studio's lumped flow (the response/verify view):

- **Engine (`engine.rs`):** in the lumped path (`synthesize_lumped` → `ladder_s21`), add a finite-Q sweep
  `ladder_s21_lossy(&ladder, f, q_unloaded)` over the same grid, plus a `q_unloaded` field on the lumped
  view (default `Q_u = 100`, a realistic chip-component value). Expose the realized midband insertion loss
  (the headline number a VNA measures).
- **UI (`stages.rs`):** a `Q_u` control (numeric input / slider, e.g. 20–300) on the lumped synthesis (or
  verify) stage, and the realistic finite-Q curve overlaid on the ideal lossless curve in the response
  plot (reuse the App.2.6 `overlay_curves` / `response_overlay` infrastructure). The ideal stays the
  reference; the finite-Q curve is the realistic companion. Label the midband insertion loss.
- Distributed (coupling-matrix `ideal_response`) flow is unchanged — `ladder_s21_lossy` is the lumped
  ladder's response; finite-Q surfacing is the lumped track's feature.

## Validation

- **Engine unit test (`yee-studio-web`, fast, non-`#[ignore]`'d):** for the 3-pole 0.5 dB Cheb at
  `Q_u = 100`, the engine's finite-Q midband insertion loss matches Cohn's `4.343·Σg/(Q_u·FBW)` (≈ 1.86 dB)
  within the narrowband tolerance (≤ 15 %), AND the ideal sweep stays ≈ 0 dB midband. Non-circular (Cohn
  from Σg/Q_u/FBW; the sweep from `ladder_s21_lossy`). Mirrors `lumped-q-001` at the studio-engine layer.
- The existing **`wasm-build` CI job** must stay green (the studio + the new UI compile to WASM;
  `ladder_s21_lossy` is pure-math, no new non-WASM dep). clippy/fmt clean.

## Consequences

- The deployed studio shows the realistic insertion-loss response on the lumped track — a concrete
  fidelity upgrade for the deliverable, consuming the validated ADR-0160 work.
- Scope: `crates/yee-studio-web/src/{engine.rs,stages.rs}` (+ `svg.rs` if the overlay needs a tweak).
  Small App-track feature; this ADR is the design record.
- **Not in scope:** distinct `Q_L`/`Q_C` (ADR-0160 single `Q_u`), the distributed-flow response (lumped
  feature), a finite-Q full-wave EM curve (a different track).

## References
- Realistic response: `yee_filter::ladder_s21_lossy` (ADR-0160); Cohn `4.343·Σg/(Q_u·FBW)`.
- Studio: `crates/yee-studio-web/src/{engine.rs,stages.rs,svg.rs}`; overlay infra from ADR-0143 (App.2.6).
