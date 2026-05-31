# ADR-0139: App.2.2 — Low-pass stepped-impedance flow in the studio

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0137 (F1.2.3 `dimension_stepped_impedance`, the low-pass engine being
surfaced), ADR-0136 (the recommender recommends `SteppedImpedance` — now routable),
ADR-0138 (App.2.1 hairpin lighting — the band-pass lighting precedent), ADR-0120 (the
lumped parallel flow — the response-path template), the product vision
(`docs/superpowers/specs/2026-05-31-ideal-filter-design-app-vision.md` §5),
[[lumped-lc-and-studio-redesign]]

---

## Context

The studio is band-pass-only (edge-coupled, hairpin, lumped). The stepped-impedance
**low-pass** dimensioner shipped (F1.2.3, Pozar-§8.6-gated) but is un-surfaced, and the
App.2.0 recommender's `SteppedImpedance` recommendation can only route to a band-pass
stand-in. The single biggest remaining capability gap is a **low-pass response class**
end-to-end in the visible app.

Both engine pieces already exist and are validated: the closed-form low-pass magnitude
response (`yee_filter::lowpass_s21_squared`, used by the band-pass `ideal_response` —
for low-pass, evaluate at `Ω = f/f_c`) and `dimension_stepped_impedance`. The lumped
flow (ADR-0120) is a proven parallel-response-path pattern in the studio.

## Decision

Surface a **low-pass stepped-impedance flow** in the studio, mirroring the lumped
parallel flow:

- `yee-filter`: a small public `ideal_response_lowpass(approx, order, cutoff_hz, freqs)`
  (the low-pass analogue of `ideal_response`, reusing `lowpass_s21_squared` at `Ω=f/f_c`).
- `yee-studio-web`: `Topology::SteppedImpedance` + a `SteppedLowpassDesigned` /
  `design_stepped_from` (mirror `LumpedDesigned` / `design_lumped_from`) carrying the
  stepped sections (`dimension_stepped_impedance`), the swept low-pass `|S21|` vs a
  low-pass mask + PASS/FAIL, and the board layout; `stepped_synthesis_stage` +
  `stepped_layout_stage`; a stepped rail; `StageCanvas` `stepped_flow` routing; the Spec
  form made low-pass-aware (Cutoff label, no FBW, `Response::Lowpass`) for the stepped
  technique; the gallery card lit and the recommender mapped to `Live`.

## Consequences

**Ships:** the first **low-pass** capability end-to-end in the visible app —
SteppedImpedance is a live, routable technique driven by the real F1.2.3 dimensioner +
the low-pass response; the recommender's low-pass recommendation now routes to a real
flow; Gerber/KiCad export of the stepped-Z board. No new physics — integration of two
validated engines via the proven lumped pattern.

**Gates:** (1) `yee-filter` `ideal_response_lowpass` — strong, non-vacuous, textbook:
Butterworth `|S21(f_c)|` = −3.01 dB, monotone rolloff, the `−20·N·log10` stopband
asymptote; Chebyshev equiripple edge. (2) Studio — `dx build` EXIT 0 + a non-vacuous
host test (`design_stepped_from` → real stepped sections + a low-pass `|S21|` ≈ −3 dB at
cutoff, proving the card routes to the real engine, not a stub) + band-pass / lumped
flows unregressed.

**Larger increment:** a parallel flow + a low-pass Spec mode. The plan's escape hatch
guarantees a clean partial — the `yee-filter` API + gate and the `SteppedLowpassDesigned`
engine + its non-vacuous test land even if the stage UI must be deferred. No gate
weakened; the EM-verify wall (ADR-0133) is untouched.

**Not in scope:** combline / interdigital; high-pass / band-stop; elliptic low-pass; EM
verify.

---

## References
- `crates/yee-filter/src/lib.rs` (`ideal_response`, `lowpass_s21_squared`,
  `dimension_stepped_impedance`); `crates/yee-studio-web/src/{engine.rs, stages.rs, main.rs}`
  (the lumped parallel-flow template).
- `docs/superpowers/specs/2026-05-31-app-2-2-studio-lowpass-stepped-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-2-studio-lowpass-stepped.md`.
