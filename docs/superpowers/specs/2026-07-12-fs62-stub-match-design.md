# FS.6.2 — single-stub matching: synthesis (a) + full-wave verification (b)

**Date:** 2026-07-12 · **Track:** FS.6 (network algebra / matching)
**Plan:** `docs/superpowers/plans/2026-07-12-fs62-stub-match.md`

## FS.6.2a — synthesis (SHIPPED, ADR-0219)

`yee_layout::single_stub_match(gamma, beta) -> StubMatch { d_m, l_open_m, b }`
— the classic Smith-chart construction in closed form: rotate Γ toward
the generator to the `g = 1` circle (`Re[Γ(d)] = −|Γ|²`), cancel the
residual `+jb` with an open stub (`tan(βl) = −b`); both lengths reduced
mod λ_g/2, the smaller-`d` crossing returned. Normalized admittance
throughout, so Z₀ cancels.

Gates (instant, GREEN first run): `pozar_example_5_2_stub_position`
(d = 0.1104 λ vs the published 0.110 λ) and the machine contract —
96 loads over the passive Γ-disk (|Γ| 0.1–0.8, 12 phases) each null the
combined reflection below 1e-9.

## FS.6.2b — full-wave verification (queued; the FS.6 roadmap gate)

"A match synthesized from a **measured** antenna Γ improves its
**measured** S11." Design sketch for the next session:

- DUT: edge-fed patch at 2.45 GHz (A.0 topology — genuinely mismatched,
  |Γ| ≈ 0.5–0.7 at resonance), but with a **long feed** (~55 mm): the
  synthesized stub lands up to λ_g/2 ≈ 33.6 mm behind the Γ reference
  plane, which must stay on-line and clear of the port.
- Measure complex Γ at f₀ with a 3-probe triple
  (`sparams::complex_reflection`; Γ referenced at probe 0 — the
  `fit_standing_wave` convention), β from the same fit.
- Synthesize (d, l_open); regenerate the layout with the shunt open stub
  at `x_ref − d` (the S.6 stub topology); re-measure.
- Gate `match-em-001`: unmatched |S11(f₀)| ≥ −6 dB (the mismatch is
  real), matched |S11(f₀)| improves by ≥ 6 dB (pin from the first green
  run — the stub also radiates/couples, so expect degradation vs the TL
  ideal). Use `snap_edges` (ADR-0218): d and l are continuous outputs.

## Non-goals

Lumped-element matching (needs engine lumped components — the
`fdtd_lumped_001` path is yee-fdtd only), double-stub, wideband match.
