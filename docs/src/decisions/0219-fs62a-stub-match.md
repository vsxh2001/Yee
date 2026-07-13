# ADR-0219: FS.6.2a — single-stub matching synthesis

**Date:** 2026-07-12 · **Status:** accepted · **Track:** FS.6
**Spec:** `docs/superpowers/specs/2026-07-12-fs62-stub-match-design.md`

## Decision

`yee_layout::single_stub_match(gamma, beta)` — the Smith-chart
single-shunt-open-stub construction in closed form (rotate to the g = 1
circle: `Re[Γ(d)] = −|Γ|²`; cancel `+jb` with `tan(βl) = −b`; smallest
non-negative lengths mod λ_g/2; smaller-d crossing returned). Lives in
yee-layout because its consumers are layout generators (the stub becomes
trace geometry); it consumes the Γ that `yee_engine::sparams::
complex_reflection` measures (referenced at probe 0 per
`fit_standing_wave`).

## Verification

Gate `stub-match-001` (instant, GREEN first run): Pozar Example 5.2
position d = 0.1104 λ (published 0.110 λ); and the **machine contract**
— for 96 loads across the passive Γ-disk (|Γ| ∈ [0.1, 0.8], 12 phases),
the synthesized (d, l_open) pair nulls the combined reflection below
1e-9. The contract gate was chosen over pinning textbook stub *lengths*:
the published tables mix open/short-stub branches (memory of Example 5.2
produced the shorted-stub length), while the null contract is
self-verifying and branch-independent.

## FS.6.2b — the full-wave loop (SHIPPED; gate `match-em-001` GREEN)

The FS.6 roadmap gate: an edge-fed 2.45 GHz patch (measured
|Γ| ≈ 0.6, −4.50 dB) is matched by a stub synthesized from its
**measured** Γ — matched |Γ| = 0.282 (−11.00 dB), **improvement
6.49 dB** (bar ≥ 6 dB), 667 s release, in the antenna CI job.

It took three runs, and the failure chain is the real lesson
(**measurement-plane hygiene near resonant radiators**):

1. Run 1 (+4.30 dB, short of bar): judging plane P at 3 mm read
   |Γ_P| = 0.636 vs the plane-invariant 0.464 — the aperture port's
   evanescent near-zone. Also the free-β standing-wave fit at plane A
   (12 mm from the patch edge) read **β = 107.4 rad/m = the
   bulk-substrate velocity** (107.7), not the line's 93.7: the patch's
   resonant near-field dominates there.
2. Run 2 (−0.27 dB, *worse*): "fixing" only β (HJ closed form) while
   keeping the near-field-corrupted Γ_A broke the partial error
   cancellation of run 1 — a clean β with a dirty Γ is worse than a
   consistently dirty pair.
3. Run 3 (GREEN): `sparams::fit_standing_wave_known_beta` — the wave
   split with β **known** (overdetermined least squares; its residual
   flags non-β contamination; the unit gate caught a Cramer conjugate
   swap on the first cut) — plus Γ measured at plane P (46 mm from the
   patch; constrained-fit residuals 0.004–0.01 at both planes) and
   λ_g/2-periodic stub placement into the feasibility window.

Standing rule extracted: **synthesize from a Γ measured on a plane far
from resonant structures, with the wave split constrained to the line's
known β, and check the fit residual** — a free-β fit near a radiator
locks onto whatever wave dominates locally and silently corrupts both
outputs.
