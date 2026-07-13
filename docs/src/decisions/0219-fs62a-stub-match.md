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

## Queued

FS.6.2b — the full-wave loop (measured Γ → synthesized stub → measured
improvement), designed in the spec; the FS.6 roadmap gate.
