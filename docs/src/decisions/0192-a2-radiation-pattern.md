# ADR-0192: A.2 — per-face CPML + far field over the protocol; the patch beams broadside

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0186 (per-axis CPML — this generalizes it), ADR-0177 (E.5a NTFF host
adapter — the pattern reused), ADR-0191 (A.1, whose match diagnosis motivated the open
top).

## Decision

1. **Per-face CPML** (`yee-compute`): `CpmlConfig` gains `faces: [[bool; 2]; 3]`
   (`with_faces`; `with_axes` sets both faces, so every existing config is bit-identical
   — all CPML gates stayed green untouched). A disabled face stays PEC. The antenna
   boundary is **absorbing side walls + open top, PEC ground**:
   `faces [[t,t],[t,t],[f,t]]`. GPU rejects face-asymmetric configs
   (`ComputeError::Unsupported`) until the WGSL mask grows face bits.
   Protocol: `BoundarySpec::Cpml` gains serde-defaulted `faces: Option<…>`.
2. **Far field over the protocol**: `JobSpec.ntff: Option<NtffSpec>` (frequency, box
   margin, optional grounded-antenna bottom face `k_min`, direction list) →
   `JobResult.far_field: Option<Vec<f64>>`. The transform is the validated
   `yee_fdtd::NtffState` (new additive constructor `with_bounds` for asymmetric boxes),
   sampled per step through the E.5a host-side grid adapter; `yee-engine` now depends on
   `yee-fdtd` for it. CPU-only (explicit `gpu` errors, `auto` falls back).

## Measured — gate `engine-antenna-003`

The A.1 inset patch under the open boundary, E-plane cut at its measured 2.425 GHz
resonance (one ~11-minute release solve with per-step NTFF):

| θ | φ=0° | φ=180° |
|---|---|---|
| 0° (broadside) | 0.0 dB | — |
| 20° | −3.3 dB | **+1.4 dB** |
| 40° | −6.9 dB | +0.5 dB |
| 60° | −10.1 dB | −2.9 dB |
| 80° | −18.2 dB | −11.7 dB |

A genuine upper-hemisphere patch beam, including the **textbook feed-side beam squint**
(the inset feed radiates too, tilting the beam a few degrees away from the feed — the
+1.4 dB at θ = 20°, φ = 180°). Asserts: broadside radiation captured; broadside beats
every θ ≥ 60° direction both sides of the cut. Documented approximation: the NTFF
bottom face (`k_min = 1`) hugs the ground and crosses the substrate.

## Consequences

The engine measures antennas end-to-end: geometry → S11 → radiation pattern, all over
the job protocol from any client. A.3 (the antenna design loop on the inset depth)
closes the track.
