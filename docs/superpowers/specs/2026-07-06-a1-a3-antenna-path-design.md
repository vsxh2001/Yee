# A.1–A.3 — The antenna path: matching, radiation, design loop

**Date:** 2026-07-06
**Phases:** A.1 / A.2 / A.3 (ENGINE-STUDIO-ROADMAP Part 3), continuing A.0 (ADR-0190).
**Plan:** `docs/superpowers/plans/2026-07-06-a1-a3-antenna-path.md`

## A.1 — Inset-fed matched patch + directional S11

- **Synthesis** (`yee-layout`): `inset_fed_patch(f0, substrate, z0)`. Edge resistance
  from the Balanis slot-conductance model (`G₁` piecewise in `W/λ₀`; walking skeleton
  ignores the mutual term `G₁₂`, documented), inset depth
  `x₀ = (L/π)·acos(√(Z₀/R_edge))` from the `cos²` current profile. Geometry: metal as a
  union of four rectangles (two outer patch bands, the centre band beyond the inset,
  the feed running through the notch), notch gap = feed width.
- **Measurement upgrade** (`sparams::directional_reflection_db`): |S11| = |bwd|/|fwd|
  from ONE run's three-probe standing-wave fit — no reference run, no subtraction
  artifacts (the A.0 caveat), half the solve cost.
- **Gate engine-antenna-002**: dip within ±10 % of f₀ AND return loss ≤ −7 dB at the
  dip (the matched-ness A.0 could not assert; closed-form insets typically land
  −10…−20 dB — final number recorded on measurement).

## A.2 — Radiation pattern over the protocol

- **Per-face CPML** (`yee-compute`): `CpmlConfig` gains face-level enables
  (`with_faces([[x−,x+],[y−,y+],[z−,z+]])`); the antenna boundary is side walls +
  **open top** absorbing, PEC ground at z-min. Protocol: `BoundarySpec::Cpml` gains
  serde-defaulted `faces: Option<[[bool;2];3]>` (None → axes semantics).
- **NTFF on the protocol**: `JobSpec.ntff: Option<NtffSpec>` (surface margin,
  frequency, list of (θ, φ) directions) → `JobResult.far_field: Option<Vec<f64>>`
  (|E| per direction), wired through the E.5a host-adapter path on CPU.
- **Gate engine-antenna-003**: the patch's upper-hemisphere E-plane cut — broadside
  within a few dB of the pattern max and strongly above the horizon directions.
  Known approximation, documented: the NTFF box crosses the substrate (equivalence
  surface not fully in homogeneous space) — standard practice, qualitative asserts.

## A.3 — Antenna design loop

Seed from the **crude** closed form (no ΔL fringing correction → resonance a few %
high — a genuine model error, not an artificial detune), then the S.11/S.12 secant on
the synthesis frequency against the measured dip. Gate engine-antenna-004: final
resonance error ≤ 2 % and ≤ half the seed error.

## Non-goals

Gain/efficiency numbers; circular polarization; arrays; probe-fed (coax) patches;
GPU NTFF-protocol path (CPU first, GPU follows the aperture-port GPU work).
