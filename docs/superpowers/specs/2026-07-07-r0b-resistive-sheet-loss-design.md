# R.0b — conductor loss: resistive-sheet boundary on masked trace edges

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.0b
**Related:** ADR-0194 (R.0 dielectric loss — the same exact-at-f_ref
single-frequency pattern), fdtd-205 (`ohmic_skin_depth.rs` — the validated
volumetric Ohmic update this deliberately does NOT reuse at copper σ).

## Problem

Skin depth in copper at GHz (~1 µm) is unresolvable at board cell sizes
(0.1–0.3 mm): a volumetric σ_Cu layer self-shields wrongly by construction
(the roadmap's standing warning). Without conductor loss, insertion-loss/Q
predictions are dielectric-only.

## Approach — walking skeleton: a resistive sheet at f_ref

The masked planar trace (pec_mask_ex/ey) currently enforces `E_tan = 0`.
The resistive-sheet boundary replaces that with the Leontovich sheet
relation at the design frequency:

`E_tan = R_s · K`,  `K = ẑ × (H_above − H_below)`

i.e. per Yee edge, after the H half-step:
`E_x[i,j,k] = R_s·(H_y[i,j,k−1] − H_y[i,j,k])`,
`E_y[i,j,k] = R_s·(H_x[i,j,k] − H_x[i,j,k−1])`,
with `R_s = √(π f_ref μ₀ / σ)` (surface resistance; the microstrip current
is bottom-face-dominant, so one R_s — not R_s/2 — matches the Pozar
convention). Local, explicit, dissipative; `R_s = 0` degenerates to the PEC
mask bit-exactly. Frequency-flat like R.0's dielectric map: exact at f_ref,
first-order elsewhere; documented.

**Scope limits (documented):** planar z-normal sheets only (the ex/ey mask
use in every board flow); Ez masks (vias) stay PEC; the ground plane (a
boundary face, not a mask) stays lossless — its share of α_c is a follow-on
(boundary-face sheet). GPU: `Unsupported` rejection (Auto falls back), the
R.3-era pattern; the WGSL port is queued.

## Gates

- **compute-017** (unit, yee-compute): with a masked sheet in a driven
  cavity, field energy decays strictly below the PEC-mask run; `R_s = 0`
  reproduces the PEC run bit-exactly. A wrong sheet sign self-amplifies —
  this gate catches it structurally.
- **engine-closs-001** (release, yee-engine): the R.0 loss-shaped fixture
  (driven line, two directional probe triples, α from the forward-wave
  ratio) with tan δ = 0 and the sheet ON, vs **Pozar's strip conductor
  loss** `α_c = R_s/(Z₀·W)`. The gate runs at an engineered
  σ = 5.8e4 S/m (R_s ≈ 0.58 Ω) because real copper's α_c over this span
  (~0.1 dB) sits below the fixture's measured ±0.24 dB ripple — the sheet
  mechanics and the √(f/σ) scaling are what the gate pins; real-copper
  validation needs a high-Q resonator scenario (follow-on). Tolerance set
  from the first honest run (edge crowding makes FDTD read above the
  uniform-current formula).

## Consequences

The board path gains its last first-order loss term at walking-skeleton
fidelity; `yee_voxel::surface_resistance_ohm` gives designers the f/σ map.
Follow-ons: ground-plane sheet, GPU kernel, frequency-dependent SIBC
(vector-fitted convolution) if a use case demands broadband loss accuracy.
