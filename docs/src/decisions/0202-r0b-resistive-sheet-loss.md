# ADR-0202: R.0b — conductor loss: the resistive-sheet trace boundary

**Status:** Accepted
**Date:** 2026-07-07
**Related:** ADR-0194 (R.0 — the exact-at-f_ref single-frequency pattern this
follows), fdtd-205 (the volumetric Ohmic reference this deliberately does not
reuse at copper σ), RF-TOOL-ROADMAP R.0b.
**Spec:** `docs/superpowers/specs/2026-07-07-r0b-resistive-sheet-loss-design.md`

## Decision

Masked planar traces can now carry the **Leontovich resistive-sheet
relation** instead of the PEC clamp: `Materials::sheet_r_ohm` /
`MaterialsSpec::sheet_r_ohm` set `E_tan = R_s·K`, `K = ẑ×(H_above−H_below)`,
on every masked `E_x`/`E_y` edge after the E half-step —
`R_s = √(π f_ref μ₀/σ)` (`yee_voxel::surface_resistance_ohm`). Local,
explicit, dissipative; `R_s = 0`/`None` is the PEC clamp **bit-exactly**
(gate compute-017, which also pins dissipation with the correct sign via a
matched-absorber septum: a wrong sign self-amplifies). Ez masks (vias) stay
PEC. GPU rejects with `Unsupported` (Auto falls back) — WGSL port queued.

## Gate `engine-closs-001` (two release solves) — what it pins and why

The R.0 loss-shaped fixture (directional forward-wave ratio, reflections
drop out) at an **engineered σ = 5.8e4 S/m** (real copper's ~0.1 dB over
this span sits below the fixture's ±0.24 dB ripple):

1. **Linearity in R_s**: α(R_s)/α(R_s/2) = **2.124** (gate 2 ±10 %) — the
   sheet mechanics, immune to closed-form ambiguity.
2. **Absolute band**: α_meas/α_Pozar = **0.415** in the documented
   strip-only band [0.30, 0.60]. Decomposition: Pozar's `R_s/(Z₀W)` is the
   **strip + ground** first-order total (ground treated strip-width, ~half
   each); this sheet losses the strip only (the ground is a boundary face,
   not a mask — scope limit), and a single zero-thickness sheet dissipates
   on the **net** current `|H_b − H_a|²` where a thick two-faced conductor
   dissipates `|H_b|² + |H_a|²` — a further modest undercount. Both effects
   are physics of the declared model, recorded, not hidden.

## Consequences and queued follow-ons

Insertion-loss/Q predictions gain their first-order conductor term
(strip-side, at f_ref). Queued: ground-plane sheet (boundary-face variant —
recovers most of the remaining 0.5 share), the WGSL kernel, a two-faced /
frequency-dependent SIBC if broadband loss accuracy is ever the binding
constraint, and a real-copper validation on a high-Q resonator scenario.
