# ADR-0194: R.0+R.1 — dielectric loss (0.1 % vs Pozar) and vias on the board path

**Status:** Accepted
**Date:** 2026-07-07
**Related:** RF-TOOL-ROADMAP (track opened per the 2026-07-07 gap assessment),
ADR-0182 (the sigma_cells protocol field), ADR-0176 (the E.1 lossy CA/CB update),
ADR-0189 (the directional observable reused for the loss measurement).
**Spec:** `docs/superpowers/specs/2026-07-07-r0-r2-losses-vias-sparams-design.md`

## R.0 — Dielectric loss

`yee_voxel::substrate_sigma_cells(&model, tan_d, f_ref)` maps the substrate loss
tangent to per-cell conductivity, `σ = 2π f_ref ε₀ ε_r tan δ` (FDTD σ is
frequency-flat — exact at the design frequency, the standard single-frequency
approximation elsewhere; documented). It rides the existing S.5 `sigma_cells`
protocol field into the E.1 lossy update — no engine change.

**Gate `engine-loss-001`** (one release solve): 6 λ_g FR-4 line at tan δ = 0.05, the
attenuation measured **loss-shaped** — two directional probe triples ~3 λ_g apart,
`α = ln(|fwd_A|/|fwd_B|)/d`, so reflections and the backward wave drop out entirely.
Measured **α = 4.328 Np/m vs the Pozar §3.199 closed form 4.323 Np/m — 0.1 % error**
(3.71 dB over the span). Gate ±5 %.

**Conductor loss is deliberately excluded** (R.0b): skin depth at GHz (~2 µm) is
unresolvable at dx = 0.3 mm; a volumetric σ_Cu would be wrong by construction. The
honest implementation is a surface-impedance boundary on masked cells — its own phase.

## R.1 — Vias

`yee_voxel::with_via_at_cell(&mut model, i, j, k_top)` — a vertical PEC column of
`E_z` edges ground→trace, attaching `pec_mask_ez` (the protocol already carries the
mask; the engine already clamps it). Cell-indexed placement because the voxel model
does not store its board-frame origin; a first-class `Layout.vias` schema field waits
until the export writers learn vias too, so the schema changes once.

**Gate `engine-via-001`** (differential, three release solves on one grid): the S.6
λ/4 open stub notches at **−45.1 dB @ 4.875 GHz** (control — deeper than the S.6-era
−36.8 dB, itself evidence of the S.9/S.10 fidelity stack); with a via shorting the
stub's far end the same frequency reads **−1.3 dB** — the notch annihilated, exactly
the transmission-line prediction (a shorted λ/4 stub is an open circuit at its input).

## Consequences

Real board physics — substrate loss and ground vias — now runs over the job protocol.
Insertion loss and Q predictions become possible for lossy designs; grounded stubs,
shorted resonators (combline!), and via-fenced structures become expressible. Next:
R.2 complex S-parameters + Touchstone export.
