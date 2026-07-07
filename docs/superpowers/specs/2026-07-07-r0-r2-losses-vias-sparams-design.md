# R.0–R.2 — Losses, vias, complex S-parameters (RF-TOOL-ROADMAP)

**Date:** 2026-07-07
**Plan:** `docs/superpowers/plans/2026-07-07-r0-r2-losses-vias-sparams.md`

## R.0 — Dielectric loss

FDTD σ is frequency-flat; a `tan δ` material is mapped at a reference frequency:
`σ = 2π f_ref ε₀ ε_r tan δ` (the standard single-frequency approximation, exact at
`f_ref`, documented). New `yee_voxel::substrate_sigma_cells(&model, tan_d, f_ref)` maps
every substrate cell (ε_r > 1) — additive helper, no signature breaks; the S.5
`MaterialsSpec.sigma_cells` field carries it, and the E.1 lossy CA/CB update consumes it.

**Gate engine-loss-001**: one run of the 6 λ_g line at `tan δ = 0.05` (large enough for
~3.7 dB over the measurement span — robust SNR), two directional probe triples ~3 λ_g
apart; `α_meas = ln(|fwd_A|/|fwd_B|)/d` at f₀ vs Pozar §3.199
`α_d = k₀ ε_r (ε_eff−1) tan δ / (2 √ε_eff (ε_r−1))`, walking-skeleton ±20 %. A lossless
control assert (α ≈ 0) is free from the same machinery on the S.5 gate's scenario — the
gate asserts the measured α is dominated by the dielectric term.

Conductor loss is **deliberately excluded** (R.0b): skin depth at GHz is ~2 µm versus
the 0.3 mm cell — a volumetric σ_Cu would be wrong by construction; the honest
implementation is a surface-impedance boundary condition on masked cells.

## R.1 — Vias

`yee_voxel::with_via(&mut model, x_m, y_m)`: a vertical PEC column of `E_z` edges
`k = 0..k_top` at the nearest cell column — attaches/extends `pec_mask_ez` (the protocol
already carries it; the engine already clamps it). No `Layout` schema change (additive
model-level helper; a first-class `Layout.vias` field waits until the export writers
also learn vias, so the schema changes once).

**Gate engine-via-001** (differential, 3 runs on one grid): the S.6 λ/4 open-stub
scenario, plus a variant with a via at the stub's far end. A shorted λ/4 stub is an
open circuit at its input — the 5 GHz notch must **vanish** (|S21| recovers) while the
no-via control keeps it deep. Same reference run serves both.

## R.2 — Complex S-parameters + Touchstone

The standing-wave fit already yields complex forward/backward phasors and β per
frequency. R.2 adds: complex `s21(f)`/`s11(f)` with **de-embedded reference planes**
(rotate phasors by `e^{±jβΔ}` from the probe plane to the port plane), and a writer
path: engine measurement → `yee_io` Touchstone `.s2p`.

**Gates**: (a) through-line complex S21 phase vs `−β·l` with β from Hammerstad–Jensen
ε_eff across the band (±5 % on unwrapped phase slope); (b) `.s2p` write→read round-trip
of an engine-measured two-port equals the in-memory values (the existing `yee-io`
fidelity gate extended to engine data).

## Non-goals (this slice)

R.0b surface impedance; multilayer stackups; via pads/antipads in export writers;
renormalization to non-50 Ω references.
