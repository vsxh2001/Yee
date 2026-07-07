# RF Tool Roadmap (R.*) — from design demonstrations to a full RF design tool

**Date opened:** 2026-07-07
**Companion to:** `ENGINE-STUDIO-ROADMAP.md` (the engine + design loops this builds on),
`FILTER-DESIGN-ROADMAP.md`, `ROADMAP.md`.

The engine now *designs* — filters (S.8–S.12: synthesized → verified → refined to 1.0 %
cutoff error) and antennas (A.0–A.3: matched to −25.7 dB by measurement). This track
closes the gaps between "demonstrates design" and "a tool people design real hardware
with", in the order assessed on 2026-07-07: **(1) losses + vias, (2) complex S-matrix +
Touchstone export, (3) GPU parity for the design flows, (4) BPF end-to-end with
surrogate BO, (5) the studio spec→loop→export flow.**

Conventions as everywhere in this repo: every phase ships behind a machine-checkable
validation gate against a strong reference; walking-skeleton first; ADRs for decisions.

| Phase | Scope | Gate | Status |
|-------|-------|------|--------|
| **R.0** | **Dielectric loss on the board path**: substrate `tan δ` → per-cell σ (`σ = 2π f_ref ε₀ ε_r tan δ`, the standard single-frequency loss map — FDTD σ is frequency-flat, documented) riding the existing S.5 `sigma_cells` protocol field; `yee_voxel::substrate_sigma_cells` helper | `engine-loss-001`: measured dielectric attenuation α (forward-wave ratio between two directional probe triples ~3 λ_g apart) vs the **Pozar closed form** — **α = 4.328 vs 4.323 Np/m, 0.1 %** (±5 % gate) | **SHIPPED** (ADR-0194) |
| **R.0b** | Conductor loss: skin depth (~1–2 µm) is unresolvable at dx = 0.3 mm, so a volumetric σ_Cu is wrong by construction — needs a **surface-impedance boundary** on masked cells (own multi-step phase; do NOT fake it with bulk σ) | Ohmic attenuation of a line vs the closed-form α_c; the yee-fdtd skin-depth machinery is the reference | queued (after R.2) |
| **R.1** | **Vias**: vertical PEC columns ground→trace (`yee_voxel::with_via`, extends/attaches `pec_mask_ez`; protocol already carries the mask) | `engine-via-001` (differential, 3 runs): control notch **−45.1 dB @ 4.875 GHz**; with the via the same frequency reads **−1.3 dB** — the notch annihilated per TL theory | **SHIPPED** (ADR-0194) |
| **R.2** | **Complex 2-port S-matrix + Touchstone export**: `sparams` returns complex S11/S21 with **de-embedded reference planes** (phase-shift by the fitted β from the standing-wave fit — already measured per frequency); `yee-io` Touchstone writer wired to engine measurements | phase gate: through-line S21 phase vs β·l from Hammerstad–Jensen ε_eff across the band; `.s2p` write→read round-trip of an engine-measured response | queued |
| **R.3** | **GPU parity for the design flows**: aperture ports (per-port column-reduction kernel), per-face CPML mask bits, and (optionally) protocol NTFF in WGSL; certified differentially vs CPU on llvmpipe, perf on the GPU nightly | differential gates ≤ drift-class error vs the bit-exact CPU path; nightly perf numbers recorded (the 20× dGPU target) | queued |
| **R.4** | **BPF end-to-end + multi-knob BO**: per-section hairpin gaps + qe→tap (the deferred F1.2.1 core), then `yee-surrogate` BO driving engine jobs over 2+ knobs | a synthesized BPF verified full-wave against its coupling-matrix response; BO closes centre-frequency + bandwidth to spec | queued |
| **R.5** | **Studio spec→loop→export**: the Tauri studio drives spec entry → design loop (live progress over the existing job events) → response/pattern plots → Touchstone/Gerber export | studio e2e test drives a scripted loop; export artifacts byte-checked | queued |

*Last updated: 2026-07-07 — R.0 SHIPPED (dielectric loss, 0.1 % vs Pozar) and R.1 SHIPPED (vias; the λ/4 notch annihilated by a shorting via) — ADR-0194. Next: R.2 complex S-matrix + Touchstone export.*
