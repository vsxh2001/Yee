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
| **R.2** | **Complex 2-port S-matrix + Touchstone export**: `sparams::forward_transfer` (complex plane-to-plane `T = fwd_B/fwd_A`, backward wave AND the port-to-port multi-bounce resonance factor cancel in the ratio) + `sparams::complex_reflection` (complex Γ = bwd/fwd — the DUT's reflection for a **one-port** DUT, per the A.1/A.3 usage); measured through-line exported `.s2p` via `yee_io::touchstone` (S11 = S22 = 0, exact for a uniform line in its own reference impedance) with **per-frequency passivity enforcement** (σ_max closed-form `max\|a±b\|` for the symmetric fill; correction asserted ≤ 0.5 dB) | `engine-sparams-002`: unwrapped-phase slope vs HJ closed form `−2π d √ε_eff/c` — **2.91 %** (±5 % gate); worst \|T\| **+0.24 dB** (±1 dB); `.s2p` write→read round-trip at writer precision | **SHIPPED** (ADR-0195) |
| **R.2b** | **Measured through-line S11**: de-embedding/calibration so the exported S11 is the DUT's, not the fixture's (the ADR-0195 negative result: on a THRU, plane-A Γ is the load-port's reflection, and window-truncated multi-bounce pushes raw \|Γ\| to 2.4 — needs ring-down-complete windows or shorter lines + port calibration) | complex Γ of a known one-port (via-shorted line: Γ = −1 rotated by 2βl) vs TL theory, magnitude AND phase | queued |
| **R.3** | **GPU parity for the design flows**: aperture ports (per-port column-reduction kernel), per-face CPML mask bits, and (optionally) protocol NTFF in WGSL; certified differentially vs CPU on llvmpipe, perf on the GPU nightly | differential gates ≤ drift-class error vs the bit-exact CPU path; nightly perf numbers recorded (the 20× dGPU target) | queued |
| **R.4** | **BPF end-to-end + multi-knob BO**: per-section hairpin gaps + qe→tap (the deferred F1.2.1 core), then `yee-surrogate` BO driving engine jobs over 2+ knobs | a synthesized BPF verified full-wave against its coupling-matrix response; BO closes centre-frequency + bandwidth to spec | queued |
| **R.5** | **Studio spec→loop→export**: the Tauri studio drives spec entry → design loop (live progress over the existing job events) → response/pattern plots → Touchstone/Gerber export | studio e2e test drives a scripted loop; export artifacts byte-checked | queued |

*Last updated: 2026-07-07 — R.2 SHIPPED (complex T/Γ observables; engine-measured `.s2p` with passivity enforcement at the export boundary; phase slope 2.91 % vs HJ; negative result: THRU-line Γ is the fixture's, not the DUT's → R.2b queued) — ADR-0195, after R.0/R.1 (ADR-0194). Next: R.3 GPU parity for the design flows.*
