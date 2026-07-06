# S.6 — S-parameters on the engine (walking skeleton): two-run transmission over the job protocol

**Date:** 2026-07-06
**Phase:** S.6 (ENGINE-STUDIO-ROADMAP), building directly on S.5 (ADR-0182).
**Plan:** `docs/superpowers/plans/2026-07-06-s6-engine-sparams.md`

## Problem

S.5 put voxelized layouts on the job protocol; the filter verify (F1.3: simulated response
vs spec mask) additionally needs |S21|(f) — a scalar transmission spectrum — out of engine
jobs. Nothing in the protocol or the engine computes spectra today; the studio's TS
`dftMagnitude` is client-side display math, not a reusable Rust API.

## Design

### Method: two-run transmission ratio (reference / DUT)

The classic FDTD S-parameter recipe (Sheen et al. 1990, adapted to lumped ports):

- **Reference run:** the bare feed line, drive port (50 Ω, modulated-Gaussian EMF) at one
  end, *passive* 50 Ω port (`v0 = 0` — the resistive-port update with zero EMF is exactly
  a lumped resistor load) at the other, one E_z probe under the trace near the load.
- **DUT run:** identical everything + the device (here: a λ/4 open stub) inserted mid-line.
- `|S21|(f) = |DFT(probe_dut)(f)| / |DFT(probe_ref)(f)|` — feed-line loss, launch
  discontinuity, probe coupling, and (to first order) termination mismatch divide out.

Both runs are ordinary S.5 `JobSpec`s (materials + dt over the protocol, CPU backend);
no engine-core (`yee-compute`) changes at all.

### New API: `yee_engine::sparams` (pure post-processing)

- `single_bin_dft(series, dt_s, f_hz) -> (re, im)` — the same single-bin DFT the ε_eff
  gates use, as a reusable function.
- `transmission_db(dut, reference, dt_s, freqs_hz) -> Vec<f64>` — 20·log₁₀ magnitude
  ratio per frequency.

Pure functions on `JobResult.probes` series: usable identically by the studio (via a
future Tauri command), Python, WS clients, and tests. No new protocol fields.

## Validation gates

- **engine-sparams-001** (`yee-engine/tests/sparams_stub_notch.rs`, `#[ignore]`, release
  CI): FR-4 microstrip line (the S.5-certified stack: W = 3 mm, h = 1.6 mm, ε_r = 4.4,
  dx = 0.3 mm) with a **quarter-wave open-circuited stub** — the textbook bandstop
  (Pozar, *Microwave Engineering*, transmission-line theory): |S21| notches where the
  stub is λ/4. Stub length is sized `L_s = λ_g/4 − ΔL` with the Hammerstad open-end
  correction ΔL, predicting the notch at f₀ = 5 GHz from closed forms only. Assert the
  measured |S21|(f) minimum sits within **±15 %** of the prediction (the pipeline's
  walking-skeleton band), the notch is ≥ 8 dB deep, and the passband edges stay shallow
  (a genuine dip, not broadband loss).
- **Fast, non-ignored** unit tests on `sparams`: `single_bin_dft` recovers the amplitude
  and phase of a known sinusoid over integer periods; `transmission_db` of a half-scaled
  copy is −6.02 dB across the band.

## Non-goals

S11 / full 2-port extraction (needs incident/reflected separation at the drive plane —
next slice); spec-mask overlay + pass/fail (F1.3 proper); on-GPU DFT of probe series
(host-side is microseconds); protocol changes (none needed).
