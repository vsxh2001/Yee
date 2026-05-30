# Phase 2.fdtd.6.10 — CW capacitor steady-state diagnostic — Design Spec

**ADR:** ADR-0127 · **Date:** 2026-05-31 · **Status:** Accepted (shipped)

## Problem

ADR-0126: F2.3 on the aperture port loads the line but no band-pass forms, and a
longer record makes the in-band loss WORSE — the shunt-tank capacitor looks like a
deepening near-short. Decide: a *measurement* limit (pulse DFT of an unsettled
transient on a short standing-wave line) vs a *cap-update windup bug*.

## Goal

A CW single-frequency steady-state diagnostic that settles the question by
measuring the aperture cap's own steady-state behaviour, isolated from line-de-embed
noise — and a recorded verdict (measurement-limit → F2.3 needs CW; or cap-bug → fix).

## Method

`tests/cap_cw_001.rs` (`#[ignore]`'d, release): drive the aperture capacitor — and
a shunt L‖C tank — with a CW sinusoid at f0 (Hann-ramped) for ~200 cycles. Two
probes:
- **Probe 1 (asserted):** the isolated cap arm in the same closed back-action loop
  the driver uses, no transmission line in the feedback — measures the arm's own
  `Z = V_T/I` (sliding window) and the `V_C` envelope.
- **Probe 2 (recorded):** the same port on a short PEC guide (the `aperture_port_001`
  harness) — expected to scatter (standing-wave de-embed noise).

Assert on Probe 1: reactance `−jX` every settled window; `|Z|` drift `< 0.10`;
`V_C` growth `< 0.10` (bounded, no windup); realized `|Z| ≈ β + 1/(jωC)` analytic
`< 0.10`; the L‖C tank resonates (`|Z| > single-arm |X|`).

## Changes (`crates/yee-fdtd/**` ONLY)

- New `tests/cap_cw_001.rs`. No `src` change (this is a diagnostic; the verdict
  determines whether a later `src` fix is needed).

## DoD

1. fmt + clippy -D warnings exit 0.
2. No regression (`aperture_port_001`, `lumped_rlc_twoway_001`, lumped/cpml gates).
3. `cap_cw_001` GREEN with genuine, tight assertions; the verdict recorded in
   ADR-0127.

## Outcome

VERDICT MEASUREMENT-LIMIT — the cap arm is correct (`Z = 97.9 − j100.1 Ω`, drift
0.003, `V_C` bounded, tank resonates). No `src` change. F2.3 needs a CW
per-frequency drive (F2.3-d). See ADR-0127.

## Out of scope

The F2.3 CW driver (F2.3-d); a CI gate job (on-demand diagnostic).
