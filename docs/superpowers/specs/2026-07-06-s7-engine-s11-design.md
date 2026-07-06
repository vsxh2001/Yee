# S.7 — S11 on the engine: incident/reflected separation (walking skeleton)

**Date:** 2026-07-06
**Phase:** S.7 (ENGINE-STUDIO-ROADMAP), the follow-on ADR-0183 named.
**Plan:** `docs/superpowers/plans/2026-07-06-s7-engine-s11.md`

## Problem

S.6 delivers |S21|(f); a filter verify also needs |S11|(f) (return loss is half of every
filter spec). At a single probe near the drive, incident and reflected waves superpose —
they must be separated.

## Design

The separation falls out of the two-run method S.6 already uses, at zero extra solve
cost: probe the port-1 reference plane in BOTH runs. The runs share launch, line, and
grid, so the reference run's P1 series **is** the incident wave at that plane; the DUT
run's P1 series is incident + device-reflected. Then

`reflected(t) = dut_p1(t) − ref_p1(t)`,  `|S11|(f) = |DFT(reflected)| / |DFT(ref_p1)|`.

New pure function `sparams::reflection_db(dut_p1, ref_p1, dt_s, freqs)`; series lengths
must match (they are two probes on the same jobs). Residual error: the device-reflected
wave re-reflects off the imperfectly matched drive port and passes P1 again — a
second-order term, accepted at walking-skeleton tolerance.

## Validation gates

- **engine-sparams-001 extended** (same two release solves): the λ/4 open stub at
  resonance reflects nearly everything, so at the measured notch frequency assert
  **|S11| ≥ −4 dB** (strong reflection) and the lossless-DUT energy sanity
  **|S11|² + |S21|² ∈ [0.5, 1.3]** at the notch (the band that survives the known
  band-edge ripple; the measured value is recorded in the ADR).
- **Fast, non-ignored**: `reflection_db` of a synthetic `dut = ref + 0.25·ref` is
  −12.04 dB across the band.

## Non-goals

Full complex 2-port de-embedding (phase reference planes, port impedance
renormalization); multi-mode ports; spec-mask overlay (F1.3, next).
