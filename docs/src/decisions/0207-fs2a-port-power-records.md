# ADR-0207: FS.2a — aperture-port power records: the gain denominator

**Status:** Accepted
**Date:** 2026-07-08
**Related:** ADR-0206 (FS.1b — flagged that every pattern gate is
relative-|E| only), S.10/ADR-0187 (the aperture port whose state this
records), R.3/ADR-0196 (the GPU `Unsupported` idiom reused).
**Spec:** `docs/superpowers/specs/2026-07-08-fs2-farfield-products-design.md`

## Decision

`AperturePortSpec::record` (serde default `false`) → per-step
`(v_src, v_terminal, i_branch)` triples in `JobResult::port_records`. The
port loop already computes all three every step (the LumpedRlcPort
correction), so recording is free of new field probes and modal-
normalization guesswork. CPU-only; the GPU backend rejects recording
ports with `Unsupported` (Auto falls back) until a readback buffer lands.

## Measured lesson: account on the circuit side, not the aperture side

The first gate run summed the naive aperture-side `v_term·i` and read a
non-physical B/A energy ratio of **1.596**: the port's implicit β term
(`β = dt·h/2ε₀A` ≈ **14.5 Ω** on this fixture — comparable to the 50 Ω
branch!) sits between the terminal voltage and the resistor, so `v_term·i`
mixes real transfer with reversible discretization storage. The honest
identity uses **circuit-side quantities only**: EMF supply `Σ −v_src·i·dt`
vs branch-resistor dissipation `Σ i²R·dt`.

## Gate `engine-power-001` (release, 1 solve ≈ 55 s) — GREEN

Lossless matched through line, both ports recording:

- **closure = 0.9917** (pinned [0.95, 1.0]): E_emf = 5.557e-13 J supplied;
  A's resistor 2.708e-13 + B's resistor 2.802e-13 dissipated; 0.8 % left
  in CPML leakage + residual ring.
- accepted-by-field `E_emf − E_R(A)` = **51.3 %** of the supply — the
  textbook matched-source halving, and the FS.2b gain denominator.
- passive-branch sample-wise sanity `v_term·i ≥ 0` (a load only
  dissipates).

## Queued

FS.2b: audit the NtffState |E| normalization, `farfield::gain_dbi`,
patch 5–8 dBi textbook window + the 2×1-vs-single array-gain
differential. FS.2c: efficiency + full-sphere export. GPU port-record
readback.
