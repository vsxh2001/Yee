# ADR-0195: R.2 — complex S-parameters + Touchstone export from engine measurements

**Status:** Accepted
**Date:** 2026-07-07
**Related:** RF-TOOL-ROADMAP R.2, ADR-0189 (the 3-probe standing-wave fit both new
observables build on), ADR-0194 (R.0/R.1, the prior track increments), Touchstone
round-trip gate in `yee-io` (Phase 0).
**Spec:** `docs/superpowers/specs/2026-07-07-r0-r2-losses-vias-sparams-design.md`

## Decision

Two complex-valued observables land in `yee_engine::sparams`, and the engine's
measurements now flow into the project's primary external interface, `.sNp` via
`yee_io::touchstone`.

- **`forward_transfer(triple_a, triple_b, dt, spacing, freqs) → Vec<(f64, f64)>`** —
  the complex plane-to-plane transfer `T(f) = fwd_B(f) / fwd_A(f)`. Both triples run
  the ADR-0189 standing-wave decomposition, so the backward wave drops out of both
  numerator and denominator; for a uniform line `T = e^{−jβ(f)d}` with d the
  plane-to-plane distance.
- **`complex_reflection(triple, dt, spacing, freqs) → Vec<(f64, f64)>`** — the complex
  `Γ(f) = bwd(f) / fwd(f)` at the triple's first probe (the reference plane). The
  existing `directional_reflection_db` is now its magnitude.

Complex values stay `(re, im)` tuples on the sparams API (no num-complex dep in the
engine crate); the gate converts to `Complex64` at the yee-io boundary.

## Gate `engine-sparams-002` (one release solve)

The S.5-certified lossless 6 λ_g FR-4 line on the S.9/S.10 stack (CPML-xy, two
aperture ports), two directional probe triples d ≈ 3 λ_g apart:

1. **Phase**: least-squares slope of the unwrapped `arg T(f)` over 4–6 GHz vs the
   closed form `dφ/df = −2π d √ε_eff / c` (Hammerstad–Jensen). Measured
   **−3.8816e−9 vs −3.7720e−9 rad/Hz → 2.91 %** (gate ±5 %).
2. **Magnitude**: |T| within ±1 dB of lossless unity across the band — worst
   deviation **+0.24 dB**.
3. **Touchstone**: the measured two-port of the plane-A→B line segment — S21 =
   S12 = T measured; S11 = S22 = 0 (a uniform line in its own reference impedance
   has zero reflection; see the negative result below) — written to `.s2p` and
   read back, data asserted to round-trip at writer precision.

## Negative result: Γ on a THRU is the fixture's, not the DUT's

The first gate fill used S11 = S22 = the measured complex Γ at plane A. It produced
a **read-back rejection**: `yee_io::touchstone::read` enforces λ_max(S†S) ≤ 1 + 1e−9
and the assembled matrix violated it by up to **+10.75 dB** (σ_max = 3.45). Probe-dump
forensics (41 bins, fit vs HJ-pinned-β fit) exonerated the standing-wave fit — fitted
β tracks Hammerstad–Jensen to < 1 %, fit residual ~1e−3, and pinning β to the closed
form reproduces the same Γ. The physics:

- |fwd_A| ripples ×5 with period ≈ v_p/2L (the full line length): **multi-bounce
  between the two imperfectly matched aperture ports**. `forward_transfer` is immune —
  the resonance factor is common to both planes and cancels in fwd_B/fwd_A (why |T|
  stays at +0.24 dB) — but the 9000-step window truncates the ring-down, breaking the
  steady-state |Γ| ≤ 1 identity exactly at the ripple minima (raw |Γ| up to 2.4).
- More fundamentally, on a **through** measurement Γ at plane A is the *measurement
  fixture's* load-port reflection, not the DUT's S11: the uniform line's
  own-reference S11 is 0. Stuffing fixture Γ into the DUT S11 slot builds a
  fictitious, rightly-rejected non-passive hybrid.

So the export fills S11 = S22 = 0 (exact for the DUT), `complex_reflection` carries an
interpretation caveat in its docs (for a **one-port** DUT — antenna, shorted stub — Γ
IS the DUT reflection, which is how A.1/A.3 use it), and **R.2b** is queued: measured
through-line S11 via de-embedding/calibration plus ring-down-complete windows.

## Passivity enforcement at the export boundary

With the correct fill the only over-unity content is |T|'s +0.24 dB measurement
ripple. The standard treatment of raw measured data applies: per-frequency
**passivity enforcement** — for the symmetric reciprocal fill `S = [[a, b], [b, a]]`
the singular values are `|a ± b|` in closed form, so any sample with σ_max > 1 is
scaled by 1/σ_max. The gate logs the worst correction and **asserts it ≤ 0.5 dB**,
so the export can never silently reshape a bad measurement into a plausible file.
The physics asserts (1)–(2) run on the *raw* measurement; only the exported file is
touched.

Rejected alternative: weakening the yee-io passivity check. That check is the
Touchstone reader's defense against garbage S-data from any source; measurement
noise at the writer is the caller's problem, exactly as the yee-io docs state.

## Consequences

The engine now produces industry-standard deliverables: any measured response can
leave the tool as a valid, passive `.s2p` for consumption by ADS/AWR/scikit-rf.
Complex Γ and T also unlock de-embedded input impedance (`Z_in = Z₀(1+Γ)/(1−Γ)`)
for matching work in later phases. R.2b (measured through S11) queued behind R.3
GPU parity.
