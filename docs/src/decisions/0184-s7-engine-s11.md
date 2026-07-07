# ADR-0184: S.7 S11 on the engine — incident/reflected separation

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0183 (S.6 |S21| + the two-run method this extends; named this follow-on).
**Spec:** `docs/superpowers/specs/2026-07-06-s7-engine-s11-design.md`

## Context

Return loss is half of every filter spec; S.6 delivered only |S21|. At a single probe
near the drive, incident and reflected waves superpose and must be separated.

## Decision

The separation falls out of the two-run method at **zero extra solve cost**: probe the
port-1 reference plane in both runs. The runs share launch, feed line, and grid, so the
reference run's P1 series *is* the incident wave at that plane;
`reflected(t) = dut_p1(t) − ref_p1(t)` isolates the device-caused reflection, and
`|S11|(f) = |DFT(reflected)| / |DFT(ref_p1)|`. New pure function
`sparams::reflection_db` (panics on length mismatch — the series must come from the
same job pair). Known second-order residual, accepted at walking-skeleton tolerance:
the reflected wave re-reflects off the imperfectly matched drive port and crosses the
plane again.

## Gate (extension of engine-sparams-001; same two release solves)

On the λ/4 open-stub scenario, measured over the job protocol:
**|S11| at the notch = −0.93 dB** (|S11| ≈ 0.90 — the stub reflects nearly everything
at resonance, gate ≥ −4 dB) and **|S11|² + |S21|² = 0.807** at the notch (gate: the
physical band [0.5, 1.3]; the ~19 % energy deficit is the un-subtracted second-order
re-reflection plus band-edge ripple, consistent with the +8.7 dB standing-wave ripple
already documented in ADR-0183). Fast unit gate: a synthetic `dut = 1.25·incident`
yields −12.04 dB across the band.

## Consequences

Both halves of a filter response — |S21|(f) and |S11|(f) — now come out of two engine
jobs through pure protocol-side post-processing. The F1.3 spec-mask verify (insertion
loss + return loss vs mask) has everything it needs. Follow-ons: complex S-parameters
with de-embedded reference planes, Touchstone export of engine-measured responses
(`yee-io` already writes `.s2p`), and driving this from the studio/Python surfaces.
