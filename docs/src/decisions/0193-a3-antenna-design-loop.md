# ADR-0193: A.3 — the antenna design loop delivers a −25.7 dB match

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0191 (the model gap this closes), ADR-0192 (the open boundary),
ADR-0189 (the loop discipline from the filter track).

## Decision and measured result

Gate `engine-antenna-004` (`antenna_patch_match_loop.rs`, five release solves, ~22 min,
own CI job with the other antenna gates): coarse scan of the inset depth + one
neighbour-midpoint refinement, each evaluation one engine job measuring the single-run
directional |S11| dip. The measured map:

| inset depth | dip | return loss |
|---|---|---|
| 0.10·L | 2.450 GHz | −6.5 dB |
| 0.20·L | 2.475 GHz | −13.2 dB |
| **0.25·L (refined)** | **2.475 GHz** | **−25.7 dB** |
| 0.30·L | 2.475 GHz | −9.3 dB |
| 0.40·L (≈ closed-form seed) | 2.425 GHz | −0.9 dB |

The map is cleanly unimodal with its match point at **0.25·L** — not the G₁-model's
0.396·L — and the loop turns the closed-form seed's unusable −0.9 dB into a
**−25.7 dB match** (better than a typical −15 dB commercial spec), with the matched
resonance 1.0 % from the 2.45 GHz design. Asserts bind the quality with ~10 dB
headroom: best ≤ −15 dB, improvement over seed ≥ 10 dB, resonance ±10 %.

## What this closes

The antenna track A.0–A.3 is complete: closed-form synthesis (A.0/A.1), S-parameters
and radiation pattern over the job protocol under a physically-open boundary (A.2), and
a design loop that measurably beats the closed forms (A.3) — the engine designs
antennas end-to-end, the same loop discipline the filter track proved in S.11/S.12.
Follow-ons: multi-knob antenna optimization via `yee-surrogate`, gain/efficiency
numbers, probe-fed and array geometries, GPU support for the CPU-only pieces
(aperture ports, per-face CPML, protocol NTFF) on the nightly track.
