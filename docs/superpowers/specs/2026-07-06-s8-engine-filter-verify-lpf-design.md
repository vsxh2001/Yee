# S.8 / F1.3.0 — First engine-verified synthesized filter (stepped-impedance LPF)

**Date:** 2026-07-06
**Phase:** S.8 (ENGINE-STUDIO-ROADMAP) = F1.3.0 (FILTER-DESIGN-ROADMAP stage 6 walking
skeleton). Builds on S.5–S.7 (ADR-0182..0184).
**Plan:** `docs/superpowers/plans/2026-07-06-s8-engine-filter-verify-lpf.md`

## Problem

Everything before this verified the engine against *scenarios* (a line, a stub). The
project goal is designing filters: the missing demonstration is a filter that was
**synthesized by the pipeline** (spec → prototype → dimensions → layout) being
**verified by the engine** against its own design intent. That closed loop is stage 6
of the filter pipeline (F1.3), never yet exercised with full-wave EM on a synthesized
multi-section device.

## Design

**DUT choice: N = 5 Butterworth stepped-impedance low-pass, f_c = 2 GHz, FR-4**
(`yee_synth::prototype` → `yee_filter::dimension_stepped_impedance_layout`, shipped in
F1.2.3). Deliberately the easiest real filter for a first full-wave verify:

- Non-resonant (low Q) → short ring-down, the S.6 run length suffices.
- In-line straight sections → robust voxelization (no coupled-gap sensitivity).
- Loose gates are *meaningful*: cutoff position, passband flatness, stopband rejection
  are all closed-form-predicted by `ideal_response_lowpass`.

A hairpin/edge-coupled BPF is the *wrong* first target: its tap position (`qe`→tap) and
per-section gaps are explicitly deferred to F1.2.1, so a detuned narrowband passband
would tell us nothing about the verify machinery.

**Measurement**: exactly the S.6/S.7 two-run method. Reference = a straight Z₀ through
line spanning the same port-to-port extent on the same bbox/grid; DUT = the synthesized
layout. |S21| from the transmission ratio, |S11| from incident/reflected separation.
Known measurement limits (accepted, documented): dx = 0.3 mm quantizes the ~0.55 mm
high-Z sections to ~2 cells (impedance error → cutoff shift); feeds/junctions are not
de-embedded; band-edge ripple as in ADR-0183.

**Home**: `crates/yee-filter/tests/engine_lpf_verify.rs` — F1.3 is a filter-pipeline
stage, so the gate lives with the filter crate (dev-deps: `yee-engine`, `yee-voxel`).
CI: one new step in the `compute-engine-gates` job.

## Validation gate — engine-filter-verify-001

Over 0.8–4.2 GHz (drive band sized to cover it):

1. **Cutoff**: the measured −3 dB crossing of |S21| lies within **±20 %** of the
   designed f_c = 2 GHz (walking-skeleton band; the dominant error is the staircased
   high-Z width).
2. **Passband**: |S21| at 1 GHz (Ω = 0.5, ideal −0.013 dB) ≥ **−3 dB**.
3. **Stopband**: |S21| at 4 GHz (Ω = 2, ideal −30.1 dB) ≤ **−15 dB**.
4. The measured-vs-`ideal_response_lowpass` table is printed for the record.

## Non-goals

BPF verify (needs F1.2.1 qe→tap + per-section gaps first); spec-mask data structure +
pass/fail API (next slice, once this proves the loop); EM-in-the-loop optimization
(F1.2.1, consumes this machinery); Touchstone export of measured response.
