# ADR-0185: S.8 / F1.3.0 — first engine-verified synthesized filter (stepped-impedance LPF)

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0182..0184 (the verify machinery), ADR-0108 (validated microstrip FDTD
stack), FILTER-DESIGN-ROADMAP stage 6 (F1.3).
**Spec:** `docs/superpowers/specs/2026-07-06-s8-engine-filter-verify-lpf-design.md`

## Context

Everything before this verified the engine against *scenarios* (a line, a stub). The
project goal is designing filters: the missing demonstration was a filter **synthesized
by the pipeline** being **verified by the engine** against its own design intent — the
closed loop. DUT: N = 5 Butterworth stepped-impedance LPF, f_c = 2 GHz, FR-4
(`yee_synth::prototype` → `dimension_stepped_impedance_layout`, all shipped closed
forms) — deliberately the easiest real filter (non-resonant, straight sections);
a hairpin/edge-coupled BPF would be detuned by the F1.2.1 gaps (`qe`→tap, per-section
gaps) and prove nothing about the verify machinery.

## Measurement and what was measured

S.6/S.7 two-run method; the reference is a Z₀ through line on the DUT's bbox (identical
grid). Gate `engine-filter-verify-001` (`yee-filter/tests/engine_lpf_verify.rs`,
`#[ignore]`, its own CI step; dev-deps `yee-engine`/`yee-voxel` keep the library
WASM-safe). Measured, PEC box, 2 × ~80 s release solves:

- **Cutoff: −3 dB at 1.900 GHz vs the designed 2.0 GHz — 5.0 % error** (gate ±20 %;
  found by a sustained-crossing scan so a ripple dip cannot fake it).
- **Rejection: passband(1 GHz) − stopband(4 GHz) = 30.6 dB** (+3.42 dB vs −27.19 dB;
  ideal Butterworth 30.1 dB; gate ≥ 20 dB). Stopband at 4 GHz: −28.5 dB vs ideal −30.1.
- **Ripple, bounded not hidden:** single-probe standing-wave ripple in the PEC box
  reaches +17.8 dB at the 0.8 GHz band edge and corrupts the transition band
  (+6.8 dB at 2.5 GHz vs ideal −10.1). The gate asserts *relative* quantities (cutoff
  position, rejection) that survive ripple, plus an explicit ±6 dB bound on the
  passband mean (measured +3.42 dB).

## Boundary finding (recorded, deferred)

Both boundaries were tried. **CPML (npml = 10, all faces) collapsed |S21| below −3 dB
across the entire band including the passband — non-physical for this DUT** — while the
PEC box produced the correct filter shape. The microstrip-into-CPML interaction on this
voxel stack (dielectric slab + ground mask meeting the absorber) needs its own
investigation before CPML is used for board-level verify; do not silently switch this
gate to CPML. Cleaning up the measurement (matched terminations / de-embedding /
directional separation at the transmission plane) is the highest-value follow-on for
verify fidelity.

## Consequences

The design→verify loop is closed end-to-end for the first time: spec → g-values →
dimensions → layout → voxelize → engine jobs → measured response vs design targets.
F1.3 proper (spec-mask data structure + pass/fail API + BPF verify after F1.2.1) and
F1.2.1 (EM-in-the-loop dimension refinement) can now both consume this machinery.
