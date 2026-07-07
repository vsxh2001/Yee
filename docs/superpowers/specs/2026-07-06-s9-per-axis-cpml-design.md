# S.9 — Per-axis CPML on the job protocol (root-causing the ADR-0185 collapse)

**Date:** 2026-07-06
**Phase:** S.9 (ENGINE-STUDIO-ROADMAP). Resolves the boundary finding recorded in
ADR-0185.
**Plan:** `docs/superpowers/plans/2026-07-06-s9-per-axis-cpml.md`

## Problem and hypothesis

ADR-0185 recorded that all-face CPML (npml = 10) collapsed the LPF's |S21| below −3 dB
across the whole band. The geometry supplies the hypothesis before any experiment: on
the dx = 0.3 mm voxel stack the substrate is ~5 cells tall and the trace sits at
k_top ≈ 5 — **entirely inside the 10-layer z-min CPML region**. The quasi-TEM mode then
propagates ~76 mm *through an absorber*, attenuating everything including the passband.
Nothing is wrong with the CPML; the boundary condition was wrong for the scenario.

The correct board-level open boundary is **side walls absorbing, z faces PEC** (the
ground plane at k = 0 is PEC anyway; the lid stays a lid). `yee_compute::CpmlConfig`
has supported per-axis enables (`with_axes`) since E.1 — the job protocol just never
exposed them: `BoundarySpec::Cpml { npml }` hard-codes all axes.

## Design

- `BoundarySpec::Cpml` gains `#[serde(default = all-true)] axes: [bool; 3]` — old JSON
  (no `axes` key) deserializes to today's all-face behaviour; `run_job` passes it via
  `CpmlConfig::for_spec(..).with_axes(axes)`. One field, no other protocol change.
- **Experiment** (the actual root-cause test): re-run `engine-filter-verify-001` with
  `axes: [true, true, false]`. Predicted outcome if the hypothesis holds: passband
  recovers to ≈ 0 dB, ripple shrinks versus the PEC box (side-wall cavity modes gone),
  cutoff/stopband unchanged in position.
- If confirmed **and** cleaner than the PEC box, the LPF gate switches to CPML-xy with
  tightened asserts (absolute passband bound), and the ADR-0185 "do not use CPML"
  caveat is narrowed to "do not put the board inside the absorber".

## Gates

- Fast: serde — `Cpml { npml }` JSON without `axes` still parses (backward compat);
  round-trip with explicit axes.
- Heavy: the re-measured `engine-filter-verify-001` numbers (recorded in ADR-0186)
  under whichever boundary wins.

## Non-goals

Per-*face* CPML (open top lid + PEC ground on the same axis) — the right shape for
antenna radiation, but it needs a `yee-compute` change; deferred to the antenna track.
