# ADR-0186: S.9 per-axis CPML on the job protocol — the ADR-0185 collapse root-caused

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0185 (recorded the anomaly), ADR-0176 (E.1 CPML + `with_axes` on
`yee-compute`), ADR-0182 (the protocol this extends).
**Spec:** `docs/superpowers/specs/2026-07-06-s9-per-axis-cpml-design.md`

## Root cause

The ADR-0185 finding — all-face CPML collapsing the LPF's |S21| below −3 dB across the
whole band — was a **scenario error, not a CPML defect**. On the dx = 0.3 mm voxel stack
the substrate is ~5 cells tall and the trace sits at k_top ≈ 5, while npml = 10: the
entire microstrip line sat **inside the 10-layer z-min absorber** and propagated ~76 mm
through it. The absorber did exactly its job, on the wrong cells.

## Decision

1. `BoundarySpec::Cpml` gains `axes: [bool; 3]` (`#[serde(default)]` = all-on, so the
   pre-S.9 wire format still parses), passed through as
   `CpmlConfig::for_spec(..).with_axes(axes)` — the per-axis support has existed in
   `yee-compute` since E.1; the protocol just never exposed it.
2. **Board-level open boundary = `[true, true, false]`**: absorbing side walls, PEC
   ground and lid. `engine-filter-verify-001` adopts it.
3. Scenario gates certified under the PEC box (engine-verify-001, engine-sparams-001)
   keep their boundary — their recorded numbers stay reproducible.

## Measurements (same LPF scenario, all else identical)

| boundary | passband mean @1 GHz | stopband @4 GHz | rejection | cutoff |
|---|---|---|---|---|
| PEC box | +3.42 dB | −27.19 dB | 30.6 dB | 1.900 GHz |
| CPML all faces | (collapsed: < −3 dB everywhere) | — | — | — |
| **CPML x/y, PEC z** | **+1.32 dB** | **−32.91 dB** | **34.2 dB** | **1.900 GHz** |

CPML-xy wins on every aggregate. Honest residual: band-edge ripple remains (+12.4 dB at
0.8 GHz, transition band still corrupted) — with the side-wall cavity modes gone, the
remaining standing waves come from the **lumped-port mismatch** (a single-cell resistor
across a ~5-cell substrate is a poor 50 Ω match). Matched terminations / port
de-embedding is now clearly the next fidelity lever, and it is a measurement problem,
not a boundary problem.

## Consequences

Any protocol client can now request per-axis absorption. The antenna track's radiation
boundary (open *top*, PEC *ground* — a per-**face** split on the z axis) still needs a
`yee-compute` extension; deferred to that track.
