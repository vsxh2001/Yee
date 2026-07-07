# S.10 — Distributed (substrate-column) ports for board-level verify

**Date:** 2026-07-06
**Phase:** S.10 (ENGINE-STUDIO-ROADMAP). The "next fidelity lever" named in ADR-0186.
**Plan:** `docs/superpowers/plans/2026-07-06-s10-aperture-ports.md`

## Problem

A single-cell resistive port bridges one E_z edge of a ~5-cell substrate: the other
four cells of the trace→ground column are plain field, so the port's effective source /
load impedance seen by the quasi-TEM mode is far from its nominal R. ADR-0186 attributed
the LPF gate's residual ±12 dB band-edge ripple to exactly this mismatch. The same
lesson was learned in the lumped-LC track: `yee-voxel::lumped_sim` moved to aperture
ports spanning the full substrate height (ADR-0125/0126).

## Design

**No protocol change.** `Drive.ports` (and therefore `JobSpec.ports`) already accepts
any number of ports. A distributed port is a **series stack** of single-cell resistive
ports on the E_z column under the trace, `k = 0 .. k_top`: each cell carries `R/N` and
EMF `V₀/N`, so the column totals R and V₀ while forcing the port voltage across the
full substrate — the openEMS-style lumped port, expressed in existing primitives.

New helper `yee_engine::column_port_specs(i, j, k_top, r_ohm, v0, f0_hz, bw_hz,
t0_steps) -> Vec<PortSpec>` so clients don't hand-roll the split; unit-tested for the
series-sum property (Σ 1/`resistance` per-cell relations, EMF split, cell coverage).

## Validation

- Fast: helper unit test (N cells, R/N + V₀/N each, k = 0..k_top coverage).
- Heavy: `engine-filter-verify-001` re-measured with column drive + column load —
  adopt if better on the recorded aggregates (passband mean closer to 0 dB, band-edge
  ripple down, rejection preserved). The measured table goes in ADR-0187 either way.

## Non-goals

Full (y, z) aperture (j-band × column) — the next step if the column alone is not
enough; port impedance renormalization; protocol changes.
