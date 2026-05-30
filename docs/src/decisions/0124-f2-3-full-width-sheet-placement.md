# ADR-0124: Filter Phase F2.3-b — full-width-sheet lumped-element placement

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0115 (F2.3 lumped FDTD EM sim — the gate this targets), ADR-0119/
0121/0123 (the reactive-port benches: the port is ≈0.37-accurate per-element on a
*full-width sheet*, but F2.3 places elements on a *single cell*), ADR-0116 (the
two-way port), the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

F2.3's `fdtd_lumped_001` gives a **flat** |S21| ≈ 1.0 (no selectivity). The
reactive-port research (2.fdtd.6.5/6.6/6.8) established two *separate* facts:

1. The per-element reactive port is **≈0.37-accurate** (single-cell ε_eff limit) —
   imperfect but **not** N×-wrong (that was a measurement artifact, ADR-0121).
2. The de-embed **benches use a full-width transverse sheet** of lumped ports, and
   on a sheet the element **does couple** to the line. F2.3, by contrast, places
   each element on a **single `E_z` edge** (`cell_for(cx,cy,k_elem)`) under a
   multi-cell-wide microstrip trace — so the element is a tiny fraction of the
   line's admittance and is **≈ inert**, which is the dominant cause of the flat
   response (a placement/geometry issue, *separate* from the ≈0.37 port accuracy).

A multi-cell aperture *port reformulation* (the deferred multi-week brick) would
fix both, but the **cheaper test first** is to fix F2.3's placement: distribute
each lumped element across a full-width sheet so it actually loads the line, then
see whether the band-pass shape emerges within F2.3's **loose** tolerance — the
≈0.37 port accuracy may well be good enough for a loose-tol cross-validation.

## Decision

Change `yee_voxel::simulate_lumped_board` (the F2.3 driver) to place each ladder
element as a **value-distributed full-width sheet** of `LumpedRlcPort`s spanning
the trace cross-section at the element's x-location, instead of a single cell:

- count `N` = the transverse `E_z` edges across the trace (or the relevant
  port-face span) at that location;
- distribute so the parallel/series combination equals the intended element: a
  **shunt** branch's `C → C/N`, `L → N·L` per cell (N parallel cells sum to the
  element); a **series** branch's element split consistently across its sheet;
- keep `.with_two_way()` on each cell (the stable two-way port) and the existing
  drive/load ports.

Re-run `fdtd_lumped_001` (unchanged, loose tol) in the bounded container:

- **If the FDTD |S21| now reproduces the band-pass shape within the gate's loose
  tol** (in-band ≈ 0 dB within a few dB, ≥ ~20 dB stopband rejection) → the EM-sim
  component **ships** (F2.3 merges; lumped-LC goal 5/6). The ≈0.37 port accuracy is
  sufficient at loose tol; the multi-cell aperture port is **not needed**.
- **If it still does not meet the loose tol** → the ≈0.37 single-cell accuracy
  genuinely gates it → the multi-cell aperture port (the deferred brick) is
  required. Record the FDTD |S21| achieved; do **not** weaken `fdtd_lumped_001`.

## Consequences

**Ships (if it passes):** the goal's EM-simulation component — full-wave FDTD of
the lumped-LC filter board cross-validated against the analytic ladder, at loose
tol. With F2.0/F2.1/F2.2/F2.4 that completes the lumped engine (5/6; only the
maintainer-gated polished-UI merge remains).

**Gate:** `fdtd_lumped_001` GREEN (unchanged, loose tol) on the F2.3 branch before
merge; the existing lumped/CPML gates non-regressed. Never weakened.

**Not in scope:** the multi-cell aperture *port* (the FDTD-core fallback, only if
sheet placement is insufficient); a tight-tol EM match; SRF/ESR parasitics.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-3-full-width-sheet-placement-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-3-full-width-sheet-placement.md`.
- ADR-0123 Outcome (single-cell vs sheet); ADR-0121 (≈0.37 port accuracy);
  ADR-0115 (the F2.3 gate + driver).
