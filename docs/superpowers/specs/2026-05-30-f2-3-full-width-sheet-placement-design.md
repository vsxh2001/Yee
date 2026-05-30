# Filter Phase F2.3-b — full-width-sheet lumped-element placement — Design Spec

**ADR:** ADR-0124 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

F2.3's `fdtd_lumped_001` is flat (no selectivity). Root cause (ADR-0123 Outcome):
the driver places each lumped element on a **single `E_z` edge** under a
multi-cell-wide microstrip trace → the element is ≈ inert. The reactive port
itself is ≈0.37-accurate per-element (not N×-wrong) and **does couple when placed
as a full-width sheet** (the de-embed benches show this). Fix the placement, not
(yet) the port.

## Goal

Make F2.3's lumped elements actually load the line by distributing each across a
full-width sheet, and determine — at F2.3's existing loose tolerance — whether the
band-pass shape emerges (→ EM-sim ships) or the ≈0.37 port accuracy genuinely
gates it (→ the deferred multi-cell aperture port is needed). Do not weaken the gate.

## Method

In `yee_voxel::simulate_lumped_board` (the F2.3 driver, `crates/yee-voxel/src/
lumped_sim.rs`), where each ladder element is currently one
`LumpedRlcPort`(s) at `cell_for(cx, cy, k_elem)`:

- Determine the transverse span `N` of the trace at the element's x-location (the
  `E_z` edges across the trace width in y at `k_elem`, from the voxel model /
  layout trace width). Use the model's known trace-width-in-cells.
- Emit `N` parallel `LumpedRlcPort`s (one per transverse cell `j` across the trace)
  with **value-distributed** parameters so the sheet equals the intended element:
  - shunt **capacitor** (`l=0`): `C/N` per cell (N parallel → C);
  - shunt **inductor** (`c=∞`): `N·L` per cell (N parallel → L);
  - **series** R-L-C branch: split so the sheet presents the intended series Z (a
    series element spanning the broken-line gap across the width → per-cell
    `N·L`, `C/N`, `N·R` so the parallel sheet sums to the intended series arm —
    derive + document; if a series sheet is awkward, a single centred series cell
    may remain, but the dominant shunt resonators must be sheets).
  - keep `.with_two_way()` and `SERIES_ESR_OHM` on each cell.
- The drive/load matched-resistor ports are unchanged.

Re-run `fdtd_lumped_001` (unchanged) in the bounded container; read the |S21|
sweep.

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `simulate_lumped_board` (and any helper) emit the value-distributed sheet per
  element. Update the module doc (the series/shunt decomposition section) to
  describe the sheet.
- The branch must first be brought up to current `main` (it predates the canonical
  port / per-axis CPML) — merge `main` in; `Cargo.lock` conflicts → take `--theirs`,
  `cargo check`, commit (CLAUDE.md §5). Do NOT edit `yee-fdtd`/`yee-filter`.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-voxel --all-targets -- -D
   warnings` exit 0.
2. `fdtd_lumped_001` (unchanged, loose tol) run in the container — REPORT the
   |S21| sweep. Two acceptable outcomes, both honest:
   - **GREEN** (band-pass within loose tol: in-band ≈ 0 dB ±few, ≥ ~20 dB stopband)
     → EM-sim ships; proceed to merge F2.3.
   - **still failing** → record the achieved |S21| (how close), conclude the
     ≈0.37 port accuracy gates it → the multi-cell aperture port is required.
     Do NOT weaken the gate to force GREEN.
3. No regression to other yee-voxel gates (`fdtd_line_eeff_001`, voxel gates) —
   `--include-ignored`, release, as feasible in the container.

## Out of scope

The multi-cell aperture port (FDTD-core fallback); tight-tol EM; SRF/ESR.

## Why

It is the cheapest, most direct test of whether EM-sim can ship now: the port is
≈0.37-accurate and a sheet couples — so a value-distributed sheet may give a
loose-tol band-pass without a multi-week port rewrite. Either result is decisive
and honest.
