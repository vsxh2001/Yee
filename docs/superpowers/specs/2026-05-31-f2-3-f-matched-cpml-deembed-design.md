# Filter Phase F2.3-f — matched-CPML board de-embed — Design Spec

**ADR:** ADR-0131 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

F2.3's short-board lumped-resistor DUT/thru de-embed does not converge (ADR-0129:
over-unity passband, notch collapses at finer dx). The dominant wall is the
*measurement*, not the grid or the port (proven correct in isolation). Per-axis
CPML (ADR-0122) now enables an x-matched microstrip; the CW drive (ADR-0128)
sidesteps the CPML≠matched-at-DC dead end (ADR-0123).

## Goal

A matched-CPML board de-embed that yields a **physical, dx-stable** |S21| (no
over-unity), isolating whether the residual is the de-embed (fixed here) or the
aperture-port accuracy (the sub-cell correction next). Keep the strict 20 dB gate.

## Method

In the F2.3 driver (`crates/yee-voxel/src/lumped_sim.rs`, `run_board_solve`):

- Terminate the microstrip with **x-only CPML** at the output end (and input end,
  behind the source) via `CpmlParams::for_grid(&grid, npml).with_axes([true,false,
  false])` + the transverse-PEC clamp (the `cpml_per_axis_001` idiom), replacing
  the matched-`Z0` lumped-resistor load. The transmitted wave is absorbed after one
  pass — no backward reflection from the output end.
- Lengthen the board (more line between source / elements / output reference plane)
  so the reference planes sit clear of the element + port discontinuities, and the
  CPML has room.
- Under the CW per-frequency drive (settle to steady state), measure the
  steady-state transmitted-wave amplitude `|V_out,ss(f)|` at the output reference
  plane (with the matched end, this is the clean transmitted wave). DUT/thru:
  `S21(f) = |V_out,ss| / |V_thru,ss|`.
- Verify the result is **physical** (|S21| ≤ ~1, no over-unity) and **dx-stable**
  (re-check at dx = 0.4 & 0.2 mm — the de-embed should now converge, unlike F2.3-e).

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `run_board_solve`: x-CPML matched termination + the transmitted-wave reference-
  plane measurement (replace the lumped-resistor load + load-cell voltage). Keep
  the aperture-port element placement + the CW drive. Lengthen the board geometry
  as needed. Module doc updated. Bring the branch up to `main` first. Do NOT edit
  yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. fmt + clippy -D warnings on yee-voxel → exit 0.
2. The matched-CPML de-embed run + REPORTED: |S21| at the gate frequencies (2.0,
   2.4 GHz) — is it **physical** (no over-unity) and **dx-stable** (0.4 vs 0.2 mm)?
3. Two honest outcomes:
   - notch reaches the strict 20 dB (physical + dx-stable) → `fdtd_lumped_001`
     GREEN at 20 dB → EM-sim ships.
   - physical + dx-stable but notch still shallow → record the |S21| (how close) →
     the residual is cleanly the aperture-port accuracy → next = the sub-cell
     reactance correction. Do NOT weaken the gate.
4. No regression to other yee-voxel gates as feasible.

## Out of scope

The sub-cell port reactance correction (next, if the clean de-embed still caps the
notch); the studio Verify stage.

## Why

It is the dominant-symptom fix on the maintainer-chosen path, enabled by shipped
per-axis CPML + the CW drive: a physical, dx-stable de-embed either ships EM-sim
(if the proven port + clean measurement reach 20 dB) or cleanly isolates the
remaining port-accuracy residual for the sub-cell correction.
