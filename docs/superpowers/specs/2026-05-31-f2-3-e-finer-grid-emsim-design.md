# Filter Phase F2.3-e — finer-grid lumped EM sim — Design Spec

**ADR:** ADR-0129 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

The coarse-grid (dx = 0.4 mm) lumped-FDTD saturates at ~5 dB stopband rejection
(gate needs 20). The maintainer chose to invest in a finer grid. The aperture port
is dx-stable (ADR-0125), so finer dx should sharpen the resonance (the single-cell
port's "finer-was-worse" no longer applies).

## Goal

Determine whether finer dx (with the dx-stable aperture port + CW drive) climbs the
2.4 GHz rejection toward the strict 20 dB gate, and if so find the dx that meets it
(→ EM-sim ships). Keep the strict gate; never weaken.

## Method

The F2.3 driver already has the aperture port + CW drive (branch `6373bc7`).
`LumpedSimConfig.dx_m` is a config field — refining dx is a config change (the
voxelizer + aperture-spec scale with dx).

- **dx-refinement sweep:** run `simulate_lumped_board` at dx = 0.4 / 0.2 mm (and
  0.1 mm if container runtime allows) at the gate frequencies (2.0, 2.4 GHz), with
  the CW drive (settle scaled so the physical settle time is held constant as dt
  shrinks). Record the 2.4 GHz rejection + the 2.0 GHz |S21| (watch the over-unity
  passband artifact) vs dx.
- **Decide:** rejection climbs toward 20 dB ⇒ finer dx is the path (find the dx
  that meets the gate, set it as the F2.3 default, re-run `fdtd_lumped_001`). Caps
  shallow ⇒ a higher-accuracy aperture port + a matched/longer-board de-embed is
  the next sub-increment (separate ADR).
- **Cost control:** finer dx is N⁴; if dx = 0.1 mm DUT/thru is too slow, run 0.2 mm
  + extrapolate, and report the runtime per dx.

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- A scratch/exploratory dx-sweep first (report the trend), then — if finer dx
  reaches 20 dB — set the F2.3 default dx + re-run `fdtd_lumped_001`. The driver
  already supports `dx_m`; the settle scaling may need a small adjustment. Bring
  the branch up to current `main` first. Do NOT edit yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. fmt + clippy -D warnings on yee-voxel → exit 0.
2. The dx-refinement sweep run + REPORTED: 2.4 GHz rejection + 2.0 GHz |S21| vs
   dx (0.4 / 0.2 [/0.1]) + the runtime per dx. The verdict: climbs (→ the dx that
   meets 20 dB) or caps (→ higher-accuracy port next).
3. If a dx meets the strict gate: `fdtd_lumped_001` GREEN at 20 dB at that dx (no
   regression). Else: a precise "rejection vs dx" trend + the extrapolated dx (or
   "caps at X dB → port accuracy"). Do NOT weaken the gate.

## Out of scope

The higher-accuracy port (next sub-increment if needed); the studio UI (Track B).

## Why

It is the maintainer-chosen path's decisive first step: with the dx-stable aperture
port, finer dx should sharpen the resonance — find whether (and at what dx/cost) it
reaches the strict 20 dB cross-validation gate.
