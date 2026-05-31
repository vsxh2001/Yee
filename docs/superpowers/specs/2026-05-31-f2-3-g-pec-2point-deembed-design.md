# Filter Phase F2.3-g — PEC-box 2-point standing-wave CW de-embed — Design Spec

**ADR:** ADR-0132 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

Every F2.3 de-embed so far (short-board, finer-grid, matched-CPML) gives a
monotone + over-unity response, never a band-pass; matched-CPML also hit the
ADR-0108 CPML-into-substrate instability. Over-unity is a bad-de-embed signature
(no forward/backward wave separation). A proper 2-point standing-wave de-embed on a
stable PEC box (the `run_line_eeff` pattern, no CPML) is the maintainer-chosen fix.

## Goal

A physical (no over-unity) S21 from a PEC-box 2-point standing-wave CW de-embed,
and a clear verdict: band-pass forms (ship if ≥20 dB), shallow band-pass (→ port
accuracy), or still no band-pass (→ board-integration finding). Keep the strict gate.

## Method (standard CW S-param extraction via standing-wave probing)

In the F2.3 driver (`crates/yee-voxel/src/lumped_sim.rs`):

- **PEC box** (NO CPML): the voxelized microstrip in a PEC-bounded grid, line run
  long enough past each port that the elements clear the ends and a steady standing
  wave forms. (Mirror `run_line_eeff`'s PEC + CW-on-a-line stability, ADR-0108.)
- **CW drive** (Hann-ramped, settle to steady state) at the input.
- **2-point (or 3-point) standing-wave probe** at each port reference region:
  sample the steady-state line voltage phasor at ≥2 points of known spacing `d`.
  With `V(x) = a·e^{−jβx} + b·e^{+jβx}`, two phasors `V(x₀), V(x₀+d)` solve for the
  forward `a` and backward `b` travelling-wave amplitudes (use β from a thru-line
  ε_eff calibration run, or fit β from a 3rd point).
- **S21** = `b₂ / a₁` (the forward wave launched into port-2's region over the
  incident forward wave at port 1), thru-normalized: `S21(f) = (b₂/a₁)_dut /
  (b₂/a₁)_thru`. (The thru run de-embeds the launch + line.)
- Frequency set = the gate points (2.0, 2.4 GHz) + a few for the shape (bounded).

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `run_board_solve`: PEC box (drop the CPML termination from F2.3-f), the 2-point
  standing-wave forward/backward extraction, S21 = b₂/a₁ thru-normalized. Keep the
  aperture-port placement + the CW drive. Module doc updated. Bring the branch up
  to `main` first. Do NOT edit yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. fmt + clippy -D warnings on yee-voxel → exit 0.
2. The 2-point de-embed run + REPORTED |S21| at the gate freqs: is it PHYSICAL (no
   over-unity)? Does a band-pass form (peak at 2.0 GHz, notch at 2.4 GHz)?
3. Three honest outcomes:
   - band-pass + notch ≥ 20 dB + passband ≈ 0 dB → `fdtd_lumped_001` GREEN → ships.
   - physical band-pass but shallow notch → record |S21| → next = the sub-cell port
     correction.
   - still monotone / no band-pass → record it → the board integration doesn't
     resonate (deeper finding to surface). Do NOT weaken the gate in any case.
4. No regression to other yee-voxel gates as feasible.

## Out of scope

The sub-cell port correction (next, if the clean de-embed shows a shallow band-pass);
the studio Verify stage.

## Why

It is the maintainer-chosen non-CPML de-embed: a stable PEC box + proper
forward/backward wave separation removes both the over-unity artifact and the CPML
instability, giving the first trustworthy F2.3 S21 — which either reveals the
band-pass (ship / port-accuracy next) or proves the board-integration wall is real.
