# Filter Phase F2.3-c — wire F2.3 onto the aperture lumped port — Design Spec

**ADR:** ADR-0126 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

F2.3's `simulate_lumped_board` places each element as a single-edge full-width
sheet (ADR-0124) — loads the line but can't resonate (the single-cell port's
O(dx²) collapse, ADR-0124 Outcome). The aperture port (ADR-0125,
`LumpedRlcPort::aperture`) fixes that (dx-stable reactance). Wire F2.3 onto it and
re-run the gate — the decisive end-to-end EM-sim test.

## Goal

F2.3's lumped filter resonates and `fdtd_lumped_001` passes within its loose tol
(→ EM-sim ships), OR a precise "how close + CW drive needed" finding. Never weaken
the gate.

## Method

In `crates/yee-voxel/src/lumped_sim.rs`:

- Replace the per-branch element placement (the ADR-0124 `C/N`,`N·L` single-edge
  sheet) with **one `LumpedRlcPort::aperture(spec, R, L, C, src)` per branch**,
  where `spec` is the `(y,z)` port-face aperture (trace width `w` × substrate
  height `h`) at the element's x-column, and the branch carries the **aggregate**
  `R/L/C` (the aperture port handles the modal V + area-A back-action internally;
  no manual value-splitting). Series branch = aperture series-RLC; shunt branch =
  aperture pure-L ‖ aperture pure-C at the same aperture.
- Keep the air-gap-fixed line-band detection (the trace `[j_lo,j_hi)` from the
  top-metal PEC mask) and extend it to the substrate height (`k=0..n_sub`) for the
  aperture `(y,z)` extent.
- **Window:** raise `n_steps` (or the config default) enough that the capacitor's
  slow integrator tail is captured (the band-pass needs the steady-state reactance;
  ADR-0125 flagged a single short pulse reads the cap as a near-short). Try a
  generous record; if still insufficient, that's the documented CW-drive follow-on.

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `simulate_lumped_board` uses `LumpedRlcPort::aperture` per branch; module doc
  updated (aperture placement). Bring the F2.3 branch up to current `main` first
  (it needs the aperture port + the 6.x work) — `Cargo.lock` `--theirs`, keep all
  CI gate jobs (CLAUDE.md §5). Do NOT edit yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-voxel --all-targets -- -D
   warnings` exit 0.
2. `fdtd_lumped_001` (unchanged, loose tol) re-run in the container — REPORT the
   full |S21| sweep. Two honest outcomes:
   - **GREEN** (band-pass within loose tol: in-band ≈ 0 dB ±few, ≥ ~20 dB stopband)
     → EM-sim ships; branch ready for review + merge.
   - **still short** → record the achieved |S21| (e.g. "peak at 2.0 GHz, stopband
     12 dB") → the CW single-frequency drive is the next increment. Do NOT weaken
     the gate to force GREEN.
3. No regression to other yee-voxel gates (`fdtd_line_eeff_001`) as feasible.

## Out of scope

A CW drive (only if pulse + long window insufficient); tight-tol EM; the UI.

## Why

It is the decisive end-to-end test of the shipped aperture port against the goal's
EM-sim gate: the inductor is dx-stable now, so the L‖C tanks should resonate — re-run
F2.3 and find out, honestly.
