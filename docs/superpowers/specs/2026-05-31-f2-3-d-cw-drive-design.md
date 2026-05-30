# Filter Phase F2.3-d — CW per-frequency drive — Design Spec

**ADR:** ADR-0128 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

The aperture lumped port is proven correct (ADR-0127: the cap presents `1/(jωC)`,
the L‖C tank resonates under CW), but F2.3's modulated-Gaussian pulse + DFT
measures an unsettled transient on a short standing-wave line, so the Q≈10 tanks
never reach steady state and no band-pass forms. F2.3 needs a CW steady-state
measurement.

## Goal

F2.3's lumped filter resonates and `fdtd_lumped_001` passes within its loose tol
(→ EM-sim ships), via a CW per-frequency steady-state drive. Never weaken the gate.

## Method (CW per-frequency steady-state)

In `crates/yee-voxel/src/lumped_sim.rs`, replace `run_board_solve`'s single
modulated-Gaussian + broadband DFT with a per-frequency loop:

- For each measured `f`: drive a CW sinusoid `sin(2πf·n·dt)` at the source sheet,
  **Hann/linear-ramped** over the first ~M cycles to suppress the turn-on
  transient; run `n_settle` cycles (enough for the highest-Q tank to reach steady
  state — Q≈10 ⇒ ~10–30 cycles — plus the source→load line transit); then measure
  the **steady-state** load-voltage amplitude over the last few cycles (peak
  envelope, or a single-bin DFT over the settled window only).
- Run the SAME for the DUT (elements present) and the thru (elements removed);
  `S21(f) = |V_dut,ss| / |V_thru,ss|` (DUT/thru divides out the line standing wave).
- **Frequency set:** the gate-check points (2.0 GHz passband, 2.4 GHz stopband) +
  a handful for the sweep shape — NOT a fine sweep (per-freq CW = 2 solves each;
  keep it bounded). A small `Vec<f64>` of CW frequencies.

Keep the aperture-port placement (one aggregate-R/L/C `LumpedRlcPort::aperture`
per branch, air-gap-fixed band) from ADR-0124/0126.

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `run_board_solve` (or a new CW variant) + `simulate_lumped_board`: CW
  per-frequency steady-state drive + DUT/thru steady-state amplitude ratio. Module
  doc updated. Bring the F2.3 branch up to current `main` first (Cargo.lock
  `--theirs`, keep all CI gate jobs). Do NOT edit yee-fdtd/yee-filter.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-voxel --all-targets -- -D
   warnings` exit 0.
2. `fdtd_lumped_001` (unchanged, loose tol) re-run in the container — REPORT the
   CW |S21| at the gate frequencies. Two honest outcomes:
   - **GREEN** (band-pass within loose tol: in-band ≈ 0 dB ±few, ≥ ~20 dB stopband)
     → EM-sim ships; branch ready for review + merge.
   - **still short** → record the achieved CW |S21| (how close). Do NOT weaken the
     gate.
3. Runtime bounded (CW only at the gate points + a few); no regression to other
   yee-voxel gates as feasible.

## Out of scope

A fine CW sweep; tight-tol EM; the UI.

## Why

The physics is proven (ADR-0127); this is the measurement change that lets the
already-correct resonant tanks show their band-pass — the loop-closing increment
for the goal's EM-simulation component.
