# Phase 2.fdtd.6.4 — canonical per-element Taflove lumped L/C updates — Design Spec

**ADR:** ADR-0118 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

ADR-0117 proved (derivation-first, replica-confirmed) that the
RLC-in-one-implicit-`K` two-way update (ADR-0116) loads the line by the
*instantaneous* `K = R + L/dt + dt/(2C)`: a shunt inductor presents `L/dt ≈ 7.6 kΩ`
(transparent, |Γ|≈0.01) and a shunt capacitor presents `dt/(2C) ≈ 16 Ω`
(near-short, |Γ|≈1.0) — wrong in **opposite directions**, so no single coefficient
fixes both. The per-frequency branch impedance is correct; the single-step
coupling magnitude is not. F2.3's lumped filter |S21| is flat ≈ 1.0.

## Goal

Implement the **canonical Taflove–Hagness per-element lumped updates** so a shunt
inductor presents `jωL`, a shunt capacitor presents `1/(jωC)`, and a series RLC
presents `R + jωL + 1/(jωC)` — to the *line*, not just per-frequency. Keep the
validated resistor path and the public API. Then assert the reactive |Γ|.

## Method (canonical lumped-element FDTD)

The field coupling already validated for the resistor is reused: a lumped element
on an `E_z` edge has terminal voltage `V = E_z·dz` and injects a current `I` into
the `E_z` update as `E_z^{n+1} = E_z^* − (dt/(ε₀·dA))·I`, where `E_z^*` is the
ordinary Yee-updated value, `dz` the edge length, `dA` the transverse cell area.
Only the **constitutive relation between `I` and `V`** changes per element type:

- **Capacitor** (`C`): the lumped C augments the cell's natural capacitance. The
  standard result is an effective permittivity `ε_eff = ε₀ + C·dz/dA` at that
  edge: `E_z^{n+1} = E_z^n + (dt/ε_eff)(curlH/dz-terms)`. Equivalently the
  displacement current `I_C = C·dV/dt = C·dz·dE_z/dt` enters the update. Presents
  `Z_C = 1/(jωC)`. Stable (capacitance only increases).
- **Inductor** (`L`): an auxiliary state current accumulates the voltage,
  `I_L^{n+1/2} = I_L^{n−1/2} + (dt·dz/L)·E_z^n`, and the `E_z` update subtracts
  `(dt/(ε₀·dA))·I_L^{n+1/2}`. Presents `Z_L = jωL`. Stable (no CFL penalty).
- **Resistor** (`R`): unchanged — the validated `pure_resistor` update.
- **Series R-L-C**: the canonical combined update carrying both the inductor
  accumulator and the capacitor-voltage state in series with R, all sharing the
  `dz/dA` coupling (Taflove lumped-RLC `E` update).

Cross-check each against the analytic shunt reflection `Γ = −Z₀/(2Z_L+Z₀)` with
the gate's calibrated z0_eff.

## Changes (`crates/yee-fdtd/**` ONLY)

- Rework the reactive arms of `LumpedRlcPort::correct_e` (and its state:
  `inductor_current`, `capacitor_voltage`) to the canonical per-element updates
  above. Keep `series_rlc` / `pure_resistor` / `with_two_way` signatures and the
  exact resistor limit. The `l=0` (pure-C) and `c=∞` (pure-L) limits must hit the
  pure-capacitor / pure-inductor updates respectively.
- `tests/lumped_rlc_twoway_001.rs`: turn the reactive |Γ| **prints** into
  **asserts** — shunt-L, shunt-C, series-RLC |Γ| within Δ|Γ| ≤ 0.15 of the
  analytic `−Z₀/(2Z_L+Z₀)` shunt law at 4/6/9 GHz, after the test's existing
  scalar calibration. Keep resistor-exact + stability (no-NaN) asserts.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: resistor Γ still exact; `lumped_lc_resonance` + `lumped_resistor`
   green (`--include-ignored`).
3. `lumped_rlc_twoway_001` (release, `#[ignore]`'d) GREEN with the reactive
   asserts: shunt-L, shunt-C within Δ|Γ| ≤ 0.15 at 4/6/9 GHz. Series-RLC within
   the same tol **or** (escape hatch) shunt cases green + series-RLC deferred with
   a clear in-test note — not weakened to a no-op.
4. Iterated in the bounded container; GREEN before merge.

## Out of scope

The F2.3 board sim (rides on this once it ships); SRF/ESR parasitics; multi-port.

## Why

The canonical method is the textbook, validated way to load an FDTD line with a
reactive lumped element; it is the precise unblocker for the goal's "EM
simulation" (F2.3) and is bounded (not the ill-posedness that defers the MoM
port). With it, F2.3's `fdtd_lumped_001` acquires its band-pass shape.
