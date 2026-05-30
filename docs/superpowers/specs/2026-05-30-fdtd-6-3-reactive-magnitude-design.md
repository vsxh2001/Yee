# Phase 2.fdtd.6.3 — reactive-magnitude correctness of the two-way lumped port — Design Spec

**ADR:** ADR-0117 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

`LumpedRlcPort::correct_e` (two-way path, ADR-0116) is stable and reproduces a
**resistive** load Γ exactly, but its **reactive** terms are quantitatively
wrong (measured with the gate's own single-load sweep, after the existing scalar
calibration):

- pure **inductor** (L = 13.2 nH): |Γ|_fdtd ≈ 0.013–0.019, analytic 0.32–0.60 →
  the inductor is nearly **transparent** (under-couples);
- pure **capacitor** (C = 53.4 fF): |Γ|_fdtd ≈ 0.94–1.08, analytic 0.32–0.60 →
  the capacitor acts like a **near-open** (over-couples);
- **series R-L-C**: |Γ|_fdtd ≈ 0.02 (transparent).

This blocks F2.3 (ADR-0115): the lumped-LC filter FDTD |S21| is flat ≈ 1.0
(no selectivity) because the L/C elements never load the line correctly.

## Goal

Correct the reactive (L, C) coefficients of the two-way semi-implicit update so
the discrete impedance the port realizes matches `Z(ω) = R + jωL + 1/(jωC)`,
**keeping the exact resistor limit** and the stability. Then make the gate
assert the reactive magnitude (it currently only prints it).

## Method

The current two-way update (ADR-0116) is, per `E_z` lumped cell:

```
K   = R + L/dt + dt/(2C)
β   = dt·dz / (2·ε₀·dA)
I^{n+1/2} = [ (E*+E0)·dz/2 − V_src − V_C + (L/dt)·I_old ] / (K + β)
E_z^{n+1} = E* − (dt/(ε₀·dA))·I^{n+1/2}
V_C^{n+1} = V_C + (dt/C)·I^{n+1/2}
```

The resistor limit (`L = 0`, `C = ∞`) is exact, so `β`, `V = E·dz`, and the
`(dt/(ε₀·dA))·I` field back-action are individually correct. The bug is therefore
in how the **L** (`L/dt`, `(L/dt)·I_old`) and/or **C** (`dt/(2C)` in `K`,
`dt/C` in `V_C`) coefficients scale relative to those.

**Approach (derivation-first, not guess-and-check):**
1. Take the z-transform of the discrete branch update (the `I`, `V_C`, `E_z`
   recurrences) to get the discrete impedance `Z_d(z)` the port presents, then
   evaluate at `z = e^{jωdt}` for `Z_d(ω)`.
2. Compare `Z_d(ω)` term-by-term to `R + jωL + 1/(jωC)` in the low-`ωdt` limit.
   The resistor term will match; identify the mis-scaled factor on the L and/or
   C term (a missing/extra `dz`, `dA`, `ε₀`, factor-of-2, or `dt`).
3. Correct it; re-derive to confirm `Z_d(ω) → R + jωL + 1/(jωC)`; keep `K + β > 0`
   (unconditional stability) intact.

A likely shape of the bug: the lumped branch current `I` is in **amperes** but
the field back-action and the `V = E·dz` source use a `dz`/`dA` cell-geometry
conversion that must appear **identically** on the L and C reactances as it does
on R; if the reactive coefficients omit (or double) that conversion, `jωL` and
`1/(jωC)` come out scaled by `dz/dA` (or its inverse), which over-couples C and
under-couples L exactly as observed.

## Changes (`crates/yee-fdtd/**` ONLY)

- Fix the reactive coefficients in `LumpedRlcPort::correct_e` two-way path
  (`src/lumped.rs`). Keep `series_rlc` / `pure_resistor` / `with_two_way`
  signatures and the exact resistor behaviour.
- `tests/lumped_rlc_twoway_001.rs`: convert the reactive |Γ| **print** into an
  **assert** — pure-L, pure-C, and series-R-L-C |Γ| within a loose tol
  (Δ|Γ| ≤ 0.15 suggested, after the existing scalar amplitude calibration) at the
  three sweep frequencies. Keep the resistor-exact + stability asserts.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. Resistor Γ still matches analytic (no regression); fdtd-206 LC ring-down +
   `lumped_resistor` still green (`--include-ignored`).
3. `lumped_rlc_twoway_001` (release, `#[ignore]`'d) GREEN with the new reactive
   asserts: pure-L, pure-C, series-RLC |Γ| within the loose tol at 4/6/9 GHz.
4. Iterated in the bounded container; GREEN before merge; gate not weakened.

## Out of scope

The F2.3 board sim (rides on this once it ships); SRF/ESR parasitics; multi-port.

## Why

It is the precise unblocker for the goal's "EM simulation": with reactive |Γ|
correct, F2.3's `fdtd_lumped_001` passes and the lumped filter resonates. It also
completes ADR-0116's deferred reactive DoD with an honest, asserted gate.
