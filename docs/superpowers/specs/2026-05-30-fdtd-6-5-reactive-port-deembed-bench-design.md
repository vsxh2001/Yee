# Phase 2.fdtd.6.5 — reactive lumped-port V+I de-embedding bench — Design Spec

**ADR:** ADR-0119 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The reactive lumped-port research track (maintainer-approved) must first resolve a
contradiction (ADR-0117/0118): the **port-local** proxy says a canonical inductor
presents ≈`jωL` (correct), but the **line-reflection** measurement (gated DFT +
scalar `A`/`z0_eff` calibration) says it is transparent. We cannot pick a
multi-week reformulation until we know whether the *port* or the *measurement* is
wrong.

## Goal

A clean **V+I de-embedding bench** that extracts a lumped load's `Z_L(ω)` directly
from measured voltage and current at a reference plane on the TEM line, and a
recorded verdict: port-correct (→ measurement/placement fix) vs port-wrong (→
reformulation).

## Method (VNA-style 1-port de-embed)

On the existing parallel-plate TEM line (`tests/lumped_rlc_twoway_001.rs` harness):

1. **Line `Z₀(ω)`**: from a *matched / open* run, take the incident travelling
   wave; `Z₀(ω) = V_inc(ω)/I_inc(ω)`, where `V = Σ_z E_z·dz` across the plate gap
   at the reference plane and `I = ∮ H·dl` (Ampère loop) around the conductor at
   the same plane. This is a measured property of the discrete line — no fitting.
2. **Load run**: place a single canonical lumped element (pure-L `c=∞`, pure-C
   `l=0`, series-RLC) as a full-width shunt sheet at the load plane. Single-bin
   DFT `V(ω)`, `I(ω)` at the reference plane.
3. **De-embed**: `Z_in(ω) = V(ω)/I(ω)`; for a shunt load on a continuing line,
   `Z_L = (Z_in·Z₀)/(Z₀ − Z_in)` (or terminate the line and use `Z_L = Z_in`
   directly — pick the cleaner topology). Form `Γ = (Z_in−Z₀)/(Z_in+Z₀)`.
4. **Compare** `Z_L(ω)` to `R + jωL + 1/(jωC)` at 3+ frequencies.

Use the canonical per-element updates from branch `feature/fdtd-6-4-canonical-lc`
(`021bed2`) — bring that code onto this branch (it is per-edge-verified and the
right foundation), so the bench measures the canonical port.

## Changes (`crates/yee-fdtd/**` ONLY)

- Adopt the canonical per-element `correct_e` updates (from `021bed2`).
- New `tests/reactive_deembed_001.rs`: the V+I bench above. ASSERT the resistor
  case (`Z_L(ω) → R` within a loose tol — the known-good anchor that proves the
  bench is honest). PRINT/record the reactive `Z_L(ω)` table and assert whichever
  verdict the data supports (if reactive matches, assert it; if not, assert the
  resistor anchor + a documented `// VERDICT:` note with the measured numbers).
  Never weaken the resistor anchor; never fake.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: existing lumped gates (`lumped_lc_resonance`, `lumped_resistor`,
   `lumped_rlc_twoway_001`) green (`--include-ignored`).
3. `reactive_deembed_001` GREEN: resistor `Z_L → R` asserted; the reactive
   `Z_L(ω)` table recorded with an explicit verdict (port-correct vs port-wrong).
4. Iterated in the bounded container; the verdict written back into ADR-0119.

## Out of scope

The port reformulation (increment 2, gated on the verdict); F2.3 placement; the
studio UI (Track B).

## Why

It de-risks the whole research track: a few hours of measurement work decides
whether EM-sim is a measurement/placement fix (close) or a multi-week port rewrite
(far) — and either way leaves a trustworthy reactive-load bench behind.
