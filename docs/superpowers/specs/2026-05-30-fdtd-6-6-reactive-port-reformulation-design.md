# Phase 2.fdtd.6.6 — reactive lumped-port reformulation (sheet→mode coupling) — Design Spec

**ADR:** ADR-0121 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The de-embed bench (ADR-0119, `reactive_deembed_001`) proved PORT-WRONG: the
canonical shunt capacitor presents `Z_in ≈ 94 Ω` vs the expected `κ/(ωC) ≈ 3175 Ω`
(over-coupled ~N×), while the per-element constitutive is verified correct
per-edge (`L → +488j`, `C → −496j`). The bug is the **sheet → guide-mode
coupling**: a full-width sheet of `N` transverse cells each carrying the full
element value sums to the wrong total admittance.

## Goal

Make the de-embedded `Z_L(ω)` match `R + jωL + 1/(jωC)` within a loose tol on the
bench, by fixing the sheet→mode coupling — not the constitutive law. Then assert
the reactive arms in `reactive_deembed_001` (shunt-C first; the others as the
de-embed conditioning allows).

## Method (bench-iterable hypotheses, in order)

The bench is the harness: after each change, read `Z_L(ω)` for shunt-C (the
well-conditioned, load-bearing case), shunt-L, series-RLC, and the resistor anchor.

1. **Sheet value-normalization.** Count `N` = the number of interior transverse
   `E_z` edges the sheet spans. Distribute the lumped value so the parallel
   combination equals the intended element: shunt `C → C/N` per cell, shunt
   `L → N·L` per cell (series-RLC: dual — `R/N`? derive from the topology). Verify
   the well-conditioned capacitor `Z_L` jumps from ~94 Ω toward ~3175 Ω.
2. **Modal coupling factor.** If (1) is insufficient, the per-cell back-action
   `E_z -= (dt/(ε₀·dA))·I` and/or the source `V = E_z·dz` must carry the modal
   normalization that makes `V = ∫E·dz` and `I = ∮H·dl` consistent with a single
   lumped `Z_L`. The resistor's measured `κ` (the bench already extracts it)
   encodes this factor — use it to derive the correct reactive coupling, then
   confirm the resistor anchor still holds exactly.
3. **Modal lumped port.** If still wrong, enforce `V = Z_L·I` directly on the
   measured modal `V`/`I` (the bench's quantities), distributing the injected
   current over the port face ∝ the mode. Larger change; only if (1)/(2) fail.

Keep the public API (`series_rlc`/`pure_resistor`/`with_two_way`); the resistor
limit must stay exact; `K + β > 0` stability preserved; the one-way / fdtd-206
paths untouched.

## Changes (`crates/yee-fdtd/**` ONLY)

- The sheet→mode coupling in `LumpedRlcPort` (`src/lumped.rs`) and/or how the
  driver places the sheet (if the per-cell value split lives at the call site,
  do it inside the port or expose a sheet-aware constructor — keep API stable).
- `tests/reactive_deembed_001.rs`: once the reactive arms match, turn their
  recorded prints into **assertions** (shunt-C within a loose `Δ` of `1/(jωC)` at
  the sweep frequencies; shunt-L / series-RLC as conditioning allows). Keep the
  resistor anchor + the `// VERDICT` accounting. Flip the module verdict to
  PORT-CORRECT only when the data supports it.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: resistor anchor still exact; `lumped_lc_resonance`,
   `lumped_resistor`, `lumped_rlc_twoway_001` green (`--include-ignored`).
3. `reactive_deembed_001` GREEN with the **well-conditioned shunt-capacitor** arm
   asserted to match `1/(jωC)` within a loose tol (the minimum bar); shunt-L /
   series-RLC asserted if the de-embed conditioning permits, else recorded with a
   precise residual + a documented note (NOT a no-op, NOT faked).
4. Iterated in the bounded container; GREEN before merge.

## Out of scope

F2.3's board sim + element placement (follow-on, once the port is correct);
SRF/ESR parasitics; the studio UI (Track B).

## Why

It is the research track's increment 2, now bench-iterable. Fixing the
well-conditioned capacitor (the cleanest signal) is the decisive step; it converts
"multi-week unknown" into a measured, validated correctness fix and is the direct
unblocker for F2.3.
