# Phase 2.fdtd.6.8 — matched-line reactive de-embed bench — Design Spec

**ADR:** ADR-0123 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The PEC-source de-embed bench (ADR-0119/0121) can't give a long-window +
echo-free + clean-anchor reactive measurement — its `Z₀`/κ calibration is tied to
the reflecting source. So the reactive port's ≈0.37 residual is only *bracketed*,
and "PORT-CORRECT vs single-cell limit" is unresolved. Per-axis CPML (ADR-0122,
shipped) now lets us build a matched line that removes the echoes and the hack.

## Goal

A matched-line bench that **pins** the reactive port's `Z_L(ω)` definitively, with
an asserted resistor anchor, and decides: PORT-CORRECT (capacitor within tol → the
blocker was measurement, re-run F2.3) vs single-cell-limit-confirmed (→ brick 3,
the multi-cell port).

## Method (matched-line, standard incident/reflected de-embed)

Reuse the geometry idiom of `reactive_deembed_001.rs` (parallel-plate guide,
full-width source/load sheets, V=∫E·dz, modal I=∮H·dl) but change the boundaries
and the de-embed:

- **Boundaries:** `CpmlParams::for_grid(&grid, npml).with_axes([true, false,
  false])` → x-only CPML at both x-ends (absorbing); PEC on the transverse y/z
  faces (the guide). Drive a soft Gaussian `E_z` source; place a full-width lumped
  load at the load plane; a reference plane sits between source and load.
- **Reference (incident) run:** no load (or a matched/transparent load). At the
  reference plane single-bin-DFT `V_inc(ω) = ∫E·dz`, `I_inc(ω) = ∮H·dl`. Measured
  line `Z₀(ω) = V_inc/I_inc` (no fitting; matched line → clean travelling wave).
- **Load run:** total `V(ω)`, `I(ω)` at the reference plane. Because both x-ends
  absorb, there are **no multiple bounces** — the reflected wave passes the
  reference plane once and is absorbed at the source end. `Z_in(ω) = V/I`;
  `Γ = (Z_in−Z₀)/(Z_in+Z₀)`; de-embed `Z_L(ω)` (shunt: `Z_L = Z_in·Z₀/(Z₀−Z_in)`,
  or terminate so `Z_L = Z_in` — pick the cleaner topology, document it).
- **Window:** long enough to capture the full dispersive reactive tail (no echo to
  truncate against now — the only limit is the run length, set generously).
- **Honesty anchor (asserted):** a known resistor de-embeds to `Z_L → R` within a
  loose tol. If the anchor fails, the bench is wrong — fix the bench, do not
  proceed to a verdict.

## Changes (`crates/yee-fdtd/**` ONLY)

- `tests/reactive_deembed_matched_001.rs` (NEW, `#[ignore]`'d, release): the bench
  above. Assert the resistor anchor; measure + record the reactive `Z_L(ω)` table;
  assert the reactive arms to the pinned result (PORT-CORRECT if the capacitor is
  within `react_tol` — flip the verdict honestly; else pin the confirmed residual
  with a `// VERDICT:` note). Never weaken the anchor; never fake.
- A `ci.yml` release-gate job (mirror `fdtd-per-axis-cpml-gate`).
- `src/cpml.rs` only if a tiny accessor is genuinely needed (prefer not — the
  `with_axes` API + existing public surface should suffice).

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: `cpml_per_axis_001`, `cpml_reflection`, `reactive_deembed_001`,
   `lumped_lc_resonance`, `lumped_resistor`, `lumped_rlc_twoway_001` green
   (`--include-ignored`, release).
3. `reactive_deembed_matched_001` GREEN: resistor anchor asserted (`Z_L → R`);
   the pinned reactive `Z_L(ω)` recorded; the capacitor arm asserted to its
   verdict (PORT-CORRECT to `1/(jωC)` if within tol, else the confirmed residual
   pinned — not a no-op, not faked).
4. Iterated in the bounded container; GREEN before merge.

## Out of scope

The multi-cell aperture port (brick 3, only if this confirms the limit); the F2.3
re-run (follow-on); the studio UI.

## Why

It converts the ADR-0121 *bracketed* ≈0.37 into a *pinned* verdict on a
trustworthy matched line, deciding the cheapest remaining path: re-run F2.3
(if PORT-CORRECT) vs build the multi-cell port (if the limit is real).
