# Phase 2.fdtd.6.9 — multi-cell aperture lumped port — Design Spec

**ADR:** ADR-0125 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The single-cell `LumpedRlcPort` reactance is fundamentally too coarse for a sharp
L‖C resonance (ADR-0124 investigation): the two-way inductor back-action collapses
as **O(dx²)** (it references the bare single-cell `dA = dx²`) while the capacitor
saturates at a fixed per-cell short, so finer dx converges to a transparent line
and no dx meets F2.3's 20 dB gate. The fix must reference the field coupling to the
**modal port-face area**, not one Yee cell.

## Goal

A reactive lumped FDTD port whose de-embedded `Z_L(ω)` is **accurate AND dx-stable**
(the `O(dx²)` collapse gone), validated on a de-embed bench, keeping the
resistor-exact path. This unblocks F2.3's resonance.

## Method (the formulation the investigation pins)

A lumped element placed across a **port-face aperture** `(w × h)` = trace width ×
substrate height, on a microstrip/parallel-plate guide:

- **Branch voltage = modal `V = ∫E_z·dz`** over the full substrate height (all
  `n_sub` vertical `E_z` edges in the column), not one `E_z` edge.
- **Branch current → field back-action referenced to the aperture area `A = w·h`**
  (the modal cross-section the current threads), not the single-cell `dA = dx²`.
  Concretely the two-way update's `(dt/(ε₀·dA))·I` back-action and the `K`/`β`
  coefficients must use the aperture `A` (and the modal `V`), so the realized
  `Z_L` is **independent of the cell count / dx** — the root fix for the `O(dx²)`
  inductor collapse and the dx-frozen capacitor short.
- **Distribute the lumped element over the `(y, z)` aperture cells** with a value
  scaling that holds the **aggregate `Z_L` fixed** (the per-cell values follow from
  the aperture normalization above, not the ad-hoc `C/N`,`N·L` of ADR-0124).
- Keep `K + β > 0` (unconditional stability); the pure-R limit must reduce to the
  validated resistor exactly.

Expose additively: a new `LumpedRlcPort::aperture(...)` constructor / an aperture
spec (cells spanning `(y,z)` + the modal `V`/`A`), leaving `series_rlc` /
`pure_resistor` / `with_two_way` untouched.

## Changes (`crates/yee-fdtd/**` ONLY)

- `src/lumped.rs`: the aperture-port coupling (modal `V`, aperture-`A` back-action,
  `(y,z)` distribution). Additive API; resistor-exact preserved; stability kept.
- `tests/aperture_port_001.rs` (NEW, `#[ignore]`'d, release): de-embed a pure-L,
  pure-C, series-RLC **aperture** port on a guide at **two dx values** (e.g. 0.4 &
  0.2 mm) and assert (a) the resistor anchor `Z_L → R` (honest), (b) the reactive
  `Z_L(ω)` matches `R + jωL + 1/(jωC)` within a loose tol, and (c) **dx-stability**:
  the reactive `Z_L` at the two dx values agree within a loose tol (the `O(dx²)`
  collapse is gone — the decisive new assertion). Reuse the `reactive_deembed_001`
  PEC-source harness (the trustworthy one; NOT the matched-line, ADR-0123) +
  per-axis CPML if helpful. Never weaken the anchor; never fake.

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: `reactive_deembed_001`, `lumped_lc_resonance`, `lumped_resistor`,
   `lumped_rlc_twoway_001`, `cpml_per_axis_001`, `cpml_reflection` green
   (`--include-ignored`, release).
3. `aperture_port_001` GREEN: resistor anchor asserted; reactive `Z_L` within loose
   tol AND **dx-stable across two dx** (the `O(dx²)` collapse gone). If the reactive
   tol is not yet met but dx-stability IS achieved, that is a real partial — record
   it precisely (the collapse is the thing to kill first).
4. Iterated in the bounded container; GREEN before merge.

## Out of scope

F2.3's re-run on the aperture port (the next sub-increment, once the bench shows
accurate + dx-stable reactance); tight-tol EM; SRF/ESR; the studio UI.

## Why

It is the de-risked, mechanism-specified EM-sim unblocker: kill the `O(dx²)`
inductor collapse via aperture-`A` back-action + modal `V`, prove it on the bench
(accurate + dx-stable), then F2.3's tanks can resonate. The dx-stability assertion
is the concrete, decisive success signal the investigation handed us.
