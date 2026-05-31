# Filter Phase F2.3-i — re-scope the lumped EM-sim gate — Design Spec

**ADR:** ADR-0134 · **Date:** 2026-05-31 · **Status:** Accepted

## Problem

`fdtd_lumped_001`'s ≥20 dB full-board band-pass cross-validation is a fundamental
FDTD-measurement wall (ADR-0133: CW high-Q microstrip in a stable box is
cavity-dominated; matched CPML unstable into substrate). The physics is validated
(port in isolation; circuit `ladder_s21`); the full-board measurement is
intractable. Maintainer chose to re-scope the gate to an achievable bar + ship.

## Goal

`fdtd_lumped_001` GREEN at a **principled, achievable, NON-VACUOUS** bar that
validates what the full-wave FDTD EM-sim genuinely demonstrates, with the sharp
response delegated to `ladder_s21` and the wall documented. Then F2.3 merges →
EM-sim ships (6/6). NEVER fake.

## Method

In `crates/yee-voxel/tests/fdtd_lumped_001.rs` (and minimal `lumped_sim.rs` if a
helper is needed):

1. Keep the EM-sim pipeline (synthesize_lumped → lumped_board → voxelize → aperture
   ports → FDTD solve → S21 sweep). Use the cleanest of the F2.3 de-embeds (the
   PEC-box 2-point physical de-embed, F2.3-g, or the simplest stable variant) at
   dx=0.4 mm, bounded.
2. Replace the ≥20 dB band-pass assertions with a **principled achievable** set,
   chosen from what the actual F2.3 board data RELIABLY shows:
   - **pipeline finite/non-trivial:** the sweep is finite (no NaN/Inf), non-empty.
   - **elements LOAD the line (the real EM contribution, MUST be non-vacuous):**
     the loaded-DUT response differs MEANINGFULLY from the bare-thru by a real
     margin and/or is meaningfully frequency-dependent — i.e. it would FAIL for the
     inert single-cell-placement flat≈1 response (ADR-0124) or a broken sim. Pick
     the metric + threshold from the data (e.g. "max |20·log10(S21_dut/S21_thru)|
     across the band ≥ A dB", with A set well above the inert-noise floor but below
     what the cavity-limited measurement reliably delivers).
3. Docstring: delegate the SHARP-response cross-validation to `ladder_s21` (F2.0,
   vs Pozar) + the per-element reactance to `aperture_port_001`/`cap_cw_001`, and
   document the cavity wall (ADR-0133) — why full-board ≥20 dB is not FDTD-achievable.

## Changes (`crates/yee-voxel/**` ONLY, on the F2.3 branch)

- `tests/fdtd_lumped_001.rs`: the re-scoped assertions + the docstring (delegation +
  wall). `src/lumped_sim.rs` only if a small helper/metric is needed. Bring the
  branch up to `main` first. Do NOT edit yee-fdtd/yee-filter. Do NOT touch the
  isolation gates (aperture_port_001/cap_cw_001) or ladder_s21.

## DoD (machine-checkable; container-iterated)

1. fmt + clippy -D warnings on yee-voxel → exit 0.
2. `fdtd_lumped_001` GREEN at the re-scoped achievable bar, in the container.
3. **The re-scoped assertion is NON-VACUOUS** — demonstrate it would FAIL for an
   inert/flat response (e.g. note the inert single-cell margin ≈0 dB vs the asserted
   threshold; or a brief check). The reviewer confirms this.
4. No regression to aperture_port_001 / cap_cw_001 / lumped_lc_resonance /
   lumped_resistor / voxel gates.

## Out of scope

Chasing ≥20 dB full-board (the wall); the studio Verify wiring; a stable-non-CPML
absorber research track.

## Why

It ships the goal's EM-sim component honestly at the bar the method achieves —
full-wave FDTD board + components loading the line (real EM) + circuit-level sharp
cross-validation — per the maintainer's explicit re-scope choice, without faking
the unreachable ≥20 dB full-board measurement.
