# Phase 2.fdtd.6.7 — per-axis CPML face selection — Design Spec

**ADR:** ADR-0122 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

The reactive-port research track's increment 3 (ADR-0121) needs a matched-line
de-embed bench: a parallel-plate guide absorbing at both x-ends but PEC on the
transverse (y, z) faces. The repo CPML (`cpml.rs`) is symmetric on all six faces
and destroys the guide mode when applied; 2.fdtd.6.6 confirmed every source-end
absorber broke the bench. A matched line is unbuildable without per-axis CPML.

## Goal

Generalize the CPML to absorb on a caller-selected subset of axes (x-only for the
matched line), default unchanged (all three axes), with a validation gate. Enable
the matched-line bench (the next sub-increment) without regressing any existing
FDTD gate.

## Method

`cpml.rs` currently applies CPML stretching on all three axes symmetrically
(`CpmlState::update_e`/`update_h` loop over the x/y/z PML regions; `pml_depth(i,n)`
returns the in-PML depth + side per axis). Add an axis mask:

- `CpmlParams` gains `axes: [bool; 3]` (x, y, z), **default `[true; 3]`** — so
  every existing caller and gate is byte-identical. A `with_axes([bool;3])`
  builder (or a `CpmlParams::for_grid_axes`) sets it.
- `CpmlState` stores the mask; `update_e`/`update_h` skip the PML
  contribution for any axis whose flag is `false` (that axis's faces then carry no
  CPML — the caller applies `apply_pec` there as usual). `pml_depth` returns
  `None` for a disabled axis.
- Net effect for `axes = [true,false,false]`: CPML absorbs only at x=0 and x=nx;
  y and z faces are left for PEC. The guide mode (E_z between PEC z-plates,
  bounded by PEC y-walls) propagates in x and is absorbed at both x-ends.

Keep the public API additive; the existing `CpmlParams::for_grid` / `CpmlState::new`
paths unchanged in behaviour (mask defaults all-true).

## Changes (`crates/yee-fdtd/**` ONLY)

- `src/cpml.rs`: the `axes` mask on `CpmlParams` + `CpmlState`, honoured in
  `update_e`/`update_h`/`pml_depth`; a `with_axes`-style builder. Default all-true.
- `tests/cpml_per_axis_001.rs` (NEW, `#[ignore]`'d, release): x-only CPML on a
  guide with PEC transverse walls; an x-travelling Gaussian pulse; measure E-field
  reflection at an interior probe vs an all-PEC control → **≥30 dB reduction**
  (mirror `cpml_reflection.rs`). PLUS assert the transverse PEC walls are intact
  (tangential E on the y/z faces ≈ 0 — the guide is preserved, not absorbed).

## DoD (machine-checkable; container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0.
2. No regression: `cpml_reflection` (all-faces, default mask) green; the FDTD
   line/coupling gates that use CPML green (`--include-ignored`, release).
3. `cpml_per_axis_001` GREEN: x-only CPML ≥30 dB reduction vs PEC control + the
   transverse PEC walls verified intact.
4. Iterated in the bounded container; GREEN before merge; gate not weakened.

## Out of scope

The matched-line de-embed bench (next sub-increment, rides on this); the multi-cell
aperture port; F2.3; asymmetric per-*face* (vs per-axis) slabs.

## Why

It is the enabling first brick of increment 3 (the maintainer-approved multi-week
EM-sim path), independently useful (per-axis CPML = any waveguide / matched-line
FDTD), and bounded + validatable on its own — concrete progress, not deferral.
