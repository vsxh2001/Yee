# FEM numerical-eigenmode microstrip port — design

**Date:** 2026-06-03
**ADR:** [ADR-0154](../../src/decisions/0154-fem-numerical-eigenmode-microstrip-port.md)
**Forks from:** ADR-0153 (B4 GO point / B7 v1-port ceiling)

## Problem

ADR-0153's FEM driven-sweep track hit the ε_eff milestone (B4: 0.61 % of HJ, real physics) but the
3-pole filter S21 (B7) floors at ~−42 dB. Root cause, reviewer-confirmed: the v1 analytic flat-E_z
port `modal_e_t` only ~9 % overlaps the true FEM eigenmode → ~−21 dB modal-overlap loss per port,
squared across two ports. The fix is a **higher-fidelity modal shape** at the port face.

A decisive de-risk probe (verified independently in the box) proved the fix: feeding `yee-mom`'s
shipped quasi-TEM cross-section eigenmode as the port's `modal_e_t` lifts **|S21| 0.089 → 0.778**
(8.74×) and collapses **|S11| 0.573 → 0.087** (a genuinely matched port), with ε_eff phase
unchanged at 0.61 % (only the shape moved). The maintainer picked "Full Option-1 now."

## Goal

Productionize the numerical-eigenmode port (promote the validated probe code from test into
`yee-fem/src`), gate it on the straight line with a real |S21| lower-bound tripwire, then re-grade
the 3-pole filter S21 against the strict Chebyshev mask.

## Non-goals

- Reopening the ADR-0064 planar-MoM port, ADR-0133 FDTD cavity wall, or `fem-eig-006` eigen port.
- The full Sommerfeld tail or `bicgstab` scaling (ADR-0153 B5b deferred — N1–N3 fit the per-ω LU).
- Taking β from the eigensolve: N1 keeps analytic-HJ β (validated 0.61 %); the numerical mode
  supplies **only the shape** (the ADR-0153 B4 non-circularity probe showed driven ε_eff is robust
  to a mistuned port β, and the eigensolve's own ε_eff has an ~8 % box-truncation gap).

## Architecture

Three ordered bricks (N1 → N2 → N3); N1+N2 ship as one de-risked increment, N3 follows.

### N1 — `microstrip_port_numerical` production API (`yee-fem/src`)

- **Cargo:** promote `yee-mom = { workspace = true }` from `[dev-dependencies]` to
  `[dependencies]` in `crates/yee-fem/Cargo.toml`. Acyclic (yee-mom → core/mesh/io only).
- **New module** `crates/yee-fem/src/microstrip_port_numerical.rs` (or extend `microstrip_port.rs`),
  promoting the three probe helpers from `tests/microstrip_eeff.rs` into public `src` items:
  - `fn microstrip_cross_section(box_w, box_h, sub_h, trace_w, t_strip, nx, ny) -> TriMesh2D` —
    strip-as-hole builder (mirrors `yee-mom/tests/eigensolver_microstrip_quasi_tem.rs`), FR-4 tag 1
    below `sub_h`, air tag 0 above, signal strip a rectangular hole (its border edges = PEC).
  - a one-shot `NumericalCrossSection::with_quasi_tem` solve, **Arc-shared** across the two port
    closures (cheap ~ms solve, ran once).
  - `fn microstrip_port_numerical(geom: &MicrostripPortGeom, f_hz: f64) -> Result<PortDefinition>` —
    returns a `PortDefinition::single_mode` with analytic-HJ β + the numerical `modal_e_t` closure.
    **Frame map** (validated by the probe): yee-mom `(coord0=x_width, coord1=substrate-normal)` →
    FEM `(x, z)`; sample `mode.e_tangential_at(p.x, p.z)` → `[ex, e_normal]` → emit
    `Vector3::new(ex, 0.0, e_normal)` (substrate-normal on FEM ẑ; propagation-ŷ component 0).
  - a small `MicrostripPortGeom` struct (trace_w, sub_h, eps_r, box_w, box_h) so callers don't pass
    7 positional args; the cross-section density (≈20×12) is an internal validated default.
- **DoD gate:** clippy `-D warnings` + fmt clean; a fast non-ignored unit test asserts the port's
  `modal_e_t` is finite + nonzero + **E_z-dominant in the gap and decaying into the air** (mirror
  the existing `modal_e_t_is_ez_dominant_in_gap_and_decays_in_air` test for the v1 port), and β
  matches `yee_layout::eps_eff` to the existing tolerance.

### N2 — validated straight-line |S21| gate (`fem_line_eeff_numerical_001`)

- In `crates/yee-fem/tests/microstrip_eeff.rs`, promote the probe MEASUREMENT into a PASS/FAIL gate
  that calls the new `src` API (delete the probe's inline helpers — they now live in `src`).
- Same mesh/box/interior-PEC/coupled-Whitney/two-length extraction as `fem_line_eeff_001`.
- **DoD gate:** `|S21|(L2) ≥ 0.6` (lower-bound tripwire — the probe measured 0.778, so 0.6 is a
  safe floor that still catches a re-flooring regression) **AND** `|S11|(L2) ≤ 0.2` **AND**
  `ε_eff` within 5 % of HJ. `#[ignore]`'d + run in the existing `fem-eigen` `--release` gate job.

### N3 — filter S21 re-grade (`microstrip_filter_s21.rs`)

- Swap the filter's port construction to `microstrip_port_numerical`; re-run the driven sweep.
- **DoD gate (HONEST):** record the measured in-band peak, band-edge depths, and mask margin; assert
  the **measured lift** over the v1 −42 dB floor + that the **asymmetry discriminator still fires**
  (depth(1.6 GHz) > depth(2.4 GHz)). If |S21| clears the `oracle_grade` strict mask, assert that;
  if it lifts dramatically but stops short, assert the measured margin honestly (a real bandpass
  with a quantified gap — not a fake pass). Never assert a match-by-construction.

## Data flow

`MicrostripPortGeom + f` → (N1) build (x,z) cross-section TriMesh2D → `NumericalCrossSection::
with_quasi_tem().solve(f)` (once, Arc) → `PortDefinition{ beta: analytic-HJ, modal_e_t: sample
e_tangential_at }` → (N2/N3) `OpenBoundarySolver::sweep_matrix` per-ω complex LU → S-params.

## Testing

- N1: fast unit test (non-ignored, debug-safe): modal-shape + β assertions.
- N2: `#[ignore]`'d release gate (`fem-eigen` job): the |S21|≥0.6 / |S11|≤0.2 / ε_eff≤5 % tripwire.
- N3: `#[ignore]`'d release gate: honest mask-margin + asymmetry assertions.
- All heavy runs boxed (`scripts/yee-box.sh`, ≤14 g / 3 cpu).

## Risks

- **N3 strict-mask clearance (moderate):** the line is proven; the filter adds resonator coupling +
  gap-mesh sensitivity. Mitigation: honest gate — a dramatic-lift-but-short result is documented,
  not failed; the follow-on is finer mesh, not a port wall.
- **Cross-section convergence (low):** the probe used the validated ≈20×12 density (NOT the 8×8 that
  failed in `mom_002_numerical_waveport.rs`). Keep that density as the internal default.
