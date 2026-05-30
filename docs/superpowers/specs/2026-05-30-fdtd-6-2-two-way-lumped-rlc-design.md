# Phase 2.fdtd.6.2 — stable two-way lumped RLC port — Design Spec

**ADR:** ADR-0116 · **Date:** 2026-05-30 · **Status:** Accepted

## Problem

`yee_fdtd::LumpedRlcPort` cannot model a low-loss reactive S-parameter element
(blocks F2.3 lumped-LC EM sim, ADR-0115). Confirmed via container experiments:

- the `l > 0` (inductor / series-RLC) branch is **one-way** (circuit current is
  driven by the field, but the lumped current does **not** feed back into the
  `E_z` update) ⇒ a source-free inductor is **inert**, and a shunt L‖C never
  resonates;
- the only two-way arm (pure capacitor, `l = 0`) is **unstable below ~196 Ω**
  ESR (≈ η₀/√3, grid-independent) — it diverges to NaN for any low-loss value.

So no parameter choice is simultaneously two-way, stable, and low-loss.

## Goal

Implement the **stable, two-way** lumped-element `E_z` update for a series
**R–L–C** at a cell — the Piket-May / Taflove–Hagness semi-implicit formulation
(unconditionally stable; the lumped branch current updates implicitly with the
field so energy flows both ways). This is the canonical FDTD lumped-element
method; it replaces the one-way inductor path and the unstable capacitor path
with one stable two-way update covering R, L, C (and the `l=0` / `c=∞` limits).

## Method (Piket-May–Taflove lumped-element FDTD)

At the lumped cell, the standard `E_z^{n+1}` Yee update gains a lumped-current
term. For a series R–L–C between the `E_z` node, write the branch current with a
semi-implicit (Crank-Nicolson-style) discretisation of all three of R, L, C and
solve the coupled `E_z^{n+1}` / branch-state equations **together** (not the
field first then a one-way correction). The resulting update is unconditionally
stable for any non-negative R, L and positive C (Taflove–Hagness §… lumped
elements). Reduce cleanly to: pure R (`pure_resistor` — keep current behaviour),
pure C (`l=0`), pure L (`c=∞`), and a Thévenin source (`SourceWaveform`).

## Changes (`crates/yee-fdtd/**` ONLY)

- Rework `LumpedRlcPort::correct_e` (and its state: `e_z_prev`,
  `inductor_current`, `capacitor_voltage`) to the stable two-way semi-implicit
  update for the general series R–L–C. Keep the public constructors
  (`series_rlc`, `pure_resistor`) and signature stable; this is an internal
  correctness fix to the update math.
- Keep `pure_resistor`'s existing validated behaviour (the resistor gate must
  stay green).

## DoD (machine-checkable; CI-gated, container-iterated)

1. `cargo fmt --check --all` + `cargo clippy -p yee-fdtd --all-targets -- -D
   warnings` exit 0 (local, light).
2. Existing yee-fdtd lumped tests (the resistor/series-RLC validation referenced
   in `lumped.rs` docs) **still pass** — no regression.
3. New `#[ignore]`'d gate `lumped_rlc_twoway_001` (CI `--release`):
   - **Stability:** a low-loss reactive element (e.g. a shunt C, or a series L–C
     at the F2.0 values, ESR ≈ 1e-3 Ω) runs the full record with **no NaN/Inf**
     and bounded fields (the ≥196 Ω instability is gone).
   - **Two-way correctness:** a single lumped load terminating the line reflects
     with `|S11|` (or the equivalent ratio) matching the analytic lumped-load
     Γ = `(Z_L − Z0)/(Z_L + Z0)`, `Z_L = R + jωL + 1/(jωC)`, within a loose tol
     at a few frequencies (incl. one where the reactance dominates) — i.e. a
     source-free inductor/capacitor is **not inert** and resonates correctly.
4. Iterate in the bounded container; the gate is GREEN in CI on the branch before
   merge (CLAUDE.md §4); not weakened to a no-op.

## Out of scope

Parallel-RLC as a single primitive (parallel is composed from two series ports,
per F2.3); shunt-vs-series wiring (that's the F2.3 driver's job); non-`E_z`
orientations; the F2.3 board sim itself (rides on this once it ships).

## Why

It is the precise, literature-bounded unblocker for the goal's "EM simulation"
(F2.3): F2.3's driver + gate are already correct and pass **unchanged** once this
primitive is two-way and stable. It also upgrades a core FDTD capability (proper
lumped elements) used beyond filters.
