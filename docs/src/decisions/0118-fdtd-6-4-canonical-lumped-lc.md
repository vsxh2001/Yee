# ADR-0118: Phase 2.fdtd.6.4 — canonical per-element Taflove lumped L/C updates

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0117 (2.fdtd.6.3 — found the reactive defect is structural, not a
coefficient), ADR-0116 (2.fdtd.6.2 — the stable two-way *resistor* port stands),
ADR-0115 (F2.3 lumped FDTD — blocked on this), the lumped-LC → PCB goal,
[[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0117 proved (derivation-first) that the RLC-in-one-implicit-`K` two-way update
loads the line by the **instantaneous** `K = R + L/dt + dt/(2C)`, so a shunt
inductor presents `L/dt ≈ 7.6 kΩ` (transparent) and a shunt capacitor presents
`dt/(2C) ≈ 16 Ω` (near-short) — wrong in **opposite directions**, uncorrectable by
any single coefficient. The per-frequency branch impedance is right; the
single-step *coupling magnitude* is not. F2.3's lumped filter therefore shows a
flat |S21| ≈ 1.0 (no selectivity).

## Decision

Replace the reactive coupling with the **canonical Taflove–Hagness per-element
lumped updates**, which couple to the field correctly per timestep:

- **Lumped capacitor** at an `E_z` edge → augment the cell's update with the lumped
  capacitance, equivalent to a local effective permittivity
  `ε_eff = ε₀ + C·dz/dA`: the standard lumped-C FDTD update
  `E^{n+1} = E^n + (dt/ε_eff)(∇×H − …)`. Presents `1/(jωC)`. Unconditionally
  stable (it only *raises* the cell capacitance).
- **Lumped inductor** at an `E_z` edge → an auxiliary branch current that
  *accumulates* the voltage: `I_L^{n+1/2} = I_L^{n−1/2} + (dt·dz/L)·E_z^n`, with
  the `E_z` update gaining `−(dt/(ε₀·dA))·I_L^{n+1/2}`. Presents `jωL`. Stable
  (an inductor adds no CFL constraint).
- **Lumped resistor** → unchanged (the validated `pure_resistor` path).
- **Series R-L-C** (the F2.3 series branch) → the canonical combined series-RLC
  `E_z` update (resistor + accumulated inductor current + capacitor-voltage
  state), the Taflove lumped-RLC formula, sharing the same `dz/dA` field coupling.

Expose via the existing builder surface (a `with_two_way()`-style opt-in or a new
per-element constructor) without breaking `series_rlc` / `pure_resistor` /
`with_two_way` callers. ADR-0116's stable two-way *resistor* behaviour is retained.

`lumped_rlc_twoway_001` gains **asserted** reactive |Γ| (shunt-L, shunt-C,
series-RLC within a loose tol of the analytic `−Z₀/(2Z_L+Z₀)` shunt law at
4/6/9 GHz), keeping the resistor-exact + stability asserts. Container-iterated;
GREEN before merge; never weakened.

## Consequences

**Ships:** physically-correct reactive lumped-element FDTD loading — the canonical
method, validated. **Unblocks F2.3**: the L‖C resonators load the line, the FDTD
|S21| acquires its band-pass shape, `fdtd_lumped_001` passes → the goal's "EM
simulation" component completes (5/6).

**Gate:** `lumped_rlc_twoway_001` GREEN with the asserted reactive |Γ|; F2.3's
`fdtd_lumped_001` GREEN on top (separate follow-up merge).

**Escape hatch:** if the series-RLC combined case resists the loose tol after the
shunt-L / shunt-C cases pass, ship the shunt cases (which give F2.3 its dominant
selectivity) and defer the combined series-RLC as a documented follow-on — do not
weaken or fake. The shunt L/C alone is the increment's floor.

**Not in scope:** SRF/ESR vendor parasitics (F2.1b); the F2.3 board sim itself
(rides on this); multi-port S-params beyond S21.

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-4-canonical-lumped-lc-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-4-canonical-lumped-lc.md`.
- Taflove & Hagness, *Computational Electrodynamics: The FDTD Method*, 3rd ed.,
  Ch. on lumped circuit elements (Piket-May, Taflove & Baron 1994): the canonical
  lumped resistor / capacitor (effective-permittivity) / inductor (accumulated
  current) / series-RLC `E` updates.
