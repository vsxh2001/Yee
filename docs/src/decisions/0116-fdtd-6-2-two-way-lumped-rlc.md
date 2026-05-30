# ADR-0116: Phase 2.fdtd.6.2 — stable two-way lumped RLC port

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0017 (`LumpedRlcPort`), ADR-0115 (F2.3 lumped FDTD — blocked on
this), the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

F2.3 (lumped-LC FDTD EM sim, ADR-0115) is blocked: `LumpedRlcPort::series_rlc`
cannot model a low-loss reactive S-parameter element. Container experiments
confirmed the `l>0` branch is **one-way** (lumped current never feeds back into
`E_z` ⇒ a source-free inductor is inert; a shunt L‖C never resonates) and the
only two-way arm (pure C, `l=0`) is **unstable below ~196 Ω** ESR (≈ η₀/√3).
The two-way RLC port was always the un-shipped Phase 2.fdtd.6.2.

## Decision

Implement the **stable, two-way** lumped-element `E_z` update (Piket-May /
Taflove–Hagness semi-implicit formulation) for a series R-L-C in
`LumpedRlcPort::correct_e`: solve the `E_z^{n+1}` and lumped branch state
implicitly together so the lumped current couples back into the field
(two-way) and the update is unconditionally stable for any R≥0, L≥0, C>0. Keep
the public constructors/signatures (`series_rlc`, `pure_resistor`) and the
resistor's validated behaviour; this is an internal correctness fix.

Gate `lumped_rlc_twoway_001` (`#[ignore]`'d, CI `--release`): stability of a
low-loss reactive element (no NaN; the ≥196 Ω limit gone) + two-way correctness
(single lumped load reflects with Γ = `(Z_L−Z0)/(Z_L+Z0)` matching analytic).
Iterated in the bounded container; GREEN in CI on the branch before merge.

## Consequences

**Ships:** a correct, stable, two-way lumped-element FDTD port — a core capability
beyond filters. **Unblocks F2.3**: its driver + gate (`311a796`,
`feature/filter-f2-3-lumped-fdtd`) pass **unchanged** once this lands, completing
the lumped-LC "EM simulation" goal component.

**Gate:** `lumped_rlc_twoway_001` GREEN in CI before merge; existing lumped tests
non-regressed.

**Not in scope:** parallel-RLC primitive (composed from two series ports in F2.3);
the F2.3 board sim (rides on this).

---

## References
- `docs/superpowers/specs/2026-05-30-fdtd-6-2-two-way-lumped-rlc-design.md`;
  `docs/superpowers/plans/2026-05-30-fdtd-6-2-two-way-lumped-rlc.md`.
- Piket-May, Taflove & Baron (1994); Taflove & Hagness, *Computational
  Electrodynamics*, lumped-element FDTD chapter.
