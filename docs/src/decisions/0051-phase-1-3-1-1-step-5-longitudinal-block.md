# ADR-0051 — Phase 1.3.1.1 step 5: mixed (E_t, E_z) longitudinal block wire-in

**Status:** Accepted
**Date:** 2026-05-23
**Context Phase:** 1.3.1.1 step 5 (quasi-TEM microstrip wave-ports)

## Context

`NumericalCrossSection::solve` ships a **transverse-only** Nedelec
eigensolve. That is exact for homogeneous (air-filled) guides but wrong
for the **inhomogeneous** cross-sections the cross-section eigensolver
exists to handle (microstrip, partial dielectric fill, CPW), where the
longitudinal field `E_z` couples through the dielectric interface. The
mixed-formulation longitudinal element matrices (`local_a_zz`,
`local_b_zz`, `local_b_ze`) are already staged and unit-tested in
`eigensolver/assembly.rs` but unused.

Three design questions:

1. **Formulation** — full mixed `(E_t, E_z)` block eigenproblem
   (Lee-Sun-Cendes 1991) vs a scalar quasi-static / transverse-only
   approximation with an ε_eff correction.
2. **Z_w definition** — keep the TE-mode `η₀k₀/β` approximation vs a
   numerical line-integral extraction off the eigenvector.
3. **Inhomogeneous validation reference** — published transcendental
   dielectric-loaded-guide root vs an internal physics inequality.

## Decision

1. **Full mixed `(E_t, E_z)` block formulation** — assemble the staged
   longitudinal blocks into `A x = β² B x` with `x = [E_t; E_z]`,
   reusing the existing element matrices. No new approximation layer.
2. **Numerical Z_w** — voltage line-integral + power, reducing to the
   TE form on the homogeneous guide as a regression guard.
3. **Validation:** published transcendental dielectric-slab-loaded-guide
   β as the primary gate (loose ≤5% tolerance per the
   placeholder-tolerance policy), with a physics-inequality +
   regression fallback (`k₀ < β < k₀√ε_r,max`) mirroring the original
   spec's septum Case B if the transcendental reference trips the
   implementation escape-hatch.

## Rationale

(1) The mixed formulation is the *reason the (E_t, E_z) element
matrices were staged* — steps 2-3 deliberately froze the transverse
path first (walking-skeleton discipline) and deferred the longitudinal
block to here. Using a scalar ε_eff approximation instead would
duplicate Hammerstad-Jensen curve-fits the spec explicitly rejects
(they break for thick conductors / arbitrary stack-ups).

(2) The TE-mode `η₀k₀/β` Z_w is only correct for the air-filled guide;
on a dielectric stack-up the wave impedance is genuinely different and
a line-integral extraction is the standard quasi-TEM definition (Jin
§8.4). Keeping the regression-reduction to the TE form guards the
formula.

(3) CLAUDE.md §4 requires a published benchmark; the dielectric-loaded
guide has one (Pozar §3 / Collin). The inequality fallback keeps the
gate complete and shippable even if the transcendental reference solver
is deferred to step-5.1 — it is the same gate class as the existing
septum Case B.

## Consequences

* New `assemble_mixed` + `AssembledMixed` (block pencil) and
  `solve_dense_mixed` in `eigensolver/`; `assemble_transverse` /
  `solve_dense` retained for the homogeneous path + tests.
* New `mode_profile_ez` field on `NumericalCrossSection`; the
  `e_tangential_at` + `Numerical2D` wave-port RHS contract is preserved
  (the transverse field form is unchanged; its value shifts on
  inhomogeneous guides — the intended effect).
* Highest implementation risk is the block sign/placement convention of
  the staged element matrices; the homogeneous-guide β regression
  (DoD-V1) is the canary.
* Sparse mixed solve, CPW multi-conductor Z₀ matrix, and the yee-py E_z
  binding are out of scope (later steps).

## References

* Lee, Sun, Cendes, IEEE MTT 39(8), 1991.
* Jin, *FEM in EM* 3rd ed. §8.4.
* Pozar, *Microwave Engineering* 4th ed. §3.
* Step 5 spec + plan
  `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-5-longitudinal-block-design.md`,
  `docs/superpowers/plans/2026-05-23-phase-1-3-1-1-step-5-longitudinal-block.md`.
* ADR-0022 / ADR-0023 — Phase 1.3.1.1 eigensolver spec + TriMesh2D stub.
