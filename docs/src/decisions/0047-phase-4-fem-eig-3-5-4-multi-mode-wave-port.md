# ADR-0047 — Phase 4.fem.eig.3.5.4 multi-mode wave-port API extension

**Status:** Accepted
**Date:** 2026-05-20
**Context Phase:** 4.fem.eig.3.5.4 (post-v3.5.3 single-mode wave-port,
ADR-0046)

## Context

ADR-0046 §Decision (5) blueprinted the v3.5.4 multi-mode wave-port
extension as the binding follow-on to v3.5.3's single-mode W1
attempt. The v3.5.3 measurement (`|S_{11}|(30 GHz) = 0.925644`)
matched the v3.5.2 CFS-PML floor (0.926) to four decimals — both
were one-dimensional modal projections onto the TE_{10} basis
vector, and the residual reflection is exactly the field content
orthogonal to that basis (TE_{20} at `f_c = 3.0 GHz` and TE_{01} at
`f_c = 15.0 GHz`, both propagating at 30 GHz on the +x face
cross-section).

The API question is **how to express the multi-mode basis in the
Rust surface**:

* **Option A** — extend `PortDefinition` to carry
  `Vec<PortMode>`. Every wave-port is structurally a modal-basis
  projection; the single-mode case is the degenerate one-element
  vector.
* **Option B** — introduce a parallel `MultiModeWavePort` type that
  coexists with `PortDefinition`. Single-mode call sites unchanged;
  multi-mode call sites use the new type. `FaceKind::WavePort(p)`
  resolves to either by table lookup.
* **Option C** — encode the modal basis as a free-function
  callback `Box<dyn Fn(world_point) -> Vec<(Vector3<f64>, Complex64)>>`
  returning `(e_t_mode_m, a_inc_mode_m)` tuples. No struct change;
  the caller flattens the basis into one closure.

## Decision

**Option A**: extend `PortDefinition` to
`{ modes: Vec<PortMode> }`. Add a `single_mode` constructor that
collapses the v3.5.3 / fem-eig-004 / fem-eig-005 call shape to one
line.

## Rationale

(1) **Structural truth.** Every modal wave-port is a basis
projection. The single-mode case is not a different *kind* of
port — it is the same kind with a trivial basis. Carrying that in
the type system avoids two parallel surfaces drifting (Option B
risks bug-fix asymmetry).

(2) **Forward compatibility.** Future absorbing-mode ports (Phase
4.fem.eig.3.5.6, Lee-Mittra 1997 §IV) plug into the same
`Vec<PortMode>` slot — an absorbing-mode termination is structurally
a modal basis with a tower of evanescent modes. Option C bundles
the modal structure into an opaque closure that hides the
per-mode introspection points (`beta_mode(ω)`, normalisation
inner products) downstream code needs for S-parameter
post-processing.

(3) **Migration cost.** `PortDefinition::single_mode(beta, e_t)`
mechanically rewrites the three v3.5.3 call sites; no semantic
change at any caller that doesn't want multi-mode behaviour.
Option B forces every consumer of `FaceKind::WavePort` to branch
on the lookup-table result type.

## Consequences

* Public API change on `crates/yee-fem`: `PortDefinition` field
  rename from `{ beta_mode, modal_e_t }` to `{ modes: Vec<PortMode> }`.
  Pre-v3.5.4 external consumers (none documented; the API is
  Phase-4 internal) need to adopt `PortDefinition::single_mode`.
* Assembly path (`crates/yee-fem/src/element/port.rs`) loops over
  modes; v3.5.3 numerics preserved because the single-mode loop
  body matches the v3.5.3 expression.
* S-parameter extraction picks the driving mode (`a_inc != 0`) for
  `S_{p, p}`. Multi-driving-mode dual-feed is reserved for
  v3.5.5; v3.5.4 errors out with `MultipleDrivingModes` if more
  than one mode carries `a_inc != 0`.
* yee-py multi-mode kwarg shape is **out of scope** for v3.5.4;
  the Python binding continues to construct single-mode
  `PortDefinition` via `single_mode` until v3.5.4.1 lands a
  `modes: [...]` kwarg.
* If `|S_{11}|(30 GHz) < 0.1` is not reached at the three-mode
  basis (TE_{10}, TE_{20}, TE_{01}), v3.5.5 adds TE_{11}; v3.5.6
  replaces the modal-basis projection with a Lee-Mittra
  absorbing-mode wave-port. Both extensions inherit the
  `Vec<PortMode>` slot.

## References

* ADR-0046 §Decision (5)
  `docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md`.
* Phase 4.fem.eig.3.5.4 spec
  `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-4-design.md`.
* Phase 4.fem.eig.3.5.4 plan
  `docs/superpowers/plans/2026-05-20-phase-4-fem-eig-3-5-4.md`.
* Jin, *FEM in EM*, 3rd ed., §10.6 multi-mode wave-port.
* Pozar, *Microwave Engineering*, 4th ed., §3.3 TE_{mn} field patterns.
