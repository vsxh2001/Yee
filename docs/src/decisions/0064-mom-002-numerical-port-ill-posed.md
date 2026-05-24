# ADR-0064 — mom-002 numerical wave-port frame-mapping is ill-posed for planar MoM (stop)

**Status:** Accepted (negative decision — do NOT pursue)
**Date:** 2026-05-24
**Context Phase:** MoM beachhead follow-on (closes the ADR-0061 frame-mapping question)

## Context

ADR-0061 (Phase 1.3.1.2-B) shipped a numerical quasi-TEM microstrip
wave-port into the mom-002 line, reached `|Z_in| = 378 Ω` (vs 51 Ω HJ,
674 Ω delta-gap), and localized the residual to a coordinate-frame issue,
recommending a "2-D-cross-section → 2.5-D-RWG frame mapping" follow-on. A
read-only feasibility scoping was run before committing effort to that
follow-on (the standing rule: approach grind-risky tracks as bounded
experiments / verify feasibility first).

## Decision

**Do NOT pursue the frame-mapping follow-on. The mom-002 numerical
wave-port (the `ModalDistribution::Numerical2D` arm) is the wrong vehicle
for a microstrip port; its residual is partly intrinsic to the planar
formulation.** Stop chasing this port mechanism; rotate for breadth.

## Rationale (verified against source)

(1) **The dominant microstrip mode is orthogonal to the port basis.** The
mom-002 port is 16 y-aligned RWG edges at a single x = L/2, z = 0
(`crates/yee-validation/src/lib.rs:761`; port-edge rule
`crates/yee-mom/src/basis.rs:174`). The physical transverse plane there is
(y, z); the quasi-TEM mode's *defining* field is the substrate-normal
`E_z`. The Numerical2D RHS projects the modal field onto the port-edge
**tangent**, which is pure MoM-y (`ports.rs:983`) — so `E_z` projects to
~zero. A correct (y, z) relabel of the cross-section therefore drives the
RHS *toward zero*, making it worse. The current `378 Ω` arises only
because the diagnostic builds the cross-section in the *wrong* (x, y)
longitudinal frame, which happens to expose an in-plane component
(`ε_eff = 2.105`, not the line's 3.33 — it is not even the physical mode).

(2) **Fundamental dimensionality mismatch.** The planar MoM unknown is an
in-plane surface current `J_s` (x/y components only; RWG `eval`,
`basis.rs:241`). A microstrip quasi-TEM port aperture is a 2-D (y, z) face
whose mode lives out-of-plane; the planar "aperture" is a 1-D line of
in-plane current edges. There is **no DoF** that can carry or be excited by
the substrate-normal field. The frame identification at `ports.rs:979`
(`e_tangential_at(mid_x, mid_y)`, z dropped) cannot be fixed by a relabel.

(3) **Contrast with the validated WR-90 case.** The Numerical2D arm passes
<1% for WR-90 (`tests/wave_port_numerical_te10.rs`) because there the port
face *is* the 2-D cross-section and the TE₁₀ mode is fully in-plane —
both preconditions the microstrip geometry violates.

(4) **The kernel is not the bottleneck.** The Greens kernel is exonerated
to 1.83% of HJ (ADR-0061). So the residual is the port, and the port is
intrinsically limited here — not closable cheaply.

## Consequences

* The `Numerical2D` arm stays as-is (correct + validated for waveguide
  ports); it is **not** extended for microstrip. The `#[ignore]`'d
  `mom_002_numerical_waveport.rs` diagnostic stands as the documented
  finding.
* mom-002 stays at loose tolerance with the TEM-smoothed delta-gap as its
  production port; this is now understood as a planar-formulation limit,
  not a fixable wiring gap.
* The principled path to a true microstrip Z₀ (a new port formulation —
  aperture/frill reciprocity excitation, or TL-based Z₀ de-embedding from
  the solved line currents) is a **major multi-week track**, deferred; it
  is not opened now (the autonomous loop favors bounded, dispatchable
  tracks, and mom-002 is a known quagmire with an exonerated kernel).
* Near-term effort rotates to breadth instead.

## References

* The feasibility scoping (2026-05-24, read-only) — `ports.rs:933-994`
  (Numerical2D arm; frame ID `:979`), `:689-759` (`e_tangential_at`);
  `basis.rs:174,241`; `crates/yee-validation/src/lib.rs:708-784,829`;
  `tests/wave_port_numerical_te10.rs`; `tests/mom_002_numerical_waveport.rs:143`.
* ADR-0061 (the Phase-B result this closes), ADR-0036/0037 (the kernel
  exoneration), ADR-0059/0060.
