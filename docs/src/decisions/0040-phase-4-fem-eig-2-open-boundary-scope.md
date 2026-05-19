# ADR-0040: Phase 4.fem.eig.2 scope — open-boundary FEM (ABC + wave ports)

## Status

Accepted — 2026-05-19 (spec + plan; implementation deferred to
follow-up tracks).

## Context

Phase 4.fem.eig.0 shipped a real-valued, lossless, free-space
(`ε_r = μ_r = 1`) eigensolver on tetrahedral Whitney-1 Nedelec edge
elements (ADR-0029, ADR-0032). Phase 4.fem.eig.1 lifted it to lossy
single-pole dispersive `ε(ω)` via a Newton-Raphson `ω`-tracker
(ADR-0039) — `fem-eig-002` clears the lossy-SiO₂-filled cavity at
1.3e-3 Re(f) and 3e-3 Im(f). Both solvers enforce PEC tangential-`E`
zero on the *entire* exterior boundary via Dirichlet row/column
elimination — the only boundary the FEM stack knows how to apply.

Real microwave / antenna / accelerator workloads need open
boundaries: waveguide-fed `S_{11}` analysis, ABC-terminated radiation
domains, coaxial / SMA ports, and driven Q-extraction on iris-coupled
cavity filters. `CLAUDE.md` §4 forbids shipping a solver feature
without a published-benchmark validation case; the closed-cavity
gates `fem-eig-001` and `fem-eig-002` do not exercise an open
boundary or a wave-port at all.

Phase 1.3.1.1 already ships the **2-D Nedelec cross-section
eigensolver** in `crates/yee-mom/src/eigensolver/`, with a public
`NumericalCrossSection::e_tangential_at(x, y)` accessor returning the
dominant-mode tangential profile. Reusing it on the FEM port face is
the natural Phase 4.fem.eig.2 deliverable.

The radiation-condition surrogate has two textbook options:

1. **1st-order Engquist–Majda ABC** (Engquist & Majda, *Math. Comp.*
   1977) — a local boundary-term contribution `+ j k₀ (n̂ × E)`
   added per face to the curl-curl bilinear form. Reflection floor
   `~ −40 dB` at normal incidence; weaker off-normal.
2. **PML / UPML / CFS-PML** — a stretched-coordinate or anisotropic
   absorbing layer that suppresses reflection to numerical-floor
   levels. Requires extending the FEM domain with absorbing tetrahedra
   and re-doing the Nedelec basis on stretched coordinates.

This ADR records the v2 choice between them.

## Decision

Phase 4.fem.eig.2 ships **1st-order Engquist–Majda ABC + single-mode
modal wave-ports** only. PML / 2nd-order ABC / Higdon ABC are deferred
to Phase 4.fem.eig.2.5 if and when a published case shows up that
1st-order Engquist–Majda cannot meet at the required tolerance.

Six load-bearing decisions:

1. **1st-order Engquist–Majda ABC per spec §4.2.** Per ABC-tagged
   face, add `+ j k₀ (1/μ_r) · area · (n̂ × N_i) · (n̂ × N_j)` to the
   global complex stiffness matrix. Promotes the eigenproblem from
   real to complex-symmetric (not Hermitian) even with real `ε_r`;
   this is the same mechanism by which radiation loss appears as a
   negative imaginary eigenvalue. Reflection floor `~ −40 dB` is
   accepted as the v0 physics floor.
2. **PML / 2nd-order ABC deferred.** PML requires re-doing the
   curl-conforming basis on stretched coordinates; 2nd-order
   Engquist–Majda adds an auxiliary unknown per ABC face. Neither
   carries its weight on the v0 gate (`fem-eig-003` is a uniform-mode
   stub where the 1st-order floor is well within `−35 dB`). Phase
   4.fem.eig.2.5 reserved.
3. **Modal wave-port via Phase 1.3.1.1
   `NumericalCrossSection::e_tangential_at`.** Single source of
   modal-profile truth across MoM and FEM. v0 samples `e_mode` per
   Gauss point on the FEM port face (nearest sample on the 2-D
   cross-section eigensolver mesh); cubic / barycentric
   interpolation lands in Phase 4.fem.eig.2.0.1 if `fem-eig-003`
   shows interpolation-floor error.
4. **Single dominant mode per port.** TE_{10} for rectangular
   waveguide, TEM for coax. Higher-order modes are *captured* in the
   reflection spectrum but not *driven*. Multi-mode incident
   excitation is Phase 4.fem.eig.2.0.2.
5. **Existing closed-cavity API unchanged.** An `OpenBoundarySolver`
   with empty `abc_faces` and empty `ports` is ill-posed (no
   excitation) and is rejected at construction; the v0
   `FemEigenAssembly` and v1 `DispersiveSolver` surfaces are
   strictly unchanged. Every `fem-eig-001` and `fem-eig-002` caller
   compiles and passes unmodified.
6. **`fem-eig-003` is the production target.** WR-90 air-filled stub
   (22.86 × 10.16 × 30 mm), one end PEC, other end ABC, TE_{10}
   wave-port at the open end. `|S_{11}(f)|` across 50 sweep points
   in 8–12 GHz matches Pozar §3.3 closed-form within ±0.5 dB on
   the *sweep shape*; absolute floor falls in the `[−45, −35] dB`
   window per the 1st-order Engquist–Majda physics floor.

CPU-only, single-threaded, FP64 complex. No GPU. Single incident
mode per port. Scalar isotropic real `ε_r`, `μ_r` on the driven
sweep (combining v2 with v1's dispersive Newton tracker is Phase
4.fem.eig.2.1). Complex-valued eigenproblem stays the default for
any cavity with an ABC face (radiation loss is structurally complex).

## Consequences

- **Driven open-region workloads become reachable.** S-parameter
  analysis of waveguide stubs, irises, coax-fed dipoles inside FEM
  boxes, and slot-antenna feed problems all land on the same
  `OpenBoundarySolver` surface.
- **The existing `fem-eig-001` and `fem-eig-002` gates stay green.**
  v0/v1 callers do not see any API drift; `FemEigenAssembly` and
  `DispersiveSolver` are unchanged. New consumers explicitly
  construct an `OpenBoundarySolver`.
- **The complex-valued eigenproblem stays the default for any
  cavity with an ABC face.** Even with real `ε_r`, an ABC face
  promotes `K(ω)` to complex-symmetric. This matches the v1
  dispersive-cavity behaviour and reuses the Phase 4.fem.eig.1
  complex sparse LU surface unchanged.
- **1st-order Engquist–Majda reflection floor is the load-bearing
  correctness risk.** The `fem-eig-003` gate window
  `[−45, −35] dB` accepts the documented floor; tighter floors
  require Phase 4.fem.eig.2.5 (2nd-order ABC or PML). The gate
  asserts the *sweep shape* against Pozar §3.3 within ±0.5 dB,
  decoupling the absolute floor from the shape canary.
- **Cross-lane consumption of `yee-mom::NumericalCrossSection`
  becomes load-bearing for FEM.** The FEM port reads `e_mode` via
  the existing public `e_tangential_at` API; no changes to MoM.
  Any future MoM-side API change must keep that accessor stable
  or break this gate.
- **PML / 2nd-order ABC stays available** as a Phase 4.fem.eig.2.5
  upgrade path if 1st-order Engquist–Majda cannot meet some future
  validation tolerance. No code freeze on that option.

## References

- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
- `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-2-open-boundary.md`
- ADR-0029 — Phase 4.fem.eig.0 scope lock (this ADR's grandparent).
- ADR-0032 — Phase 4.fem.eig.0 implementation plan.
- ADR-0039 — Phase 4.fem.eig.1 dispersive scope (this ADR's parent).
- B. Engquist and A. Majda, "Absorbing boundary conditions for the
  numerical simulation of waves", *Math. Comp.* 31 (1977),
  pp. 629–651.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd
  ed., Wiley 2014, Ch. 10 (driven FEM analysis; §10.4 ABC,
  §10.5 wave-port modal decomposition, §10.7 S-parameters).
- D. M. Pozar, *Microwave Engineering*, 4th ed., 2012, §3.3
  (waveguide TE/TM modal characterisation; closed-form
  propagation constants).
- W.-J. Beyn, *Linear Algebra Appl.* 436 (2012) — Phase
  4.fem.eig.1.5 reserved.
- J.-P. Berenger, *J. Comput. Phys.* 114 (1994) — PML reference,
  Phase 4.fem.eig.2.5 reserved.
- CLAUDE.md §3, §4.
