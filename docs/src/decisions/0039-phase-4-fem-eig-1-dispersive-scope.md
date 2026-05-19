# ADR-0039: Phase 4.fem.eig.1 scope — dispersive `ε_r(ω)` on tet FEM eigensolver

## Status

Accepted — 2026-05-19 (spec + plan; implementation deferred to
follow-up tracks).

## Context

Phase 4.fem.eig.0 shipped a real-valued, lossless, free-space
(`ε_r = μ_r = 1`) eigensolver on tetrahedral Whitney-1 Nedelec
edge elements (ADR-0029, ADR-0032). Validation gate
`fem-eig-001` clears TE_{101} at 9.660 GHz on a WR-90-based
cavity to 0.09 % rel.err. That solver is correct for empty
metallic enclosures and only for empty metallic enclosures.

Real cavity-filter, dielectric-resonator-antenna, and
accelerator workloads need lossy, frequency-dependent media
inside the cavity. `CLAUDE.md` §4 forbids shipping a solver
feature without a published-benchmark validation case; lossy
materials add a complex eigenvalue path which `fem-eig-001`
does not exercise at all.

Phase 2.fdtd.3 (ADR — not yet recorded here, see
`crates/yee-fdtd/src/material.rs`) already ships the
single-pole **Drude / Lorentz / Debye ADE `Material` enum**
with a `permittivity(omega) -> Complex64` accessor. Reusing it
on the FEM side is the natural Phase 4.fem.eig.1 deliverable.

The mathematical challenge is that
`K(ω) e = (ω/c)² M(ω) e` is **nonlinear in ω** — the matrix
entries depend on the unknown angular frequency. Two textbook
solver-side options: Newton-Raphson tracking along a frequency
sweep, or Beyn 2012 / Sakurai–Sugiura contour-integral
nonlinear eigensolve. This ADR records the v1 choice between
them.

## Decision

Phase 4.fem.eig.1 ships **Newton-Raphson `ω`-tracking only**.
Beyn 2012 is deferred to Phase 4.fem.eig.1.5 if and when a
Newton-with-bisection-fallback failure mode shows up that the
documented escape hatches cannot recover from.

Six load-bearing decisions:

1. **Newton tracker per the spec's §4.2 pseudocode.** Outer
   loop assembles `K(ω)`, `M(ω)`, solves the *linearised*
   complex generalised eigenproblem at trial ω, takes a
   Hellmann–Feynman analytic Newton step in ω until the
   fixed-point condition `θ = (ω/c)²` holds. Converges in 3–5
   iterations from a real-valued warm-start for the published
   gate. Bisection fallback wraps the case where
   `|F'(ω)| → 0`.
2. **Beyn 2012 deferred.** Contour-integral nonlinear
   eigensolve needs a complex-contour-integration framework
   Yee does not have and produces all modes inside a contour
   regardless of warm-start. The Newton tracker handles one
   mode at a time given a warm-start, which is what
   `fem-eig-002` needs.
3. **`Material` enum reused verbatim from Phase 2.fdtd.3.**
   Plan step D3 relocates the type from `yee-fdtd::material`
   to `yee-core::material` with a `pub use` re-export so all
   existing yee-fdtd callers compile unchanged. Single source
   of truth for Drude / Lorentz / Debye across the workspace.
4. **Existing free-space yee-fem path stays as-is.** The
   complex-coefficient lift of `assemble_tet_element` accepts
   `Complex64` ε / μ; v0 free-space callers do an explicit
   `Complex64::from(f64)` lift at the boundary. Phase
   4.fem.eig.0's `fem-eig-001` gate stays green unmodified.
5. **Complex inverse-iter is a search-and-replace lift.**
   `ComplexInverseIterEigen` sits next to `InverseIterEigen`
   behind a new `SparseEigenComplex` trait — a peer of the v0
   `SparseEigen<f64>` rather than a parametric unification.
   Two traits is cleaner than one parametric trait given the
   complex symmetric (not Hermitian) inner-product
   conventions and the pivoting differences in complex LU.
6. **`fem-eig-002` is the production target.** Lossy
   single-pole-Drude-SiO₂-filled cavity (a = 10 mm, b = 5 mm,
   d = 20 mm), ±0.5 % on Re(f), ±5 % on Im(f). Tighter Re(f)
   bound than fem-eig-001's ±0.3 % because v0 has 0.09 %
   headroom; ±5 % on Im(f) is consistent with Pozar's
   published wall-loss Q tolerances.

CPU-only, single-threaded, FP64 complex. No GPU.
Single-pole dispersion only. Scalar isotropic complex `ε(ω)`,
real `μ_r` (magnetic dispersion is Phase 4.fem.eig.1.2). PEC
closed cavity only.

## Consequences

- **Lossy resonator workloads become reachable.** Q-factor
  extraction, lossy-filter cavity tuning, dispersive DRA modes
  (after Phase 4.fem.eig.3) all land on the same surface.
- **Existing `fem-eig-001` gate stays green.** v0 callers do
  the `Complex64::from(f64)` lift at the boundary; no v0
  semantics change.
- **Newton convergence is the load-bearing correctness
  risk.** Bisection fallback inside `track_mode` is the
  documented escape hatch; the fem-eig-002 gate asserts the
  fallback is **not** triggered on the published case, which
  makes the gate the canary for Newton-basin failures.
- **`Material` relocation to `yee-core` cleans up a latent
  bug.** Until now `yee-fdtd::material::Material` was the
  single source of truth across crates that all secretly
  depended on each other; the relocation makes that explicit.
- **Beyn 2012 stays available** as a Phase 4.fem.eig.1.5
  upgrade path if Newton-with-bisection cannot converge on
  some future validation case. No code freeze on that option.

## References

- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
- `docs/superpowers/plans/2026-05-19-phase-4-fem-eig-1-dispersive.md`
- ADR-0029 — Phase 4.fem.eig.0 scope lock (this ADR's parent).
- ADR-0032 — Phase 4.fem.eig.0 implementation plan.
- J.-M. Jin, *The Finite Element Method in Electromagnetics*,
  3rd ed., Wiley 2014, §9.5 (lossy-material FEM eigenvalue
  problems and Hellmann–Feynman differentiation).
- D. M. Pozar, *Microwave Engineering*, 4th ed., 2012, §3.1
  (lossy-waveguide propagation constants), §6.3 (wall-loss Q).
- W.-J. Beyn, *Linear Algebra Appl.* 436 (2012) — Phase
  4.fem.eig.1.5 reserved.
- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed.,
  Artech 2005, Ch. 9 — Phase 2.fdtd.3 ADE reference.
- `crates/yee-fdtd/src/material.rs` — `Material` enum reused
  verbatim.
- CLAUDE.md §3, §4.
