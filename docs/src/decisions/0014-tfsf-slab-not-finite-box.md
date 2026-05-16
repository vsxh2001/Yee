# ADR-0014: Phase 2.fdtd.5 ships a TF/SF *slab*, not a finite TF/SF box

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Phase 2.fdtd.5 introduces a **Total-Field / Scattered-Field (TF/SF)**
source region to `yee-fdtd`. TF/SF is the standard FDTD plane-wave
injection technique (Taflove & Hagness, *Computational Electrodynamics*
3rd ed., §6.2 "Three-dimensional total-field/scattered-field
formulation"): the simulation volume is partitioned into an inner
"total-field" region where Maxwell's equations are integrated with the
plane wave present, and an outer "scattered-field" region where only
the scattered (target-radiated) field lives. Consistency is maintained
by correction terms applied on the **TF/SF boundary surface** — for the
canonical Taflove formulation, a six-face rectangular box around the
inner region.

The original Phase 2.fdtd.5 scope as written into `ROADMAP.md` was
"full 6-face TF/SF box", on the basis that Taflove §6 derives the
algorithm in 3D for a closed box and the full box is the textbook
deliverable.

When the Phase 2.fdtd.5 brief was being decomposed, two scope
constraints became clear:

- **The brief explicitly limited Phase 2.fdtd.5 to a single
  propagation direction (`+x`, `Ez`-polarised).** This was the
  walking-skeleton choice from `CLAUDE.md §3`: ship the minimum that
  exercises the contract before generalising. A `+x` plane wave at
  normal incidence is the canonical first TF/SF test case in every
  FDTD textbook and is sufficient to validate the consistency
  bookkeeping.
- **The 6-face box adds non-trivial implementation surface that is
  unused at `+x` normal incidence.** The four "side" faces (the
  `±y` and `±z` faces of the box) carry corrections only when the
  injected wave has a non-zero component into / out of those faces.
  At `+x` propagation with a `Ez` field, the `±y` and `±z` corrections
  vanish identically — they contribute zero update on every time step.
  Only the `±x` faces (call them `i0` and `i1`) carry non-zero
  corrections.

Two structural responses present themselves:

1. **Implement the full 6-face box.** Write the corrections on all
   six faces, even though four of them evaluate to zero. Future
   oblique-incidence work then needs no new structure, only new
   1-D incident-field bookkeeping.
2. **Implement a TF/SF *slab*.** Restrict corrections to the two
   faces normal to the propagation direction (`i0` and `i1`), and
   let the TF region span the **full `j` and `k` extent of the
   grid**. The corrections on the side faces are not absent but
   *not present*: the TF region in `y` and `z` is the whole grid,
   so there are no `±y` or `±z` faces to correct against.

The two formulations are **physically equivalent at normal incidence**:
a slab and a box with zero side-face corrections produce bit-identical
fields. They are **not equivalent at oblique incidence** (a slab cannot
model a wave whose wavevector has non-zero `y` or `z` components,
because the wave would clip the lateral PML; a finite box can).

The choice is between **shipping a slab now** (smaller agent brief,
faster to validate) **and deferring the finite-box generalisation to
Phase 2.fdtd.5.1**, vs. shipping the full box up front.

The walking-skeleton-first principle (`CLAUDE.md §3`) settles it: ship
the slab first. Side-face corrections would roughly double the agent
brief — separate 1-D incident-field grids per face, separate sign
conventions for `H × n̂` on each face, separate index permutations for
the four side faces' `J_surface` arrays — and **none of that code
runs** at `+x` normal incidence. The walking skeleton stays minimal
and the finite-box generalisation becomes a focused follow-up.

The validation gate is **TF/SF contrast**: inject a plane wave in the
TF region and measure the ratio of TF energy to SF energy. For a
correctly-implemented TF/SF source, the SF region should contain only
**numerical leakage** (typically a few parts in `10⁴` of the TF
amplitude), giving a TF/SF energy ratio of `≥ 10×` (10 dB) at the
loose end and `≥ 1000×` (30 dB) at the textbook-quality end. The
gate is set at **≥ 10× (10 dB)** to leave headroom for grid-resolution
choices that future tutorials may pick.

The shipped slab implementation measured a TF/SF energy ratio of
**2676×** on the validation grid — equivalent to **68.5 dB** — well
past the 10× gate and within an order of magnitude of the Taflove
"textbook quality" reference.

## Decision

`yee-fdtd` Phase 2.fdtd.5 ships TF/SF as a **slab** geometry, not as
a finite box.

Concretely:

- The TF region is a slab `i ∈ [i0, i1]` spanning the full `j` and
  `k` grid extents. The SF region is `i < i0` and `i > i1`.
- **Corrections are applied only on the `i0` and `i1` faces.**
  - On `i0` (the back face of the slab from the wave's perspective):
    `E_y(i0, j, k)` and `E_z(i0, j, k)` updates get an additional
    `H × n̂` term from the incident `H_y^{inc}`, `H_z^{inc}` 1-D
    auxiliary grid.
  - Symmetric correction on `i1`.
- The `±y` and `±z` faces of the conceptual TF/SF box **do not
  exist** in this implementation. The slab abuts the lateral PMLs
  directly.
- The incident field is propagated on a **1-D auxiliary grid** along
  the propagation axis (`+x`), with its own CPML termination at both
  ends. This is the standard Taflove §6.2 1-D-incident-field
  technique, restricted to a single axis.
- The wave is hard-coded to **`+x` propagation, `Ez` polarisation**
  in Phase 2.fdtd.5. The slab `i0, i1` extents and the source
  waveform (Gaussian-modulated sinusoid by default) are
  configurable; the propagation axis and polarisation are not.

**Validation.** `crates/yee-fdtd/tests/tfsf_slab.rs` injects a
plane-wave pulse, integrates for a window long enough for the wave to
cross the slab, and measures
`E_TF² / E_SF²` integrated over the TF and SF regions. Gate:
`≥ 10×` (10 dB). Phase 2.fdtd.5 shipped at **2676× (68.5 dB)**.

**Phase 2.fdtd.5.1 (deferred follow-up).** Generalisation to a
finite TF/SF box with side-face corrections, parameterised
propagation direction `k̂`, parameterised polarisation
`ê ⊥ k̂`, and the four side-face `J_surface` correction arrays.
The slab implementation in Phase 2.fdtd.5 is a strict subset and
will continue to work as the `k̂ = +x̂, ê = ẑ` special case of
the Phase 2.fdtd.5.1 generalisation.

## Consequences

**What becomes easier:**

- The Phase 2.fdtd.5 walking skeleton ships with a clean validation
  gate (10× contrast, achieved 2676×) and a small, reviewable agent
  brief. The implementation is roughly half the code surface of a
  full 6-face box.
- The slab geometry is the right answer for the entire set of TF/SF
  tutorial cases that Phase 2 surfaces: a plane-wave pulse hitting
  a finite scatterer (sphere, dipole, microstrip patch), measured
  via scattered-field probes or NTFF-extrapolated far-field
  patterns. None of these need a finite TF/SF box; the lateral
  PMLs absorb the unscattered plane wave wings without artefact.
- The 1-D auxiliary incident-field grid (with its own CPML) is the
  exact same infrastructure that a future finite-box implementation
  needs along the propagation axis; nothing in Phase 2.fdtd.5 has
  to be rewritten to support 5.1.

**What becomes harder:**

- **Users cannot model a small TF box embedded inside a larger SF
  region.** The canonical case the slab cannot do is: a scatterer
  (say, a finite metallic object) placed more than ~100 cells from
  the PML, with a TF region tightly wrapping the scatterer so the
  scattered field can be probed in a clean SF region all the way
  out to the PML. With a slab, the TF region is the full lateral
  grid, so there is no SF region in the `y` or `z` directions —
  scattered probes in those directions see total field, not
  scattered field. This case requires Phase 2.fdtd.5.1.
- **Oblique-incidence and non-normal-polarisation use cases are
  blocked.** A wave with `k̂ = (cos θ, sin θ, 0)` would clip the
  lateral PMLs in the slab geometry; the 1-D auxiliary grid
  technique only generalises cleanly when there are matching side-
  face corrections to absorb the lateral components of the
  incident field. Users wanting RCS / bistatic scattering from a
  finite object must wait for 5.1.
- The "TF/SF slab" terminology is non-standard; Taflove uses
  "TF/SF region" or "TF/SF box" exclusively. The `yee-fdtd`
  documentation has to explicitly call out the slab restriction
  and link to the Phase 2.fdtd.5.1 plan.

**What's now closed off:**

- Implementing the four side-face corrections inside Phase 2.fdtd.5.
  They are out of scope by the brief and would have doubled the
  agent surface; they are a separate sub-project (5.1).
- A configurable propagation axis or polarisation in Phase
  2.fdtd.5. The hard-coded `+x` / `Ez` choice is the minimum that
  validates the slab plumbing; making it configurable without
  side-face corrections would invite a class of bugs where the
  configuration value is changed but the implementation
  silently keeps assuming `+x`.

## References

- `crates/yee-fdtd/src/tfsf/` — the slab implementation.
- `crates/yee-fdtd/src/tfsf/slab.rs` — `i0` / `i1` face
  corrections and 1-D auxiliary incident-field integration.
- `crates/yee-fdtd/tests/tfsf_slab.rs` — TF/SF contrast gate
  (≥ 10×; shipped value 2676× / 68.5 dB).
- `docs/src/theory/fdtd-details.md` — TF/SF chapter, including
  the slab vs finite-box distinction.
- A. Taflove and S. C. Hagness, *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method*, 3rd ed., Artech
  House, 2005, §6.2 (3-D TF/SF formulation) and §6.5 (1-D
  auxiliary incident-field source).
- D. M. Sullivan, *Electromagnetic Simulation Using the FDTD
  Method*, 2nd ed., IEEE Press, 2013, §3.2 (TF/SF in 1-D and 2-D,
  the slab special case under a different name).
- ADR-0008 — validation aggregator; `tfsf_slab.rs` reports
  through this.
- Phase 2.fdtd.5.1 (queued, not yet specced) — finite TF/SF box,
  oblique incidence, arbitrary polarisation.
