# ADR-0021: Phase 2.fdtd.5.2 extends TF/SF corrections to all six box faces

**Status:** Accepted (supersedes part of ADR-0014's "finite-box deferred" plan)
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0014 shipped TF/SF as a *slab* (corrections on the `±x` faces
only, lateral `±y` / `±z` faces collapsed into the grid PML)
because at `+x` normal incidence the four lateral face corrections
vanish identically. The slab achieved 2676× (68.5 dB) TF/SF contrast
on the validation grid and was sufficient for plane-wave scattering
problems where lateral scattering can be probed by NTFF or far-field
extrapolation rather than by a clean SF region around the scatterer.

ADR-0014 also flagged the scope gap: a *finite TF/SF box* — TF
region tightly wrapping a scatterer with a clean SF region in all
three directions — needs the four lateral-face corrections. That
work was queued as "Phase 2.fdtd.5.1, deferred." Phase 2.fdtd.5.1
landed the i-face-only finite-box variant (TF region of finite
extent in `x`, full grid in `y` and `z`). That step ran but produced
only 6× (~15.6 dB) contrast on the finite-box test — well below the
≥1000× textbook quality bar and below even the 10× minimum gate.

The 5.1 contrast was bounded by the missing lateral-face corrections.
At `+x` propagation with `E_z` polarisation, four of the twelve
Yee-stencil curl components straddle the TF/SF boundary:

- `∂E_z/∂x` at the i-face — applied in 5.1
- `∂H_y/∂x` at the i-face — applied in 5.1
- `∂E_z/∂y` at the j-face — **missing** in 5.1, applied in 5.2
- `∂H_y/∂z` at the k-face — **missing** in 5.1, applied in 5.2

Every other curl component pairs an `E_z` or `H_y` field with a
zero-incident neighbour and needs no correction at this polarisation.
With four of the four straddled-stencil corrections in place, the
finite-box test should reach the FP64 floor rather than the 15.6 dB
soft limit that 5.1 produced.

A second issue surfaced during the 5.2 implementation: an
**off-by-one in the 5.1 i-face k-loop** (`k0..k_hi` exclusive,
dropping the upper TF slice). This was a no-op on slab geometry
(where `k1 == nz` so the inclusive and exclusive bounds coincide)
and therefore undetected by the Phase 2.fdtd.5 slab regression
test. It became visible — and fixable — only when the j/k face
corrections gave the contrast metric enough dynamic range to expose
boundary anomalies.

The structural question for **this** ADR was whether to ship 5.2 as
a focused fix (extend the kernel from i-face to all six faces, fold
in the off-by-one) or to wait for Phase 2.fdtd.5.3 (arbitrary
oblique incidence, ADR-0026) and ship a unified oblique-capable
kernel. The walking-skeleton-first principle (CLAUDE.md §3) selects
the focused fix: a finite-box normal-incidence kernel at
textbook-quality contrast is a strict prerequisite for the oblique
kernel (oblique reduces to normal-incidence as `θ → π/2, φ → 0`,
and the normal-incidence regression has to be bit-clean before
oblique scaffolding can build on it).

## Decision

`yee-fdtd` Phase 2.fdtd.5.2 extends the TF/SF correction kernel
from the i-face-only finite-box variant (5.1) to **all six box
faces**, and fixes the 5.1 i-face k-loop off-by-one. After the
change, the finite-box test reaches the FP64 floor.

**Six-face correction map** (for `+x` propagation, `E_z`
polarisation):

| Face | Update term      | Module function              | Phase |
|------|-----------------|------------------------------|-------|
| i0   | `H_y`: `∂E_z/∂x` | (existing)                   | 5.1   |
| i1   | `E_z`: `∂H_y/∂x` | (existing)                   | 5.1   |
| j0/j1| `H_x`: `∂E_z/∂y` | `correct_h_jface_plus_x`     | 5.2   |
| k0/k1| `E_x`: `∂H_y/∂z` | `correct_e_kface_plus_x`     | 5.2   |

The k-face `E_x` correction's i-range is `[i0, i1−1]` rather than
`[i0, i1]`: `H_y` on the back i-boundary (`i = i1`) is SF, so its
`∂H_y/∂z` stencil never straddles the k-face. The j-face `H_x`
correction has no analogous narrowing because both `H_x` and `E_z`
live on integer-x. These index ranges are derived inline in the
design-notes block at the top of `sources.rs` (the lane is
`crates/yee-fdtd/**`, so there is no separate spec under
`docs/superpowers/specs/`).

**Validation gates tightened in lockstep.**

- Finite-box `plane_wave_finite_box.rs`: contrast gate ≥ 100× (was
  ≥ 5× to accommodate the i-face-only 5.1 value of 6×). Empirical
  contrast on the 80³ vacuum grid is **~7.5×10¹⁴ (~298 dB)** —
  effectively roundoff. The 100× threshold is a guardrail against
  regressions in either the new j/k kernels or the underlying
  i-face kernel.
- Slab `plane_wave_propagation.rs`: contrast gate ≥ 1000× (was
  ≥ 10×). Empirical value remains 2676× — unchanged by the j/k face
  work, because the j/k corrections are no-ops on slab geometry
  (the SF rows they target sit outside the `H_x` / `E_x` extent).
  The tightened gate guards against future regressions in the
  i-face kernel or in lateral CPML.

## Alternatives considered

1. **Stick with slab-only TF/SF (cancel 5.1 / 5.2).** Rejected.
   CLAUDE.md §10 had documented the slab-only limitation as a
   workaround, but the limitation forecloses an entire class of
   tutorials (RCS / bistatic scattering from a finite object with a
   clean SF region around it). The fix is one bounded sub-project,
   not a permanent gap.
2. **Berenger-style absorbing layers around the TF region instead
   of a finite TF/SF box.** Rejected. This is a different
   architecture — replacing the TF/SF boundary with an absorbing
   region — and not a Phase 2.fdtd.5.x extension. CPML already
   provides the absorbing role on the grid boundary; layering an
   absorber around an internal TF region duplicates that role
   without a clear win. If a future use case wants it, it lives in
   a different ADR.
3. **Ship 5.2 as part of Phase 2.fdtd.5.3 oblique-incidence (skip
   the focused fix).** Rejected. Oblique-incidence (ADR-0026)
   reduces to normal-incidence in its limiting case, and the
   oblique kernel's normal-incidence regression has to be
   bit-clean. Shipping 5.2 first gives the oblique work a stable
   regression target; bundling them would conflate two distinct
   correctness issues.

## Consequences

**What becomes easier:**

- **Finite-box TF/SF is now textbook-quality.** A scatterer placed
  inside a finite TF region with a clean SF region on all sides
  can be probed without lateral-PML interference. The
  Taflove-canonical 3-D TF/SF formulation is now the *default*
  finite-box behaviour in `yee-fdtd`.
- **The slab geometry of ADR-0014 is preserved bit-for-bit.** No
  slab regression test or slab tutorial changes behaviour.
- **Phase 2.fdtd.5.3 oblique incidence (ADR-0026) has a clean
  regression target.** Its normal-incidence dispatch
  (`legacy_plus_x` path) compares against the 5.2 finite-box
  result.

**What becomes harder:**

- **Future TF/SF kernel work has to be tested against the tighter
  gates.** The 5× → 100× finite-box gate and the 10× → 1000× slab
  gate are a deliberately small fraction of the empirical headroom
  (~7×10¹² and ~2.7× respectively); regressions on either kernel
  are now visible to CI.
- **The i-face k-loop off-by-one is part of the public history.**
  The 5.1 i-face-only finite-box behaviour was technically incorrect
  on the upper k slice; anyone relying on 5.1 output (none, as 5.1
  was a single-commit intermediate) would see a behaviour change at
  5.2.

**What's now closed off:**

- A "five-face-and-a-half" intermediate variant. The six faces are
  all-or-nothing; partial correction kernels were rejected in 5.1
  on the same grounds they are rejected now.
- Quietly leaving the off-by-one in place. The 5.2 fix is the
  authoritative version of the i-face k-loop bounds.

## References

- `crates/yee-fdtd/src/tfsf/sources.rs` — six-face correction
  kernel; inline design-notes block at the top documents the
  straddled-stencil derivation and the j/k-face / k-face index
  ranges.
- `crates/yee-fdtd/tests/plane_wave_finite_box.rs` — 100× gate;
  empirical ~7.5×10¹⁴.
- `crates/yee-fdtd/tests/plane_wave_propagation.rs` — 1000× slab
  gate; empirical 2676×.
- Commits ef19cda (j/k-face corrections), 87caca5 (tightened gates),
  00ffe4a (Track PPPP merge).
- A. Taflove and S. C. Hagness, *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method*, 3rd ed., Artech
  House, 2005, §5.10 (3-D TF/SF formulation).
- ADR-0014 — TF/SF slab; this ADR closes the "finite-box deferred"
  scope gap that ADR-0014 explicitly queued.
- ADR-0026 — Phase 2.fdtd.5.3 oblique incidence; consumes this
  ADR's normal-incidence regression as its bisection reference.
- CLAUDE.md §10 — TF/SF status; the slab-only caveat is now
  superseded by this ADR.
