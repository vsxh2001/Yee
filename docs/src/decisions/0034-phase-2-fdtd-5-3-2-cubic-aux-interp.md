# ADR-0034: Phase 2.fdtd.5.3.2 — cubic Lagrange aux-grid interpolation clears the oblique TF/SF >1000× DoD

## Status

Accepted — 2026-05-18.

## Context

ADR-0026 shipped the Phase 2.fdtd.5.3 oblique-incidence TF/SF
kernel with two known caveats: (a) the 1-D auxiliary Yee grid
dispersion mismatched the 3-D Yee grid at oblique angles, and
(b) the 30°/45° `E ∥ ê_φ` test cleared a placeholder `>10×`
gate at ~14.5× contrast — well below the original `>1000×` DoD.
Phase 2.fdtd.5.3.1 (dispersion-matched aux step) was queued as
the suspected dominant remaining error term.

Track KKKKK (merge `f878bdd`) ran the face-stencil audit that
Phase 2.fdtd.5.3.2 was scoped to do, expecting to find a sign
or cross-section bug in `correct_h_oblique` / `correct_e_oblique`.
The audit found something else: the H_z k-range off-by-one fix
on the i/j-face stencil (commit d9f0a5a) was a necessary precursor
but turned out **not** to be the dominant residual SF leakage
contributor. The dominant contributor was the
**linear interpolation of the 1-D auxiliary incident-field grid**
in `sample_inc_e` / `sample_inc_h`. At `dx/λ = 0.05` the
per-sample linear-interp residual was ~0.1%, which compounds
across the 12-stencil correction on each box face and floors the
oblique contrast around 15×.

The fix is a 4-point cubic Lagrange interpolation via a new
`sample_aux_cubic()` helper (commit f9b9cef): per-sample
residual drops from ~0.1% to ~1e-4, ~60 dB of headroom. The
30°/45° `E ∥ ê_φ` contrast progression is:

- Phase 2.fdtd.5.3 ship (linear interp): **14.5×**
- Phase 2.fdtd.5.3.1 dispersion-matched aux step: **15.6×**
- Phase 2.fdtd.5.3.2 cubic interp: **1027× ≈ 60.2 dB → clears DoD**

## Decision

The Phase 2.fdtd.5.3.2 implementation locks four load-bearing
choices:

1. **Two-stage fix in this phase's diff,** not a single
   single-commit fix: (a) the H_z k-range off-by-one fix on the
   i/j-face stencil (commit d9f0a5a) is the necessary precursor
   — without it the cubic interpolation cannot help — and
   (b) the cubic Lagrange aux-grid interpolator (commit f9b9cef)
   is the dominant contributor, ~60 dB headroom by itself. The
   ADR documents both because re-deriving the priority order is
   the kind of audit cost we don't want to pay twice.
2. **4-point cubic Lagrange with a linear fallback at aux-grid
   boundaries.** Cubic needs two samples on each side of the
   query point; near the aux-grid ends fewer are available. The
   boundary pad cells already provide a buffer so the linear
   fallback path is rarely hit in practice; the fallback stays
   in place as cheap insurance.
3. **The 30°/45° `E ∥ ê_φ` gate test is tightened from `>10×`
   back to the original `>1000×` DoD,** and the ADR-0026
   in-test comment that flagged "Phase 2.fdtd.5.3 DoD of >1000×
   is gated on Phase 2.fdtd.5.3.2" is now obsolete and removed.
   The 1027× empirical contrast leaves ~3% headroom above the
   gate, which is the standard validation margin for this crate.
4. **The `no_match` baseline regression test is re-floored at
   `>50×`,** up from the pre-cubic 14.5× baseline. Cubic
   interpolation benefits both `dispersion_match=true` and
   `dispersion_match=false` paths; the old 14.5× contrast
   floor is no longer the right invariant to assert against.
   The 50× floor is conservatively below the empirical no-match
   contrast post-cubic so the gate measures the regression we
   actually want to catch.

## Consequences

The Phase 2.fdtd.5.3 DoD now clears at >1000× empirically, and
oblique TF/SF in `yee-fdtd` is production-quality. The
Phase 2.fdtd.5.x sub-ledger closes: 5.0 slab → 5.1 normal
finite-box → 5.2 j/k-face → 5.3 oblique + sign → 5.3.1
dispersion-matched aux step → 5.3.2 cubic interp. The
"contrast-limited at ~14×" caveat from ADR-0026's Consequences
section is now historical context, not a live limitation.

Forward-looking, the active FDTD horizons become **Phase
2.fdtd.6 (lumped RLC port Γ-against-analytic validation, deferred
from ADR-0017)** and **Phase 2.fdtd.7 (subgridding, spec/plan in
ADR-0027/ADR-0030)**. No further TF/SF-adjacent work is queued
on the roadmap.

## References

- ADR-0026 — Phase 2.fdtd.5.3 oblique TF/SF kernel + sign
  convention (this ADR's predecessor; the `>10×` placeholder
  gate it documented is now retired).
- Track KKKKK merge SHA `f878bdd`.
- Commit d9f0a5a — H_z k-range off-by-one fix on the i/j-face
  stencil (necessary precursor).
- Commit f9b9cef — `sample_aux_cubic()` 4-point cubic
  Lagrange aux-grid interpolator (dominant contributor).
- `crates/yee-fdtd/src/sources.rs` — `sample_inc_e` /
  `sample_inc_h` / `sample_aux_cubic` aux-grid sampling.
- `crates/yee-fdtd/tests/plane_wave_oblique.rs` — 30°/45°
  `E ∥ ê_φ` gate retightened to `>1000×`; `no_match`
  regression re-floored at `>50×`.
- A. Taflove and S. C. Hagness, *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method*, 3rd ed., Artech
  House, 2005, §5.10.5 (1-D aux source, dispersion considerations,
  interpolation order for the projection onto the 3-D grid).
- CLAUDE.md §4 — no-feature-without-a-gate; the original
  `>1000×` DoD is now the live gate.
