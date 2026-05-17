# ADR-0027: Phase 2.fdtd.7.0 subgridding walking-skeleton scope

## Status

Accepted — 2026-05-18 (spec only; implementation deferred to
follow-up tracks — see ADR-0030).

## Context

`YeeGrid` today is uniform: one `dx = dy = dz` set at construction
and applied across the entire volume. Fine for the resonant-cavity
and free-space-dipole gates we ship, but it fails on Phase 2
production targets: a 2.4 GHz patch on 0.508 mm RO4003C demands
`dx ≤ 50 µm` to resolve the substrate as ≥ 10 cells across,
forcing a `2000³ ≈ 8·10⁹` cell global domain. Slot antennas
(`λ/200` apertures) and conductor-corner singularities have the
same shape. Subgridding is the textbook fix (Taflove §13; Chevalier
1997; Berenger 2006) and complements conformal cell-cut techniques,
not replaces them. Track MMMMM (merge `d0d7af0`) lands the spec.

## Decision

Phase 2.fdtd.7.0 ships the minimum end-to-end nested pipeline:

1. **Single axis-aligned cuboidal nest, 2× isotropic refinement.**
   `dx_fine = dx_coarse / 2` in all three axes. Multiple nests,
   nest-in-nest, anisotropic and >2 ratios deferred (7.1, 7.4).
2. **2× time-subcycling, linear spatial + temporal interpolation**
   (Chevalier 1997 §III). Coarse → fine `E_t` linearly
   interpolated between bracketing coarse edges; blended in time
   at `t = n·dt_coarse + {¼,¾}·dt_coarse`. Fine → coarse `H_t`
   uses **area-averaging of the four fine cells covering each
   coarse face** (Chevalier 1997 §IV) — the closure keeping the
   discrete energy balance second-order in `dx_coarse`.
3. **EM-only inside the fine region.** Non-dispersive isotropic
   scalar `ε_r` / `μ_r`; dispersive ADE in-nest deferred to 7.2.
4. **No CPML or TF/SF face inside the fine region.** Co-location
   is a documented **runtime error** in 7.0; lands in 7.3.

Validation gate **fdtd-007 — Maloney-Smith 1993 dielectric-loaded
thin slot** (`w = 0.5 mm`, `L = 30 mm`, `ε_r = 2.2`, `h = 1.524
mm`): resonance ±2%, `|S_11|` ±1 dB, plus 0.3% / 0.3 dB sanity
check vs uniform-fine reference. Stability companion: 10 000-step
round-trip energy-drift gate (≤ 0.5%) — the canary for the
classical asymmetric-coupling failure mode (Berenger 2003 §IV).

CPU-only, single-threaded, FP64. Lane: `crates/yee-fdtd/**`. No
CLI / Python / GUI exposure in 7.0.

## Consequences

- **Enables sub-cell features without globally refining.** Thin
  substrates, slot antennas, fillets become tractable inside the
  nest while the coarse grid carries the far field. Phase
  2.fdtd.7.5 (RO4003C patch vs openEMS) becomes reachable.
- **Defers the hard composability questions** — ADE-in-nest
  (7.2), CPML / TF-SF co-location (7.3), higher-order spatial
  interpolation (7.4) — each as a one-paper follow-up, not an
  FDTD rewrite.
- **Late-time interface instability is the load-bearing risk.**
  The 10 000-step gate is the canary; documented fallback is
  Berenger 2006's Huygens-surface variant.
- **No CLI / Python / GUI in 7.0.** Direct Rust API only;
  consumer wiring is a 7.0.1 follow-up.

## References

- `docs/superpowers/specs/2026-05-18-phase-2-fdtd-7-subgridding-design.md`
- Track MMMMM merge commit `d0d7af0`.
- ADR-0030 — Phase 2.fdtd.7 implementation plan (companion).
- M. W. Chevalier et al., *IEEE Trans. Antennas Propag.* 45(3),
  1997 — interpolation + area-average closure.
- J.-P. Berenger, *IEEE Trans. Antennas Propag.* 54(12), 2006 —
  stability analysis and documented fallback.
- J. G. Maloney, G. S. Smith, *IEEE Trans. Antennas Propag.*
  41(5), 1993, Fig. 9 — fdtd-007 reference.
- CLAUDE.md §3, §4.
