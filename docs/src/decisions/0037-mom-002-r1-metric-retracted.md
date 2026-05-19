# ADR-0037: MMMMMMM R1 ε_eff-biasing metric retracted

**Status:** Accepted (2026-05-19, Track NNNNNNN investigation).

## Context

Track MMMMMMM (commit `c7b001e`) ran a three-probe diagnostic on the
post-IIIIIII reframed mom-002 (L = 82 mm centered uniform on FR-4 at
1 GHz) and reported:

> **R1 ε_eff biasing — DETECTED.** Linear-fit of `arg(G_φ(ρ, 0, 0))`
> against ρ ∈ [40, 124] mm gives `k_eff/k_0 = 1.168`, so
> `ε_eff_solver = 1.364` vs Hammerstad-Jensen analytic `3.32`. The
> solver **under-estimates** dielectric loading by **−59 %**.

Track NNNNNNN investigated the proposed fix — adding a missing
`1/ε_r` divisor inside `MultilayerGreens::scalar_scalar` — and proved
this fix cannot move the metric. Empirical confirmation: dividing
`G_φ(ρ, 0, 0)` by `ε_r = 4.4` or by `ε_eff_HJ = 3.32` yields exactly
the same `k_eff/k_0 = 1.168` linear-fit slope. Any global real
divisor `c · G_φ` rescales the complex value uniformly, and the linear
fit on `arg(G_φ)` is mathematically invariant under such a transformation
(`arg(c · z) = arg(z)` for positive real `c`).

NNNNNNN's further per-term decomposition (free-space-only, image-sum-
only, surface-wave-only, TM₀ pole) confirmed no kernel component
carries a phase velocity near `k_0 · √3.32 ≈ 1.82 · k_0`:

| Kernel term      | `k_eff / k_0` | `ε_eff` |
|------------------|---------------|---------|
| Full kernel      | 1.168         | 1.364   |
| Free-space only  | 1.000         | 1.000   |
| Image sum only   | 0.001         | ≈ 0     |
| Surface-wave     | 1.037         | 1.076   |
| TM₀ pole         | 1.0003        | 1.0006  |

## Decision

MMMMMMM's R1 verdict is a **measurement artifact**, not a kernel bug.
`arg(G_φ(ρ; 0, 0))` measures the **point-source scalar-potential
phase decay** between two field/source points on the air side of the
slab. The Hammerstad-Jensen `ε_eff = 3.32` is a **strip-eigenmode**
property — the propagation constant of the integral-equation
solution `Z · I = V` for a finite microstrip on the slab — and is
not directly extractable from `G_φ(ρ; 0, 0)` at any radial sampling
of the kernel.

The R1 probe is retracted. The diagnostic file
`crates/yee-mom/tests/mom_002_13x_residual_diagnostic.rs` is kept for
the R2 (frequency sweep) and R3 (width sweep) probes which remain
valuable independent measurements, but its R1 verdict block should
not be used as evidence for an `ε_r` weighting bug.

## Consequences

* **No code change in this session.** No `1/ε_r` factor is added to
  `scalar_scalar`; the kernel matches the Aksun 1996 / Michalski-Mosig
  1997 published forms as audited by Track NNNNNNN.
* **The `|Im(Z)| = 674 Ω` residual at 1 GHz on the L = 82 mm reframed
  mom-002 remains unexplained.** Track NNNNNNN's recommendation for
  the next diagnostic:
  1. Extract the strip eigenmode directly from the assembled impedance
     matrix `Z`. The smallest-singular-value right eigenvector gives
     the dominant current distribution; its phase-vs-x slope along the
     strip's longitudinal axis is the propagation constant `β`.
  2. Compare `β / k_0` against `√ε_eff_HJ ≈ 1.82`. If they agree, the
     kernel is fine and the `|Im(Z)|` residual is a port-excitation or
     edge-singularity discretisation effect. If they disagree, the
     kernel's **strip eigenmode** physics is genuinely off — the
     `K^A` vector-potential image train (TE channel, inductive part of
     the line) is the most likely site.
* **Methodological note for future diagnostic tracks:** probes must
  measure quantities physically equivalent to the property under
  test. A point-source kernel-phase metric does not measure a
  strip-eigenmode property, even when both have units of "effective
  dielectric constant."

## References

* Track MMMMMMM diagnostic commit (`c7b001e`, 2026-05-19).
* Track NNNNNNN read-only audit + per-term decomposition (no commit,
  2026-05-19).
* Michalski, K. A. and Mosig, J. R., "Multilayered Media Green's
  Functions in Integral Equation Formulations," IEEE Trans. Antennas
  Propag., vol. 45, no. 3, pp. 508–519, Mar 1997, eqs. 6–12.
* Aksun, M. I., "A Robust Approach for the Derivation of Closed-Form
  Green's Functions," IEEE Trans. Microwave Theory Tech., vol. 44,
  no. 5, pp. 651–658, May 1996, eq. 9.
* Hammerstad, E. and Jensen, Ø., "Accurate Models for Microstrip
  Computer-Aided Design," MTT-S Digest, 1980.
