# Multilayer Green's Function — DCIM and Surface-Wave Poles

This page is the theory-of-operation reference for the multilayer
mixed-potential Green's function in `yee-mom` — the structural
analog of the [planar MoM](./planar-mom.md) and
[RWG / MPIE](./mom-rwg-mpie.md) chapters for the case where the
free-space `exp(-jk₀R) / (4πR)` is no longer the right kernel
because a finite-thickness dielectric slab over a PEC ground plane
breaks translation invariance in `z`. Audience: an engineer reading
`crates/yee-mom/src/multilayer.rs` with a textbook open. Equations
are plain-text ASCII — mdBook has no math preprocessor (see
`docs/book.toml`).

## 1. Introduction

The planar-MoM kernel `Z_{mn}` of [planar-mom.md §4](./planar-mom.md)
carries a Green's function `G(r, r')` inside two surface integrals.
In free space the kernel is `exp(-jk₀R) / (4πR)`. On a **grounded
dielectric slab** — PEC ground at `z = -h`, dielectric `ε_r` of
thickness `h` in `-h < z < 0`, free space above hosting the metal
traces — `G` carries contributions from infinitely many image
reflections off the ground plane *and* from a discrete set of
**surface-wave bound modes** propagating radially along the
substrate. Closed-form evaluation is no longer possible; the kernel
is most naturally written as a Sommerfeld integral in the spectral
domain, and the central problem of multilayer MoM is reducing that
integral to something tractable.

Phase 1.1.0 shipped a one-image DCIM placeholder; Phase 1.1.1.0
(Track OOOO) replaced it with N-image DCIM via GPOF; Phase 1.1.1.2
(DDDDD spec; JJJJJ in flight) adds analytic surface-wave pole
extraction to recover the radial Hankel tail pure DCIM cannot
represent. Validation target: **mom-002** — 50 Ω microstrip on
FR-4 (`ε_r = 4.4`, `h = 1.6 mm`) at 1 GHz, `|Z_in|` in the
Hammerstad-Jensen `[35, 75] Ω` corridor (Pozar §3.8). Gate stays
loose until Phase 1.1.1.2 lands — see §11, `CLAUDE.md` §4.

## 2. Spectral-domain Green's function

The grounded slab is translationally symmetric in `(x, y)` with
layering only in `z`. Fourier transform in the transverse plane
decouples each transverse wavenumber `k_ρ = sqrt(k_x² + k_y²)`.
The vector and scalar potentials take the **Sommerfeld-integral**
form

```text
G(ρ, z, z')  =  (1 / (4π)) · ∫_C  G̃(k_ρ; z, z')  H_0^{(2)}(k_ρ ρ)  k_ρ dk_ρ,
```

with `C` the standard real-axis contour with radiation-condition
indentations (Mosig 1989; Michalski-Mosig 1997). The spectral
kernel `G̃` is a closed-form transverse-line equivalent: for the
TE / TM split (Pozar §3.7, eqs. 3.196-3.199) the channel
denominators have the form

```text
D_TE(k_ρ)  =  k_z0  +  k_zd · cot(k_zd · h),
D_TM(k_ρ)  =  ε_r · k_z0  +  j · k_zd · tan(k_zd · h),
k_z0       =  sqrt( k_0²       -  k_ρ² ),    Im(k_z0) ≤ 0  (outgoing),
k_zd       =  sqrt( ε_r · k_0² -  k_ρ² ),    principal branch.
```

(Exact tan / cot signs are source-position dependent; DDDDD plan
Step 0 re-derives them from PEC short at `z = -h` and radiation
at `z → +∞`.) Two obstacles make the integral hard:

- **Highly oscillatory integrand.** For `k_ρ ρ → ∞`,
  `H_0^{(2)}(k_ρ ρ)` is a slowly decaying complex exponential;
  conditional convergence only — the **Sommerfeld tail**.
- **Singularities near the real axis.** Branch points at
  `k_ρ = ±k_0` and a finite set of **surface-wave poles** — zeros
  of `D_TE`, `D_TM` on the proper Riemann sheet — sit between `k_0`
  and `sqrt(ε_r) k_0`. TM₀ is always present for `ε_r > 1`; higher
  poles appear for thick or high-`ε_r` substrates.

Brute-force real-axis quadrature is dead on arrival. Production
solutions: (a) **DCIM**, deforming the contour off the real axis
and fitting the smooth remainder by complex exponentials that map
to free-space-like image sources, and (b) **analytic pole
extraction**, subtracting surface-wave residues in closed form
before the fit.

## 3. The DCIM image-method approach

Aksun 1996 observes that on the **deformed contour**
`k_z0(t) = k_0 (1 - j t)` for `t ∈ [0, T_max]`, `G̃` reparameterised
in `k_z0` is smooth, no longer oscillatory, and well-approximated
by a short sum of complex exponentials:

```text
G̃(k_ρ)      ≈  Σ_{n=1}^N  α_n · exp( -j k_z0(k_ρ) · a_n ),
G(ρ, z, z') ≈  Σ_n  b_n · exp(-j k_0 R_n) / (4π R_n),
R_n²        =  ρ² + (z + z' - a_n)².
```

Each spectral exponential maps to a **complex image source** via
the **Sommerfeld identity** `exp(-j k_z0 a) / k_z0 ↔ exp(-j k_0 R)
/ R` with `R = sqrt(ρ² + a²)`, `a` a complex offset. `(b_n, a_n)`
are recovered from the `(α_n, β_n)` fit (`a_n = -β_n / k_0`,
`b_n = α_n · exp(j k_0 a_n)`; see `multilayer.rs` §1). Vector and
scalar potentials carry independent TE / TM fits — the Phase 1.1.0
channel collapse is resolved in Phase 1.1.1.0 (ADR-0020). Every
`Z_{mn}` integral then becomes a sum of free-space-style integrals
handed off to the Khayat-Wilton singular machinery of mom-001
(`planar-mom.md` §6).

## 4. GPOF — generalised pencil of function

The fitting step: given uniformly-sampled
`y_m = G̃(k_ρ(m · Δt))` of `y(t) = Σ_n α_n · exp(β_n · t)`,
recover the `N` pairs `(α_n, β_n)`. This is the **Hua-Sarkar
matrix-pencil** problem (Hua & Sarkar 1989). The closed-form
pipeline:

1. **Sample.** `M` uniform points on the Aksun contour
   (`M ∈ [30, 100]`, `M ≥ 2N` floor).
2. **Hankel matrices.** `Y_1[i,j] = y_{i+j}`,
   `Y_2[i,j] = y_{i+j+1}` of shape `(M - L) × L`, `L = M - N`.
3. **SVD-truncate.** `Y_1 = U Σ V^H`; keep top `N` components.
4. **Reduced-pencil eigenproblem.**
   `Z = Σ_N^{-1} · U_N^H · Y_2 · V_N`. The eigenvalues `z_n`
   satisfy `z_n = exp(β_n · Δt)`, recovering the exponents.
5. **Vandermonde least-squares.** Solve `V α = y` for amplitudes
   (`V_{m,n} = z_n^m`) via SVD-backed `solve_lstsq`
   (`crates/yee-mom/src/gpof.rs`).

Pure linear algebra, no iteration, no convergence tolerance. The
**noise-floor vs truncation-order trade-off**: more `N` resolves
more complex spectra at the cost of pushing `Σ_N`'s smallest
singular value toward the floating-point floor, where overfitting
injects spurious near-cancelling image pairs. The sweet spot for
grounded-slab DCIM is `N ∈ [5, 10]`; `n_images = 5` is the mom-002
production choice (ADR-0020). GPOF is numerically equivalent to
Prony in exact arithmetic but materially more stable in floating
point — SVD truncation is a built-in regulariser.

## 5. Why pure DCIM fails near the surface-wave pole

The exponential basis `exp(-j k_z0 a)` cannot fit a function with a
**pole**. If `G̃ ~ R_p / (k_ρ - k_p)` near a pole, no finite sum of
exponentials in `k_z0` reproduces that singularity. GPOF converges
in coefficient norm but never resolves the residue. Space-domain
consequence: the image sum decays as `1/R` (free-space-like), but
the true Green's function at a surface-wave pole decays only as
`ρ^{-1/2}` — the Hankel tail.

Track AAAAA's mom-002 signature: `|Z_in|` flooring at ~2 kΩ and
`Im(Z_in) ≈ -2.1 kΩ` across mesh refinements `nz = 8 → 32`.
Refinement does nothing because **the residual is not mesh-bound**
— it is a missing analytic term. At FR-4 / 1.6 mm / 1 GHz the TM₀
pole sits at `k_ρ / k_0 ≈ 1.6`, between `k_0` and
`sqrt(ε_r) k_0 ≈ 2.1 k_0`. Microstrip geometry (`L = 30 mm ≫ h`)
amplifies the contribution because `ρ ~ L` is firmly in the
Hankel-tail regime. This is the mom-002 failure mode that drove
Phase 1.1.1.2.

## 6. Surface-wave pole extraction

Surface-wave poles are zeros of `D_TE(k_ρ)` and `D_TM(k_ρ)` on the
proper Riemann sheet, found by **Newton-Raphson** in the complex
`k_ρ` plane:

```text
k_ρ^{(n+1)}  =  k_ρ^{(n)}  -  D(k_ρ^{(n)}) / D'(k_ρ^{(n)}),
k_{ρ,0}      =  k_0 · sqrt( (ε_r + 1) / 2 )    (quasi-static guess).
```

Tolerance `|D(k_p)| < 1e-12`, max 50 iterations (typical 5-10). The
quasi-static guess places the TM₀ seed at the effective-permittivity
wavenumber — the quasi-TEM-microstrip approximation. For FR-4 at
1 GHz this gives `k_{ρ,0} / k_0 ≈ 1.64`, well inside Newton's basin
for the true pole at `k_p / k_0 ≈ 1.6`. The Jacobian is closed-form
via implicit differentiation: `∂k_z0/∂k_ρ = -k_ρ/k_z0`,
`∂k_zd/∂k_ρ = -k_ρ/k_zd`, chained through `d/dx[cot x] = -(1 + cot² x)`.

At each converged pole `k_p`, the **residue** of `G̃` is
`R_p = N(k_p) / D'(k_p)`, where `N` is the TE / TM channel
numerator (Pozar §3.7; Michalski-Mosig 1997 eqs. 16-19). Lossless
slab → purely real residue; lossy substrates → complex. The
**degenerate-pole** branch triggers when `|D'(k_p)| < 1e-10` (near
an adjacent mode's cutoff); the DDDDD escape hatch falls back to
finite-differenced `D'`.

## 7. Pole-subtracted GPOF + Hankel reconstruction

With the pole list `{(k_{p,j}, R_{p,j})}` in hand, the spectral
Green's function splits as

```text
G̃(k_ρ)  =  G̃_pole(k_ρ)  +  G̃_residual(k_ρ),
G̃_pole      =  Σ_j  R_{p,j}  /  (k_ρ  -  k_{p,j}).
```

`G̃_residual` is **smooth everywhere on the Aksun contour** — the
singularities that defeated pure GPOF are analytically removed.
The OOOO GPOF machinery then fits `G̃_residual` cleanly with
`N ≈ 5` images. TE and TM are pole-subtracted **independently**
using their own zero sets. The space-domain reconstruction is

```text
G(ρ, z, z')  ≈   Σ_n  b_n · exp(-j k_0 R_n) / (4π R_n)
              +  Σ_j  (-j / 4) · R_{p,j} · H_0^{(2)}(k_{p,j} · ρ)
                       · ψ_{p,j}(z) · ψ_{p,j}(z').
```

The first sum is the image train of §3. The second is the
**analytic surface-wave contribution**, with `H_0^{(2)}` the
zeroth-order Hankel function of the second kind; its large-arg
asymptotic `H_0^{(2)}(z) ≈ sqrt(2/(π z)) · exp(-j(z - π/4))` for
`|z| ≥ 8` gives the canonical `ρ^{-1/2}` decay no exponential
basis can reproduce. `ψ_{p,j}(z)` is the modal `z`-profile:
sinusoid inside the slab tied to PEC ground
(`ψ ∝ cos(k_zd (z + h)) / cos(k_zd h)`), exponential decay above
(`ψ ∝ exp(-α_0 z)`, `α_0 = sqrt(k_p² - k_0²)`). For mom-002 strips
on the slab top, `z = z' = 0` collapses the modal product to unity.
The Phase 1.1.1.2 kernel is a finite image sum plus a finite
Hankel-pole sum, slotted into the RWG matrix fill of
`planar-mom.md` §4-§6. Construction cost: one Newton solve per
channel per frequency plus one GPOF — milliseconds, negligible vs.
the per-frequency `O(N²)` matrix fill.

## 8. Higher-order poles, lossy substrates, leaky modes

Three extensions matter; only the first two are realistic for
Phase 1.1.1.x.

**Higher-order poles.** Thick or high-`ε_r` substrates support
multiple TM_n / TE_n modes above their cutoffs. One Newton solve
per pole, seeded from the next mode's cutoff (`k_0 sqrt(ε_r)` for
TE₁), accepting a converged pole only if it differs from earlier
poles by more than `0.01 k_0`. For FR-4 / 1.6 mm only TM₀ is bound
up to ~10 GHz; the second `n_surface_wave_poles = 2` slot
gracefully no-ops when Newton finds no second mode.

**Lossy substrates.** `ε_r = ε_r' - j ε_r''` moves the dispersion
equation into the complex domain; poles slide off the real axis
into the lower half-plane. Newton remains well-posed with a complex
initial guess; `R_p` becomes complex, contributing an extra
radial-decay envelope through `Im(k_p) ρ` in the Hankel asymptotic.

**Leaky modes.** Improper Riemann-sheet poles, radiating energy out
of the substrate; mathematically genuine but require the contour
to cross a branch cut. Out of scope for Phase 1.1.1.x — needs the
Phase 1.1.1.3 full-Sommerfeld contour-deformation machinery.

## 9. Validation: mom-002 microstrip Z₀

The canonical Phase 1.1.x case is **mom-002**: a 50 Ω microstrip on
FR-4 (`ε_r = 4.4`, `h = 1.6 mm`, `W ≈ h`) at 1 GHz, delta-gap fed
(`planar-mom.md` §7) with a matched load via sufficient line
length. `|Z_in|` is expected in the Hammerstad-Jensen `[35, 75] Ω`
corridor (Pozar §3.8 — `Z_0 ≈ 50 Ω ± 50%` covers manufacturing
tolerance on `W/h ≈ 1`). Four reference points trace convergence:

- **Phase 1.1.0 (one-image placeholder).** `|Z_in| ≈ 14 kΩ`; loose
  `[1, 100 kΩ]` non-degeneracy gate only.
- **Phase 1.1.1.0 (N-image DCIM, OOOO).** `|Z_in| ≈ 2.7 kΩ`. ~5×
  better, still outside Hammerstad-Jensen. ADR-0020 holds the
  loose gate because the residual is dominated by the missing
  pole, not GPOF resolution.
- **Phase 1.1.1.1 (edge-clustered mesh, AAAAA).** Refinement
  `nz = 8 → 32` does not move `|Z_in|` — residual is not
  mesh-bound, the diagnostic that drove DDDDD.
- **Phase 1.1.1.2 (pole-subtracted GPOF + Hankel, JJJJJ — in
  flight).** Hammerstad-Jensen gate tightened to hard pass.

mom-001 (NEC-4 87 + j41 Ω) is unaffected: the surface-wave term is
gated on the `MicrostripSommerfeld` `GreensSpec` variant;
`FreeSpace` / `Microstrip*` execute the OOOO codepath verbatim and
mom-001's NEC-4 gate stays green (DDDDD DoD item 3; `CLAUDE.md` §4).

## 10. Numerical considerations

**Contour choice.** Aksun 1996's single-segment
`k_z0(t) = k_0 (1 - j t)` suffices because the pole is removed
analytically before sampling. The T-shaped two-segment contour
(Aksun-1996 §III) is a refinement for kernels passing close to
poles or branch points. `T_max ≈ 10` and `M ∈ [30, 100]` samples
are the production envelope; outside this the GPOF condition
number degrades.

**GPOF truncation order.** `n_images = 5` is production for mom-002
(ADR-0020). Below `N = 3` the fit cannot resolve `G̃_residual`'s
structure; above `N = 10` the smallest singular value of `Σ_N`
drops below `1e-12` and spurious near-cancelling image pairs
appear. `n_images` is exposed on the constructor.

**Newton convergence radius.** The quasi-static guess lies inside
Newton's basin for TM₀ on every FR-4 / RO4003C / Rogers-class
substrate in the design corridor
(`ε_r ∈ [2, 12], h ∈ [0.1, 3] mm`). The DDDDD escape hatch falls
back to **Müller's method** (three-point parabolic extrapolation;
quadratic convergence near a simple pole) if Newton overshoots.

**Hankel-function evaluation.** The mom-002 radial range gives
`|k_p ρ| ∈ [3.3e-3, 3.3]` at 1 GHz — straddling the small-arg /
large-arg switchover. The implementation uses the small-argument
series for `|z| < 8` (with `ln(z/2) + γ` handling) and the
asymptotic expansion for `|z| ≥ 8`; cross-checking at the seam to
`1e-10` relative is the sanity test (DDDDD plan Step 3).

## 11. Limitations and roadmap

The scope is intentionally narrow; each restriction maps to a
deferred sub-phase.

- **Single pole per channel (TM₀ + optional TE₁).** Phase 1.1.1.2
  ships `n_surface_wave_poles ≤ 2`; multi-pole search with
  separation guards lands in Phase 1.1.1.3.
- **Phase 1.1.1.2 in flight.** DDDDD spec + ADR-0025 are on `main`;
  the JJJJJ track is mid-flight at this chapter's base SHA. mom-002
  stays at the OOOO loose gate until JJJJJ merges (`CLAUDE.md`
  §4 / §10).
- **Single dielectric layer.** Multi-layer stacks are Phase 1.1.2;
  the transverse-line-equivalent kernel generalises cleanly
  (Michalski-Mosig 1997 §II) — GPOF and pole-search reuse.
- **Isotropic dielectric.** Anisotropic substrates (ferrite,
  LCP-laminate, engineered-tensor oxide stacks) are Phase 1.1.3+;
  the kernel loses TE / TM decoupling, pole search runs in a
  coupled 4×4 transfer-matrix framework.
- **Full Sommerfeld integration with pole-locus contour
  deformation.** Phase 1.1.1.3 — strict-tolerance alternative.

## 12. References

- Sommerfeld, A. "Über die Ausbreitung der Wellen in der
  drahtlosen Telegraphie." *Ann. Phys.* 28 (1909), pp. 665–736. —
  Original Sommerfeld integral; surface-wave decomposition.
- Hua, Y., and Sarkar, T. K. "Generalized pencil-of-function
  method." *IEEE Trans. Antennas Propag.* 37.2 (Feb 1989),
  pp. 229–234. — Foundational matrix-pencil method (§4).
- Aksun, M. I. "A robust approach for the derivation of closed-form
  Green's functions." *IEEE Trans. Microw. Theory Tech.* 44.5
  (May 1996), pp. 651–658. — Canonical DCIM reference; introduces
  the deformed contour and explicitly calls for pole extraction
  **before** GPOF — the step Phase 1.1.1.2 implements.
- Michalski, K. A., and Mosig, J. R. "Multilayered media Green's
  functions in integral equation formulations." *IEEE Trans.
  Antennas Propag.* 45.3 (Mar 1997), pp. 508–519. — Modern review;
  eqs. 16-19 carry the residue formula §7 follows.
- Mosig, J. R. "Integral equation technique for planar
  geometries." In *Numerical Techniques for Microwave and
  Millimeter-Wave Passive Structures*, T. Itoh (ed.), Wiley, 1989,
  Ch. 3. — MPIE on a grounded slab.
- Yang, J. J., and Aksun, M. I. *IEEE Trans. Antennas Propag.*
  39.7 (Jul 1991), pp. 1042–1046. — Pre-DCIM spectral-domain
  formulation referenced in §3.
- Felsen, L. B., and Marcuvitz, N. *Radiation and Scattering of
  Waves.* IEEE Press, 1994. Ch. 5 — Hankel-asymptotic derivation
  of the surface-wave residue term in §7.
- Pozar, D. M. *Microwave Engineering.* 4th ed. Wiley, 2011. §3.7
  (grounded-slab dispersion, §2), §3.8 (Hammerstad-Jensen `Z_0`,
  mom-002 reference in §9).
- Chow, Y. L., et al. "A closed-form spatial Green's function for
  the thick microstrip substrate." *IEEE Trans. Microw. Theory
  Tech.* 39.3 (Mar 1991), pp. 588–592. — Pre-Aksun closed-form
  image method.
- ADR-0020 (multi-image DCIM via GPOF) and ADR-0025 (Phase
  1.1.1.2 Sommerfeld pole extraction spec) — locked-in decisions
  for §4 and §6-§7.
