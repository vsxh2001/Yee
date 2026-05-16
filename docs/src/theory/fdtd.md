# Finite-Difference Time-Domain — Theory of Operation

This page is the theory-of-operation reference for Yee's
Finite-Difference Time-Domain (FDTD) solver, implemented in the
`yee-fdtd` crate. Same audience as the planar-MoM page (an engineer
reading source code with a textbook open), same conventions
(plain-text math, inline citations, source-file references in inline
code).

## 1. Overview

FDTD discretises Maxwell's curl equations directly on a 3D lattice and
marches them forward in time with a second-order leapfrog scheme. No
matrix to invert, no integral kernel: every cell update is a handful
of floating-point operations on its six nearest field neighbours. The
price is volumetric meshing — every cubic wavelength consumes memory
and compute — and an explicit Courant time-step constraint.

K. S. Yee introduced the staggered-grid updates in 1966 ("Numerical
Solution of Initial Boundary Value Problems Involving Maxwell's
Equations in Isotropic Media", *IEEE Trans. Antennas Propag.* 14.3,
pp. 302–307). Through the 1980s Taflove and collaborators extended
the scheme to anisotropic and dispersive media, absorbing boundaries,
and the TF/SF source formulation, culminating in Taflove and Hagness,
*Computational Electrodynamics: The Finite-Difference Time-Domain
Method*, 3rd ed., Artech House (2005) — still the standard reference.
CPML (Roden & Gedney 2000) and dispersive ADE / PLRC (Luebbers &
Hunsberger 1992) closed the remaining production gaps. Today FDTD
underlies most commercial volumetric transient EM solvers — CST
T-solver, Remcom XFdtd, Lumerical FDTD, openEMS.

FDTD and planar MoM are complementary. MoM excels on conductive
surfaces in stratified media (PCBs, microstrip filters, wire
antennas). FDTD excels in volume, on broadband transients (one run
yields the entire pulse spectrum after one FFT), and on dispersive
or open-region radiation problems.

## 2. Maxwell's curl equations on a staggered grid

The starting point is the time-domain curl pair (SI units, no
impressed sources):

```text
∂B/∂t  =  -∇ × E,          B = μ₀ μ_r H
∂D/∂t  =   ∇ × H - J,      D = ε₀ ε_r E
```

In a linear isotropic non-dispersive medium this expands into six
coupled scalar equations. Yee's insight: they naturally live on a
lattice staggered in both space and time. Each `E` component sits on
the *edges* of a primary cubic cell aligned with its own axis; each
`H` component sits on the *faces* of that cell, normal to its own
axis. With cell origin at `(i, j, k) · (dx, dy, dz)`:

| field | grid offset                | array shape (Yee impl)  |
|-------|----------------------------|-------------------------|
| E_x   | `(i+1/2, j,    k   )`      | `[nx,   ny+1, nz+1]`    |
| E_y   | `(i,    j+1/2, k   )`      | `[nx+1, ny,   nz+1]`    |
| E_z   | `(i,    j,    k+1/2)`      | `[nx+1, ny+1, nz  ]`    |
| H_x   | `(i,    j+1/2, k+1/2)`     | `[nx+1, ny,   nz  ]`    |
| H_y   | `(i+1/2, j,    k+1/2)`     | `[nx,   ny+1, nz  ]`    |
| H_z   | `(i+1/2, j+1/2, k   )`     | `[nx,   ny,   nz+1]`    |

ASCII sketch of a single cell, looking down `+x` (so `y` horizontal, `z`
vertical):

```text
              E_z
              ↑
        +-----|-----+
        |           |
   E_y ←|   • H_x   |→ E_y
        |           |
        +-----|-----+
              ↑
              E_z
```

The face carries one `H` component normal to it; the edges carry `E`
components along them. Cyclically permute `(x, y, z)` for the other
two face-normal `H` components. In time, `E` is sampled at integer
levels `n · dt`, `H` at half-integer levels `(n + 1/2) · dt`. The
implementation lives in `crates/yee-fdtd/src/grid.rs` as six
`Array3<f64>` arrays of the shapes above.

Every successful volumetric finite-difference scheme for Maxwell uses
some variant of Yee's staggering, for geometric reasons. The curl
operator maps `E`-edges to `H`-faces, so each discrete curl is a
perfect local circulation: four edge values around one face. The
identities `∇·B = 0` and `∇·D = ρ` are preserved exactly in vacuum
(round-off only); the scheme is second-order accurate in space and
time on uniform grids. The staggered grid is a discrete realisation
of the de Rham complex on a cubic lattice.

## 3. Leapfrog time stepping

Replace each space derivative with a centred difference across the
half-cell gap and each time derivative with a centred difference
across the half-step gap. The textbook Yee update for `H_x` (to fix
the indexing convention):

```text
H_x(i,j,k) [n+1/2] = H_x(i,j,k) [n-1/2]
                    + (dt / (μ₀ μ_r)) · ( (E_y(i,j,k+1) - E_y(i,j,k)) / dz
                                        - (E_z(i,j+1,k) - E_z(i,j,k)) / dy )
```

The five symmetric counterparts follow by cyclic permutation of
`(x, y, z)` and the sign change between curl pairs;
`crates/yee-fdtd/src/update.rs` lists all six in its module doc and
walks them with triply-nested loops in `update_h` and `update_e`. The
half-step pair is the only temporal coupling: `H` at `n+1/2` depends
on `E` at `n` and `H` at `n-1/2`; `E` at `n+1` depends on `H` at
`n+1/2` and `E` at `n`. That decoupling makes FDTD trivially explicit
and trivially parallel on a GPU.

The per-step order in `WalkingSkeletonSolver::step` is
`update_h → apply_pec → update_e → apply_pec`; the clamp is applied
twice so both the half-step `H` and full-step `E` read consistent
outer-face values. Taflove and Hagness §3 (3rd ed., 2005, eqns.
3.27–3.32) derive these updates in this exact form.

## 4. Courant–Friedrichs–Lewy stability

Leapfrog time stepping is conditionally stable. A plane-wave ansatz
`exp(jk·r - jωt)` in the update equations, with the amplification
factor required to sit on the unit circle for all spatial Fourier
modes, yields the 3D Courant condition

```text
dt  ≤  1 / ( c · sqrt(1/dx² + 1/dy² + 1/dz²) )
```

(Taflove & Hagness §4.7, eqn. 4.60). Physically: the time step must
be shorter than a light-speed wave needs to cross the diagonal of the
smallest cell. Violating the limit produces exponential blow-up
within a few dozen steps.

The implementation exposes the limit as `YeeGrid::courant_limit` in
`crates/yee-fdtd/src/grid.rs` and pins `dt = 0.9 · courant_limit()`
in the `vacuum` constructor. The 0.9 safety factor (~10% headroom) is
standard FDTD practice: it absorbs stability erosion from non-vacuum
materials, round-off, and the stencil modifications introduced by
CPML and dispersive-material updates. Running closer to the limit
buys ~10% wall-clock and loses robustness margin; running well below
wastes time without improving accuracy (numerical phase error is
minimised at the Courant limit, not at half of it).

## 5. Sources

A **soft point source** adds an analytic time function to a field
component each step (`E_z(i,j,k) += g(t)`). The cell continues to
respond normally, so the source point does not become a reflector.
`gaussian_pulse_ez` in `crates/yee-fdtd/src/sources.rs` injects
`exp(-((t - t0) / sigma)²)` into the chosen cell. The Gaussian's
closed-form Fourier transform (another Gaussian, width `1/sigma`)
makes it a natural broadband stimulus for transfer-function
extraction.

A **hard source** overwrites the field (`E_z(i,j,k) = g(t)`). Simple,
but it reflects every wave returning to the source cell; soft sources
are preferred unless a Dirichlet boundary is wanted.

**Plane-wave injection** uses a surface of soft-source cells emitting
a coherent uniform-amplitude wave. Naïvely this pollutes the entire
domain. The clean fix is the **total-field / scattered-field (TF/SF)**
decomposition (Taflove & Hagness §5): a closed surface separates a
total-field interior (incident plus scattered) from a scattered-field
exterior, with injection only on that surface and consistency-
correcting updates suppressing leakage.

Phase 2 walking skeleton ships only the Gaussian-on-`E_z` soft
source. Hard sources, plane-wave injection, TF/SF, modal port
sources, and lumped-element sources are Phase 2.1+ work.

## 6. Absorbing boundaries (CPML)

Untreated walls reflect every outgoing wave. For closed-cavity
problems that *is* the physics and a PEC clamp on tangential `E` is
correct. For radiation, scattering, and antenna problems the walls are
fictitious — the real boundary is at infinity — and the only
acceptable behaviour is broadband absorption with negligible
reflection.

The modern solution is the **convolutional perfectly matched layer
(CPML)** of Roden and Gedney, "Convolution PML (CPML): An efficient
FDTD implementation of the CFS-PML for arbitrary media", *Microwave
Opt. Technol. Lett.* 27.5 (2000), pp. 334–339. CPML descends from
Berenger's 1994 split-field PML through a complex-frequency-shifted
(CFS) stretched-coordinate reformulation, recasting frequency-domain
stretching as a time-domain auxiliary state on a thin shell of cells.

In sketch: replace coordinate `x` in the curl equations with a
stretched coordinate `x̃` whose differential is `s_x(ω) dx`, where

```text
s_x(ω)  =  κ_x  +  σ_x / (α_x + jωε₀)
```

ramps polynomially from `(κ_x, σ_x, α_x) = (1, 0, free)` at the inner
PML face to large values at the outer face. The CFS term `α_x`
absorbs the low-frequency and evanescent components Berenger's
original handled poorly. Discretised on the Yee grid the stretching
becomes a recursive convolution: each update inside the PML adds an
auxiliary state vector (`ψ_E`, `ψ_H`) integrating field history with
an exponential kernel determined by `σ_x` and `α_x`. Polynomial
grading is `(σ, κ)(d) = (σ_max, κ_max - 1) · (d / L)^m + (0, 1)` with
`m = 3` or `4` and `L = 8–12` cells.

Phase 2.0 uses a hard PEC clamp on all six outer faces via
`crates/yee-fdtd/src/boundary.rs::apply_pec`. Phase 2.1 introduces
CPML in `crates/yee-fdtd/src/cpml.rs`.

## 7. Dispersive materials (forward-looking)

In real media — silicon, gold at optical frequencies, biological
tissue, loss-tangent dielectrics — permittivity depends on frequency
and a constant `ε_r` is the wrong model. FDTD accommodates dispersion
by promoting `D = ε E` to a differential relation integrated
alongside Ampère's law each step.

Two equivalent recipes dominate. **Auxiliary Differential Equation
(ADE)** (Luebbers and Hunsberger, "FDTD for Nth-order dispersive
media", *IEEE Trans. Antennas Propag.* 40.11 (1992), pp. 1297–1301)
writes each pole as an ODE in a polarisation auxiliary variable and
time-steps it on the same leapfrog stencil as the fields. **Piecewise
Linear Recursive Convolution (PLRC)** expresses
`D(t) = ε₀ ε_∞ E(t) + ∫ χ(τ) E(t-τ) dτ`, exploiting the exponential
kernel per pole so the integral collapses to one multiply-add per
pole per cell. ADE and PLRC are mathematically equivalent for the
standard pole models.

Standard pole models: **Debye**
`ε(ω) = ε_∞ + (ε_s - ε_∞) / (1 + jωτ)` (one real pole; water, tissues);
**Lorentz** (one complex-conjugate pair; narrow-band resonances);
**Drude** `ε(ω) = ε_∞ - ω_p² / (ω² + jωγ)` (one pole at zero; metals
below their plasma frequency); and **multi-pole** sums for broadband
fits.

Phase 2 walking skeleton is vacuum-only; dispersive materials enter
in Phase 2.2 alongside CPML.

## 8. Near-to-far-field (NTFF) transformation (forward-looking)

Antenna patterns live in the far field but FDTD only stores the near
field. The standard bridge is the NTFF transformation: enclose the
radiator in a virtual integration surface `S` inside the PML, record
tangential `E` and `H` on `S` throughout the run, and post-process
into surface equivalent currents

```text
J_s  =  n̂ × H,        M_s  =  -n̂ × E         (on S)
```

By the surface equivalence principle these currents radiate the same
far-field as the original sources. The far-field integral

```text
E_far(r̂)  ∝   r̂ × ( r̂ × ∫_S J_s(r') exp(jk · r') dS'
                   - η ∫_S M_s(r') exp(jk · r') dS' )
```

(with `k = ω/c · r̂`, `η = sqrt(μ₀/ε₀)`) can be evaluated
frequency-by-frequency or directly in time. Taflove and Hagness §8
give the full derivation and the discrete surface weights consistent
with Yee staggering.

NTFF is Phase 2.3, on top of CPML and stable broadband sources. The
implementation will record tangential E/H on a configurable inner box
and emit `(theta, phi) → gain, directivity, axial ratio` pattern
files for round-trip validation against planar MoM.

## 9. Material handling

`YeeGrid` carries `eps_r` and `mu_r` as scalar `f64` constants applied
uniformly (Phase 2.0; `crates/yee-fdtd/src/grid.rs`). Update kernels
read them through `coeff = dt / (EPS0 * grid.eps_r)`, so per-cell
heterogeneity is a structural promotion to `Array3<f64>`: each cell
carries its own `eps_r`, `mu_r`, conductivity `σ_e`, and magnetic
loss `σ_m`. The update coefficient then takes the standard two-term
form `C_E = (1 - σ_e dt / (2ε)) / (1 + σ_e dt / (2ε))` and
`D_E = (dt / ε) / (1 + σ_e dt / (2ε))`.

Phase 2.1 promotes material storage to per-cell arrays and adds the
lossy-dielectric two-term update. Arbitrary curved geometry cutting
across the cubic-cell axes is **conformal FDTD** (Phase 2.3+):
**Dey-Mittra** (partially-filled-cell modification of curl integration
weights; *S. Dey and R. Mittra*, *IEEE MGWL* 7.9, 1997) or the
simpler staircased fallback. Conformal techniques typically recover
one to two orders of magnitude in geometric accuracy at fixed cell
size.

## 10. Multi-GPU domain decomposition (forward-looking)

Single-GPU FDTD reaches the memory-bandwidth wall quickly: arithmetic
intensity is around one flop per byte. Beyond a single-GPU working
set (24–80 GB) the only path forward is domain decomposition: split
the grid across GPUs, each owning a slab, exchanging a one-cell-deep
tangential field halo each leapfrog step.

The Yee plan (`ROADMAP.md` Phase 2 and Phase 4) is one rank per GPU,
NCCL for halo exchange (cudarc ships safe bindings), with the boundary
swap overlapped against interior kernels: interior first, boundary —
consuming the freshly arrived halo — second on the same CUDA stream.
Bandwidth scales with slab cross-section, so overhead falls as the
per-GPU domain grows. Phase 4 also handles load-balancing and CFS-PML
corner / edge bookkeeping across decomposition planes.

Phase 2.0 is deliberately single-CPU, single-threaded, scalar,
vacuum-only so the walking skeleton can integrate with the rest of
the workspace before kernels are tuned.

## 11. Limitations and known gaps (Phase 2.0)

The walking skeleton is a correctness-first reference, not
performance-tuned production code. Known gaps, in roadmap order:

- **Hard PEC outer boundary.** All six outer faces clamp tangential
  `E` to zero (`crates/yee-fdtd/src/boundary.rs::apply_pec`). Fine
  for closed-cavity tests, wrong for radiation. CPML lands in
  Phase 2.1 (`ROADMAP.md` Phase 2; planned
  `crates/yee-fdtd/src/cpml.rs`).
- **Vacuum-only materials.** `YeeGrid::vacuum` is the only public
  constructor. Per-cell materials and the two-term lossy update land
  in Phase 2.1; dispersive media via ADE or PLRC in Phase 2.2
  (`ROADMAP.md` Phase 2).
- **Single-CPU scalar (no SIMD, no rayon, no GPU).** Kernels in
  `crates/yee-fdtd/src/update.rs` are triply-nested `f64` loops.
  CUDA E/H kernels are core Phase 2; multi-GPU NCCL decomposition is
  Phase 4.
- **No subgridding, no conformal cells.** Geometry is staircased.
  Dey-Mittra cells are Phase 2.3; subgridding is post-Phase 2.
- **No NTFF, no lumped ports, no waveguide ports.** Pattern
  extraction and S-parameter ports are Phase 2.2–2.3.

## 12. References

- Yee, K. S. "Numerical Solution of Initial Boundary Value Problems
  Involving Maxwell's Equations in Isotropic Media." *IEEE Trans.
  Antennas Propag.* 14.3 (1966), pp. 302–307.
- Taflove, A., and Hagness, S. C. *Computational Electrodynamics: The
  Finite-Difference Time-Domain Method.* 3rd ed. Artech House, 2005.
  (§3 Yee update; §4.7 Courant; §5 TF/SF; §7 CPML; §8 NTFF;
  §9 dispersive media.)
- Gedney, S. D. *Introduction to the Finite-Difference Time-Domain
  (FDTD) Method for Electromagnetics.* Morgan & Claypool, 2011.
- Roden, J. A., and Gedney, S. D. "Convolution PML (CPML): An
  Efficient FDTD Implementation of the CFS-PML for Arbitrary Media."
  *Microwave Opt. Technol. Lett.* 27.5 (2000), pp. 334–339.
- Namiki, T. "A New FDTD Algorithm Based on Alternating-Direction
  Implicit Method." *IEEE Trans. Microwave Theory Tech.* 47.10
  (1999), pp. 2003–2007. (ADI-FDTD context.)
- Luebbers, R., and Hunsberger, F. "FDTD for Nth-Order Dispersive
  Media." *IEEE Trans. Antennas Propag.* 40.11 (1992), pp. 1297–1301.
- Berenger, J.-P. "A Perfectly Matched Layer for the Absorption of
  Electromagnetic Waves." *J. Comput. Phys.* 114.2 (1994),
  pp. 185–200.
- Dey, S., and Mittra, R. "A Locally Conformal Finite-Difference
  Time-Domain (FDTD) Algorithm for Modeling Three-Dimensional
  Perfectly Conducting Objects." *IEEE Microwave Guided Wave Lett.*
  7.9 (1997), pp. 273–275.
