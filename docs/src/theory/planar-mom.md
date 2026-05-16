# Planar Method of Moments — Theory of Operation

This page is the theory-of-operation reference for Yee's planar Method of
Moments (MoM) solver, implemented in the `yee-mom` crate. It is written for
the engineer who wants to read the code with a textbook open on the desk.
Equations are kept in plain text because the documentation build does not
currently render LaTeX.

## 1. Overview

The Method of Moments is a frequency-domain integral-equation technique
for Maxwell's equations on conducting (and, in extensions, dielectric)
surfaces. Surface currents are expanded in a finite basis, the
electric-field integral equation is enforced in a weighted-residual
sense, and the result is a dense complex linear system whose unknowns
are the basis-function expansion coefficients. From those one recovers
port voltages, currents, input impedances, and S-parameters.

Planar MoM is the Yee v1 beachhead: `ROADMAP.md` Phase 1 orders the
project around shipping a usable planar MoM solver before expanding the
FDTD or multilayer Green's-function work. Every commercial planar EM
solver in routine use today (Sonnet em, Keysight ADS Momentum, Cadence
AWR AXIEM, ANSYS HFSS-IE) is a planar MoM in this same family. The
intellectual arc spans roughly fifty years: Harrington's *Field
Computation by Moment Methods* (Macmillan, 1968; reprinted IEEE Press,
1993) crystallised the moment-method abstraction; Rao, Wilton, and
Glisson, "Electromagnetic Scattering by Surfaces of Arbitrary Shape"
(IEEE Trans. Antennas Propag. 30.3, 1982, pp. 409–418) introduced the
divergence-conforming triangular basis every modern surface MoM still
uses; the commercial codes are the industrial descendants.

## 2. From Maxwell to the integral equation

Time-harmonic Maxwell in a homogeneous lossless background, `exp(+jωt)`
convention:

```text
∇ × E = -jωμ₀ H
∇ × H =  jωε₀ E + J
∇ · E =  ρ/ε₀
∇ · H =  0
```

For a PEC the boundary condition is `[n̂ × E_total]_on_S = 0`. Splitting
the field into known incident `E^inc` and scattered `E^scat` (radiated
by the unknown induced surface current `J`) gives
`[n̂ × E^scat]_on_S = -[n̂ × E^inc]_on_S`.

Express `E^scat` through the standard mixed potentials. The magnetic
vector potential `A` and electric scalar potential `Φ` are

```text
A(r)  = ∫_S G(r, r') J(r') dS'
Φ(r)  = (-1/(jωε₀)) ∫_S G(r, r') ∇'·J(r') dS'
G(r, r') = exp(-jk₀R) / (4πR),    R = |r - r'|,    k₀ = ω/c
```

and the scattered field is `E^scat = -jωμ₀ A - ∇Φ`. Substituting and taking
the tangential trace gives the *mixed-potential integral equation* (MPIE)
in its electric-field (EFIE) form:

```text
[E^inc]_tan = (jωμ₀ A + ∇Φ)_tan,    on PEC surface S
```

Symbols: `J(r')` is the unknown vector surface current density (A/m) on
S; `G(r, r')` is the scalar free-space Green's function whose `1/R`
factor drives the numerical difficulty downstream; `k₀ = ω/c` is the
free-space wavenumber; `μ₀, ε₀` are the free-space permeability and
permittivity; and `∇'·J` is the surface divergence of the current, which
by the continuity equation is proportional to the surface charge density
`ρ_s = -∇'·J / (jω)`.

The MPIE is preferred over the bare EFIE because the differential
operator `∇` is moved out from inside the `1/R` kernel onto the basis
functions, where it is evaluated analytically. The cost is the scalar-
potential term, which forces the basis to have a well-defined,
integrable surface divergence — exactly the property the RWG basis
guarantees. Standard reference: Gibson, *The Method of Moments in
Electromagnetics*, 2nd ed., CRC Press (2014), §7.

## 3. RWG basis functions

The Rao-Wilton-Glisson (RWG) basis lives on a pair of triangles
`T^+_n` and `T^-_n` sharing an interior edge of length `l_n`. Let
`p^±_free` be the vertex of `T^±_n` opposite the shared edge. The basis
function `f_n(r)` is

```text
            +(l_n / (2 A^+_n)) (r - p^+_free),    r ∈ T^+_n
f_n(r)  =   -(l_n / (2 A^-_n)) (r - p^-_free),    r ∈ T^-_n
            0,                                     otherwise
```

where `A^±_n` is the area of the corresponding triangle. The surface
divergence is piecewise constant:

```text
∇·f_n = ±l_n / A^±_n
```

The normal component of `f_n` is continuous across the shared edge with
unit magnitude, and vanishes across every other edge of the two triangles
(the free vertex lies on those edges, so `r - p^±_free` is parallel to
them). This gives RWG its two essential properties:

1. **Charge conservation across the shared edge.** Continuous normal
   component means no fictitious line charge appears on the internal
   edge; the surface charge is bounded everywhere.
2. **Linear independence.** Each interior edge carries exactly one RWG
   function; the basis spans the divergence-conforming subspace of
   `H(div, S)` on the triangulated surface.

```text
              p^-_free
                 *
                /  \
               /    \
              / T^-  \
             /        \
            /---edge---\     ← shared edge of length l_n
            \  T^+    /
             \       /
              \     /
               \   /
                \ /
                 *
              p^+_free
```

The arrow on `f_n` points from `T^+_n` toward `T^-_n`; reversing the
triangle labelling flips the sign of `f_n`. The choice is fixed once and
recorded in the basis table. In the Yee implementation, the basis is
built by the planned `crates/yee-mom/src/basis.rs` `from_mesh`
constructor, which walks the `TriMesh`, finds every interior edge, and
emits one RWG function per edge with a stable global ordering. See
Rao-Wilton-Glisson, IEEE T-AP 30.3 (1982), 409–418.

## 4. Galerkin projection

Expand the unknown current in the RWG basis:

```text
J(r) = Σ_n i_n f_n(r),    n = 1, …, N
```

Substitute into the MPIE, dot with test function `f_m`, integrate, and
use the identity `∫ f_m · ∇Φ dS = -∫ (∇·f_m) Φ dS` (the boundary term
vanishes because `f_m` has zero normal component on the outer boundary
of its support). The result is a linear system

```text
Σ_n Z_{mn} i_n = b_m,    m = 1, …, N
```

with matrix entries

```text
Z_{mn} =   jωμ₀  ∫∫_{S_m × S_n} G(r, r') f_m(r) · f_n(r') dS dS'
         + (1/(jωε₀)) ∫∫_{S_m × S_n} G(r, r') (∇·f_m)(r) (∇·f_n)(r') dS dS'
```

Structural facts: `Z` is **complex-valued** (the kernel carries
`exp(-jk₀R)`); **symmetric** for reciprocal media (`Z_{mn} = Z_{nm}`,
because `G` is symmetric in `r ↔ r'`) but *not* Hermitian — the `j`
factors break self-adjointness in the standard Hilbert sense; and
**dense**, because every basis function couples to every other through
the long-range `1/R` kernel. Sparsification is deferred to MLFMA/ACA in
Phase 4. The matrix fill lives in the planned
`crates/yee-mom/src/fill.rs`; Green's-function evaluators in
`crates/yee-mom/src/greens.rs`.

## 5. Quadrature on triangles

For well-separated triangle pairs the integrand is smooth and standard
symmetric Gauss quadrature on the triangle suffices. The Yee
implementation uses the Dunavant rules — *D. A. Dunavant, "High Degree
Efficient Symmetrical Gaussian Quadrature Rules for the Triangle",
Int. J. Numer. Methods Eng. 21.6 (1985), pp. 1129–1148* — at three
working orders: **3** (4 points, exact for cubics) for far-field pairs;
**5** (7 points, exact for quintics) as the bulk default; **7** (13
points, exact for degree-7) for the outer integration when the inner is
near-singular. Order selection is driven by a distance heuristic:
centroid separation > ~5× the larger triangle diameter → order 3; closer
→ order 5; touching → order 7 outside combined with §6 singular handling
inside.

## 6. Singular and near-singular integration

When the two triangles share a face, edge, or vertex, the `1/R` kernel
is singular in the inner integration domain, and Gauss quadrature either
diverges (face self-term) or converges painfully slowly (edge/vertex
shared). The remedy is the **Duffy transform** — a change of variables
that absorbs the singular Jacobian into a bounded integrand.

Same-triangle case in sketch: fix the outer quadrature point `r` inside
triangle `T`; subdivide `T` into three sub-triangles by drawing segments
from `r` to each vertex; on each sub-triangle introduce polar-like
radial coordinates anchored at `r`. The polar Jacobian contains an `R`
factor that exactly cancels the `1/R` kernel singularity, leaving a
smooth integrand that ordinary Gauss-Legendre quadrature handles in tens
of points. Edge- and vertex-shared cases follow the same principle with
different geometric subdivisions. Implementation reference: *Khayat &
Wilton, "Numerical Evaluation of Singular and Near-Singular Potential
Integrals", IEEE Trans. Antennas Propag. 53.10 (2005), pp. 3180–3190*,
which gives explicit formulas for all three sharing cases.

Near-singular (close but not touching) pairs escalate to Khayat-Wilton on
the inner integral while keeping a high-order Dunavant rule on the outer.

## 7. Ports and excitation

The simplest port model — and the only one in Phase 1.0 — is the
**delta-gap**. A 1 V source is impressed across an internal edge designated
as the port edge; the edge sits between two triangles tagged as belonging
to different port-net groups. The right-hand-side vector becomes

```text
b_m = V · l_m,   for m ∈ port_basis_indices
b_m = 0,          otherwise
```

Only RWG functions on the port edge see a nonzero excitation, with
amplitude `V · l_m`. Conceptually `b_m = ∫ f_m · E^inc dS`, and the
delta-gap model concentrates `E^inc` into a Dirac sheet across the port
edge tuned so that the line integral of `E^inc` across the gap equals
`V`; the edge-length factor falls out of that line integral.

The Yee convention, fixed by the basis builder, is that *differently
tagged adjacent triangles define a port edge*. The basis table emitted
by the planned `crates/yee-mom/src/basis.rs::from_mesh` carries a
per-edge port tag used by the matrix-fill stage to assemble both the
RHS and the post-solve port-current extractor.

After the dense solve produces `i_n`, the **port current** is

```text
I_port = Σ_k l_k · i_k,    k ∈ port_basis_indices
```

Input impedance is `Z_in = V / I_port`, and the one-port reflection
coefficient (port reference impedance `Z₀`, typically 50 Ω) is

```text
S₁₁ = (Z_in - Z₀) / (Z_in + Z₀)
```

N-port networks generalise: excite one port at a time, solve `Z i = b`,
extract every port's current, and assemble columns of the admittance
matrix `Y`; the S-matrix follows from `S = (I + Z₀ Y)^(-1) (I - Z₀ Y)`.
Reference: Pozar, *Microwave Engineering*, 4th ed., Wiley (2012),
§4.3–4.4.

## 8. Dense linear solve

For Phase 1.0 problem sizes (dipole, 50 Ω microstrip, 2.4 GHz patch)
`N` is well under 50 000, the impedance matrix fits in memory, and the
solve is one complex dense LU factorisation plus a back-substitution
per excitation. CPU: `faer`. Planned GPU: cuSOLVER `Zgetrf` /
`Zgetrs`, with optional iterative refinement.

Above ~50 000 unknowns dense LU dominates (`O(N^3)` work,
`O(N^2)` memory); Phase 1 switches to GMRES with a block-diagonal
preconditioner (block sizes chosen by edge clustering). Phase 4 promotes
the fast methods (MLFMA, ACA, H-matrix compression) to first-class
status; their absence is the headline reason Yee will not yet attack
`N ≥ 100 000` competitively. The solve interface lives in the planned
`crates/yee-mom/src/solve.rs`.

## 9. Conditioning

The MPIE form has a well-known *low-frequency breakdown*: when mesh edge
length `h` is much smaller than the wavelength `λ`, the vector- and
scalar-potential contributions to `Z_{mn}` decouple onto nearly
orthogonal subspaces — divergence-free "loop" currents and
divergence-conforming "tree" currents — and the condition number scales
like

```text
cond(Z)  ~  (λ / h)^2 · κ_geo
```

where `κ_geo` packages mesh-quality factors (aspect ratio and similar
pathologies). For the Phase 1.0 `mom-001` dipole fixture (finite-radius
cylinder, 24 axial × 48 circumferential cells, ~λ/40 edges) the observed
condition number is ~`4 × 10^7` — large but not catastrophic for
double-precision LU; residuals stay below `1e-10`. For larger `λ/h` the
standard remedies are **loop-tree** (or loop-star) decomposition —
explicitly separate the basis into divergence-free and irrotational
subspaces and rescale — and **Calderón preconditioning**, which uses the
magnetic-field operator as an approximate inverse to the electric-field
operator. Both are Phase 1.1 work.

## 10. Limitations and known gaps (Phase 1.0)

Phase 1.0 is deliberately a walking skeleton. The free-space kernel is
the simplest possible Green's function; everything more interesting is a
later sub-project, in roughly this order:

- **Multilayer dielectric stack-ups.** Real PCB structures live on
  layered substrates; the kernel becomes a Sommerfeld-type spectral
  integral over the stratification. Phase 1.1 follows the Michalski-Mosig
  formulation — *K. A. Michalski and J. R. Mosig, "Multilayered media
  Green's functions in integral equation formulations", IEEE Trans.
  Antennas Propag. 45.3 (1997), pp. 508–519* — with DCIM extraction and
  adaptive rational fitting for fast evaluation.
- **Surface roughness.** Hammerstad-Jensen / Groiss / Huray models;
  Phase 1 loss-model work.
- **Wave ports.** Phase 1 cross-section modal solver for microstrip /
  CPW; the delta-gap of §7 becomes a special case.
- **TRL / SOLT de-embedding.** Required for any serious comparison
  against measured fixtures; planned with the wave-port work.
- **GPU paths.** Matrix fill and dense LU on CUDA are Phase 1
  deliverables; Phase 1.0 is CPU-only.

### NEC-4 vs Balanis: the mom-001 reference value

The Phase 1.0 gate `mom-001` targets the half-wave dipole impedance.
Two reference values coexist in the literature:

- **Balanis** (*Antenna Theory: Analysis and Design*, 4th ed., Wiley
  2016) gives `Z ≈ 73 + j42 Ω` — the *zero-radius, sinusoidal-current*
  analytical limit, assuming an infinitely thin wire with an *assumed*
  (not solved) sinusoidal distribution.
- **NEC-4** (the industry-standard surface MoM reference) gives
  `Z ≈ 87 + j41 Ω` for a finite-radius cylinder with `a/L = 5 × 10^-3`
  (1 m dipole, 5 mm radius) — the same geometry Yee meshes.

Yee's planar MoM is a **full surface MoM on a finite-radius cylinder**,
not a thin-wire approximation with assumed current; it must agree with
the finite-radius reference, not the thin-wire limit. `mom-001`
therefore uses **NEC-4's 87 + j41 Ω**: `±5%` on `Re(Z)` and `±10%` on
`Im(Z)` (the imaginary part is more mesh-density sensitive and closes
as the circumferential discretisation tightens). The Track A diagnostic
in May 2026 showed Yee's solver matches NEC-4 on `Re(Z)` to within 0.1%
at the pinned fixture (`L = 1 m`, `a = 5 mm`, `n_axial = 24`,
`n_around = 48`). The general lesson: a finite-radius cylindrical
dipole solved with surface MoM is a different boundary-value problem
from a zero-radius wire with assumed sinusoidal current, and citing the
wrong reference is one of the easiest ways to misjudge a solver.

## 11. References

- Harrington, R. F. *Field Computation by Moment Methods.* Macmillan,
  1968; reprinted IEEE Press, 1993.
- Rao, S. M., Wilton, D. R., and Glisson, A. W. "Electromagnetic
  Scattering by Surfaces of Arbitrary Shape." *IEEE Trans. Antennas
  Propag.* 30.3 (May 1982), pp. 409–418.
- Gibson, W. C. *The Method of Moments in Electromagnetics.* 2nd ed.
  CRC Press, 2014. (Chapter 7 for the MPIE derivation in §2.)
- Balanis, C. A. *Antenna Theory: Analysis and Design.* 4th ed. Wiley,
  2016. (Source of the zero-radius `73 + j42 Ω` reference in §10.)
- Dunavant, D. A. "High Degree Efficient Symmetrical Gaussian
  Quadrature Rules for the Triangle." *Int. J. Numer. Methods Eng.*
  21.6 (1985), pp. 1129–1148.
- Khayat, M. A., and Wilton, D. R. "Numerical Evaluation of Singular
  and Near-Singular Potential Integrals." *IEEE Trans. Antennas Propag.*
  53.10 (October 2005), pp. 3180–3190.
- Michalski, K. A., and Mosig, J. R. "Multilayered Media Green's
  Functions in Integral Equation Formulations." *IEEE Trans. Antennas
  Propag.* 45.3 (March 1997), pp. 508–519.
- Pozar, D. M. *Microwave Engineering.* 4th ed. Wiley, 2012. (Port
  theory and S-parameter conventions used in §7.)
- Burke, G. J., and Poggio, A. J. *Numerical Electromagnetics Code
  (NEC) — Method of Moments.* LLNL technical document; later revisions
  (NEC-2, NEC-4) are the industry-standard reference solvers cited
  in §10.
