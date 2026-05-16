# RWG Basis, MPIE, and Phase 1.0 Quadrature — Theory of Operation

This page is the theory-of-operation reference for the discretisation
choices that landed in Yee's Phase 1.0 free-space MoM solver. The
companion chapter [Planar Method of Moments](./planar-mom.md) covers
the integral-equation setup at the textbook level — Maxwell to MPIE,
Galerkin projection, dense LU, port theory. This chapter is the deep
dive into what `crates/yee-mom/` actually computes: the Rao–Wilton–
Glisson (RWG) basis, the mixed-potential integral equation (MPIE)
matrix entries, the Dunavant Gauss quadrature schemes used for
well-separated triangle pairs, the Duffy-transform singular-pair
handling, and the delta-gap port excitation.

Same audience as the GP / FDTD chapters: an engineer reading source
with a textbook open on the desk. Equations are written in KaTeX so the
inline notation can stay close to the Rust source.

## 1. Introduction

The Phase 1.0 deliverable was a free-space surface MoM that reproduces
NEC-4's published half-wave-dipole input impedance. The reference
fixture is a 1 m cylinder of radius $a = 5\,\text{mm}$ excited by a
unit delta-gap at the centre; the NEC-4 result is
$Z \approx 87 + j41\,\Omega$ at resonance. Yee's solver currently
returns

$$
Z \approx 87 + j50\,\Omega
$$

on a 24-axial × 176-circumferential cylinder mesh — inside the
`mom-001` validation envelope ($\pm 5\%$ on $\operatorname{Re}(Z)$,
$\pm 10\%$ on $\operatorname{Im}(Z)$, with the imaginary tolerance
acknowledging mesh-density sensitivity that closes as the
circumferential cell count tightens).

Hitting that envelope from a cold start required four pieces working
together:

- a **divergence-conforming surface-current basis** (RWG, §2) so the
  charge integral in the MPIE has a well-defined surface divergence;
- the **MPIE matrix entries** (§3) — a symmetric vector- and
  scalar-potential pair, not the raw EFIE, so the $\nabla$ stays out
  of the singular kernel;
- **Dunavant symmetric Gauss quadrature** (§4) for the well-separated
  triangle pairs that dominate the matrix-fill cost;
- **singularity-subtraction plus a Duffy-transform path** (§5–§6) for
  the same-triangle and shared-vertex pairs where the $1/R$ kernel
  goes singular and ordinary Gauss quadrature diverges.

Section 7 covers the delta-gap port treatment that turns the dense
solve into a one-port S-parameter. Section 8 lists the pieces that are
**out of scope** for this chapter and which other Phase 1.x
sub-projects own them.

## 2. RWG basis functions

The Rao–Wilton–Glisson basis (Rao, Wilton, Glisson 1982) is the de
facto standard for surface-MoM on PEC objects. Each interior edge of
the triangulation carries exactly one basis function, defined on the
pair of triangles sharing that edge.

Let edge $e_n$ of length $\ell_n$ be shared by triangles $T_n^+$ and
$T_n^-$ with areas $A_n^+, A_n^-$, and let $\mathbf{r}_n^\pm$ be the
vertex of $T_n^\pm$ opposite the shared edge (the "free" vertex). The
RWG basis function is

$$
\mathbf{f}_n(\mathbf{r}) =
\begin{cases}
+\dfrac{\ell_n}{2 A_n^+}\bigl(\mathbf{r} - \mathbf{r}_n^+\bigr), & \mathbf{r} \in T_n^+, \\[4pt]
-\dfrac{\ell_n}{2 A_n^-}\bigl(\mathbf{r} - \mathbf{r}_n^-\bigr), & \mathbf{r} \in T_n^-, \\[4pt]
\mathbf{0}, & \text{elsewhere.}
\end{cases}
$$

Its surface divergence is **piecewise constant**:

$$
\nabla_s \cdot \mathbf{f}_n =
\begin{cases}
+\ell_n / A_n^+, & \mathbf{r} \in T_n^+, \\
-\ell_n / A_n^-, & \mathbf{r} \in T_n^-, \\
0, & \text{elsewhere,}
\end{cases}
$$

so the scalar-potential integral in the MPIE reduces to two constants
times surface integrals of the Green's function. The normal component
of $\mathbf{f}_n$ at the shared edge evaluates to

$$
\hat{\mathbf{n}} \cdot \mathbf{f}_n \bigr|_{\text{shared edge}} = \ell_n / \ell_n = 1,
$$

i.e. unit current crosses the edge with continuous normal trace and
no fictitious line charge.

Two implementation-side conventions that the rest of the chapter
relies on:

**Sign convention.** The $+$ / $-$ labelling of the two adjacent
triangles is whichever order the mesh enumerator visits them; it
carries no geometric meaning. Reversing the labelling flips
$\mathbf{f}_n$'s sign, which the linear system absorbs into a sign on
the expansion coefficient $i_n$. The Galerkin matrix and the
post-solve port-current extraction are unaffected.

**Port-edge tagging.** Yee tags each triangle with a `port_tag : u32`.
An edge between two triangles whose `port_tag` values are **different
and both non-zero** is a port edge — it sits on the boundary between
two driven nets. Edges between same-tag or zero-tag triangles are
ordinary interior edges. Boundary edges (touched by one triangle
only) are silently dropped: they carry no RWG basis function. This
convention is enforced inside `RwgBasis::from_mesh`.

A typical caller builds the basis directly from a `TriMesh`:

```rust,ignore
use yee_mom::__internal::{build_basis, RwgBasis};
use yee_mesh::TriMesh;

let mesh: TriMesh = /* ... build vertices, triangles, port tags ... */;
let basis: RwgBasis = build_basis(&mesh)?;
// basis.edges.len() is the number of RWG basis functions
// (one per non-boundary edge); each basis.edges[k] carries
// per-edge length, port_tag, and the +/- triangle and free-vertex
// indices used by the matrix-fill stage.
```

`RwgBasis::from_mesh` rejects two degenerate cases: any non-positive
or NaN triangle area (`Error::Invalid`), and any edge shared by three
or more triangles (non-manifold mesh, also `Error::Invalid`).

## 3. The MPIE matrix entries

The companion chapter derives the mixed-potential integral equation
from Stratton's EFIE; we restate only the discretised matrix entry
here. With the current expanded in RWG functions
$\mathbf{J}(\mathbf{r}) = \sum_n i_n \mathbf{f}_n(\mathbf{r})$ and a
Galerkin test
$\mathbf{f}_m$, the impedance-matrix entry is

$$
Z_{mn}
= j\omega\mu_0 \!\!\iint_{S_m \times S_n}\!\! G(\mathbf{r}, \mathbf{r}')\,\mathbf{f}_m(\mathbf{r}) \cdot \mathbf{f}_n(\mathbf{r}')\, dS\, dS'
+ \frac{1}{j\omega\epsilon_0} \!\!\iint_{S_m \times S_n}\!\! G(\mathbf{r}, \mathbf{r}')\bigl(\nabla_s \cdot \mathbf{f}_m\bigr)(\mathbf{r}) \bigl(\nabla_s \cdot \mathbf{f}_n\bigr)(\mathbf{r}')\, dS\, dS'.
$$

The kernel is the **free-space scalar Green's function**

$$
G(\mathbf{r}, \mathbf{r}') = \frac{e^{-j k_0 R}}{4\pi R}, \qquad R = \lVert \mathbf{r} - \mathbf{r}' \rVert, \qquad k_0 = \omega / c.
$$

Two structural facts follow from the MPIE form:

- The **divergence factors are piecewise constant** (§2), so the
  scalar-potential term reduces to a sum of four
  scalar-Green's-on-triangle integrals weighted by
  $\pm \ell_m / A_m^\pm \cdot \pm \ell_n / A_n^\pm$. No spatial
  derivative of $G$ appears anywhere — the $\nabla$ moves onto the
  basis functions where it is evaluated analytically.
- The **vector-potential factor** $\mathbf{f}_m \cdot \mathbf{f}_n$ is
  polynomial in $\mathbf{r}, \mathbf{r}'$ (linear in each), so its
  contribution to the smooth integrand reduces to a low-degree
  polynomial against the scalar kernel.

Both terms share the same scalar $G$, and the matrix is therefore
**symmetric but not Hermitian** — the $j$ factors break the standard
inner-product self-adjointness. Phase 1.0 ships only the free-space
kernel; the [`Greens`] trait abstraction (in `crates/yee-mom/src/greens.rs`)
exists as the swap point for multilayer kernels and is consumed by
the multilayer placeholder.

## 4. Dunavant Gauss quadrature

For triangle pairs that are well separated relative to the larger
triangle's diameter, both integrands are smooth functions of
$(\mathbf{r}, \mathbf{r}')$, and the obvious thing — symmetric Gauss
quadrature on the triangle — works. Yee uses the **Dunavant rules**
(Dunavant 1985): families of barycentric-symmetric quadrature points
$\{(\xi_i, \eta_i, w_i)\}$ on the master triangle $T_0$ with vertices
$(0,0), (1,0), (0,1)$, exact for polynomials up to a stated degree.

For a triangle $T$ with vertices $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$
and area $A$, the affine map

$$
\mathbf{r}(\xi, \eta) = (1 - \xi - \eta)\mathbf{v}_0 + \xi \mathbf{v}_1 + \eta \mathbf{v}_2
$$

sends $T_0 \to T$ with Jacobian $2A$. A scalar integral transforms as

$$
\int_T g(\mathbf{r})\, dS = 2 A \int_{T_0} g(\mathbf{r}(\xi, \eta))\, d\xi\, d\eta \approx 2 A \sum_{i=1}^{N_q} w_i\, g\bigl(\mathbf{r}(\xi_i, \eta_i)\bigr).
$$

Yee implements three Dunavant orders, selected per pair by a
distance heuristic on the centroid separation $d$ relative to the
larger triangle diameter $h$:

| Order | $N_q$ | Polynomial exact to | Use site |
|-------|-------|---------------------|----------|
| 3     | 4     | cubics              | far-field pairs, $d \gtrsim 5h$ |
| 5     | 7     | quintics            | bulk default, near-pairs outer integral |
| 7     | 13    | degree-7            | near-singular outer when inner uses Duffy |

The bulk default is order 5; the matrix-fill module reaches for order
7 on the outer integral whenever the inner integral has switched to
the Duffy-transform path (§6) to keep the outer truncation error from
becoming the bottleneck.

## 5. Singularity handling: subtraction

When the inner triangle $T_n$ overlaps the outer triangle $T_m$ or
shares an edge or a vertex with it, $R \to 0$ within the inner
integration domain, the kernel $1/R$ diverges, and tensor-product
Gauss quadrature on the master triangle fails — either by literally
sampling at the singularity (face self-term) or by converging
painfully slowly (shared edge / shared vertex). Yee handles this in
two complementary ways: **singularity subtraction** in the smooth path
of the [`Greens`] trait (§5), and a **Duffy-transform reparametrisation**
on the singular-pair fill path (§6).

The subtraction split is

$$
G(\mathbf{r}, \mathbf{r}') = \underbrace{G(\mathbf{r}, \mathbf{r}') - \frac{1}{4\pi R}}_{G_{\text{smooth}}} + \underbrace{\frac{1}{4\pi R}}_{G_{\text{sing}}}.
$$

$G_{\text{smooth}}$ is bounded as $R \to 0$ — its limit is

$$
\lim_{R \to 0} G_{\text{smooth}}(\mathbf{r}, \mathbf{r}') = -\frac{j k_0}{4\pi},
$$

so an ordinary Gauss quadrature handles it (Wilton et al. 1984). The
singular part $1/(4\pi R)$ has a closed-form analytic integral over a
flat triangle — the standard Wilton–Rao–Glisson–Schaubert closed-form
expressions, parameterised by the projection of the source vertex onto
the target triangle's plane plus the three signed in-plane distances
to the triangle's edges. The smooth and singular integrals are
computed separately and summed.

In `yee-mom`, `Greens::scalar_vector_smooth` and
`Greens::scalar_scalar_smooth` return the subtracted kernel directly,
including the $R \to 0$ limit so that callers do not have to handle
the coincident-point branch themselves.

## 6. Singular pairs: the Duffy transform

Subtraction alone is not enough for the **shared-vertex** case. Even
after $1/R$ is removed, the original $1/R$ integral still has to be
evaluated on a pair of triangles meeting at a corner, and the corner
singularity is too sharp for standard tensor-product quadrature on the
master triangle to converge cleanly.

The remedy is the **Duffy transform** (Duffy 1982). On a triangle
$T$ with a singular vertex at $(0,0)$ in the master coordinates, the
mapping

$$
(\xi, \eta) = (u, u v), \qquad u \in [0, 1], \quad v \in [0, 1],
$$

sends the unit square $[0,1]^2$ onto $T_0$ with Jacobian $u$. The
$1/R$ singularity at the corner becomes proportional to $1/u$ in the
new variables; the Jacobian's leading $u$ factor cancels it, leaving
a bounded integrand of the form $\sin\theta \cdot (\text{bounded})$
(after a final angular reinterpretation) that ordinary Gauss–Legendre
on $[0,1]^2$ resolves in tens of points.

In `yee-mom`, the singular-pair fill path applies the Duffy transform
to the **full** Green's function $G$ — *not* `Greens::scalar_smooth`.
Using the smooth kernel inside Duffy is a subtle bug: the Duffy
Jacobian is the regulator that makes the $1/R$ part integrable, so
subtracting $1/R$ first means the Duffy path then double-subtracts
it. An early Phase 1.0 prototype hit exactly this trap; the dipole
$\operatorname{Re}(Z)$ came in roughly 25% high until the singular
fill path was switched to call `scalar` directly (with a bit-exact
$r > 0$ guard returning the analytic $-j k_0 / (4\pi)$ limit at
coincident points to keep the kernel bounded).

The corresponding source guard in `greens.rs` is

```rust,ignore
pub fn scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64 {
    let r = (r1 - r2).norm();
    if r == 0.0 {
        // Analytic limit of G - 1/(4πR) at R = 0.
        return Complex64::new(0.0, -self.k0.re / (4.0 * std::f64::consts::PI));
    }
    let g = self.scalar(r1, r2);
    g - Complex64::new(1.0 / (4.0 * std::f64::consts::PI * r), 0.0)
}
```

— the coincident branch returns the analytic limit so that any caller
who happens to evaluate `scalar_smooth` at $R = 0$ (the Gauss point
on $T_0$ that exactly maps to the singular vertex, weighted by zero
Jacobian) gets a finite number instead of `NaN`.

## 7. Delta-gap port excitation

Phase 1.0 ships exactly one port model: the **delta-gap**. A 1 V
source is impressed across the port edge — concretely, across every
edge whose `port_tag` matches the driven-port label. The right-hand-
side vector of the dense linear system $Z \mathbf{i} = \mathbf{b}$
is

$$
b_n =
\begin{cases}
V_{\text{port}} \cdot \ell_n, & n \in \mathcal{P}, \\
0, & \text{otherwise,}
\end{cases}
$$

where $\mathcal{P}$ is the set of basis indices whose RWG edge has
the matching port tag and $V_{\text{port}} = 1\,\text{V}$ for the
canonical excitation. Conceptually $b_n = \int \mathbf{f}_n \cdot
\mathbf{E}^{\text{inc}}\, dS$, and the delta-gap model concentrates
$\mathbf{E}^{\text{inc}}$ into a Dirac sheet across the port edge
tuned so that the line integral across the gap equals $V_{\text{port}}$;
the edge-length factor falls out of that line integral.

After the LU solve produces $\mathbf{i}$, the **port current** is the
Galerkin projection of the current onto the port edges,

$$
I_{\text{port}} = \sum_{n \in \mathcal{P}} \ell_n \cdot i_n,
$$

which is the same inner-product structure as the RHS by construction.
The **input impedance** at the port is

$$
Z_{\text{in}} = \frac{V_{\text{port}}}{I_{\text{port}}},
$$

and the one-port reflection coefficient against the port reference
impedance $Z_0$ (typically $50\,\Omega$) is
$S_{11} = (Z_{\text{in}} - Z_0) / (Z_{\text{in}} + Z_0)$.

The full Phase 1.0 dipole sweep is one call:

```rust,ignore
use yee_core::FreqRange;
use yee_mom::{PlanarMoM, SParameters};
use yee_mesh::TriMesh;

let mesh: TriMesh = /* dipole cylinder mesh, centre triangles tagged 1/2 */;
let freq: FreqRange = FreqRange::new(50e6, 500e6, 51)?;
let solver = PlanarMoM::default();
let s: SParameters = solver.run(&mesh, freq)?;
// s.freq_hz / s.data carry the swept S₁₁; write to Touchstone via
// SParameters::write_touchstone(&path, 50.0).
```

`PlanarMoM::default().run` wires up the canonical 1 V delta-gap on
`port_tag = 1`, the free-space Green's function, the RWG basis builder,
the dense LU, and the Touchstone wrapper. For multilayer / multi-port
configurations where the kernel needs to vary by frequency or by
substrate stack-up, the lower-level

```rust,ignore
use yee_mom::__internal::{z_in_with_greens, z_in_free_space};

let z_in = z_in_free_space(&mesh, /*port_tag=*/ 1, /*freq_hz=*/ 1.5e9)?;
```

helpers expose the basis-build → fill → solve → port-extract path
parameterised over any `Greens` implementation. The `z_in_with_greens`
generic entry point is the integration-test surface used to validate
the Phase 1.1 multilayer kernel against the free-space baseline on
the same mesh, and it is where the Track RRR multi-port extension
will hang the per-port-net excitation loop (forward-looking — RRR has
not landed in `main` yet).

## 8. What's not in this chapter

Several Phase 1.x deliverables touch the same code paths discussed
above but introduce orthogonal physics that deserves its own
treatment. They are deliberately out of scope here:

- **Multilayer Green's functions.** The Phase 1.1 `MultilayerGreens`
  type plugs into the [`Greens`] trait but currently ships a
  **one-image Discrete Complex Image Method (DCIM) placeholder**;
  proper Sommerfeld-integral / multi-image extraction is the Phase
  1.1.1 deliverable. `mom-002` (microstrip $Z_0$) and `mom-003`
  (2.4 GHz patch) run with loose tolerances until then. The companion
  [planar-mom](./planar-mom.md) §10 plus the Michalski–Mosig
  references collect the formulation references.
- **Wave-port modal extraction.** Phase 1.3.0 ships a `WavePort` API
  placeholder whose modal distribution is bit-for-bit identical to
  the delta-gap behaviour. The cross-section eigenmode solver lands
  in Phase 1.3.1; until then a microstrip wave port and a delta gap
  produce the same result, which is intentional.
- **Surface-roughness loss models.** The Hammerstad–Jensen, Groiss,
  and Huray models shipped in `yee-mom::roughness` are out of scope
  for the matrix-fill chapter; they enter only as a frequency-dependent
  scalar multiplier on the assembled $Z$ matrix.
- **MLFMA and ACA.** Fast methods for the dense $O(N^2)$ matrix and
  $O(N^3)$ LU are Phase 4 work. The Phase 1.0 numbers in §1 come from
  unaccelerated dense LU on `faer`.
- **GPU matrix fill.** The Phase 1.5 cuSOLVER work shipped the LU
  step only (`Zgetrf` / `Zgetrs`); matrix fill remains on the CPU.
  When the GPU fill lands, the quadrature scheme in §4 is the part
  that ports first because it has no data-dependent branching;
  the singular path in §6 is comparatively gnarly.

## 9. References

- Rao, S. M., Wilton, D. R., and Glisson, A. W. "Electromagnetic
  Scattering by Surfaces of Arbitrary Shape." *IEEE Trans. Antennas
  Propag.* 30.3 (May 1982), pp. 409–418. — The RWG basis paper.
- Wilton, D. R., Rao, S. M., Glisson, A. W., Schaubert, D. H.,
  Al-Bundak, O. M., and Butler, C. M. "Potential Integrals for
  Uniform and Linear Source Distributions on Polygonal and
  Polyhedral Domains." *IEEE Trans. Antennas Propag.* 32.3 (March
  1984), pp. 276–281. — Closed-form $1/R$ integrals on flat
  triangles; the singularity-subtraction reference for §5.
- Dunavant, D. A. "High Degree Efficient Symmetrical Gaussian
  Quadrature Rules for the Triangle." *Int. J. Numer. Methods Eng.*
  21.6 (1985), pp. 1129–1148. — The symmetric Gauss rules on the
  master triangle used in §4.
- Duffy, M. G. "Quadrature Over a Pyramid or Cube of Integrands with
  a Singularity at a Vertex." *SIAM J. Numer. Anal.* 19.6 (December
  1982), pp. 1260–1262. — The original Duffy-transform paper that
  §6 builds on.
- Mosig, J. R. "Integral Equation Technique." In T. Itoh, ed.,
  *Numerical Techniques for Microwave and Millimeter Wave Passive
  Structures*, Wiley, 1989, ch. 3. — Early MPIE formulation
  reference; the modern form descends from this treatment.
- Burke, G. J. "Numerical Electromagnetics Code — NEC-4: Method of
  Moments, Part I — User's Manual." LLNL UCRL-MA-109338 (January
  1992). — The finite-radius cylindrical-dipole reference value
  $Z \approx 87 + j41\,\Omega$ cited in §1 and used as the
  `mom-001` gate.
