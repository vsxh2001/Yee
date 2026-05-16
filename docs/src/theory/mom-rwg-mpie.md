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
- **Dunavant symmetric Gauss quadrature** for the well-separated
  triangle pairs that dominate the matrix-fill cost (forthcoming §4);
- **singularity-subtraction plus a Duffy-transform path** for the
  same-triangle and shared-vertex pairs where the $1/R$ kernel goes
  singular (forthcoming §§5–6).

The forthcoming sections cover the delta-gap port treatment that turns
the dense solve into a one-port S-parameter, the pieces that are
out of scope for this chapter, and the reference list.

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
  When the GPU fill lands, the quadrature scheme is the part that
  ports first because it has no data-dependent branching; the
  singular path is comparatively gnarly.

## 9. References

- Rao, S. M., Wilton, D. R., and Glisson, A. W. "Electromagnetic
  Scattering by Surfaces of Arbitrary Shape." *IEEE Trans. Antennas
  Propag.* 30.3 (May 1982), pp. 409–418. — The RWG basis paper.
- Wilton, D. R., Rao, S. M., Glisson, A. W., Schaubert, D. H.,
  Al-Bundak, O. M., and Butler, C. M. "Potential Integrals for
  Uniform and Linear Source Distributions on Polygonal and
  Polyhedral Domains." *IEEE Trans. Antennas Propag.* 32.3 (March
  1984), pp. 276–281. — Closed-form $1/R$ integrals on flat
  triangles; the singularity-subtraction reference.
- Dunavant, D. A. "High Degree Efficient Symmetrical Gaussian
  Quadrature Rules for the Triangle." *Int. J. Numer. Methods Eng.*
  21.6 (1985), pp. 1129–1148. — The symmetric Gauss rules on the
  master triangle used in the forthcoming quadrature section.
- Duffy, M. G. "Quadrature Over a Pyramid or Cube of Integrands with
  a Singularity at a Vertex." *SIAM J. Numer. Anal.* 19.6 (December
  1982), pp. 1260–1262. — The original Duffy-transform paper.
- Mosig, J. R. "Integral Equation Technique." In T. Itoh, ed.,
  *Numerical Techniques for Microwave and Millimeter Wave Passive
  Structures*, Wiley, 1989, ch. 3. — Early MPIE formulation
  reference; the modern form descends from this treatment.
- Burke, G. J. "Numerical Electromagnetics Code — NEC-4: Method of
  Moments, Part I — User's Manual." LLNL UCRL-MA-109338 (January
  1992). — The finite-radius cylindrical-dipole reference value
  $Z \approx 87 + j41\,\Omega$ cited in §1 and used as the
  `mom-001` gate.
