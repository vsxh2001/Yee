# FEM Eigenmode Solver (3-D Nedelec) — Theory of Operation

This page is the theory-of-operation reference for the 3-D tetrahedral
Nedelec edge-element eigenmode solver planned for Phase 4.fem.eig.0 of
`yee-fem`. It is the structural analog of the 2-D sibling
[Waveguide Cross-Section Eigenmode Solver](./waveguide-eigenmode.md) —
basis derivation, element-matrix integrals, signed assembly, numerical
gotchas — for a closed 3-D cavity rather than a translationally
invariant cross-section. Audience: an engineer implementing Phase
4.fem.eig.0 (`crates/yee-fem/src/{assembly,solve}.rs`) or debugging a
stray near-zero eigenvalue. Equations are plain-text ASCII — no math
preprocessor (see `docs/book.toml`). Companion spec:
`docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`.

## 1. Introduction

Phase 1.3.1.1 shipped a 2-D Nedelec eigensolver for guided waveguide
cross-sections (`E ∝ E_t(x,y)·e^{-jβz}`). Phase 4.fem.eig.0 is the 3-D
analog: given a closed PEC cavity discretised by tetrahedra, find the
resonant wavenumbers `k_n = ω_n / c` at which standing fields exist with
no excitation. The use cases — combline / iris-coupled filter design,
Q-factor for lossy cavities, DRA modal analysis, accelerator RF cells —
are the slice FDTD cannot reach efficiently (lossless ringdown never
decays) and MoM cannot reach at all (no volumetric interior). The fix
is a volumetric FEM eigensolve in `H(curl; Ω)`.

The v0 walking-skeleton scope: a tetrahedral mesh from
`yee-mesh::TetMesh3D`, first-order Nedelec (Whitney-1) edges, PEC
tangential-`E`-zero Dirichlet, sparse `K·e = k²·M·e` solved with
shift-invert to skip the gradient-kernel cluster at `k² ≈ 0`. Gate:
Pozar §6.3 cavity at `a = 22.86 mm`, `b = 10.16 mm`, `d = 30 mm`,
lowest TE_{101} at `≈ 9.660 GHz` within ±0.3% (§9, fem-eig-001).
§2–§8 derive the formulation; §9 covers the gate; §10–§11 list
numerical considerations and deferred sub-phases.

## 2. Maxwell on a closed 3-D cavity

A closed lossless `Ω ⊂ ℝ³` with PEC boundary `∂Ω` filled with isotropic
`(ε_r(x), μ_r(x))`. Time-harmonic source-free Maxwell with the
`exp(+jωt)` convention reads `∇×E = -jωμ₀μ_r H`, `∇×H = jωε₀ε_r E`.
Eliminating `H` gives the vector wave equation

```text
∇ × ( (1/μ_r) ∇ × E )  =  k₀² ε_r E      in Ω,         (★)
n̂ × E                  =  0                on ∂Ω,
```

with `k₀² = ω² ε₀ μ₀ = (ω/c)²`. Equation (★) is an **eigenproblem** in
`(k₀², E)`: `k_n` are the resonant frequencies up to
`f_n = c·k_n / (2π)`, and `E_n` is the modal field. Multiplying (★) by
`v ∈ H₀(curl; Ω)` and integrating by parts (using
`∫_Ω ∇·(A × B) dV = ∮_∂Ω (A × B)·n̂ dS` and `n̂ × v = 0`) gives the
weak form: find `(k₀², E) ∈ (ℝ₊, H₀(curl; Ω))` such that for all
`v ∈ H₀(curl; Ω)`,

```text
∫_Ω (1/μ_r) (∇×E) · (∇×v) dV  =  k₀² ∫_Ω ε_r E · v dV.
```

Define `K(E, v) := ∫_Ω (1/μ_r) (∇×E)·(∇×v) dV` (curl-curl stiffness)
and `M(E, v) := ∫_Ω ε_r E · v dV` (vector mass). Galerkin on
`V_h ⊂ H₀(curl; Ω)` yields the matrix eigenproblem
`K · e = k₀² · M · e`. Jin *FEM in EM* (3rd ed., Wiley 2014) Ch. 9–10
is the textbook treatment.

## 3. Why curl-conforming Nedelec edge elements

The 2-D argument from [waveguide-eigenmode.md](./waveguide-eigenmode.md)
§4 lifts directly, and the consequences are sharper. Nodal `H¹`
Lagrange basis applied to (★) component-wise admits a forest of
**spurious modes**: the discrete kernel of `∇×` is larger than the
continuous one (which is `∇ H¹`), and the extra null-space leaks into
the spectrum as non-physical resonances. In 3-D the spurious population
grows roughly as the cube of mesh refinement — historically the reason
FEM vector eigensolves had a bad reputation in the 1980s.

The fix is the **first-order Nedelec (Whitney-1) edge element** (Nedelec
1980; Whitney 1957): one DoF per mesh edge, tangentially continuous
across faces, curl-conforming. The de Rham complex

```text
H¹(Ω)  --∇-->  H(curl; Ω)  --∇×-->  H(div; Ω)  --∇·-->  L²(Ω)
```

places Nedelec at the `H(curl)` slot. Boffi-Brezzi-Demkowicz 2013 Ch. 7
"discrete compactness" guarantees the discrete `∇×` kernel on the
Nedelec space is **exactly** `∇ V_h^Lagrange` — discrete spurious modes
coincide with continuous gauge freedom and cluster harmlessly at
`k₀² = 0`. Shift-invert at `σ > 0` (§8) does not see them.

## 4. Tetrahedral element with 6 edge DoFs

On a tetrahedron `T` with vertices `(v_0, v_1, v_2, v_3)` and signed
volume `V_T = (1/6)(v_1−v_0)·((v_2−v_0)×(v_3−v_0)) > 0`, the four
barycentric `P_1` coordinates `λ_i(x)` satisfy `λ_i(v_j) = δ_{ij}` and
have **piecewise-constant** gradients
`∇λ_i = (1/(6 V_T))·((v_{i+2}−v_{i+1})×(v_{i+3}−v_{i+1}))` (mod 4, sign
flip per inverted permutation; Jin §9.4). The six **Whitney-1 edge
basis functions** are indexed by the unordered vertex pairs `{(0,1),
(0,2), (0,3), (1,2), (1,3), (2,3)}`; for edge `i → j` (canonical
`i < j`),

```text
W_{ij}(x)  =  λ_i(x) ∇λ_j  −  λ_j(x) ∇λ_i,
∇ × W_{ij} =  2 (∇λ_i × ∇λ_j).
```

The curl is **constant** on the tet (since `∇λ_j` is constant, the
curl of `λ_i ∇λ_j` is `∇λ_i × ∇λ_j`). Two consequences: the curl-curl
integral is **exact** without quadrature, and `W_{ij}` is **tangentially
continuous** across faces (adjacent tets see the same edge with the
same tangent up to the orientation sign in §5) with unit tangential
circulation along its own edge — the geometric statement placing it
at the `H(curl)` slot of the de Rham complex.

### 4.1 Local curl-curl stiffness and mass

Constant curl + scalar `μ_r,T` gives the closed-form 6×6 stiffness

```text
K^e_{αβ}  =  (4 V_T / μ_r,T)
             · (∇λ_{i_α} × ∇λ_{j_α}) · (∇λ_{i_β} × ∇λ_{j_β}),
```

with `α = (i_α, j_α)` the unordered pair indexing local edge `α`
(Bossavit §5; Jin §9.4). The mass integrand `W_α · W_β` is quadratic
in the `λ`'s, so 4-point Gauss quadrature on the reference tet (exact
through degree 2) integrates `M^e_{αβ} = ε_r,T ∫_T W_α · W_β dV`
exactly. The closed-form alternative invokes the simplex moment
identity
`∫_T λ_p^a λ_q^b λ_r^c λ_s^d dV = V_T · a!b!c!d! · 3! / (a+b+c+d+3)!`
and reduces every term to dot products of `∇λ`'s; Jin §9.4 tabulates
the resulting 6×6 mass with the orientation-sign convention matched to
Boffi-Brezzi §3.5. The 2-D analog is `assembly::local_a_ee_curl` in
`crates/yee-mom/src/eigensolver/assembly.rs`; the planned 3-D version
is `crates/yee-fem/src/assembly.rs::local_k_curl` (Phase 4 plan T3).

### 4.2 Worked unit-tetrahedron example

For the unit tet `v_0=(0,0,0)`, `v_1=(1,0,0)`, `v_2=(0,1,0)`,
`v_3=(0,0,1)`: `V_T = 1/6`, `∇λ_0 = (-1,-1,-1)`, `∇λ_j` the standard
basis vector for `j ≥ 1`. The `(0,j)` cross products give
`K^e_{(0,j),(0,j)} = (4/6)·2 = 4/3`. The 3-D counterpart of the unit-
right-triangle example in
[waveguide-eigenmode.md](./waveguide-eigenmode.md) §5.3; the planned
`fem::tests::local_curl_matrix_unit_tet` exercises it.

## 5. Edge orientation and signed assembly

Each global edge is keyed by its endpoint pair in canonical
**lower-vertex-index → higher-vertex-index** order. A tet's local edge
`α = (i, j)` matches the global orientation with `σ_α = +1` if
`vertex_id(i) < vertex_id(j)`, else `σ_α = −1`. The local 6×6 block
scatters with a diagonal sign matrix `D_T = diag(σ_0, …, σ_5)` on both
sides:

```text
K_global += D_T · K^e · D_T,
M_global += D_T · M^e · D_T,
```

i.e. per-pair off-diagonal `K^e_{αβ}` scatters with sign `σ_α σ_β`.
This is the **single most common cancellation bug** when porting an
FEM edge-element solver: off-diagonal signs on edges where adjacent
tets disagree on local direction come out wrong, and the dominant
eigenvalue is unreasonably small or negative. The 2-D solver has the
same gotcha; see `crates/yee-mom/src/eigensolver/mesh.rs::EdgeKey` and
the sign threading in `local_a_ee_curl`. The 3-D `EdgeKey` parallels
it with an extra index.

## 6. PEC boundary conditions

PEC walls impose `n̂ × E = 0` on `∂Ω`. For a Nedelec basis this means
**the tangential circulation of E along any boundary edge vanishes** —
the DoF on a boundary edge is identically zero. Concretely: enumerate
`boundary_edges()` (edges belonging to a boundary face) and **drop
their rows and columns** from `K` and `M` before the eigensolve. The
reduced system is `(n_interior × n_interior)` where
`n_interior = n_edges − n_boundary_edges`.

The 3-D counterpart of the 2-D Dirichlet elimination in
[waveguide-eigenmode.md](./waveguide-eigenmode.md) §6, and
**fundamentally different** from nodal Lagrange Dirichlet (which
eliminates boundary **vertices**, not **edges**). The penalty /
row-replacement alternative avoids resizing the DoF map at the cost of
one large eigenvalue per Dirichlet edge that shift-invert trivially
skips; v0 uses row elimination (spec §6 bullet 4).

## 7. The generalised sparse eigenproblem

After assembly and PEC elimination the discrete problem is
`K · e = k² · M · e`, `e ∈ ℝ^{n_interior}`, with `M` symmetric positive
definite (Gram matrix of the Nedelec basis; v0 requires `ε_r > 0` on
every tet) and `K` symmetric positive **semi-definite**. The
non-trivial kernel of `K` comes from gradients of nodal Lagrange
functions on the same mesh:

```text
ker(K)  ⊇  ∇ V_h^{Lagrange, 0}        (homogeneous on ∂Ω),
dim ker(K)  =  n_vertices_interior − 1     (simply connected Ω).
```

The `−1` accounts for the single global constant `φ ≡ const`; on a
multiply-connected domain (coaxial, torus) harmonic forms add
`dim H¹(Ω; ℝ)` ranks. v0 targets simply connected cavities.

These zero eigenvalues are the **spurious gradient modes** of §3 in
discrete form — not a numerical artifact: they are the exact Nedelec
discretisation of the continuous gauge freedom and an honest
eigensolver will return them. §8's job is to **skip them** at solve
time, not detect-and-discard.

## 8. Shift-invert Arnoldi / LOBPCG

The physical modes of interest are the lowest-`k²` eigenpairs **above**
the gradient kernel at `k² = 0`. A naive Arnoldi run on `(K, M)`
converges to extremal eigenvalues — the zero cluster on the small-`k²`
end — and is useless. The fix is **shift-invert**: pick a positive
`σ ∈ (0, k²_{smallest physical})` and target

```text
(K − σ M)^{−1} M · e  =  θ · e,     k²  =  σ + 1/θ.
```

The transformed operator has its largest-magnitude eigenvalues at the
original `k²` nearest the shift, where Arnoldi or LOBPCG converges fast
and zero modes (mapped to `θ = −1/σ < 0`) sit at the opposite end. The
primitive is **one sparse LU of `(K − σ M)`**, reused across every
inner-product evaluation.

Choice of `σ`. Too small and the zero cluster lands close to the
physical eigenvalues in transformed magnitude and contaminates the
converged subspace. Too large and `(K − σ M)` is numerically singular;
LU pivoting degrades. Spec §6 recommends `σ = 0.5 · k²_{lowest
expected}` from a cheap analytic estimate (for fem-eig-001,
`(π/a)² + (π/d)²`).

Library choice (spec §8). The 2-D solver escape-hatched to dense
`SymmetricEigen` at ≤ 500 DoFs; dense is infeasible in 3-D where
`(K − σ M)` is sparse with `n ∼ 10⁴`–`10⁶`. Phase 4.fem.eig.0 ships a
`SparseEigen` trait with a default LOBPCG backend (Knyazev 2001,
pure-Rust `lobpcg` crate) plus a `faer` sparse LU preconditioner; an
ARPACK binding is gated behind `eig-arpack`. Trait is load-bearing,
library is the swap point — same pattern as `yee_cuda::backend`.

## 9. Validation gate — fem-eig-001 rectangular cavity

The Phase 4.fem.eig.0 production gate is the canonical Pozar §6.3
example: a lossless air-filled rectangular PEC cavity with broad wall
`a = 22.86 mm`, narrow wall `b = 10.16 mm`, length `d = 30 mm` (WR-90
cross-section closed off 30 mm apart). The analytic TE_{mnp} / TM_{mnp}
resonant frequencies are

```text
f_{mnp}  =  (c / (2 · sqrt(μ_r · ε_r)))
            · sqrt( (m/a)² + (n/b)² + (p/d)² ),
```

with integer indices, not all three zero. For air-filled and TE_{101},
`f_{101} = (c/2)·sqrt((1/0.02286)² + (1/0.030)²) ≈ 9.660 GHz`.

**Gate criteria** (spec §9): (1) lowest extracted resonance satisfies
`|f_FEM − 9.660 GHz| / 9.660 GHz ≤ 0.3%`; (2) lowest ten resonances
sorted ascending agree pairwise with the analytic Pozar table within
±1% RMS with ordering matching; (3) no spurious mode appears below the
analytic TE_{101}; (4) informational, end-to-end runs in `< 60 s` on a
~30k-edge mesh `--release`.

The ±0.3% tolerance is tighter than the 2-D TE10 gate's 0.1% because
3-D discretisation error per DoF is larger; it matches Pozar's
4-significant-digit tabulation. The hand-rolled fixture decomposes the
cavity into ≤ 6 tets per cube cell on a regular grid; Gmsh `tetgen` is
the production path. Further cases (fem-eig-002 lossy Q, fem-eig-003
DRA) are deferred to 4.fem.eig.2 / 4.fem.eig.3 per §11.

## 10. Numerical considerations

**Mesh quality.** Cavity-mode accuracy is strongly sensitive to tet
aspect ratio. Sliver / needle / cap / wedge tets inflate `cond(K)` and
degrade the inner LU. Gmsh's default `tetgen` with
`MeshSizeFactor ≤ 0.05·λ` is acceptable for fem-eig-001; cube-to-six-
tets fixtures with aspect ratio near 1 are fine. Sliver tets (worst
dihedral angle `< 10°`) can produce 5%+ error on the lowest mode even
with correct assembly — document in
`crates/yee-fem/validation/README.md`.

**Memory at production scale.** A 1 M-edge mesh has ~7 nnz per row in
`K` and `M` — `~14·10⁶` non-zeros per matrix, ~0.5 GB CSR FP64. Sparse
LU of `(K − σ M)` is the dominant memory user; AMG preconditioning
would scale better but is out of v0 scope.

**Sparse format.** CSR for assembled matrices, COO as intermediate
scatter target. COO → CSR is one sort plus one prefix scan;
`faer::sparse` provides both.

**Shift-invert convergence rate.** LOBPCG convergence is controlled by
the spectral gap between targeted eigenvalues and the rest of the
post-shift spectrum. With `σ` below the lowest physical mode, the gap
to the gradient cluster (opposite end after `1/θ` mapping) is large
and ten modes typically converge in `< 50` iterations. Poor `σ`
manifests as 10× more iterations or non-convergence; diagnose by
logging per-iteration residual and re-running with `σ` halved.

## 11. Limitations and roadmap

First-order Nedelec on tets converges as `O(h²)` for eigenvalues under
quasi-uniform refinement (Babuska-Osborn 1991), `O(h)` energy-norm for
eigenvectors. Deferred sub-phases (spec §13):

- **Higher-order Nedelec** (Phase 4.fem.eig.1). Hierarchical `p ≤ 2`
  takes the eigenvalue rate to `O(h⁴)`. The curl is no longer constant
  per tet, so 4-point Gauss quadrature on the curl-curl integrand
  becomes load-bearing.
- **Lossy / complex `ε_r` Q-factor extraction** (Phase 4.fem.eig.2).
  Complex `(K, M)` with the same sparse Arnoldi machinery on complex
  `(K − σ M)`; output complex `k_n²` with `Q = Re(k²) / (2 Im(k²))`.
  fem-eig-002 (Pozar wall-loss Q, ±5%).
- **Dielectric-resonator antenna** (Phase 4.fem.eig.3). Piecewise `ε_r`
  across air / puck. fem-eig-003 (Petosa DRA Handbook Ch. 3, ±2%).
- **Open-region cavities** (Phase 4.fem.eig.x). Absorbing boundary or
  FEM-BEM hybrid; v0 is closed-cavity only.
- **GPU acceleration.** Not in v0; `SparseEigen` is the swap point.
- **Periodic / Floquet BCs.** Deferred to 4.fem.eig.4+.

## 12. References

- Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd ed.
  Wiley-IEEE, 2014. Ch. 9 (Nedelec tets, 6×6 mass / stiffness §9.4),
  Ch. 10 (cavity eigenproblems). §4 / §8 follow this book.
- Boffi, D., Brezzi, F., and Demkowicz, L. F. *Mixed Finite Element
  Methods and Applications.* Springer, 2013. Ch. 7 — curl-conforming
  spaces, discrete compactness; grounds §3.
- Nedelec, J.-C. "Mixed finite elements in ℝ³." *Numer. Math.* 35.3
  (1980), pp. 315–341.
- Whitney, H. *Geometric Integration Theory.* Princeton, 1957.
- Pozar, D. M. *Microwave Engineering*, 4th ed. Wiley, 2012. §6.3 —
  rectangular cavity; closed-form TE / TM the gate matches.
- Bossavit, A. *Computational Electromagnetism.* Academic, 1998. Ch. 5.
- Saad, Y. *Numerical Methods for Large Eigenvalue Problems*, 2nd ed.
  SIAM, 2011. Ch. 4–6 — shift-invert Arnoldi (§8).
- Knyazev, A. V. "LOBPCG." *SIAM J. Sci. Comput.* 23.2 (2001),
  pp. 517–541.
- Babuska, I., and Osborn, J. "Eigenvalue Problems." *Handbook of
  Numerical Analysis* II, Elsevier, 1991, pp. 641–787. — `O(h²)`
  convergence of first-order Nedelec.
- Demkowicz, L. *Computing with hp-Adaptive Finite Elements*, vol. 2.
  CRC, 2007. Ch. 4 — de Rham complex, hierarchical Nedelec.
- Companion spec:
  `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md`;
  2-D sibling: `docs/src/theory/waveguide-eigenmode.md`.
