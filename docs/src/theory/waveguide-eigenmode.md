# Waveguide Cross-Section Eigenmode Solver — Theory of Operation

This page is the theory-of-operation reference for the 2-D Nedelec
edge-element eigenmode solver shipped in Phase 1.3.1.1 and used by
`yee-mom::ports::NumericalCrossSection` to extract the dominant
propagation constant `β` of an arbitrary waveguide cross-section. It
is the structural analog of [RWG / MPIE](./mom-rwg-mpie.md) — basis
derivation, element-matrix integrals, assembly, numerical gotchas —
for the cross-section vector-Helmholtz eigenproblem. Audience: an
engineer porting the solver, debugging a stray near-zero eigenvalue,
or extending to higher-order Nedelec / sparse Arnoldi. Equations are
plain-text ASCII / LaTeX-style — the mdBook build has no math
preprocessor (see `docs/book.toml`).

## 1. Introduction

Phase 1.3.1.0 shipped a closed-form rectangular-waveguide TE10 wave
port (Pozar §3.3). Phase 1.3.1.1 generalises that path to
**arbitrary** cross-sections — microstrip, GCPW, slot, multi-
conductor — by computing the mode profile numerically on a 2-D
triangular mesh. The validation gate is WR-90 TE10 at 10 GHz on a
6×6 quad-diagonal mesh: analytic `β = 158.238256 rad/m`, numerical
`β = 158.150550 rad/m`, error 0.055%. Hitting that number requires
five pieces working together: a **transverse-only weak formulation**
(§2–§3) derived from source-free Maxwell on a `z`-translation-
invariant waveguide; **first-order Nedelec (Whitney-1) edge
elements** (§4–§5) — the curl-conforming basis with no spurious
modes inside the spectrum of interest; **Dirichlet PEC elimination**
(§6) on boundary edges; a **Cholesky-symmetrised standard
eigenproblem** (§7); and **propagation-constant recovery** (§8) via
`β² = k_0² − k_c²`. Sections 9–10 cover numerical considerations
and out-of-scope work.

## 2. Maxwell on a z-translation-invariant cross-section

Consider a waveguide whose cross-section `Ω` is uniform along `z`,
filled with isotropic lossless media `(ε_r(x, y), μ_r(x, y))`, with
PEC walls `∂Ω`. Time-harmonic source-free Maxwell with the
`exp(+jωt)` convention is `∇ × E = -jωμ₀ μ_r H`,
`∇ × H = jωε₀ ε_r E`, `∇·(ε_r E) = 0`, `∇·(μ_r H) = 0`. The guided-
mode ansatz is `E = [E_t(x,y) + ẑ E_z(x,y)] · exp(j(ωt − βz))` (and
likewise `H`) with `E_t` the in-plane component and `β` the unknown
propagation constant. Substitution splits the curl equations into
transverse and longitudinal blocks (Pozar §3.1). For a
**homogeneously filled** PEC waveguide the dominant TE family has
`E_z = 0`, and the eigenproblem reduces to

```text
∇_t × ( (1/μ_r) ∇_t × E_t )  −  k_0² ε_r  E_t  =  −β² E_t          (★)
```

with `k_0 = ω/c`. Defining the **cutoff wavenumber**
`k_c² := k_0² − β²`, (★) becomes
`∇_t × ((1/μ_r) ∇_t × E_t) = k_c² ε_r E_t`, the form
`eigensolver::assemble_transverse` solves. The mixed `(E_t, E_z)`
formulation (Lee-Sun-Cendes 1991) is required for **inhomogeneous**
fills like microstrip on FR-4; the assembly module already produces
the longitudinal-block matrices (`local_a_zz`, `local_b_zz`,
`local_b_ze`) so it slots in when Phase 1.3.1.2 lands.

## 3. Weak form

Multiply (★) by a test field `E_t'` and integrate over `Ω`. Using
`∇·(A × B) = (∇×A)·B − A·(∇×B)` and the divergence theorem,

```text
∫_Ω (1/μ_r) (∇_t × E_t) · (∇_t × E_t') dA
  −  ∮_∂Ω (1/μ_r) (n̂ × E_t') · (∇_t × E_t) dℓ
   =  k_c² ∫_Ω ε_r  E_t · E_t' dA.
```

The PEC condition `n̂ × E_t = 0` on `∂Ω` is imposed on the trial
space; for a matching test space the boundary integral vanishes.
Define the bilinear forms
`S(E_t, E_t') := ∫_Ω (1/μ_r) (∇_t × E_t)·(∇_t × E_t') dA`
(curl-curl stiffness) and
`T(E_t, E_t') := ∫_Ω ε_r E_t · E_t' dA` (ε_r-weighted mass), so the
weak eigenproblem is `S(E_t, E_t') = k_c² T(E_t, E_t')`. A Galerkin
discretisation on a finite subspace `V_h ⊂ V` produces
`S · e = k_c² · T · e`, `e ∈ ℝ^n`.

## 4. Why Nedelec edge elements

Expanding each Cartesian component of `E_t` in nodal Lagrange
(Whitney-0) basis is **wrong** for the vector Helmholtz operator:
the discrete null-space of `curl` no longer matches the analytic
kernel (gradients of scalar potentials), and the discretisation
injects a forest of **spurious modes** with non-zero `k_c²`
indistinguishable from physical TE / TM modes — historically the
reason FEM vector eigensolves had a bad reputation in the 1980s.
The fix is the **first-order Nedelec (Whitney-1) edge element**
(Nedelec 1980; Whitney 1957):

- **Curl-conforming**, not gradient-conforming: only the
  **tangential** component is continuous across element edges,
  matching the physical interface condition for `E_t` across a
  dielectric discontinuity.
- The discrete `curl` kernel on a Nedelec space is **exactly** the
  gradient of the scalar nodal-Lagrange space (Boffi-Brezzi-
  Demkowicz 2013, §5). Spurious gradient null-modes cluster at
  `k_c² = 0`, where a small threshold filters them.

The price is one DoF per edge (three per triangle for first-order
Nedelec versus six for vector Lagrange); the win is a spurious-mode-
free spectrum and a clean dominant-cutoff extraction without mode
sorting. Jin §9.1–9.3 is the textbook treatment.

## 5. Element-level matrices

On a triangle `T` with CCW vertices `(v_0, v_1, v_2)` and area `A`,
the three barycentric coordinates `λ_i` satisfy `λ_i(v_j) = δ_ij`
and have piecewise-constant gradients

```text
∇λ_i = (b_i, c_i) / (2 A),
b_i  = y_{(i+1) mod 3} − y_{(i+2) mod 3},
c_i  = x_{(i+2) mod 3} − x_{(i+1) mod 3}.
```

The three Nedelec basis functions on `T` are indexed by local edges
`e ∈ {0, 1, 2}` (`eigensolver::mesh` convention: local edge `e`
lies opposite local vertex `e`, traversed CCW):

```text
N_e(r)  =  ℓ_e σ_e ( λ_a ∇λ_b  −  λ_b ∇λ_a ),
```

with `(a, b)` the local-vertex endpoints of edge `e`, `ℓ_e` the
edge length, and `σ_e ∈ {+1, −1}` the per-triangle orientation
against the canonical global-edge direction (`from < to` in
vertex-index order — `eigensolver::mesh::EdgeKey`). The Whitney
form's curl is **constant** on the triangle:

```text
∇_t × N_e  =  ( ℓ_e σ_e / A ) ẑ.
```

### 5.1 Local curl-curl stiffness

For uniform `μ_r` on `T`, the constant curl gives

```text
S^e_{ij}  =  ∫_T (1/μ_r) (∇_t × N_i) · (∇_t × N_j) dA
          =  σ_i σ_j  ℓ_i ℓ_j  /  ( μ_r · A ).
```

No quadrature is needed; `assembly::local_a_ee_curl` writes this
closed form directly.

### 5.2 Local Nedelec mass

The mass integrand `N_i · N_j` is quadratic in the `λ`'s. The
linear-triangle moment identity

```text
∫_T λ_p λ_q dA  =  (A / 12)( 1 + δ_{pq} )
```

reduces every term. For local edges `i, j` with endpoints
`(a_i, b_i)` and `(a_j, b_j)`,

```text
T^e_{ij}  =  ε_r  σ_i σ_j  ℓ_i ℓ_j  · [
                 (A/12)(1 + δ_{a_i a_j}) (∇λ_{b_i} · ∇λ_{b_j})
               − (A/12)(1 + δ_{a_i b_j}) (∇λ_{b_i} · ∇λ_{a_j})
               − (A/12)(1 + δ_{b_i a_j}) (∇λ_{a_i} · ∇λ_{b_j})
               + (A/12)(1 + δ_{b_i b_j}) (∇λ_{a_i} · ∇λ_{a_j})
             ].
```

`assembly::local_b_ee_mass` is the direct transcription. Jin §9.4
is the canonical reference (with the orientation-sign convention
matched to Boffi-Brezzi-Demkowicz §3.5).

### 5.3 Worked example

For the unit right triangle `v_0 = (0,0)`, `v_1 = (1,0)`,
`v_2 = (0,1)`: `A = 1/2`, edge lengths `(√2, 1, 1)`, all
`σ_e = +1`. For `μ_r = 1`: `S^e_{00} = (√2)²/(1·1/2) = 4`,
`S^e_{11} = S^e_{22} = 1²/(1·1/2) = 2`, `S^e_{01} = √2/(1/2) = 2√2`.
This is what
`eigensolver::tests::local_curl_matrix_symmetric_on_unit_triangle`
checks. The mass block is more tedious but equally mechanical; the
Rust code is canonical, the formulas here are for cross-checking on
port.

## 6. Global assembly and Dirichlet PEC BCs

Global assembly is the standard scatter pattern: walk every
triangle, compute local `S^e` and `T^e`, scatter into globals by
edge index. Two pieces of bookkeeping are non-obvious:

**Edge incidence and orientation.** Each global edge is keyed by
the canonical pair `(from, to)` with `from < to`. A triangle's
local edge `(a, b)` matches that orientation with `σ = +1` if its
traversal direction agrees, else `σ = −1`. Without the signed map
the off-diagonal `S^e_{ij}` entries come out with the wrong sign on
every edge where adjacent triangles disagree on local direction —
roughly half the interior edges — and the dominant eigenvalue is
meaningless. `eigensolver::mesh::EdgeTable` owns this bookkeeping.

**PEC walls = boundary-edge elimination.** A boundary edge touches
**exactly one** triangle. PEC `n̂ × E_t = 0` makes the tangential
component vanish, so the Nedelec DoF is identically zero. The
implementation **drops** boundary-edge rows and columns during
scatter: only `is_boundary = false` edges get interior-DoF indices,
so the assembled `S` and `T` are `(n_interior × n_interior)`. This
is static condensation (Jin §4.6) — equivalent to row-and-column
zeroing + deletion but cheaper. The reduced eigenproblem is
`S · e_int = k_c² · T · e_int`. The 6×6 WR-90 mesh in
`tests/eigensolver_wr90.rs` has `n_interior ≈ 84` — inside the
"few hundred DoF" envelope where dense linear algebra is fine.

## 7. Cholesky-symmetrised generalised eigenproblem

The naive approach is `M = T^{-1} S` followed by `eigenvalues(M)`.
**This does not work.** `S` and `T` are both real-symmetric (Gram
matrices of symmetric bilinear forms), but `T^{-1} S` is in general
**not symmetric** — `(T^{-1} S)^T = S T^{-1}`, which equals
`T^{-1} S` only when they commute, which they do not. Feeding the
asymmetric product to `nalgebra::DMatrix::eigenvalues` triggers a
non-symmetric QR / Hessenberg path that hung in Phase 1.3.1.1 step
3 bring-up on matrices of size `≥ 50` (EEEEE gotcha — worked around
rather than reported upstream).

The fix is **Cholesky symmetrisation** (Golub & Van Loan §8.7).
Since `T` is symmetric positive-definite, factor `T = L · L^T` with
`L` lower-triangular. Substitute `e = L^{-T} u` into
`S e = k_c² T e` and left-multiply by `L^{-1}`:

```text
L^{-1} · S · L^{-T} · u  =  k_c² · u,    i.e.  M · u = k_c² u,
M  :=  L^{-1} · S · L^{-T}.
```

`M` is **symmetric** (`S` symmetric and the two triangular factors
are transposes), so the symmetric-tridiagonal QR path
(`nalgebra::SymmetricEigen`) applies and is bulletproof. The
eigenvector back-substitution is `e = L^{-T} u`; §8 only needs the
eigenvalue. `eigensolver::solve::solve_dense` performs the
symmetrisation via two triangular solves (`L · Y = S`, then
`L · Z = Y^T`, finally `M = Z^T`) and an explicit
`M ← 0.5 (M + M^T)` to suppress floating-point asymmetry drift.

## 8. Propagation constant recovery

Given the smallest strictly-positive `k_c²` from §7, the
propagation constant at angular frequency `ω` is
`β(ω) = sqrt(k_0² − k_c²)`, `k_0 = ω / c`. For `k_0² > k_c²` the
mode propagates (`β ∈ ℝ`); for `k_0² < k_c²` it is **evanescent**
(pure imaginary `β`) and `solve_dense` surfaces `Error::Numerical`
rather than returning a complex `β` — the lossless transverse
formulation has no business returning one. The cutoff frequency is
`f_c = c · sqrt(k_c²) / (2π)`.

**Spurious-mode filtering.** Gradient null-modes (§4) cluster at
`k_c² ≈ 0`. `solve_dense` rejects any
`k_c² ≤ 1e-6 · max(eigenvalues)` and returns the smallest
strictly-positive residue. The `1e-6` factor is conservative
(gradient null-modes typically land at `< 1e-10` of the dominant)
and is empirically chosen on the WR-90 mesh.

### Worked WR-90 number

WR-90 has `a = 22.86 mm`, `b = 10.16 mm`, air-filled. The analytic
TE10 cutoff is `k_c,TE10 = π / a`, `f_c = c / (2 a) ≈ 6.5570 GHz`.
At `f = 10 GHz` (`k_0 ≈ 209.585 rad/m`),
`β_analytic = sqrt(k_0² − (π/a)²) ≈ 158.238256 rad/m`. The 6×6
quad-diagonal mesh in `tests/eigensolver_wr90.rs` yields
`β_numerical ≈ 158.150550 rad/m`, error `|Δβ/β| ≈ 5.5 × 10⁻⁴`.
This 0.055% error on `n_interior ≈ 84` DoFs is the shipped accuracy
floor; refining the mesh drives error down at the expected first-
order Nedelec rate (§9). `eigensolver_wr90_te10.ipynb` (Track LLLLL)
reproduces this from Python.

## 9. Numerical considerations

**Convergence rate.** First-order Nedelec delivers `O(h)` energy-
norm convergence for smooth solutions on a quasi-uniform mesh,
translating to `O(h²)` on `k_c²` (eigenvalues converge twice as
fast as eigenvectors — Babuska-Osborn 1991). Refining WR-90 from
`n × n = 6` to `12` should drop the error by roughly `4×`;
`eigensolver_wr90_te10.ipynb` sweeps this explicitly.

**Mesh quality.** Mass-matrix conditioning `cond(T)` scales with
worst-case triangle aspect ratio. For Cholesky to succeed without
spurious indefiniteness, keep minimum angles `≥ 20°`. Gmsh's
default Delaunay output is well above this floor; hand-rolled grids
of very skewed triangles can trip the Cholesky path.

**Real-arithmetic restriction.** `solve_dense` operates on the real
parts of `S` and `T` and surfaces `Error::Unimplemented` if the
imaginary parts exceed `1e-9` of the real norm. The assembly module
stores `DMatrix<Complex64>` to keep the API future-proof for lossy
fills (Phase 1.3.1.2), but v0 mode extraction is lossless-only.
Loss enters via complex `ε_r`, `μ_r`, and complex `β` — for which
`SymmetricEigen` no longer applies and a complex-symmetric (LDLᵀ-
style) eigensolve is required.

## 10. Limitations and roadmap

- **Higher-order Nedelec.** Phase 1.3.1.x will land second-order
  edge elements (Whitney-1 plus quadratic bubble), bringing
  convergence to `O(h²)` energy / `O(h⁴)` eigenvalue. The curl is
  no longer constant per triangle, so embedded Gauss quadrature
  per element becomes load-bearing.
- **Sparse shift-and-invert eigensolve.** Dense `SymmetricEigen`
  is `O(n³)` and viable only up to a few hundred DoF. Phase
  1.3.1.1 step 4 wires in a sparse Arnoldi / LOBPCG path
  (`arpack-rs` was escape-hatched; current candidate is hand-rolled
  LOBPCG against `faer`-sparse). See spec
  `2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`.
- **Mixed `(E_t, E_z)` formulation for inhomogeneous fills.** The
  longitudinal-block element matrices are implemented and unit-
  tested; the eigensolve wires in transverse-only. Phase 1.3.1.1
  step 5 turns on the mixed formulation, required for quasi-TEM
  microstrip on a layered substrate where `E_z = 0` breaks down.
- **Anisotropic / dispersive media.** Tensor `ε_r`, `μ_r` are Phase
  1.3.3+. Dispersive (Drude/Lorentz/Debye) materials live in
  `yee-fdtd` and are not wired into the MoM port path.
- **Multi-mode extraction and tracking.** Phase 1.3.2 / 1.3.4.

## 11. References

- Pozar, D. M. *Microwave Engineering*, 4th ed. Wiley, 2011. Ch. 3
  — TE / TM mode decomposition; closed-form TE10 (§8 gate).
- Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd
  ed. Wiley-IEEE, 2014. Ch. 8 (waveguide eigenproblem), Ch. 9
  (Nedelec / Whitney edge elements). Element-matrix formulas in §5
  follow this book's §9.4.
- Boffi, D., Brezzi, F., and Demkowicz, L. F. *Mixed Finite Element
  Methods and Applications.* Springer, 2013. Ch. 5 — modern de Rham
  / cochain treatment of why curl-conforming edge elements admit no
  spurious modes inside the spectrum.
- Nedelec, J.-C. "Mixed finite elements in ℝ³." *Numer. Math.* 35.3
  (1980), pp. 315–341. — Original curl-conforming elements on
  tetrahedra and triangles.
- Whitney, H. *Geometric Integration Theory.* Princeton University
  Press, 1957. — Simplicial-complex precursor of Nedelec's elements.
- Lee, J.-F., Sun, D.-K., and Cendes, Z. J. "Full-wave analysis of
  dielectric waveguides using tangential vector finite elements."
  *IEEE Trans. Microwave Theory Tech.* 39.8 (1991), pp. 1262–1271.
  — Mixed `(E_t, E_z)` formulation Phase 1.3.1.1 step 5 will switch
  on.
- Golub, G. H., and Van Loan, C. F. *Matrix Computations*, 4th ed.
  Johns Hopkins University Press, 2013. Ch. 8 — the Cholesky-
  symmetrised generalised eigenproblem treatment §7 follows.
- Babuska, I., and Osborn, J. "Eigenvalue Problems." In *Handbook of
  Numerical Analysis*, vol. II, Elsevier, 1991, pp. 641–787. — The
  "eigenvalues converge twice as fast as eigenvectors" claim (§9).
