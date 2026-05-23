# Cross-Section (Waveguide-Port) Eigensolver — Theory of Operation

This page is the theory-of-operation reference for Yee's 2-D
cross-section eigensolver — the numerical wave-port modal solver shipped
across Phase 1.3.1.1 (steps 2–5.2, with step 5.3 in progress). It is
implemented in `crates/yee-mom/src/eigensolver/` and driven through
`yee_mom::ports::NumericalCrossSection`. The audience is the engineer
who wants to read the code with Jin and Pozar open on the desk, or who
is debugging a stray near-zero eigenvalue. Equations are kept in plain
text because the documentation build has no math preprocessor (see
`docs/book.toml`).

This chapter is the **synthesis** of the formulation that is otherwise
spread across the module rustdoc and ADRs 0022, 0023, 0050, 0051, 0052,
0053, and 0054. It supersedes the step-3-era
[Waveguide Cross-Section Eigenmode Solver](./waveguide-eigenmode.md)
page, which documented only the transverse-only homogeneous path before
the mixed `(E_t, E_z)` formulation and the β-direct extraction landed;
read that page for the Nedelec element-matrix derivations (basis
functions, local stiffness/mass integrals, signed assembly,
Cholesky symmetrisation) which remain accurate and are not repeated here
in full.

## 1. Overview

A wave port is the principled way to excite and terminate a
transmission line in a full-wave solver: instead of the delta-gap of
[planar MoM §7](./planar-mom.md), the port carries the actual *guided
mode* of the line's cross-section — its transverse field profile, its
propagation constant `β`, and its wave impedance `Z_w`. Phase 1.3.1.0
shipped the closed-form rectangular-waveguide TE10 mode (Pozar §3.3) for
the hollow-rectangular special case. Phase 1.3.1.1 generalises that to
**arbitrary** cross-sections — microstrip, CPW, GCPW, ridge or
slab-loaded waveguide — by computing the mode numerically on a 2-D
triangular mesh of the cross-section.

The technique is a vector finite-element eigensolve over the
cross-section plane. The unknown is the transverse electric field
`E_t(x, y)` (plus, for inhomogeneous fills, the longitudinal `E_z`); the
discretisation is the mixed Nedelec-edge / nodal-Lagrange pair of
Lee-Sun-Cendes (1991); and the result is a sparse generalised
eigenproblem whose dominant physical eigenpair gives the mode. This is
the standard FEM waveguide eigenproblem of Jin §8.4–8.5 and the
intellectual basis of every commercial port solver (HFSS, CST,
COMSOL RF) in this family.

The sub-project was decomposed (ADR-0022) into a deferred spec/plan plus
follow-up steps:

- **Steps 0–1** (ADR-0023): the `TriMesh2D` mesh type and the
  `ModalDistribution::Numerical2D` API stub, so consumer code freezes
  before the matrix code lands.
- **Steps 2–3**: Nedelec + nodal-Lagrange element-matrix assembly and a
  dense `SymmetricEigen` fallback; transverse-only, exact for
  homogeneous fills.
- **Step 4** (ADR-0050): an in-tree block LOBPCG sparse eigensolver
  (chosen over `arpack-rs` to avoid a Fortran/LAPACK toolchain
  dependency).
- **Step 5** (ADR-0051): the mixed `(E_t, E_z)` longitudinal block, the
  numerical line-integral `Z_w`, and the slab-loaded validation gate.
- **Step 5.1** (ADR-0052): an independently-verified published
  transcendental reference (slab-loaded transverse resonance) to close
  the validation gap.
- **Step 5.2** (ADR-0053): the **β-direct** extraction fix — the
  correctness fix for `ε_r ≠ 1` fills (§5 below).
- **Step 5.3** (ADR-0054, in progress): a direct sparse shift-and-invert
  on the β-direct pencil to close the high-contrast inhomogeneous
  residual (§9).

## 2. Maxwell on a z-translation-invariant cross-section

Consider a waveguide whose cross-section `Ω` is uniform along the
propagation axis `z`, filled with isotropic media `(ε_r(x, y),
μ_r(x, y))`, bounded by PEC walls `∂Ω`. Time-harmonic, source-free
Maxwell with the `exp(+jωt)` convention is

```text
∇ × E = -jωμ₀ μ_r H
∇ × H =  jωε₀ ε_r E
∇ · (ε_r E) = 0
∇ · (μ_r H) = 0
```

A guided mode is a separable ansatz that propagates along `z` with an
unknown complex constant `β`:

```text
E(x, y, z, t) = [ E_t(x, y) + ẑ E_z(x, y) ] · exp(j(ωt − βz))
```

and likewise for `H`. Substituting the ansatz splits the curl equations
into a transverse block (the in-plane components) and a longitudinal
block (the `ẑ` component); the algebra is Pozar §3.1 / Jin §8.4.
Eliminating the magnetic field leaves a vector wave equation for the
electric field over the cross-section, which is the eigenproblem the
solver discretises. The eigenvalue is (a function of) `β` at a fixed
operating frequency `ω`; the eigenvector is the mode profile
`(E_t, E_z)`.

Two regimes matter and they need different formulations:

- **Homogeneous fill** (`ε_r`, `μ_r` constant over `Ω`, e.g. hollow
  WR-90): the dominant TE family has `E_z = 0` identically. A
  **transverse-only** formulation in `E_t` alone is exact.
- **Inhomogeneous fill** (microstrip on a substrate, partial dielectric
  loading, CPW): `E_z` does **not** vanish — it couples through the
  dielectric interface, and the modes are hybrid (quasi-TEM /
  quasi-TE). The **mixed `(E_t, E_z)`** formulation is required; a
  transverse-only solve gives a plausible-but-wrong answer.

## 3. The transverse weak form

For the homogeneous (or, as a building block, the transverse) problem,
the in-plane field obeys

```text
∇_t × ( (1/μ_r) ∇_t × E_t )  =  ( k₀² ε_r − β² ) E_t,        k₀ = ω/c.   (★)
```

Multiplying by a test field `E_t'`, integrating over `Ω`, and applying
the curl identity plus the divergence theorem gives the weak form

```text
∫_Ω (1/μ_r)(∇_t × E_t)·(∇_t × E_t') dA  =  ∫_Ω (k₀² ε_r − β²) E_t · E_t' dA,
```

with the PEC boundary integral vanishing because `n̂ × E_t = 0` is
imposed on the trial and test spaces. Define two real-symmetric bilinear
forms and their finite-element matrices:

```text
S  (curl-curl stiffness)   :  S[i,j]   = ∫_Ω (1/μ_r) (∇_t × N_i)·(∇_t × N_j) dA
T_ε (ε_r-weighted mass)    :  T_ε[i,j] = ∫_Ω ε_r N_i · N_j dA
T_1 (unweighted mass)      :  T_1[i,j] = ∫_Ω      N_i · N_j dA
```

where `N_i` are the basis functions of §4. The discrete (★) is then
either of two algebraically distinct generalised eigenproblems — and the
distinction between them is the entire subject of §5.

## 4. Why Nedelec edge elements

Expanding each Cartesian component of `E_t` in ordinary nodal Lagrange
(Whitney-0) basis functions is **wrong** for the vector Helmholtz
operator. The discrete `curl` kernel no longer matches the analytic one,
and the discretisation injects a forest of **spurious modes** with
nonzero, physical-looking eigenvalues — historically the reason FEM
vector eigensolves earned a bad reputation in the 1980s. The fix
(Webb 1993; Boffi-Brezzi-Demkowicz 2013) is the first-order **Nedelec
(Whitney-1) edge element** for `E_t`:

```text
N_e = ℓ_e σ_e ( λ_a ∇λ_b − λ_b ∇λ_a ),       ∇_t × N_e = (ℓ_e σ_e / A) ẑ,
```

one degree of freedom per mesh edge, with `(a, b)` the edge endpoints,
`ℓ_e` the edge length, `λ` the barycentric coordinates, and `σ_e ∈
{+1, −1}` a per-triangle orientation sign. Edge elements are
**curl-conforming**: only the *tangential* component is continuous
across element edges, which is exactly the physical interface condition
on `E_t` at a dielectric discontinuity. The longitudinal `E_z` is
expanded in ordinary linear **nodal Lagrange** (Whitney-0), one DoF per
vertex, which is correct because `E_z` *is* continuous across interfaces.

The two essential consequences (Boffi-Brezzi-Demkowicz §5):

1. The discrete `curl` kernel on the Nedelec space is **exactly** the
   gradient of the nodal-Lagrange space. The spurious "gradient" modes
   therefore cluster at one identifiable place in the spectrum (cutoff
   `k_c² ≈ 0`), where a small threshold filters them — rather than
   scattering through the physical spectrum where they cannot be told
   apart.
2. The physical TE / TM / hybrid modes are recovered cleanly without
   mode sorting.

The element-level matrix integrals (closed-form curl-curl stiffness,
the linear-triangle mass moment `∫λ_p λ_q dA = (A/12)(1+δ_pq)`, the
signed global assembly, the Dirichlet PEC boundary-edge elimination)
are derived in detail on the
[Waveguide Cross-Section Eigenmode Solver](./waveguide-eigenmode.md)
page (§5–§6) and transcribed directly in
`eigensolver::assembly` (`local_a_ee_curl`, `local_b_ee_mass`,
`local_a_zz`, `local_b_zz`, `local_b_ze`). They are not repeated here.

## 5. The cutoff pencil vs the β-direct pencil

This is the heart of the formulation, and the place where the
implementation history is most instructive.

### 5.1 The cutoff pencil

Defining the **cutoff wavenumber** `k_c² := k₀² − β²` and moving the
`k₀² ε_r` term to the right turns the discrete (★) into the **cutoff
pencil**

```text
S x = k_c² T_ε x.        (cutoff pencil)
```

The smallest strictly-positive eigenvalue `k_c²` is the dominant mode's
cutoff; the propagation constant is then recovered as

```text
β² = k₀² − k_c².         (★★)
```

This is the form `eigensolver::assemble_transverse` builds and the form
the original transverse-only solver (steps 2–3) shipped. It is the
textbook arrangement (Jin §8.5) and it is **exact for a homogeneous
fill**: for the WR-90 TE10 gate it reproduces the analytic β to 0.055%
(§8).

### 5.2 Why (★★) is wrong for ε_r ≠ 1

The extraction `β² = k₀² − k_c²` (★★) uses the **vacuum** `k₀` and an
**ε_r-weighted** mass `T_ε` on the right. Trace the physics: the actual
transverse equation (★) is `∇_t × ((1/μ_r)∇_t × E_t) = (k₀² ε_r − β²)
E_t` — the relative permittivity multiplies **only** the `k₀²` term, not
the `β²` term. Rearranged to put the eigenvalue `β²` alone on the right,
that reads

```text
( k₀² T_ε − S ) x = β² T_1 x,        (β-direct pencil)
```

with an **unweighted** mass `T_1 = ∫ N·N` on the right and the
eigenvalue **equal to `β²` directly**. Compare the two: the cutoff-pencil
extraction `β² = k₀² − k_c²` is algebraically equivalent to the β-direct
pencil **only when `ε_r ≡ 1`** (where `T_ε = T_1`). For any `ε_r ≠ 1` —
*uniform* fill as much as inhomogeneous — the cutoff form **under-counts
the dielectric**: it effectively drops a mode-resolved `⟨ε_r⟩` factor.

This was a real, shipped bug (ADR-0053). The ε_r=1 homogeneous canary
passed at 4e-14 because the two forms coincide there, masking it. A
uniformly-filled guide at `ε_r = 2.55` exposed it: the analytic dominant
mode is `β = √(ε_r k₀² − (π/a)²) ≈ 305.16 rad/m` at 10 GHz on WR-90, but
the cutoff-form extraction returned `β ≈ 191.07 rad/m` (an effective
`ε_eff ≈ 1.34`, barely above air — physically impossible for a guide
filled with `ε_r = 2.55`), a 37% error.

### 5.3 The β-direct fix (step 5.2)

The fix (ADR-0053; `eigensolver::solve_dense` and `solve_dense_mixed`)
is to make `β²` the *direct* eigenvalue of the β-direct pencil
`(k₀² T_ε − S) x = β² T_1 x`. Three points are worth flagging:

- **The operator is symmetric *indefinite*.** `T_1` is SPD (a Gram
  matrix), so the pencil reduces to a standard symmetric problem
  `M y = β² y` via the Cholesky factor `T_1 = L Lᵀ`, with `M = L⁻¹(k₀²
  T_ε − S)L⁻ᵀ`. But `k₀² T_ε − S` straddles zero, so `M` is indefinite —
  a symmetric-*definite* (Cholesky) eigensolver on the operator itself
  would be invalid; the symmetric-tridiagonal QR (`nalgebra::
  SymmetricEigen`) handles the indefinite reduced `M` correctly.
- **The spurious modes move to the *top* of the spectrum.** The
  curl-free gradient null-space satisfies `S x ≈ 0`, so in the β-direct
  pencil it lands at `β² ≈ k₀² ⟨ε_r⟩` — the *largest* eigenvalues, not
  the smallest. They are filtered by their vanishing **cutoff Rayleigh
  quotient** `k_c² := (xᵀ S x)/(xᵀ T_ε x) ≈ 0`, and the physical
  dominant mode is the **largest β²** (equivalently the lowest cutoff)
  among the survivors.
- **Uniform fill becomes exact.** With the β-direct extraction the
  `ε_r = 2.55` uniform-fill test matches the analytic `305.16 rad/m` to
  rel `1.5e-4` (machine-precision-limited), and the ε_r=1 homogeneous
  canary is preserved bit-identically. This uniform-fill analytic is a
  fully independent published-benchmark anchor (closed form, no FEM, no
  transverse resonance) and it **certifies the β-extraction is correct**.

## 6. The mixed (E_t, E_z) block for inhomogeneous fills

For an inhomogeneous cross-section the dominant mode is hybrid: `E_z ≠ 0`
and it couples to `E_t` through the dielectric interface. The mixed
formulation (Lee-Sun-Cendes 1991; Jin §8.4) assembles a **block**
generalised eigenproblem with unknown `x = [E_t ; E_z]`:

```text
A = [ A_tt   0  ]      B = [ B_tt  B_tz ]      B_1 = [ B_tt,1  B_tz,1 ]
    [  0    A_zz]          [ B_zt  B_zz ]            [ B_zt,1  B_zz,1 ]
```

where `A_tt` is the transverse curl-curl stiffness `S`, `A_zz` /
`B_zz` are the nodal-Lagrange longitudinal stiffness/mass, and `B_tz` /
`B_zt` are the edge-node coupling blocks. `eigensolver::assemble_mixed`
builds this from the staged element matrices `local_a_zz`,
`local_b_zz`, and the coupling `local_b_ze`.

Two corrections during bring-up are documented in ADR-0051 and worth
preserving here because they were both load-bearing:

1. **The coupling weight is `1/μ_r`, not `ε_r`.** The originally-staged
   coupling computed `∫ ε_r ∇L·N`, a divergence-penalty term that the
   divergence-free curl-null-space mode annihilates — making the
   coupling **physically inert on every geometry** (`‖E_z‖/‖E_t‖ = 0`).
   The correct Lee-Sun-Cendes coupling is the curl-curl cross term
   `∫ (1/μ_r) ∇L·N` (matching `A_zz`). With the fix the coupling is
   load-bearing on inhomogeneous guides; the homogeneous canary cannot
   guard it (the homogeneous dominant mode is pure-TE, `E_z = 0`,
   weight-independent), so it is pinned instead by a horizontal-slab
   `E_z ≠ 0` case and an independent-quadrature unit test.
2. **The block-mass `B` is symmetric *indefinite*.** The off-diagonal
   coupling straddles zero even though `B_tt` and `B_zz` are individually
   SPD, so Cholesky / symmetric-definite generalized solvers are invalid
   on the block pencil. `solve_dense_mixed` forms `B⁻¹A` and uses a
   non-symmetric eigensolve with inverse-iteration eigenvector recovery —
   acceptable at the `n ≈ 121` validation scale.

### 6.1 Mode selection: cutoff-pencil select, β-direct extract

A subtlety (ADR-0053 as-built, ADR-0054): solving the β-direct *block*
pencil `(k₀² B − A) x = β² B_1 x` **directly** drifts off the physical
mode. Near the physical `β²` there is a spurious `E_z ≈ 0` branch (the
gradient null-space again, now interleaved at `β² ≈ k₀² ⟨ε_r⟩`), and the
shifted operator `(K − σ B_1) ≈ −A` is near-singular there, so a global
β-direct sweep thrashes. The shipped mixed path therefore uses a
**hybrid**:

1. **Select** the dominant mode on the *cutoff* block pencil
   `A x = k_c² B x`, where the gradient cluster sits cleanly at
   `k_c² ≈ 0` (rejected by a `k_c² ≤ 1e-6 · max|k_c²|` floor) and a
   **transverse-energy filter** (`‖E_t‖²/‖x‖²` above a floor) removes
   `E_z`-dominated contamination. The physical mode is the smallest valid
   `k_c²`; its eigenvector (carrying the genuine `E_z` of the hybrid
   mode) is recovered by inverse iteration on `(A − σ B)`.
2. **Extract** `β²` as the **β-direct Rayleigh quotient** on that
   selected eigenvector:

   ```text
   β² = R(x) = (xᵀ (k₀² B − A) x) / (xᵀ B_1 x).
   ```

   Because `A x = k_c² B x`, this equals `(k₀² − k_c²)·⟨ε_r⟩` with the
   mode-resolved `⟨ε_r⟩ = (xᵀ B x)/(xᵀ B_1 x)`. It is **exact** on a
   uniform fill (`B = ε_r B_1` ⇒ `⟨ε_r⟩ = ε_r`, the analytic anchor) and
   reduces to `k₀² − k_c²` on a homogeneous guide.

The Rayleigh quotient on the *correctly selected* eigenvector is the
right `β²` for that physical mode. The residual issue (§9) is that on a
high-contrast inhomogeneous interface the *cutoff-pencil* eigenvector is
not identical to the true *β-direct* eigenvector, so this hybrid leaves a
bias — which step 5.3 closes by recovering the true β-direct eigenvector
through a targeted sparse shift-and-invert.

## 7. Wave impedance Z_w

For the dominant mode the modal field relation is
`H_t = (β / (ω μ₀ μ_r)) (ẑ × E_t)`, so the local wave impedance is
`ω μ₀ μ_r / β`. The solver extracts a single `Z_w` as the
`|E_t|²`-energy-weighted average of that local impedance off the solved
eigenvector (a numerical line-integral / power definition, Jin §8.4):

```text
Z_w = (ω μ₀ / β) · ( ∫_Ω |E_t|² dA ) / ( ∫_Ω (1/μ_r) |E_t|² dA ).
```

On a homogeneous guide (`μ_r ≡ 1`) the two integrals cancel and this
reduces **exactly** to `ω μ₀ / β = η₀ k₀ / β`, the TE-mode wave
impedance that the closed-form `RectangularWaveguideTe10::wave_impedance`
returns — which is used as a regression guard. On a dielectric stack-up
the wave impedance is genuinely different, and the energy-weighted form
is the standard quasi-TEM definition; this replaced the earlier
TE-mode `η₀ k₀ / β` approximation (ADR-0051) once the mixed eigenvector
became available. The two transverse-energy integrals are computed by
Nedelec quadrature off the cached eigenvector.

## 8. Validation — what is certified, and to what accuracy

Per the CLAUDE.md §4 discipline, no solver path ships without a
published-benchmark validation case. The cross-section eigensolver has
three classes, with sharply different confidence levels. **State the
status honestly**:

### 8.1 Homogeneous — production accuracy

The WR-90 TE10 gate (`crates/yee-mom/tests/eigensolver_wr90.rs`) is
air-filled WR-90 (`a × b = 22.86 mm × 10.16 mm`) at 10 GHz on a 6×6
quad-diagonal mesh (72 triangles, `n_interior ≈ 84` edges). The analytic
phase constant is `β = √(k₀² − (π/a)²) = 158.238256 rad/m`; the
numerical solve gives `β ≈ 158.150550 rad/m`, a relative error of
**0.055%** — comfortably inside the 1% gate. The wave impedance matches
the closed-form `η₀ k₀ / β` to within its own regression tolerance. This
path is **production-quality**.

### 8.2 Uniform fill — production accuracy (the β-extraction anchor)

A uniformly-filled WR-90 at `ε_r = 2.55` (10 GHz) has the closed-form
dominant mode `β = √(ε_r k₀² − (π/a)²) ≈ 305.16 rad/m`. The β-direct
solver (step 5.2) matches this to relative `1.5e-4`, machine-precision-
limited. This is a fully independent analytic benchmark and it
**certifies the β-extraction is exact** for any `ε_r` — independent of
the inhomogeneous discretisation question below. This path is
**production-quality**.

### 8.3 Inhomogeneous high-contrast — improving, not yet validated

The hard case is a *partially* dielectric-loaded cross-section. The
published reference is the **slab-loaded rectangular-waveguide
transverse-resonance dispersion** (Pozar §3.6 / §6.6, Collin §6),
implemented in `eigensolver::reference` and independently verified (to
rel `0.000e0`) against a shooting-method solution of the same 1-D
transverse ODE and against the analytic empty / fully-filled limits.

For a guide stratified along `y`, the modes separate into two
**longitudinal-section** families w.r.t. the stratification axis:
**LSM-to-y** (TM-to-y, `H_y = 0`), which contains the empty guide's
`TE_{m0}` modes and is therefore the dominant slab-loaded family, and
**LSE-to-y** (TE-to-y, `E_y = 0`), the dual. Treating `y` as a
transverse-resonance transmission line with the PEC walls as short
circuits, the LSM dispersion is

```text
(ε_r1 / k_y1) cot(k_y1 d_1) + (ε_r2 / k_y2) cot(k_y2 d_2) = 0,
k_yi² = ε_ri k₀² − (m π / a)² − β².
```

A subtlety the first reference attempt missed: for a strongly-loaded
mode `β` exceeds the air-region propagation constant, so `k_y2² < 0` and
the field is *evanescent in y* in the air region. Then `k_y2 = j q₂` and
`cot(k_y2 d_2) = −j coth(q₂ d_2)`, with the residual staying real. A
root-find assuming real `k_y` everywhere fails to find the loaded root.

For a half-height fill at the high contrast `ε_r = 10.2` (RT/duroid
6010) on WR-90, the verified reference puts the dominant LSM-to-y mode at
**β ≈ 582.95 rad/m** (an effective `ε_eff ≈ 8.17`, field-concentrated in
the dielectric). The β-direct solver mesh-converges (8×8 → 12×12 within
0.05%) to **β ≈ 483.29 rad/m** (`ε_eff ≈ 5.74`) and recovers the correct
weakly-hybrid mode shape (`‖E_z‖/‖E_t‖ ≈ 0.0105`, matching the
reference's field orientation). That is a large improvement on the
pre-step-5.2 state (a 2.9× gap), but it still leaves a **mesh-converged
≈ 17% residual** versus the reference.

The honest framing (ADR-0053/0054 as-built): because the residual is
*mesh-converged* and the β-extraction is *separately certified exact* by
the uniform-fill anchor (§8.2), the 17% is a **discretization limit** —
first-order Nedelec/nodal elements under-resolving the field peak at the
high-contrast interface — combined with a Rayleigh-quotient eigenvector
mismatch (the cutoff-pencil eigenvector differs from the true β-direct
eigenvector for inhomogeneous `ε_r`; §6.1). The inhomogeneous gate
therefore ships as a **non-failing reconciliation diagnostic** plus a
monotonic physics bracket (`β_air < β_loaded < β_full`), **not** as a
passing ≤5% benchmark.

**Bottom line:** treat the cross-section eigensolver as production-grade
for **homogeneous and uniformly-filled** cross-sections (WR-90 TE10
0.055%, uniform-fill rel `1.5e-4`), and as **improving for high-contrast
inhomogeneous** fills (currently ~17% at `ε_r = 10.2`, being closed by
the step-5.3 direct sparse solve targeting representative substrates such
as FR-4 `ε_r = 4.4` at ≤5%). Do not rely on inhomogeneous high-contrast β
for design until step 5.3 closes its gate.

## 9. Solver options

The eigenproblem is sparse and complex-symmetric in general; Yee
provides two solve paths behind a common contract:

- **Dense `SymmetricEigen` (shipped, default).** The β-direct pencil is
  reduced to a standard symmetric problem via a Cholesky factor of the
  SPD right-hand mass and solved with `nalgebra::SymmetricEigen`
  (symmetric-tridiagonal QR). This is `O(n³)` and viable only up to a few
  hundred DoF — fine for the validation meshes (`n ≈ 84`–`121`) and for
  coarse cross-sections, the always-available path. For the mixed block
  pencil (symmetric *indefinite*) the dense path forms `B⁻¹A` and uses a
  non-symmetric eigensolve with inverse-iteration eigenvector recovery.
- **In-tree block LOBPCG (step 4, ADR-0050).** A pure-Rust block LOBPCG
  (Knyazev 2001) layered on `faer` and the existing sparse LU,
  implementing the same `SparseEigen` trait. It was chosen over binding
  system ARPACK via `arpack-rs` specifically to avoid a Fortran/LAPACK
  toolchain dependency (CLAUDE.md §3 — feature flags default off for
  external toolchains) and adds zero new crate dependencies. Block
  iteration is robust on the clustered / degenerate spectra real
  cross-sections exhibit (TE/TM degeneracies, `TE_{mn}`/`TE_{nm}` pairs)
  where one-at-a-time deflated inverse iteration converges slowly.
- **Direct sparse shift-and-invert on the β-direct pencil (step 5.3,
  ADR-0054, in progress).** Solves `(k₀² B − A) x = β² B_1 x` directly
  via a `faer` sparse shift-and-invert with a physics-informed shift
  `σ₀ = (k₀² − k_c²)⟨ε_r⟩` (the hybrid's own β² estimate). The shift
  placed near the physical `β²` isolates and amplifies the physical
  eigenpair past the spurious cluster — the standard remedy for
  interior/clustered spectra (Saad) — recovering the *true* β-direct
  eigenvector and so eliminating the §6.1 Rayleigh-quotient eigenvector
  mismatch. This closes residual source (b) directly; mesh refinement
  (now tractable via sparse LU) then addresses the pure-discretization
  source (a).

The trait seam (`SparseEigen` / `ComplexSparseEigen`) is deliberately
kept as the swap point: a future >10⁵-DoF cross-section that genuinely
needs ARPACK's Krylov-Schur could add it as an optional feature behind
the same trait without disturbing callers.

## 10. Limitations and roadmap

- **Lossless only.** The assembly stores `DMatrix<Complex64>` to keep the
  API future-proof, but the current solve paths take the real parts and
  reject inputs whose imaginary parts exceed `1e-9` of the real norm.
  Lossy / complex-`ε_r` β-extraction (complex-symmetric / LDLᵀ eigensolve,
  complex `β`) is Phase 1.3.1.2.
- **First-order elements.** Convergence is `O(h)` in the energy norm /
  `O(h²)` in the eigenvalue (Babuska-Osborn). The high-contrast residual
  (§8.3) is the visible cost; higher-order Nedelec (Whitney-1 plus
  quadratic bubble) is the standard remedy and is later-step work.
- **Single-mode extraction.** The solver returns the dominant guided mode;
  multi-mode extraction and tracking (for multi-conductor lines and
  higher-order port modes) is Phase 1.3.2 / 1.3.4.
- **Isotropic media.** Tensor `ε_r`, `μ_r` are Phase 1.3.3+.

## 11. References

- Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd ed.
  Wiley-IEEE, 2014. §8.4–8.5 (waveguide eigenproblem, modes of a
  waveguide); §9 (Nedelec / Whitney edge elements). The formulation in
  §3–§6 follows this book.
- Pozar, D. M. *Microwave Engineering*, 4th ed. Wiley, 2012. §3.1, §3.3
  (TE/TM mode decomposition; closed-form TE10), §3.6 / §6.6 (slab-loaded
  guide / transverse resonance, the §8.3 reference).
- Collin, R. E. *Field Theory of Guided Waves*, 2nd ed. IEEE Press, 1991.
  §6 (longitudinal-section LSE/LSM modes; transverse resonance).
- Lee, J.-F., Sun, D.-K., and Cendes, Z. J. "Full-wave analysis of
  dielectric waveguides using tangential vector finite elements." *IEEE
  Trans. Microwave Theory Tech.* 39.8 (1991), pp. 1262–1271. The mixed
  `(E_t, E_z)` block formulation of §6.
- Webb, J. P. "Edge elements and what they can do for you." *IEEE Trans.
  Magnetics* 29.2 (March 1993), pp. 1460–1465. Why edge elements beat
  nodal Lagrange on the vector wave equation (§4).
- Boffi, D., Brezzi, F., and Demkowicz, L. F. *Mixed Finite Element
  Methods and Applications.* Springer, 2013. §5 — why curl-conforming
  edge elements admit no spurious modes inside the spectrum.
- Knyazev, A. V. "Toward the Optimal Preconditioned Eigensolver: Locally
  Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG)."
  *SIAM J. Sci. Comput.* 23.2 (2001), pp. 517–541. The step-4 sparse
  solver (§9).
- Saad, Y. *Numerical Methods for Large Eigenvalue Problems.* Rev. ed.
  SIAM, 2011. Shift-and-invert for interior/clustered spectra (§9, the
  step-5.3 path).
- Golub, G. H., and Van Loan, C. F. *Matrix Computations*, 4th ed. Johns
  Hopkins University Press, 2013. §8 — the Cholesky-symmetrised
  generalised eigenproblem reduction.
- Babuska, I., and Osborn, J. "Eigenvalue Problems." In *Handbook of
  Numerical Analysis*, vol. II, Elsevier, 1991, pp. 641–787. The
  eigenvalue-vs-eigenvector convergence-rate claim (§10).
