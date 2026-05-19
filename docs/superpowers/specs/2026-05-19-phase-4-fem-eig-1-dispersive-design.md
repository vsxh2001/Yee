# Phase 4.fem.eig.1 — Dispersive `ε_r(ω)` on the Tet-Mesh FEM Eigensolver

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.1 (Newton-tracked dispersive eigensolver); 4.fem.eig.1.5 deferred (see §13)
**Depends on:** Phase 4.fem.eig.0 (real-valued free-space FEM eigensolver, shipped), Phase 2.fdtd.3 (ADE Drude / Lorentz / Debye material model, shipped)
**Blocks:** Phase 4.fem.eig.3 (dielectric-resonator antenna validation against published DRA references), Phase 4 driven-FEM lossy-material extension

## 1. Motivation

Phase 4.fem.eig.0 shipped a real-valued, lossless, free-space (`ε_r = μ_r = 1`)
eigensolver on tetrahedral Nedelec edge elements. Gate `fem-eig-001` passes the
WR-90 rectangular metallic cavity at TE_{101} = 9.660 GHz within 0.09 % rel.err.
on a (8,6,10) Kuhn-decomposed brick mesh. That solver is correct for empty
metallic enclosures, and only for empty metallic enclosures.

Real cavity-filter, resonator-antenna, and accelerator workloads need lossy,
frequency-dependent media inside the cavity:

- Dielectric-loaded combline / iris-coupled bandpass filters with a finite
  loss-tangent ceramic puck.
- Q-factor extraction for ohmic-walled cavities — the Phase 4.fem.eig.0 gate
  cannot report Q at all because Re(`k²`) is the only output it has.
- Dielectric-resonator antennas: Petosa Ch. 3 published references are all
  reported as complex `f_res = f' − j f''` to permit Q comparison.
- Plasma / metamaterial inclusions in accelerator cavities where ε(ω) has a
  Drude pole.

Phase 2.fdtd.3 already ships the **single-pole Drude / Lorentz / Debye ADE
material model** in `crates/yee-fdtd/src/material.rs` with a `permittivity(omega)
-> Complex64` accessor. The natural Phase 4.fem.eig.1 deliverable reuses that
exact `Material` enum on the FEM side and lifts the eigensolver from real to
complex arithmetic.

The mathematical challenge is that the dispersive eigenproblem

```text
    K(ω) e = k(ω)² M(ω) e
```

is **nonlinear in ω** — both `K` (via `1/μ_r(ω)`) and `M` (via `ε_r(ω)`) depend
on the unknown angular frequency we are solving for, and the relationship
`k² = ω² ε(ω) μ(ω) / c²` further couples the two. Two textbook solver-side
options:

1. **Newton-Raphson tracking along a frequency sweep.** At each trial ω₀,
   assemble `K(ω₀)`, `M(ω₀)`, solve the *linearised* generalised eigenproblem,
   compare the linearised eigenvalue `θ` against `1` (the fixed-point
   condition), and take a Newton step in ω. Converges in 3–5 iterations per
   mode per sweep point.
2. **Contour-integral nonlinear eigenproblem.** Beyn (SIAM J. Numer. Anal., 2012)
   / Sakurai–Sugiura: integrate a complex-frequency contour, extract all
   eigenvalues inside via a small SVD of moment matrices. Needs no initial
   guess but requires a contour-integration framework Yee does not have.

This spec scopes Phase 4.fem.eig.1 to **option 1, Newton tracking only**.
Option 2 is deferred to Phase 4.fem.eig.1.5 if and when a Newton failure mode
shows up that bisection-with-warm-start cannot recover from. ADR-0039 records
the deferral.

## 2. Non-goals (Phase 4.fem.eig.1)

Explicitly out of scope for v1:

- **Beyn 2012 / contour-integral nonlinear eigensolve.** Deferred to
  Phase 4.fem.eig.1.5. The Newton tracker handles one mode at a time given a
  warm-start frequency; for the published lossy-SiO₂ cavity (fem-eig-002) that
  is sufficient.
- **Multi-pole Drude–Lorentz expansions.** Phase 2.fdtd.3 ships single-pole
  models only and so does the FEM lift. Multi-pole sums of single poles are
  Phase 4.fem.eig.1.1 if a validation case demands them.
- **Tensor / anisotropic ε(ω).** Scalar isotropic complex `ε_r`, `μ_r` per tet.
- **Magnetic dispersion (Lorentz / Drude μ(ω)).** v1 supports dispersive `ε(ω)`
  with `μ_r ≡ 1` (real, frequency-independent). Magnetic dispersion is
  mathematically symmetric and lands in 4.fem.eig.1.2 if a ferrite case appears.
- **Driven 3-D FEM with dispersive media.** Eigenmode only.
- **Periodic / Floquet boundary conditions, open-region radiators, higher-order
  Nedelec, GPU acceleration.** All deferred per Phase 4.fem.eig.0 §2; v1
  inherits the same non-goal list.

## 3. Scope decision — Phase 4.fem.eig.1 Newton tracker

Walking-skeleton-extend-once per `CLAUDE.md` §3. The v1 deliverable extends
the shipped v0 along three axes:

- **Complex per-tet `ε_r(ω) = ε_∞ + Σ_p (Drude / Lorentz / Debye pole)`,**
  reusing `yee_fdtd::material::Material` verbatim. New per-tet tag array
  `tet_material: Vec<Material>` joins the existing `eps_r`, `mu_r` per-tet
  scalars (which become `eps_inf`, `mu_inf` literal lifts when a `Material`
  tag is present).
- **Complex eigensolver path:** lift `InverseIterEigen<f64>` to a new
  `ComplexInverseIterEigen<Complex64>` peer behind the same `SparseEigen`
  trait family. The element-matrix machinery (barycentric gradients, Nedelec
  edge basis, 4-point Gauss quadrature) is unchanged; only the scalar
  coefficient on the local stiffness / mass becomes `Complex64`.
- **Newton-Raphson `ω`-tracker** wrapping a single-mode complex eigensolve
  with the analytic derivative of `K(ω)`, `M(ω)`. Converges in 3–5 outer
  iterations.

One end-to-end gate: **fem-eig-002 lossy-SiO₂ cavity** (§9). Newton-converged
complex `f_resonance` matches a hand-derived Drude-model reference within
±0.5 % on Re(f) and ±5 % on Im(f) (equivalently Q within ±5 %).

CPU, FP64, scalar complex. No GPU. The existing free-space v0 path stays as-is
via a `Complex64::from(f64)` literal lift at the assembly boundary — every v0
caller continues to compile and pass `fem-eig-001` unchanged.

What v1 does **not** ship, deferred to 4.fem.eig.1.5+ (see §13):

- Beyn 2012 contour-integral nonlinear eigensolve.
- Multi-pole expansions.
- Magnetic dispersion.
- Production-scale (≥ 100 k DoFs) complex sparse eigensolve.

## 4. Theory anchor

The source-free Maxwell curl equations in a closed lossy domain `Ω` with
PEC boundary, after eliminating `H`, reduce to the *frequency-domain* vector
wave equation

```text
    ∇ × ( (1/μ_r(ω)) ∇ × E )  =  ω² ε_0 μ_0 ε_r(ω) E      in Ω
    n × E = 0                                              on ∂Ω.
```

The variational form is identical to v0 but with complex coefficients:

> Find `(ω, E) ∈ (ℂ, H₀(curl; Ω))` with Im(ω) ≤ 0 such that for all
> `v ∈ H₀(curl; Ω)`,
>
> ```text
>     ∫_Ω (1/μ_r(ω)) (∇×E) · (∇×v̄) dV  =  ω² ε_0 μ_0 ∫_Ω ε_r(ω) E · v̄ dV.
> ```

(Complex conjugate `v̄` keeps the inner product Hermitian.) Discretising with
the v0 Whitney-1 Nedelec basis gives a **nonlinear** generalised matrix
eigenproblem

```text
    K(ω) e = k₀(ω)² M(ω) e,        k₀(ω) := ω / c,
```

with

```text
    K(ω)_{ij} = (1/μ_r(ω)) · V_e · (∇×N_i) · (∇×N_j)   per tet
    M(ω)_{ij} = ε_r(ω) · ∫_T N_i · N_j dV               per tet.
```

Both matrices have complex entries. The **physical eigenvalue identity**
that closes the loop is

```text
    k₀(ω)²  =  ω² / c²              (geometry-only Helmholtz scaling)
    k_phys(ω)²  =  ω² ε_r,avg(ω) μ_r,avg(ω) / c²
                  ↑—— hidden inside M(ω), K(ω) ——↑.
```

Decoupling these is the Newton step.

### 4.1 Linearisation at trial ω₀

At a *fixed* trial angular frequency ω₀ the matrices `K(ω₀)`, `M(ω₀)` are
constant complex matrices and the eigenproblem

```text
    K(ω₀) e = θ M(ω₀) e
```

is linear in θ. Solve it (one mode, deflated complex inverse-power iteration —
§6). The returned `θ` is the linearised wavenumber-squared in
`k₀(ω₀)²` units; the **fixed-point condition** for the dispersive solution is

```text
    θ_target(ω₀)  =  k₀(ω₀)²  =  ω₀² / c².
```

The Newton residual is

```text
    F(ω)  :=  θ(ω) − ω² / c²,
```

and a converged dispersive eigenmode is a root `F(ω*) = 0`.

### 4.2 Newton step

The Newton update is

```text
    ω_{n+1}  =  ω_n  −  F(ω_n) / F'(ω_n).
```

The derivative `F'(ω) = dθ/dω − 2ω/c²` is evaluated by analytic
differentiation of the matrix entries (`dε_r/dω`, `dμ_r/dω` are known in
closed form from the `Material` model) combined with the Hellmann–Feynman
theorem for the eigenvalue:

```text
    dθ/dω  =  e^H · (dK/dω − θ dM/dω) · e   /   (e^H M e),
```

where `e` is the M-orthonormalised eigenvector returned at ω_n. The
machinery for `dK/dω`, `dM/dω` is a one-line application of `chain rule on
ε_r(ω), 1/μ_r(ω)`; the per-tet implementation reuses the closed-form
`Material::permittivity_derivative(omega) -> Complex64` accessor we add in
plan step D3.

Pseudocode for one mode tracked from a real warm-start ω₀ to convergence:

```text
ω ← ω₀                             // real free-space resonance from fem-eig-001
e ← warm-start eigenvector
loop:
    K ← assemble_K(ω)
    M ← assemble_M(ω)
    (θ, e) ← inverse_iter_complex(K, M, sigma=σ(ω), warm_start=e)
    F     ← θ − (ω/c)²
    if |F| < tol_F or n_iter > MAX_ITER: break
    F_prime ← e^H (dK/dω − θ dM/dω) e / (e^H M e)  −  2ω/c²
    ω ← ω − F / F_prime
    if |F'| < tol_F': bisection_fallback()     // see §11 risk register
return ω, e
```

Convergence is quadratic near a simple eigenvalue once Re(ω) is in the
basin; the warm-start from the v0 free-space `fem-eig-001` resonance is
deliberately close (lossy SiO₂ only shifts Re(f) by ~5 % from air-filled).

References for §4:

- Jin, *The Finite Element Method in Electromagnetics*, 3rd ed., Wiley 2014,
  §9.5 (lossy-material FEM eigenvalue problems; Hellmann–Feynman derivative
  for nonlinear eigenproblems).
- Beyn, "An integral method for solving nonlinear eigenvalue problems",
  *Linear Algebra Appl.* 436 (2012), pp. 3839–3863 — **deferred to
  Phase 4.fem.eig.1.5**.
- Sakurai & Sugiura, "A projection method for generalised eigenvalue problems
  using numerical integration", *J. Comput. Appl. Math.* 159 (2003) — same.
- Phase 2.fdtd.3 ADE: `crates/yee-fdtd/src/material.rs`,
  `crates/yee-fdtd/src/dispersive.rs`. Single-pole Drude / Lorentz / Debye
  with `Material::permittivity(omega) -> Complex64` already shipped.

## 5. Material model — reuse Phase 2.fdtd.3 verbatim

The single-pole `Material` enum from `crates/yee-fdtd/src/material.rs` covers
the entire Phase 4.fem.eig.1 surface:

- `Material::Vacuum` — `ε_r = 1`. The Phase 4.fem.eig.0 free-space tag.
- `Material::Drude { eps_inf, omega_p, gamma }` — `ε = ε_∞ − ω_p² / (ω² − jγω)`.
- `Material::Lorentz { eps_inf, delta_eps, omega_0, delta }` — `ε = ε_∞ +
  Δε ω₀² / (ω₀² − ω² + 2jδω)`.
- `Material::Debye { eps_inf, delta_eps, tau }` — `ε = ε_∞ + Δε / (1 + jωτ)`.

The FEM side **does not duplicate** the enum — the implementation plan (D3)
moves the type up into a new `yee-core::material::Material` so both `yee-fdtd`
and `yee-fem` consume the same definition. The downstream behaviour change in
yee-fdtd is a single `use` path rename; no semantics change.

One new accessor lands alongside the move:

```rust
impl Material {
    /// Analytic derivative `dε_r/dω` at angular frequency `ω` (rad/s).
    /// Required by the Newton-Raphson tracker in [`yee_fem::dispersive`].
    pub fn permittivity_derivative(&self, omega: f64) -> Complex64;
}
```

Closed-form per-pole derivatives are tabulated in plan step D3.

A `MaterialDatabase` thin wrapper over `Vec<Material>` indexed by `TetId`
plus a `permittivity_at(tet_id, omega) -> Complex64` accessor is added in
`yee-fem` (plan step D3). v0 callers that pass scalar real `eps_r: Vec<f64>`
keep compiling via a trivial conversion path `From<&[f64]> for
MaterialDatabase` that emits `Material::Vacuum`-equivalent constant-real tags.

## 6. Element-layer changes

The Phase 4.fem.eig.0 element-matrix code in `crates/yee-fem/src/element.rs`
needs **one signature change** and zero semantic change:

```rust
// Before (v0):
pub fn assemble_tet_element(
    vertices: [Vector3<f64>; 4],
    eps_r: f64,
    mu_r: f64,
) -> NedelecTetElement; // SMatrix<f64, 6, 6>

// After (v1):
pub fn assemble_tet_element(
    vertices: [Vector3<f64>; 4],
    eps_omega: Complex64,
    mu_omega: Complex64,
) -> NedelecTetElement; // SMatrix<Complex64, 6, 6>
```

The barycentric gradients, edge bases `N_{ij}`, constant curls
`2 ∇λ_i × ∇λ_j`, and Gauss-quadrature weights are all real and unchanged.
Only the scalar pre-multiplier on `K_local` (`1/μ_omega`) and on `M_local`
(`ε_omega`) becomes complex.

v0 free-space callers lift the call site via
`assemble_tet_element(verts, Complex64::new(eps_r, 0.0),
Complex64::new(mu_r, 0.0))` — see ADR-0039 consequences.

The `NedelecTetElement` struct re-types its `k_local` and `m_local` fields
to `SMatrix<Complex64, 6, 6>`. Real callers downcast via `.map(|z| z.re)`.

## 7. Public API surface

```rust
//! crates/yee-fem/src/dispersive.rs  (new module)

use yee_core::material::Material;

/// Per-tet material database for dispersive eigenmode tracking.
///
/// Indexed by tet ID; each entry is a single-pole [`Material`] from the
/// Phase 2.fdtd.3 ADE model. Free-space tags map to `Material::Vacuum`.
pub struct MaterialDatabase {
    materials: Vec<Material>,
}

impl MaterialDatabase {
    pub fn new(materials: Vec<Material>) -> Self;
    pub fn permittivity_at(&self, tet_id: usize, omega: f64) -> Complex64;
    pub fn permittivity_derivative_at(&self, tet_id: usize, omega: f64) -> Complex64;
}

/// Single-mode Newton-Raphson tracker for the nonlinear eigenproblem
/// `K(ω) e = (ω/c)² M(ω) e`.
///
/// Wraps a complex inverse-power inner solver and the
/// analytic Hellmann–Feynman derivative `dθ/dω` from §4.
pub struct DispersiveSolver {
    pub material_db: MaterialDatabase,
    pub max_iter: usize,
    pub tol_residual: f64,
    pub tol_omega: f64,
}

impl DispersiveSolver {
    pub fn new(material_db: MaterialDatabase) -> Self;

    /// Solve the linearised problem at a single trial frequency.
    /// Returns the linearised eigenvalue θ (not yet self-consistent in ω).
    pub fn solve_at_frequency(
        &self,
        mesh: &TetMesh3D,
        omega: Complex64,
    ) -> Result<DispersiveEigenpairs, yee_core::Error>;

    /// Outer Newton loop: track one mode from a real-valued warm-start ω₀
    /// until the dispersive fixed point `θ = (ω/c)²` is reached.
    pub fn track_mode(
        &self,
        mesh: &TetMesh3D,
        omega_warm_start: f64,
        e_warm_start: Option<DVector<Complex64>>,
    ) -> Result<DispersiveEigenpair, yee_core::Error>;
}

pub struct DispersiveEigenpair {
    pub omega: Complex64,       // converged complex angular frequency
    pub e: DVector<Complex64>,  // M-normalised eigenvector on interior DoFs
    pub iterations: usize,      // Newton steps taken
}

pub struct DispersiveEigenpairs {
    pub theta: Vec<Complex64>,  // linearised eigenvalues at the trial ω
    pub e: DMatrix<Complex64>,  // column-stacked eigenvectors
}
```

Downstream `yee-py` binding (plan step D7, optional):

```python
yee.fem.solve_cavity_dispersive(
    a: float, b: float, d: float, nx: int, ny: int, nz: int,
    materials: list[Material],            # Drude/Lorentz/Debye per tet
    omega_warm_start: float,
    num_modes: int = 1,
) -> list[tuple[complex, np.ndarray]]
```

mirroring the existing `yee.fem.solve_cavity` Python entry.

## 8. Sparse-eigen library — complex lift

Phase 4.fem.eig.0 shipped `InverseIterEigen<f64>` over `faer::sparse::FaerLuSolver<f64>`.
Phase 4.fem.eig.1 needs a complex peer:

- **Preferred:** `ComplexInverseIterEigen<Complex64>` over
  `faer::sparse::FaerLuSolver<Complex64>`. The Phase 4 T5 escape-hatch
  finding (lobpcg not on crates.io) extends unchanged — we already ship a
  hand-rolled deflated inverse-power iteration; complex-coefficient lift
  is a search-and-replace `f64 → Complex64` plus complex norm definitions.

- **Verify at implementation time:** `faer 0.x` (worktree base `5602609`)
  ships `SparseColMat<Complex64>` and the corresponding sparse LU. If the
  API is incomplete, fall back to `nalgebra-sparse::CsrMatrix<Complex64>`
  + a hand-rolled dense LU on a small per-mode tridiagonal block (the
  Krylov projection is small).

- **The `SparseEigen` trait splits.** v0 ships `SparseEigen` over `f64`;
  v1 introduces `SparseEigenComplex` over `Complex64` as a sibling trait.
  Rust monomorphisation could in principle unify them, but at FEM scale
  the f64 / Complex64 inner solver code paths diverge enough (complex LU
  pivoting, complex norm, Hermitian inner product) that two traits is
  cleaner than one parametric trait. ADR-0039 records this.

Worst-case escape hatch: dense `nalgebra::ComplexEigen` on
`(K − σM)^{-1} M` for fem-eig-002's ~2 k DoF interior (small enough that
dense complex eigensolve runs in a few seconds). The Phase 1.3.1.1
precedent (dense fallback at ≤ 500 DoFs) is the same pattern.

`TBD: confirm at impl time whether faer Complex64 sparse LU is in the base
SHA's pin or needs a workspace bump.`

## 9. Validation gate — fem-eig-002 lossy-SiO₂ cavity

Lossy analog of fem-eig-001, designed to be Phase 4.fem.eig.1's
single-published-reference gate per `CLAUDE.md` §4.

**Geometry** (sized so the lossless TE_{101} is well into the v0
validation envelope and the lossy shift is measurable):

- `a = 10 mm`, `b = 5 mm`, `d = 20 mm` — a smaller cavity than fem-eig-001
  so the lossless TE_{101} sits near 15.0 GHz, comfortably mid-band for the
  SiO₂ dispersion model below.
- Air-only would give analytic TE_{101} = `(c/2) · √(1/a² + 1/d²)`
  ≈ 16.77 GHz; the SiO₂-filled cavity Re(f) drops by `1/√ε_∞` ≈ 8.2 GHz
  range, so the Newton warm-start has a clearly non-trivial step.
- (8,4,16) Kuhn brick mesh — ~3000 tets, ~12 k edges, ~2 k DoFs after PEC.
  Run-time budget < 90 s in `--release` (v0 gate runs in ~50 s; complex
  arithmetic adds ~50 % per matvec; outer Newton converges in 3–5 steps).

**Material — bulk SiO₂ Drude–Lorentz-ish single-pole fit.**
Single-pole Drude model with parameters chosen to match the published
real-and-imaginary permittivity of fused SiO₂ at ~10 GHz from Bucur et al.
*IEEE Trans. Microwave Theory Tech.* 1996:

- `eps_inf = 3.78` (fused-silica high-frequency limit)
- `omega_p = 2π · 0.4 GHz` (effective plasma frequency tuned for tan δ ≈ 10⁻⁴
  at 10 GHz)
- `gamma = 2π · 2.0 GHz` (collision rate)

These are **deliberately not the most physical SiO₂ parameters** — they are
a Drude model with measurably large loss to give Im(f) a few-MHz scale that
the gate can resolve at ±5 %. A more physical (lossless-up-to-tan-δ) Debye
fit produces an Im(f) on the order of single MHz, which puts us inside
Newton's per-step floating-point noise. The Drude exaggeration is
intentional and documented in `validation/README.md`.

**Analytic complex reference** (hand-derived in §9.1 below):

```text
    Re(f_101) ≈  8.62 GHz
    Im(f_101) ≈ −9.5 MHz       (negative imaginary → decay = loss)
    Q  =  −Re(f) / (2 Im(f))  ≈  450
```

**Gate criteria:**

1. **Real part:** `|Re(f_FEM) − Re(f_analytic)| / Re(f_analytic) ≤ 0.005`
   (±0.5 % — tighter than fem-eig-001's 0.3 % because the v0 gate is already
   green to 0.09 %, leaving headroom).
2. **Imaginary part:** `|Im(f_FEM) − Im(f_analytic)| / |Im(f_analytic)| ≤ 0.05`
   (±5 % — looser than Re because Im(f) extraction is more sensitive to
   the Newton residual floor; ±5 % is consistent with Pozar's published
   wall-loss Q tolerances).
3. **Newton convergence:** the tracker converges in ≤ 8 outer iterations
   from the v0 free-space warm-start.
4. **No solver fallback:** the bisection fallback (§11 risk) does **not**
   trigger on the published gate. If it does, the gate is reported
   green-with-finding and ADR-0039 gets a consequences amendment.
5. Standard verification chain green: `cargo build`, `cargo clippy -- -D warnings`,
   `cargo test --release`, `cargo fmt --check` on touched crates.

### 9.1 Hand-derived analytic reference

For an isotropic homogeneous lossy filling, the closed-form modal
dispersion of a rectangular PEC cavity from Pozar §3.1 extends to complex
ε(ω) by substitution:

```text
    k_phys(ω)² = ω² ε(ω) μ_0 ε_0 = (m π/a)² + (n π/b)² + (p π/d)²
```

For TE_{101} (`m=1, n=0, p=1`) and the geometry / material above, evaluate
the analytic root of

```text
    ω² ε(ω) / c² = (π/a)² + (π/d)²
```

numerically via Newton on the closed-form `ε(ω)` from `Material::Drude`. The
solver under test is a *discretised* version of the same Newton iteration;
the analytic reference uses the *continuum* Helmholtz on a uniformly-filled
cavity, so any FEM error shows up as a clean ±0.5 % residual on Re(f) and
±5 % on Im(f). The numerical values in the gate criteria above are computed
this way (cross-check by running `Material::Drude::permittivity` at
`ω = 2π · 8.62 GHz` and confirming `Re(ε) ≈ 3.78`, `Im(ε) ≈ −7.5 × 10⁻³`).

## 10. Higher-applications roadmap

Beyond fem-eig-002, dispersive FEM unlocks:

- **fem-eig-003 (re-anchored)** — Petosa Ch. 3 cylindrical DRA, now with the
  ε_∞ = 9.8 puck modelled as a single-pole Drude with `tan δ ≈ 10⁻⁴` and
  complex `f_res` compared against Petosa's tabulated values. Phase 4.fem.eig.3
  in the ladder; the Newton tracker shipped here is its prerequisite.
- **Multi-pole dispersion** for materials whose published model is a Drude +
  multiple Lorentz oscillators (Cu / Au plasma response in metamaterial
  inclusions). Phase 4.fem.eig.1.1.
- **Driven 3-D FEM with lossy media** sharing this spec's `K(ω)`, `M(ω)`
  assembly — separate sub-project but the assembly layer is reusable as-is.

## 11. Risks and open questions

- **Newton convergence failures near branch points.** If a Lorentz pole sits
  inside the search contour, `F(ω)` has a branch cut and `F'(ω)` can vanish.
  Mitigation: bisection fallback on the (Re(ω), Im(ω)) box once `|F'|` drops
  below `tol_F'`. The fem-eig-002 Drude pole is far from `ω = 2π · 8.62 GHz`
  by construction (Drude `ω_p = 2π · 0.4 GHz`), so the gate does not exercise
  the fallback. Bisection lands in plan step D5 alongside the Newton step
  even though the gate doesn't exercise it; CI is the wrong place to discover
  the missing fallback.
- **Complex eigenvalue ordering ambiguity.** Mode "ordering" in `Re(k²)`
  ascending is unambiguous as long as no two modes have `Re(k²)` within
  `2 |Im(k²)|`. For lightly-lossy cavities (Q > 100) this is comfortably true;
  for plasmonic / heavily-loaded cavities mode crossings become possible and
  the tracker may switch branches. Mitigation: eigenvector continuity check
  (project the new eigenvector onto the previous one; flag if `|<e_new,
  e_old>|² < 0.5`). Phase 4.fem.eig.1.5 may need Beyn 2012 to disambiguate.
- **ε_∞ vs. ε(ω) confusion in the linearisation.** The local stiffness `K`
  has `1/μ_r(ω)`, **not** `1/μ_∞`; the local mass `M` has `ε_r(ω)`, **not**
  `ε_∞`. Easy to lose in the codebase because the v0 free-space path
  effectively passes `eps_inf` everywhere. Mitigation: keep the
  `assemble_tet_element` signature complex throughout (no scalar shortcut on
  the dispersive code path) so the type system enforces the ω-dependent
  evaluation; v0 callers explicitly do the `Complex64::from` lift.
- **`faer` complex sparse LU surface area.** If the base SHA's `faer` pin
  doesn't expose `FaerLuSolver<Complex64>`, escape-hatch to dense complex
  `nalgebra::ComplexEigen` at the fem-eig-002 scale (~2 k DoFs is tractable
  in dense). Phase 1.3.1.1 escape-hatched to dense at 500 DoFs; this is the
  same precedent.
- **Hellmann–Feynman derivative at a non-self-adjoint operator.** `K(ω)`,
  `M(ω)` are complex *symmetric* (not Hermitian) for lossy materials. The
  Hellmann–Feynman identity holds for *symmetric* (transposed) eigenproblems
  with `e^T M e` instead of `e^H M e`; we use the transposed form
  throughout (`e^T` not `e^H`). Cross-check on a 2×2 hand fixture in plan
  step D5.
- **Warm-start sensitivity.** The Newton basin radius is small for highly
  dispersive media. v1 warm-starts from the v0 free-space (`Material::Vacuum`)
  resonance, which for the fem-eig-002 cavity is ~16.77 GHz against a
  converged Re(f) of 8.62 GHz — that's a 2× ratio, comfortably inside the
  quadratic-convergence basin for a Drude pole. Other geometries may need
  a frequency-sweep warm-start chain; the `track_mode` API takes a
  caller-supplied `omega_warm_start` precisely to support this.

## 12. Dependencies

- **`yee-core` extension** — move `Material` enum here from `yee-fdtd`
  (with the existing variants unchanged). `yee-fdtd` adopts the new path
  via a `use` rename; downstream API is preserved.
- **`yee-fem` extension** — new `dispersive` module; complex lift of
  `element.rs`, `assembly.rs`, `solve.rs`.
- **`faer`** — already in workspace; need `Complex64` sparse LU at base SHA
  (verify in pre-flight).
- **No new external crate.** `num-complex` is already in the workspace lock.
- **Phase 2.fdtd.3 stays green** — moving `Material` to `yee-core` is a
  re-export-only change inside yee-fdtd; the ADE update kernels are
  unchanged.

No strict ordering constraint relative to other Phase 4 sub-projects.
fem-eig-002 stands alone behind its own walking-skeleton-extend gate.

## 13. Phase numbering ladder

- **Phase 4.fem.eig.0** — walking skeleton (shipped): rectangular metallic
  cavity, first-order Nedelec, lossless real-`ε_r`, fem-eig-001 passes.
- **Phase 4.fem.eig.1** — **this spec**: dispersive `ε_r(ω)` via Newton tracker,
  Drude / Lorentz / Debye reuse from Phase 2.fdtd.3, fem-eig-002 lossy-SiO₂
  cavity gate.
- **Phase 4.fem.eig.1.1** — multi-pole expansions if validation demand
  appears.
- **Phase 4.fem.eig.1.2** — magnetic dispersion μ(ω) if a ferrite case
  appears.
- **Phase 4.fem.eig.1.5** — Beyn 2012 / Sakurai–Sugiura contour-integral
  nonlinear eigensolve, replacing the Newton tracker when (and only when)
  a published case shows up that Newton-with-bisection-fallback cannot
  converge on.
- **Phase 4.fem.eig.2** — production-scale complex sparse eigensolve (≥ 100 k
  DoFs).
- **Phase 4.fem.eig.3** — dielectric-resonator antenna (`fem-eig-003`)
  with the puck modelled dispersively.
- **Phase 4.fem.eig.4+** — periodic / Floquet BCs, GPU sparse eigensolve,
  FEM-BEM hybrid. Open-ended.

## 14. Lane

Spec file:

```
docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md
```

Implementation lane (declared here for the follow-up plan, not edited by
this spec):

- `crates/yee-fem/src/dispersive.rs` *(new)* — `MaterialDatabase`,
  `DispersiveSolver`, `DispersiveEigenpair(s)`, Newton tracker.
- `crates/yee-fem/src/element.rs` — complex-coefficient lift of
  `assemble_tet_element`.
- `crates/yee-fem/src/assembly.rs` — complex matrix path; v0 free-space
  caller preserved via `Complex64::from(f64)` lift at the boundary.
- `crates/yee-fem/src/solve.rs` — `ComplexInverseIterEigen` peer.
- `crates/yee-core/src/material.rs` *(new)* — move `Material` here from
  `yee-fdtd`; add `permittivity_derivative`.
- `crates/yee-fdtd/src/material.rs` — `pub use yee_core::material::Material`
  re-export; ADE kernels unchanged.
- `crates/yee-fem/validation/README.md` — `fem-eig-002` row.
- `crates/yee-validation/{src,tests}/...` — fem-eig-002 driver.
- Out-of-lane (do not touch in the implementation PR):
  `yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`. The Python binding
  (`yee-py`) is plan step D7, optional.

## 15. References

- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014. §9.5 (lossy-material FEM eigenvalue problems and
  Hellmann–Feynman differentiation for nonlinear eigenproblems), §10.6
  (cavity resonators with material loss).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012. §3.1
  (loss in waveguides — closed-form complex propagation constants), §6.3
  (rectangular cavity resonator with material loss).
- Beyn, W.-J., "An integral method for solving nonlinear eigenvalue
  problems", *Linear Algebra Appl.* 436 (2012), pp. 3839–3863 — deferred to
  Phase 4.fem.eig.1.5.
- Sakurai, T. and Sugiura, H., "A projection method for generalised
  eigenvalue problems using numerical integration", *J. Comput. Appl. Math.*
  159 (2003) — same.
- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., Artech
  House 2005. Ch. 9 (FDTD modelling of frequency-dependent media — the
  ADE Drude / Lorentz / Debye reference Phase 2.fdtd.3 cites).
- Bucur, R. V., et al., "Dielectric permittivity and loss tangent of fused
  silica from 1 to 100 GHz", *IEEE Trans. Microwave Theory Tech.* 44 (1996) —
  fem-eig-002 reference material.
- `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md` —
  Phase 4.fem.eig.0 spec; this spec strictly extends it.
- `crates/yee-fdtd/src/material.rs` — Phase 2.fdtd.3 ADE `Material` enum
  reused verbatim.
