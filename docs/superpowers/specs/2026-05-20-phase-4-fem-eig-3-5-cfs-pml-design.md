# Phase 4.fem.eig.3.5 — CFS-PML open-boundary truncation

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5 (Complex Frequency Shifted Perfectly Matched
Layer replacing the Engquist-Majda intrinsic-reflection floor on the
open-boundary FEM solver).
**Depends on:** Phase 4.fem.eig.3 (coupled-Whitney + 2nd-order
Engquist-Majda ABC + multi-port `S_{p,q}` matrix shipped F1-F7;
NNNNNNNNN mesh refinement to `(24, 12, 36)` shipped on top).
**Blocks:** retirement of the fem-eig-003 strict `[-45, -35] dB`
absorption-floor and strict `|S_{11}| < 1` continuum-limit gates;
high-aspect-ratio cavity validation (fem-eig-006).

## 1. Goal

Phase 4.fem.eig.3 shipped the 2nd-order Engquist-Majda ABC (ADR-0042 §F3
+ F4) which moves the normal-incidence reflection floor from `~ −40 dB`
(1st-order) toward the `~ −60 dB` Jin §10.4 Table 10.1 figure for plane
waves at normal incidence. Track NNNNNNNNN (2026-05-20, commit `9cf4b3a`,
merge `fc62f82`) refined the fem-eig-003 mesh from `(16, 8, 24) = 18 432
tets` to `(24, 12, 36) = 62 208 tets` (~3.4× tets, ~24 linear samples
across the WR-90 broad wall) and measured the swept WR-90-stub band
`|S_{11}(f)| ∈ [0.9976, 0.99997]` ⇒ `s11_db ∈ [-2.22e-2, -2.86e-5] dB`
across 8-12 GHz. The refined mesh is ~2× better in dB than the
JJJJJJJJJ `(16, 8, 24)` baseline but **still ~35 dB above** the spec §8
`[-45, -35] dB` window. ADR-0042 §risks called out the band-edge /
off-normal-incidence intrinsic floor of the 2nd-order Engquist-Majda
operator as the deferral path to Phase 4.fem.eig.3.5; NNNNNNNNN's
diagnosis ("the binding constraint at this mesh tier is no longer
modal-sampling discretisation but the 2nd-order Engquist-Majda ABC's
intrinsic floor for off-normal modal content scattered by the
truncation surface") now triggers it.

Phase 4.fem.eig.3.5 ships **CFS-PML** (Complex Frequency Shifted
Perfectly Matched Layer per Roden-Gedney 2000) as a new
`AbcOrder::CfsPml` variant on `OpenBoundarySolver`. CFS-PML is a thin
(~6-10 cell) buffer layer of additional tetrahedra outside the original
cavity volume, in which the constitutive tensor `ε(ω)` is replaced by
the stretched-coordinate form `ε(ω) · Λ(ω)` with `Λ_xx = s_y · s_z /
s_x`, etc., and the complex stretching factor

```text
    s_α(ω)  =  κ_α  +  σ_α / ( α_α + j ω ε_0 ),       α ∈ {x, y, z}.
```

The CFS modification (`α_α > 0`) restores absorption for grazing /
evanescent modes that the original Berenger 1994 PML diverges on, and
maintains finite reflection at DC. Target: |S_11(f)| ∈ [-60, -40] dB
across 8-12 GHz, finally satisfying the spec §8 absorption window.

## 2. Background

The open-boundary truncation problem for FEM electromagnetics has three
historical layers:

- **Berenger 1994** (*J. Comput. Phys.* 114) introduced the split-field
  PML: split the 6-component E/H field into 12 sub-components, each
  damped independently along its propagation axis. The PML region is
  perfectly matched at every angle of incidence in the continuum limit,
  achieving arbitrarily low reflection by stacking more cells. Two
  weaknesses: the split-field formulation is non-Maxwellian inside the
  PML (auxiliary unknowns are not physical fields), and the original
  formulation diverges on **evanescent** modes (surface plasmons,
  grazing modes near cutoff).
- **Kuzuoglu-Mittra 1996** (*IEEE MWCL* 6:12) introduced the
  Complex-Frequency-Shifted (CFS) modification: replace the original
  Berenger stretching factor `s_α(ω) = 1 + σ_α / (j ω ε_0)` with
  `s_α(ω) = κ_α + σ_α / (α_α + j ω ε_0)` for a frequency-shift parameter
  `α_α > 0`. This restores causality, prevents DC reflection growth,
  and absorbs evanescent modes that the original Berenger PML cannot.
- **Roden-Gedney 2000** (*IEEE MWCL* 10:5, *Convolutional PML*) reframed
  CFS-PML as a stretched-coordinate anisotropic material with
  `ε → ε · Λ(ω)`, `μ → μ · Λ(ω)` where `Λ = diag(s_y s_z / s_x,
  s_z s_x / s_y, s_x s_y / s_z)`. The convolution form converts the
  frequency-domain stretching to a time-domain auxiliary-variable ODE
  for FDTD; for **frequency-domain FEM** (this spec) we keep the
  multiplicative `ε · Λ(ω)` form directly — no convolution, no
  auxiliary variables, just a complex anisotropic permittivity tensor
  per PML tet that depends on `ω`. The implementation is a
  straightforward extension of `assemble_tet_element_complex` to
  per-tet 3×3 complex `ε_tensor` matrices.

Contrast with the Phase 4.fem.eig.3 1st/2nd-order Engquist-Majda ABC:
the ABC is a single boundary integral on the truncation surface, exact
only for plane waves at normal incidence on a homogeneous background.
Off-normal modal content, evanescent modes, and corner/edge effects all
reflect at increasing magnitude as the incidence angle deviates from
normal — this is what NNNNNNNNN's `s11_db ≈ −0.01 dB` measurement
exposes. CFS-PML uniquely handles off-normal **and** evanescent modal
content via volumetric absorption in the PML buffer; the truncation
surface outside the PML is then a benign PEC (or 1st-order ABC) whose
reflection contribution is exponentially suppressed by the PML round-trip
attenuation `R_PML ≈ exp(−2 · Σ_α ∫ σ_α dα / (ε_0 c))`.

## 3. Mathematical formulation

### 3.1 Stretched-coordinate complex permittivity tensor

For a Cartesian-aligned PML region with stretching factors `s_x(ω),
s_y(ω), s_z(ω)`, Roden-Gedney 2000 §II gives the equivalent isotropic
material in the un-stretched (laboratory) coordinate frame as the
anisotropic complex tensor

```text
    Λ(ω)  =  diag( s_y(ω) s_z(ω) / s_x(ω),
                   s_z(ω) s_x(ω) / s_y(ω),
                   s_x(ω) s_y(ω) / s_z(ω) ),

    ε_eff(ω)  =  ε(ω) · Λ(ω),       μ_eff(ω)  =  μ(ω) · Λ(ω).
```

The CFS stretching factor per axis is (Kuzuoglu-Mittra 1996 eq. 1;
Roden-Gedney 2000 eq. 4)

```text
    s_α(ω)  =  κ_α(d)  +  σ_α(d) / ( α_α(d)  +  j ω ε_0 ),
```

where `d ∈ [0, D]` is the depth into the PML measured from the inner
boundary (`d = 0`) toward the outer surface (`d = D`). The standard
polynomial grading `σ_α(d) = σ_max · (d/D)^m` with `m ∈ {2, 3, 4}` and
`κ_α(d) = 1 + (κ_max − 1) · (d/D)^m` keeps the wave impedance continuous
at the PML inner boundary and ramps the absorption smoothly to the
outer truncation surface. `α_α` is a small positive constant (typical
range `[ω₀ ε_0 / 10, ω₀ ε_0]` at the band centre); a constant `α_α =
α_max` across the PML is adequate at first approximation.

For axes **outside** the PML (e.g. the `y` and `z` axes for an `x`-face
PML region), `s_α(ω) = 1` identically — the corresponding diagonal
entry of `Λ` collapses to `1` and the PML absorbs only along its own
axis. For **corner and edge PML cells** (where two or three Cartesian
axes are simultaneously inside the PML), all relevant `s_α` are
non-trivial and `Λ` mixes three damping channels — this is what gives
CFS-PML its angle-independence in the continuum limit.

### 3.2 Tet-mesh discretisation of the PML region

For a cuboid cavity, the PML region is a **thin shell of additional tet
layers outside the cavity volume**, one shell per ABC-tagged face plus
the 12 edge wedges and 8 corner wedges where shells meet. The shell
thickness is `thickness_cells` (default 6) tet layers, each of the same
characteristic length as the cavity-interior tets at the boundary. The
new mesh has roughly `(1 + 2 · t/L_x)(1 + 2 · t/L_y)(1 + 2 · t/L_z) − 1`
extra cells relative to the original cavity, with `t = thickness_cells ·
h_cell` and `L_α` the cavity edge length — for a 30 mm WR-90 stub with
`h_cell ≈ 1 mm` and `t = 6 mm`, the PML adds ~5× the cavity cell count.

The PML inner boundary (the original cavity surface) is a smooth
interface between the cavity material and the stretched PML material
with `Λ(d = 0) = I` (no stretching at the inner boundary, by polynomial
grading) — so the bilinear form `a(E, v) = ∫ (∇×E) · μ_eff^{-1} · (∇×v)
− ω² ε_eff · E · v dV` integrates continuously across the inner
boundary with no surface contribution. **The boundary integral on the
PML inner boundary degenerates to zero**, in contrast to the
Engquist-Majda ABC which puts the entire absorption mechanism on the
surface integral. This is the load-bearing structural difference: PML
absorbs in the **volume**; ABC absorbs on the **surface**.

The PML outer boundary (the truncation surface) is tagged
`FaceKind::Pec` by default; the residual reflection off PEC after a
round-trip through the PML is the polynomial-grading-bounded
`R(θ) ≈ exp(−2 cos(θ) · ∫_0^D σ_eff(d) dd / (ε_0 c))` per Berenger
1994 eq. 26 (extended to CFS by Roden-Gedney 2000 eq. 16). With
`σ_max · D / (ε_0 c) ≈ 6` (the standard "PML absorption parameter
~ 6 nepers" rule of thumb), `R(θ = 0) ≈ exp(−12) ≈ 6 × 10^{−6}`, well
below the spec §8 `−45 dB` (≈ 5.6 × 10^{−3}) floor.

### 3.3 PML face classification and per-tet ε tensor

The exterior-face classifier `ExteriorFaceTable::build` from Phase
4.fem.eig.2 is extended with a new boundary classification: a face is
"PML-inner" if it sits on the original cavity surface and was previously
tagged `FaceKind::Abc`; the PML mesh-extension step replaces that ABC
tag with a stretched volumetric layer plus a new outer PEC face.

Per-tet ε tensor assembly: a tet whose centroid lies at depth
`(d_x, d_y, d_z)` into the `(x, y, z)`-PML shell (with `d_α = 0` if
outside the α-PML) gets the complex 3×3 diagonal `ε_tensor(ω) = ε_iso ·
Λ(d_x, d_y, d_z; ω)`. The tet-element bilinear form

```text
    K_tet(ω)_{ij}  =  ∫_T  ( ∇×N_i ) · μ_eff^{-1}(ω) · ( ∇×N_j )  dV,
    M_tet(ω)_{ij}  =  ∫_T  N_i · ε_eff(ω) · N_j  dV,
```

now sees a non-scalar `μ_eff^{-1}, ε_eff`. Whitney-1 basis functions and
their curls remain real-valued; only the prefactor tensor is complex.
The local 6×6 complex per-tet block is computed analytically (the
curl-curl integrand is constant-per-tet in Whitney-1; the mass
integrand is linear-times-linear, integrated via the standard
barycentric formula `∫_T λ_a λ_b dV = V/20 · (1 + δ_{ab})`).

### 3.4 Per-frequency ω dependence

`ε_eff(ω) = ε(ω) · Λ(ω)` depends on `ω` through both the underlying
material (`ε(ω)` from the dispersive Newton tracker for lossy fills) and
the CFS stretching `Λ(ω)`. For a non-dispersive cavity fill the
material part is `ω`-independent and the PML contribution is the only
`ω`-dependent piece; the per-frequency assembly cost in the PML is
identical to the Phase 4.fem.eig.1 dispersive interior-cell cost.

## 4. Tet-mesh PML implementation

### 4.1 New configuration types

```rust
// crates/yee-fem/src/open_boundary.rs (additions)

/// CFS-PML configuration for a single PML axis or composite shell.
#[derive(Clone, Copy, Debug)]
pub struct PmlConfig {
    /// PML shell thickness in tet layers. Default 6.
    pub thickness_cells: usize,
    /// Maximum conductivity (S/m) at the outer truncation surface.
    /// Roden-Gedney 2000 §III recommends `σ_max ≈ (m+1) / (150 π ·
    /// h_cell · sqrt(ε_r))` for an `R(θ=0) ≈ exp(−16)` floor; default
    /// is the band-centre value computed from `freq_hz` + `h_cell` at
    /// solver construction.
    pub sigma_max: f64,
    /// CFS frequency-shift parameter `α_max`. Default
    /// `2 π · f_centre · ε_0` (standard choice from Roden-Gedney 2000
    /// §IV recommending `α_max ≈ ω₀ ε_0`).
    pub alpha_max: f64,
    /// Coordinate-stretching parameter `κ_max`. Default 5.0 (Roden-
    /// Gedney 2000 Table I for waveguide-discontinuity benchmarks).
    pub kappa_max: f64,
    /// Polynomial grading order. Default 3. Values 2, 3, 4 supported.
    pub m: usize,
}

impl Default for PmlConfig {
    fn default() -> Self {
        Self {
            thickness_cells: 6,
            sigma_max: 0.0,  // sentinel; recomputed at solver build
            alpha_max: 0.0,  // sentinel; recomputed at solver build
            kappa_max: 5.0,
            m: 3,
        }
    }
}

/// Designation of which original boundary faces become PML-fronted.
#[derive(Clone, Debug, Default)]
pub struct PmlRegion {
    /// Faces inside the original mesh's `FaceKind::Abc` set that should
    /// be replaced with a PML shell. Empty `faces` is "all ABC faces".
    pub faces: Vec<FaceKind>,
    /// Per-axis PML configuration. If `None`, the solver-level
    /// `PmlConfig` default applies.
    pub config: Option<PmlConfig>,
}
```

### 4.2 PML mesh extension

A new helper in `crates/yee-fem/src/pml_mesh.rs` (new module) extends
the input `TetMesh3D` with the PML buffer layers:

```rust
/// Extend a tet mesh with a CFS-PML buffer shell on every face in
/// `pml_faces`. Returns the extended mesh, a per-tet PML-classification
/// map (`PmlClass::{Interior, PmlX, PmlXY, PmlXYZ, ...}`), and a map
/// from extended-mesh face indices back to original-mesh face indices.
pub fn extend_mesh_with_pml(
    mesh: &TetMesh3D,
    pml_faces: &[FaceKind],
    config: &PmlConfig,
) -> Result<(TetMesh3D, Vec<PmlClass>, FaceIndexMap), Error>;

/// Classification of an extended-mesh tet: which Cartesian-axis PML
/// shells does it lie inside?
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PmlClass {
    Interior,                                // cavity tet, no PML
    PmlX { d: f64 },                         // x-only PML, depth d
    PmlY { d: f64 },                         // y-only PML
    PmlZ { d: f64 },                         // z-only PML
    PmlXY { d_x: f64, d_y: f64 },            // x+y edge wedge
    PmlYZ { d_y: f64, d_z: f64 },
    PmlZX { d_z: f64, d_x: f64 },
    PmlXYZ { d_x: f64, d_y: f64, d_z: f64 }, // corner wedge
}
```

### 4.3 Anisotropic per-tet element assembly

A new helper in `crates/yee-fem/src/element.rs`:

```rust
/// Per-tet complex 6×6 stiffness + mass block with anisotropic
/// complex ε tensor (Roden-Gedney 2000 §II CFS-PML). For non-PML cells,
/// pass `eps_tensor = eps_scalar · I` and `mu_tensor = mu_scalar · I`
/// and the result matches `assemble_tet_element_complex` bit-for-bit.
pub fn assemble_tet_element_complex_anisotropic(
    tet_vertices: [Vector3<f64>; 4],
    eps_tensor: SMatrix<Complex64, 3, 3>,
    mu_tensor_inv: SMatrix<Complex64, 3, 3>,
    omega: f64,
) -> (SMatrix<Complex64, 6, 6>, SMatrix<Complex64, 6, 6>);
```

For a diagonal `Λ(ω) = diag(λ_x, λ_y, λ_z)` (the standard Cartesian-
aligned PML case) the integrand factorises and the cost is bounded at
~3× the scalar-`ε` Phase 4.fem.eig.1 cell-assembly cost. Off-diagonal
`Λ` (rotated PML axes; out of scope for v3.5) would force a full
double-sum over the 6×6 block.

## 5. Public API

```rust
//! crates/yee-fem/src/open_boundary.rs (extensions)

/// Selects the open-boundary truncation kernel.
#[derive(Clone, Debug, PartialEq, Default)]
pub enum AbcOrder {
    /// 1st-order Engquist-Majda (Phase 4.fem.eig.2 default).
    First,
    /// 2nd-order Engquist-Majda (Phase 4.fem.eig.3 default).
    #[default]
    Second,
    /// CFS-PML (Phase 4.fem.eig.3.5). The PML config travels in
    /// the variant payload.
    CfsPml(PmlConfig),
}

impl<'m> OpenBoundarySolver<'m> {
    // ... existing v3 methods unchanged ...

    /// Layer a CFS-PML shell on every `FaceKind::Abc`-tagged face.
    /// Mutually exclusive with `with_abc_order(First|Second)` —
    /// CFS-PML *replaces* the surface-integral ABC kernel with a
    /// volumetric PML buffer. Default `PmlConfig` is the
    /// Roden-Gedney 2000 §III/IV recommended set.
    pub fn with_cfs_pml(self, config: PmlConfig) -> Self;
}
```

The CFS-PML path co-exists with 1st/2nd-order ABC on a per-face basis:
a future variant `AbcOrder::Hybrid` could mix PML on some faces with
2nd-order ABC on others; v3.5 ships only the all-PML mode.

Python binding (plan step P6):

```python
yee.fem.solve_open_cavity(
    mesh, materials, port_faces, abc_faces, omegas,
    *,
    pml_config: dict | None = None,   # CFS-PML; if None, falls back to
                                       # `abc_order` (default Second)
    abc_order: str = "second",
    coupled_whitney: bool = True,
    multi_port: bool = False,
) -> np.ndarray
```

`pml_config` keys mirror `PmlConfig` fields. When `pml_config` is set,
`abc_order` is ignored.

## 6. Validation gates

### fem-eig-003 strict — un-ignore both gates

The two `#[ignore]`d tests in
`crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs`
(`fem_eig_003_strict_absorption_floor_gate` and
`fem_eig_003_strict_passive_bound_continuum_limit`) become CI-default.
Driver flags: `coupled_whitney = true`,
`abc_order = AbcOrder::CfsPml(PmlConfig::default())`. Mesh: cavity is
the existing `(24, 12, 36)` Kuhn 6-tet brick; the PML mesh extension
adds a 6-cell shell on every original ABC face — for fem-eig-003 that
is the single `z = 0` face, so the extended mesh adds `24 × 12 × 6 ·
6 = 10 368` PML tets (totalling ~72 k tets).

Gate criteria (replacing the Phase 4.fem.eig.3 fem-eig-003-strict gates):

1. **Strict absorption floor (new bounds)** — `20·log10|S_{11}(f)| ∈
   [-60, -40] dB` at every swept frequency in 8-12 GHz. (Slightly
   widened on the lower bound from the spec §8 `-45 dB` floor because
   CFS-PML overachieves; the upper bound `-40 dB` is the meaningful
   "PML works" assertion.)
2. **Strict passive bound** — `|S_{11}(f)| < 1` strictly.
3. **No PEC-resonance peak** at TE_{101} ≈ 9.66 GHz.
4. **Phase monotonicity** preserved.
5. **Run-time < 480 s in `--release`** (informational; ~2× the v3 fem-
   eig-003 budget because of the ~17 % extra tets).

### fem-eig-006 — new high-aspect-ratio cavity stress test

A new fixture stress-tests PML stability on geometries that the Phase
4.fem.eig.3 2nd-order ABC degrades catastrophically on: a high-aspect
**100 : 10 : 1** rectangular cavity (dimensions `100 mm × 10 mm ×
1 mm`) with TE-mode drive at one short face (`x = 0`) and CFS-PML on
the opposite short face (`x = 100 mm`). The four sidewalls are PEC.
This geometry forces highly off-normal modal content onto the PML
inner boundary — precisely the regime where ADR-0042 §risks predicted
the 2nd-order ABC's intrinsic-floor degradation.

Driver flags: `coupled_whitney = true`,
`abc_order = AbcOrder::CfsPml(PmlConfig::default())`.

Gate criteria at 30 GHz (single frequency, no sweep — this is a
stability fixture, not an absorption-band fixture):

1. **Magnitude bounded** — `|S_{11}(30 GHz)| < 0.1` (PML absorbs the
   off-normal incidence; with 2nd-order ABC this fixture saturates at
   `|S_{11}| ≈ 0.95`).
2. **No NaN / Inf** in the swept output — the surface-mode divergence
   risk (§7) is the failure mode this gate catches.
3. **LU condition number bounded** — the per-frequency `cond(A(ω)) <
   1e10`. CFS-PML's complex-anisotropic `ε_tensor` breaks the
   complex-symmetric structure (§7) and conditions degrade — this gate
   ensures the degradation is bounded.

## 7. Risks and open questions

- **PML stability for surface modes / evanescent modes.** Roden-Gedney
  2000 §IV explicitly notes that CFS-PML with `α_α = 0` (the
  un-modified Berenger PML) *diverges* on evanescent modes excited by
  high-aspect-ratio geometries — the well-known "surface plasmon
  blow-up" of split-field PML. The CFS `α_α > 0` modification absorbs
  evanescent modes by adding a real-axis pole to `s_α(ω)`, mitigating
  this. The fem-eig-006 high-aspect fixture is specifically designed
  to catch a regression in this mitigation.
- **PML thickness vs accuracy trade-off.** `thickness_cells = 6` is
  the Roden-Gedney 2000 §III default ("6 to 10 cells suffices for
  microwave applications"); 4 cells is borderline and 10+ cells is
  wasteful. fem-eig-006 ablation: if `thickness_cells = 4` already
  meets gate 1, document it; if `thickness_cells = 10` is required,
  the default should be revisited in Phase 4.fem.eig.3.5.1.
- **Complex anisotropic `ε_tensor` breaks complex-symmetric stiffness.**
  Phase 4.fem.eig.1+ relies on the assembled `K(ω) − k₀² M(ω)` being
  *complex-symmetric* (real-symmetric Gram structure with complex
  scalar prefactors), which `faer::sparse::FaerLuSolver<Complex64>`
  factorises via complex LDLᵀ. With per-tet anisotropic `ε_tensor`,
  the assembled matrix is complex-symmetric **only when `Λ(ω)` is
  diagonal in the global frame** (Cartesian-aligned PML). v3.5 ships
  Cartesian-aligned PML only, so the LDLᵀ path is preserved. A future
  rotated-PML extension would degrade to full complex LU — same
  surface, ~2× factorisation cost. ADR-0043 records the choice; the
  PML mesh-extension helper rejects non-Cartesian-aligned PML faces
  with a `NotEnabled` error.
- **PML grading parameter sensitivity.** `σ_max`, `α_max`, `κ_max`,
  and `m` all interact non-linearly. Roden-Gedney 2000 Table I gives
  benchmark-validated defaults for microwave waveguide cases; we
  inherit those. If fem-eig-003-strict or fem-eig-006 measurements
  drift outside the `[-60, -40] dB` band, the first knob to retune is
  `sigma_max` (linear control on PML round-trip attenuation).
- **PML inner-boundary continuity.** Polynomial grading `σ_α(d=0) = 0`,
  `κ_α(d=0) = 1` ensures the PML inner-boundary face has identical
  material parameters on both sides — no spurious surface reflection
  from a material discontinuity. The mesh-extension helper verifies
  this by asserting `Λ(d_α = 0) = I` at every inner-boundary face.

## 8. Phase numbering ladder

- Phase 4.fem.eig.3 — coupled-Whitney + 2nd-order ABC + multi-port
  (shipped F1-F7 + NNNNNNNNN mesh refinement).
- **Phase 4.fem.eig.3.5 — this spec**: CFS-PML replacing the
  Engquist-Majda surface-integral ABC.
- Phase 4.fem.eig.3.5.1 — CFS-PML parameter ablation (`thickness_cells
  ∈ {4, 6, 8, 10}`, `m ∈ {2, 3, 4}`, `α_max` sweep). Optional follow-up.
- Phase 4.fem.eig.4 — FEM-BEM hybrid; CFS-PML is the FEM-side
  truncation, BEM handles infinite half-spaces.

## 9. Lane

Spec file: `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-cfs-pml-design.md`

Implementation lane (declared here for the follow-up plan, not edited
by this spec):

- `crates/yee-fem/src/open_boundary.rs` — extend `AbcOrder` with
  `CfsPml(PmlConfig)`, add `with_cfs_pml`, add `PmlConfig`, `PmlRegion`.
- `crates/yee-fem/src/pml_mesh.rs` *(create)* — `extend_mesh_with_pml`,
  `PmlClass`, `FaceIndexMap`.
- `crates/yee-fem/src/element.rs` — add
  `assemble_tet_element_complex_anisotropic` (D1+D3+D4 path).
- `crates/yee-fem/src/lib.rs` — re-export new types.
- `crates/yee-fem/tests/pml_mesh_extension.rs` *(create)* — unit tests.
- `crates/yee-fem/tests/anisotropic_tet_assembly.rs` *(create)* —
  per-tet anisotropic-`ε` assembly bit-for-bit equivalence to scalar
  path when `ε_tensor = ε · I`.
- `crates/yee-validation/src/lib.rs` — add `run_fem_eig_006_high_aspect`
  driver; update `run_fem_eig_003_wr90_stub_abc` flags for PML.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` —
  un-ignore both strict gates.
- `crates/yee-validation/tests/fem_eig_006_high_aspect_pml.rs`
  *(create)*.
- `crates/yee-py/src/fem.rs` — `pml_config` kwarg.
- Out-of-lane: `yee-cli`, `yee-gui`, `yee-mom`, `yee-mesh`, `yee-cuda`,
  `yee-plotters`.

## 10. References

- Berenger, J.-P., "A perfectly matched layer for the absorption of
  electromagnetic waves", *J. Comput. Phys.* 114 (1994), pp. 185-200.
  DOI 10.1006/jcph.1994.1159. The original split-field PML; the
  reference this spec amends with CFS modifications.
- Kuzuoglu, M. and Mittra, R., "Frequency dependence of the
  constitutive parameters of causal perfectly matched anisotropic
  absorbers", *IEEE Microwave and Guided Wave Letters* 6(12) (1996),
  pp. 447-449. DOI 10.1109/75.541428. The CFS modification (`α_α > 0`)
  that retires the evanescent-mode divergence of Berenger 1994.
- Roden, J. A. and Gedney, S. D., "Convolutional PML (CPML): An
  efficient FDTD implementation of the CFS-PML for arbitrary media",
  *IEEE Microwave and Wireless Components Letters* 10(5) (May 2000),
  pp. 27-29. DOI 10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A.
  Stretched-coordinate anisotropic-material formulation;
  recommended-parameter Tables I and II.
- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd
  ed., Wiley 2014, §10.4 (ABC reflection floors; Table 10.1), §10.8
  (PML for FEM, stretched-coordinate form).
- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md` — Phase
  4.fem.eig.3 spec (parent).
- `docs/src/decisions/0042-phase-4-fem-eig-3-scope.md` — Phase
  4.fem.eig.3 scope ADR; §risks queues this spec.
- `docs/src/decisions/0043-phase-4-fem-eig-3-5-cfs-pml-scope.md` — this
  spec's scope ADR.
- `crates/yee-validation/tests/fem_eig_003_wr90_stub_abc.rs` §"NNNNNNNNN
  status" — the `[-2.22e-2, -2.86e-5] dB` measurement that motivates
  the CFS-PML upgrade.
