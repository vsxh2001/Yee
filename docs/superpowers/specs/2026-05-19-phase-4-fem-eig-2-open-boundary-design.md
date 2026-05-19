# Phase 4.fem.eig.2 — Open-Boundary FEM (ABC + Wave Ports)

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.2 (1st-order Engquist–Majda ABC + modal wave-port driven analysis); 4.fem.eig.2.5 deferred (see §13)
**Depends on:** Phase 4.fem.eig.0 (real-valued closed-cavity FEM eigensolver, shipped), Phase 4.fem.eig.1 (lossy dispersive `ε(ω)` Newton tracker, shipped), Phase 1.3.1.1 (2-D Nedelec cross-section eigensolver, shipped)
**Blocks:** Phase 4.fem.eig.3 (DRA validation against published radiating cases), Phase 4 driven-FEM extension into open antennas

## 1. Motivation

Phases 4.fem.eig.0 and 4.fem.eig.1 ship a clean closed-cavity FEM eigensolver,
respectively in the lossless real-`ε_r` and lossy single-pole dispersive
`ε_r(ω)` regimes. Gate `fem-eig-001` clears TE_{101} on a WR-90-based metallic
box at 0.09 % rel.err., and `fem-eig-002` clears the lossy SiO₂-filled cavity
at 1.3e-3 on Re(f) and 3e-3 on Im(f). Both solvers enforce PEC tangential-`E`
zero on the *entire* boundary via Dirichlet row/column elimination — the only
boundary the FEM stack knows how to apply.

Real microwave / antenna / accelerator workloads need open boundaries:

- Driven `S_{11}(f)` analysis of waveguide-fed stubs, irises, and cavity
  filters — the input port is *not* a PEC wall; it carries an incident TE_{10}
  (or higher-order) modal current and an outgoing reflection.
- Radiation problems where the FEM domain truncates free space: slot
  antennas, fed dipoles enclosed in a small FEM box, dielectric-resonator
  antennas with a finite air halo. The terminating surface needs to *absorb*
  outgoing waves to whatever floor the ABC permits.
- Coaxial / SMA port modelling for board-level circuit pieces: a coax line
  exiting the FEM domain looks identical to a wave-port with a TEM modal
  profile.
- Driven Q-extraction on the same lossy cavities Phase 4.fem.eig.1 covers,
  but now with *external loading* via a coupling iris on a wave-port face
  instead of a bulk-material loss term.

Closed-cavity eigenmode is the wrong tool for any of these. The natural Phase
4.fem.eig.2 deliverable adds:

1. A **1st-order Engquist–Majda absorbing boundary condition** (ABC) — the
   minimal viable radiation-condition surrogate that keeps the curl-curl
   bilinear form well-posed on an open boundary, at the cost of a `~ −40 dB`
   reflection floor for normal incidence.
2. A **modal wave-port face**, where the tangential field is decomposed into
   the cross-section's Nedelec eigenmodes from Phase 1.3.1.1 weighted by an
   incident-amplitude vector. The forward sweep returns `S` parameters
   directly from the modal weights.

The mathematics is textbook (Jin §10 "Driven FEM analysis", Pozar §3.3
"Waveguide modes"). The engineering risk is the boundary-term bookkeeping
inside the existing `assemble_tet_element_complex` element layer and the
re-use of `NumericalCrossSection` from `yee-mom`'s Phase 1.3.1.1 stack as a
mode-profile source on FEM port faces — neither involves new physics, only
careful surface-integral bookkeeping per face.

ADR-0040 records the Engquist–Majda choice and the deferral of higher-order /
PML termination to Phase 4.fem.eig.2.5.

## 2. Non-goals (Phase 4.fem.eig.2)

Explicitly out of scope for v2:

- **2nd-order Engquist–Majda ABC and Higdon ABC.** Tighter than 1st-order at
  the cost of additional auxiliary equations on the boundary; deferred to
  Phase 4.fem.eig.2.5 if `fem-eig-003` hits the 1st-order reflection floor.
- **Berenger PML on FEM.** The split-field PML formulation does not transfer
  cleanly to Nedelec curl-conforming elements; uniaxial PML (UPML) and
  complex-coordinate stretching (CFS-PML) do, but only one of those lands
  in 4.fem.eig.2.5 — never v2.
- **Multi-mode incident excitation on a port.** v2 supports a single incident
  mode per port (typically TE_{10} for rectangular waveguide; TEM for coax).
  Higher-order modes are *captured* in the modal-decomposition reflection
  spectrum — they are not *driven*.
- **FEM-BEM hybrid for finite-aperture radiation.** The ABC is the truncation
  surface for v2; a BEM-coupled exterior is Phase 4.fem.eig.4+.
- **Adaptive mesh refinement on the ABC face.** Mesh is an input contract;
  the caller is responsible for ABC-face sizing per Jin §10.4 (`h_ABC ≤ λ /
  20` at the highest swept frequency).
- **Periodic / Floquet BCs on the ABC face.** Phase 4.fem.eig.4+.
- **Driven dispersive lossy media via the Phase 4.fem.eig.1 Newton tracker.**
  v2 supports lossless or constant-real-loss interior media on the driven
  sweep; combining wave-port drive with Newton-tracked dispersive media is a
  Phase 4.fem.eig.2.1 superposition exercise.
- **GPU acceleration.** CPU/FP64 complex scalar. The complex sparse LU
  shipped in Phase 4.fem.eig.1 is reused at each sweep frequency.
- **Higher-order Nedelec on port or ABC faces.** First-order Whitney-1 only.

## 3. Scope decision — Phase 4.fem.eig.2 open-boundary v0

Walking-skeleton-extend-once per `CLAUDE.md` §3. The v2 deliverable extends
the shipped Phase 4.fem.eig.1 dispersive eigensolver along two axes:

- **1st-order Engquist–Majda ABC on a designated set of mesh faces.** Per
  spec §4 derivation: on a face with outward normal `n̂`, the curl-curl
  bilinear form picks up a *boundary term*

  ```text
      + j k₀  ∫_face  (n̂ × N_i) · (n̂ × N_j) dS
  ```

  which is added per-face into the global stiffness matrix `K(ω)`. The
  boundary term is **complex-valued for real `ε_r`** — adding an ABC face
  promotes the closed-cavity eigenproblem from real to complex even with no
  material loss; this is the same mathematical fact that lets the ABC absorb
  outgoing waves (the imaginary part is the radiation resistance).
- **Modal wave-port faces.** Per spec §4.3: on a port face, the tangential
  electric field is parameterised as

  ```text
      E_t  =  a_inc · e_mode(x,y) · e^{−j β z}  +  Σ_n  b_n · e_mode,n(x,y) · e^{+j β_n z}
  ```

  where `e_mode(x,y)` is the dominant-mode tangential profile sourced from
  Phase 1.3.1.1's `NumericalCrossSection::e_tangential_at` (the same surface
  the MoM wave-port already consumes), and the unknowns `b_n` are the modal
  reflection amplitudes that compose the `S_{11}(f)` row of the scattering
  matrix. A right-hand side is added to the driven system encoding the
  incident wave.

Single end-to-end gate: **fem-eig-003 WR-90 stub with ABC termination** (§9).
A finite-length 22.86 × 10.16 × 30 mm air-filled rectangular waveguide
section is closed at `z = 0` with PEC and terminated at `z = d` with a 1st-
order Engquist–Majda ABC. A TE_{10} wave-port at `z = d` drives the system;
the modal decomposition extracts `S_{11}(f)` across an 8–12 GHz sweep. The
gate asserts `|S_{11}(f)|` agrees with the Pozar §3.3 closed-form reference
(equivalently, a well-terminated PEC stub at TE_{10} has `|S_{11}| ≈ −60 dB`
mid-band and rises to the ABC reflection floor at the band edges) within
**±0.5 dB across 8–12 GHz**.

CPU, FP64, scalar complex. No GPU. The Phase 4.fem.eig.0 closed-cavity and
Phase 4.fem.eig.1 dispersive paths stay green and unchanged: an
`OpenBoundarySolver` consumes a mesh with optional `port_faces` and
`abc_faces` annotations, and an empty face list degrades to the closed-
cavity eigenproblem (zero boundary contribution, real or complex depending on
the material model).

What v2 does **not** ship, deferred to 4.fem.eig.2.5+ (see §13):

- 2nd-order / Higdon / UPML / CFS-PML ABC variants.
- Multi-mode incident excitation.
- FEM-BEM hybrid for fully-radiating apertures.
- Driven sweep over dispersive `ε(ω)` with the Phase 4.fem.eig.1 Newton tracker.

## 4. Theory anchor

The starting point is the *frequency-domain* vector wave equation from
Phase 4.fem.eig.0/1 §4:

```text
    ∇ × ( (1/μ_r) ∇ × E )  −  k₀² ε_r E   =  0       in  Ω
                                      k₀   =  ω / c
```

For closed cavities the boundary `∂Ω` carried tangential-PEC (`n̂ × E = 0`),
contributing nothing to the variational form via the eigenproblem's natural
Neumann boundary. Phase 4.fem.eig.2 partitions `∂Ω` into three disjoint
pieces:

```text
    ∂Ω  =  Γ_PEC  ∪  Γ_ABC  ∪  Γ_port,
```

each with its own boundary-term contribution.

### 4.1 Variational form with surface terms

Multiplying the wave equation by a test field `v ∈ H(curl; Ω)` and
integrating by parts (Jin §10.2, eq. 10.18) gives

```text
    ∫_Ω  (1/μ_r) (∇×E) · (∇×v) dV
        −  k₀²  ∫_Ω  ε_r  E · v  dV
        +  ∮_∂Ω  (n̂ × (1/μ_r ∇×E))  ·  v  dS
        =  0.
```

The surface integral over `Γ_PEC` vanishes for `v ∈ H₀,PEC(curl; Ω)`. The
remaining surface contributions are the ABC and port terms.

### 4.2 1st-order Engquist–Majda ABC

The Engquist–Majda 1977 radiation condition on a planar surface with outward
normal `n̂` reads

```text
    n̂ × ∇×E   =   −j k₀  n̂ × (n̂ × E)         on  Γ_ABC.
```

(Engquist, B. and Majda, A., "Absorbing boundary conditions for the numerical
simulation of waves", *Math. Comp.* 31 (1977), pp. 629–651.) Substituting
into the variational form's surface integral yields the per-face stiffness
contribution

```text
    K_ABC^{e,face}_{ij}  =  + j k₀  ∫_face  (1/μ_r,face)  (n̂ × N_i) · (n̂ × N_j)  dS.
```

This is a `3 × 3` block per triangular face (one DoF per face edge), added
into the global `K(ω)` at the corresponding global edge indices with the
correct local-to-global orientation flip from the v0 assembly. The factor
`j k₀` is **purely imaginary** — adding ABC faces makes the global `K`
complex-symmetric (not Hermitian) even when `ε_r`, `μ_r` are real. This is
the canonical "radiation absorbs energy → eigenvalues acquire negative
imaginary part → physically meaningful Q" identity.

The integral is computed exactly: `N_i × n̂` is constant per face (the
Whitney-1 edge basis restricted to a face is a constant tangent vector), so
the surface integral over a flat triangular face is `face_area · (n̂ × N_i) ·
(n̂ × N_j)` with the cross products evaluated once per face.

### 4.3 Modal wave-port

On `Γ_port` the tangential field is constrained to lie in the span of a
finite set of port modes computed from the cross-section eigensolver. For a
single-mode TE_{10} port:

```text
    E_t(x,y,z_port)  =  (a_inc + b)  ·  e_mode(x,y)
    H_t(x,y,z_port)  =  (a_inc − b)  ·  h_mode(x,y),
```

where `e_mode`, `h_mode` are the dominant-mode tangential profiles
(orthonormalised so `∫_port (e_mode × h_mode^*) · ẑ dS = 1`) and `b` is the
modal reflection coefficient. The right-hand side added to the driven
system per port face is

```text
    b_port,i  =  +  2 j β_mode  ·  a_inc  ·  ∫_face  N_i · e_mode  dS,
```

and the corresponding stiffness contribution mirrors the ABC form with the
mode-dependent propagation constant `β_mode = √(k₀² ε_r μ_r − k_c²)`:

```text
    K_port^{e,face}_{ij}  =  + j β_mode  ∫_face  (n̂ × N_i) · (n̂ × N_j) dS.
```

After solving `(K(ω) + j k₀ B_ABC + j β B_port) e = b_port`, the modal
reflection coefficient is extracted via

```text
    b   =   2  ·  ⟨ E_FEM,t , e_mode ⟩_port   −   a_inc,
    S_{11}(f)  =  b / a_inc.
```

The orthonormality of `e_mode` makes the inner product reduce to a per-face
quadrature against the FEM solution's edge coefficients on the port face.

`e_mode(x,y)` is sourced from `yee_mom::eigensolver::NumericalCrossSection`
(Phase 1.3.1.1) via the existing `e_tangential_at` API. The mesh-coupling
contract: the FEM mesh's port-face triangulation does **not** need to be
geometrically identical to the cross-section eigensolver's 2-D triangulation
— `e_mode` is sampled per-face-centroid (or per-Gauss point) using the
analytic projection in `e_tangential_at`. For v0 this is `nearest`-style
interpolation; cubic or barycentric interpolation lands in 4.fem.eig.2.0.1
if the gate proves sensitive.

References for §4:

- Engquist, B. and Majda, A., "Absorbing boundary conditions for the
  numerical simulation of waves", *Math. Comp.* 31 (1977), pp. 629–651 —
  the canonical 1st- and 2nd-order ABC derivation.
- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014, Ch. 10 (driven FEM analysis), §10.4 (ABC face contributions),
  §10.5 (wave-port modal decomposition).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §3.3
  (waveguide port modal characterisation, TE/TM cutoff, propagation
  constants).
- Phase 1.3.1.1 spec: `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`
  — modal profile source.

## 5. Element-layer changes

The Phase 4.fem.eig.1 `assemble_tet_element_complex` is unchanged. The
element layer gains **one new function** for the ABC face block:

```rust
// crates/yee-fem/src/element.rs — new helper

/// Per-face Engquist–Majda ABC contribution at angular frequency ω.
///
/// Returns the `3 × 3` Hermitian-symmetric face block whose entries are
/// `+ j k₀ · (1/μ_r) · area · (n̂ × N_i) · (n̂ × N_j)` for `i, j ∈
/// {edge_0, edge_1, edge_2}` on the face. `vertices[0..3]` are the three
/// face corners in CCW order with respect to the outward normal `n̂`,
/// derived from the parent tet so the orientation is consistent with the
/// existing global edge map.
pub fn assemble_abc_face_block(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    k0: f64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3>;
```

A peer `assemble_port_face_block` returns the analogous block with `j β`
instead of `j k₀` and an associated `port_rhs_block: SMatrix<Complex64, 3, 1>`
encoding the incident modal current per edge — both consumed by the
assembly layer in §6.

The `assemble_tet_element_complex` signature is extended with an optional
`abc_faces: &[FaceId]` argument; when non-empty, the function adds the
face-block contributions to the returned local stiffness matrix at the
correct local-edge indices. v0 closed-cavity callers pass `&[]` and the
function is byte-identical to its Phase 4.fem.eig.1 form.

## 6. Public API surface

```rust
//! crates/yee-fem/src/open_boundary.rs  (new module)

use yee_core::Error;

/// Tag classifying a mesh face for open-boundary FEM.
#[derive(Clone, Copy, Debug)]
pub enum FaceKind {
    Pec,                 // tangential-E-zero Dirichlet (default for unannotated faces)
    Abc,                 // 1st-order Engquist–Majda absorbing boundary
    WavePort(PortId),    // modal wave-port with a NumericalCrossSection mode source
}

/// Wave-port descriptor.
///
/// One `WavePortFace` per physical port face. The port's modal profile is
/// sourced from a pre-computed [`NumericalCrossSection`] (Phase 1.3.1.1)
/// and sampled per-Gauss point at assembly time.
pub struct WavePortFace {
    pub face_ids: Vec<FaceId>,
    pub cross_section: NumericalCrossSection,
    pub mode_index: usize,           // typically 0 for TE_{10}
    pub incident_amplitude: Complex64,
}

/// Open-boundary FEM driven solver.
///
/// Consumes a [`TetMesh3D`] with per-face [`FaceKind`] tags and a swept
/// frequency list; returns the S-parameter matrix indexed by port.
pub struct OpenBoundarySolver<'m> {
    mesh: &'m TetMesh3D,
    material_db: MaterialDatabase,
    face_kinds: Vec<FaceKind>,
    ports: Vec<WavePortFace>,
    abc_faces: Vec<FaceId>,
}

impl<'m> OpenBoundarySolver<'m> {
    pub fn new(
        mesh: &'m TetMesh3D,
        material_db: MaterialDatabase,
        face_kinds: Vec<FaceKind>,
        ports: Vec<WavePortFace>,
    ) -> Result<Self, Error>;

    /// Solve `(K(ω) + j k₀ B_ABC + j β B_port) e = b_port` at a single
    /// frequency and return the S-parameter row (one entry per port).
    pub fn solve_at_frequency(
        &self,
        omega: f64,
    ) -> Result<SParameterRow, Error>;

    /// Frequency-sweep driven solve.  Returns the full `SParameters` matrix
    /// across the swept band, one row per swept frequency.
    pub fn sweep(
        &self,
        omegas: &[f64],
    ) -> Result<SParameters, Error>;
}

/// One row of the S-parameter matrix at a single frequency.
pub struct SParameterRow {
    pub omega: f64,
    pub s: DMatrix<Complex64>, // [n_ports × n_ports]
}

/// Full frequency-swept S-parameter matrix.
pub struct SParameters {
    pub omegas: Vec<f64>,
    pub rows: Vec<SParameterRow>,
}
```

The Phase 4.fem.eig.0/1 surface (`FemEigenAssembly`, `DispersiveSolver`) is
strictly unchanged; an `OpenBoundarySolver` constructed with an empty port
list and an empty ABC list is mathematically equivalent to a closed-cavity
driven solve with no excitation (and is rejected as ill-posed).

Optional Python binding (plan step E6):

```python
yee.fem.solve_open_cavity(
    mesh: TetMesh3D,
    materials: list[Material],
    port_faces: list[dict],         # face IDs + mode index + amplitude
    abc_faces: list[int],           # face IDs
    omegas: np.ndarray,
) -> np.ndarray                     # shape (n_omegas, n_ports, n_ports)
```

mirroring the existing `yee.fem.solve_cavity_dispersive` entry from Phase
4.fem.eig.1.

## 7. Complex sparse linear solve

Phase 4.fem.eig.1 shipped `ComplexInverseIterEigen` over `faer::sparse::
FaerLuSolver<Complex64>` for the *eigen*problem. Phase 4.fem.eig.2's driven
solve replaces the inverse-iteration outer loop with a **single complex
sparse LU back-substitution per swept frequency**:

```text
    (K(ω) + j k₀ B_ABC(ω) + Σ_p j β_p B_port,p(ω))  ·  e   =   b_port(ω).
```

The same `faer::sparse::FaerLuSolver<Complex64>` is consumed — the matrix is
complex-symmetric (not Hermitian) but `faer` already handles complex-symmetric
LU via standard pivoting. No new solver dep.

Per swept frequency, one factorisation + one back-substitution dominates the
runtime; a 50-point band sweep on the fem-eig-003 mesh (~4 k DoFs) budgets
at ~60 s in `--release` per spec §9.

## 8. Validation gate — fem-eig-003 WR-90 stub with ABC termination

Canonical driven-FEM smoke per `CLAUDE.md` §4. Air-filled WR-90 rectangular
waveguide section with:

- `a = 22.86 mm` (broad wall)
- `b = 10.16 mm` (narrow wall)
- `d = 30 mm` (axial length)

Boundary partition:

- `Γ_PEC`: the four longitudinal side walls (top, bottom, both narrow
  sidewalls).
- `Γ_port`: the face at `z = d`. TE_{10} modal source sampled from
  `NumericalCrossSection` solved on the cavity cross-section at the swept
  frequency. Incident amplitude `a_inc = 1`.
- `Γ_ABC`: the face at `z = 0`. 1st-order Engquist–Majda absorbing boundary.

This is *not* the closed-stub gate fem-eig-001 covered; the closed-stub
analytic resonance at 9.66 GHz would show up as `|S_{11}| = 0 dB` (full
reflection) if the `z = 0` face were PEC. With the ABC at `z = 0` the wave
propagates through and is absorbed, modulo the 1st-order Engquist–Majda
reflection floor (`~ −40 dB` for normal incidence at TE_{10}). The expected
`|S_{11}(f)|` across 8–12 GHz is:

```text
    |S_{11}(f)|  ≈  R_ABC(f)  ≈  −40 dB    for  f  >  f_c(TE_{10})  ≈  6.56 GHz,
```

with weak frequency dependence in the propagating band 8–12 GHz. The
Pozar §3.3 closed-form modal reference for a uniformly-terminated WR-90
section in the dominant-mode band is `|S_{11}(f)| = R_ABC(f)` (the ABC is
the only reflector). Below cutoff the gate is undefined and the sweep
excludes `f < 7 GHz`.

**Gate criteria:**

1. **`|S_{11}(f)|` across the sweep** matches the analytic Pozar §3.3
   reference within **±0.5 dB at every swept frequency in 8–12 GHz** (50
   points, 80 MHz spacing).
2. **No PEC-resonance peak** appears at the closed-stub TE_{101} = 9.66 GHz
   in the swept band. The ABC must absorb the would-be standing wave; a
   visible peak with `> 3 dB` excursion indicates the ABC face block is not
   wired correctly.
3. **Phase of `S_{11}(f)`** is monotonic across the sweep (no aliasing or
   sign flip). The Phase 1.3.1.1 modal source enforces a consistent reference
   plane at `z = d`.
4. **Run-time `< 180 s` in `--release`** on the standard fem-eig-003 mesh
   (~25 k tets after Kuhn decomposition, ~4 k interior DoFs). Informational;
   not a CI gate.
5. Standard verification chain green: `cargo build`, `cargo clippy
   -- -D warnings`, `cargo test --release`, `cargo fmt --check` on the
   touched crates.

A second cross-check fixture (E5 unit test) drives a TE_{10} wave at 10 GHz
into a PEC-closed WR-90 stub (no ABC) and confirms `|S_{11}| = 0 dB ± 0.05`
— this is the "ABC-off" version of the gate and isolates the modal
projection from the ABC face block.

## 9. Higher-applications roadmap

Beyond fem-eig-003, open-boundary FEM unlocks:

- **fem-eig-004 — coax-fed dipole inside an ABC-terminated FEM box.** A
  thin-wire half-wave dipole at 1 GHz in a small (10 cm)³ air box with all
  six exterior faces tagged `Abc`. Compare driving-point impedance against
  the mom-001 NEC-4 reference (`Z ≈ 87 + j41 Ω`); the FEM result should
  agree to within the ABC floor — likely ~5 % on Re(Z), ~10 % on Im(Z) at
  the 1st-order ABC's reflection floor.
- **Slot-antenna and iris-coupling cases** from Pozar §6.4: drive the slot
  via a TE_{10} wave-port, terminate the radiating side with an ABC face,
  extract the slot impedance and compare to Pozar's tabulated closed-form.
- **Lossy dispersive cavity filter** combining Phase 4.fem.eig.1's Newton
  tracker (interior `ε(ω)`) with v2's wave-port drive — sweep `S_{21}(f)`
  through a two-port iris-coupled bandpass and compare to Pozar §8.4.
  Combining the two paths is a Phase 4.fem.eig.2.1 superposition exercise;
  the assembly layer already supports it.

These are deferred to Phase 4.fem.eig.3+ (§13).

## 10. Risks and open questions

- **1st-order Engquist–Majda reflection floor (`~ −40 dB` at normal
  incidence; worse off-normal).** This is the fundamental v0 limit. If
  `fem-eig-003` shows the swept `|S_{11}(f)|` clipping at `−40 dB` instead of
  the expected `−45 dB` or better, the gate is reported green-with-finding
  and 4.fem.eig.2.5 (2nd-order ABC or CFS-PML) becomes the next phase. The
  ±0.5 dB tolerance window in §8 is deliberately wide enough to absorb the
  ABC floor's frequency dependence within the dominant-mode band.
- **Modal projection on FEM-vs-MoM port mesh mismatch.** Phase 1.3.1.1's
  `NumericalCrossSection` mesh and the FEM port-face triangulation are
  geometrically independent. v0 samples `e_mode` per-face-centroid (or per
  Gauss point); if `fem-eig-003` shows mode-projection error > ±0.2 dB,
  upgrade to cubic interpolation in 4.fem.eig.2.0.1. Mitigation: cross-check
  with the analytic TE_{10} profile `e_mode(x,y) = ŷ sin(π x / a)`, which
  has a closed form for the WR-90 air-filled case — discrepancy between the
  numerical and analytic modal projections quantifies the interpolation
  error.
- **Port reference plane and phase consistency.** The S-parameter extraction
  in §4.3 implicitly assumes the port face is at the reference plane
  `z = z_port`. If the FEM mesh's port face is slightly offset (sub-mm)
  from the analytic plane the `S_{11}(f)` phase rotates linearly with
  frequency — easy to mistake for a real result. Mitigation: explicit
  reference-plane offset parameter on `WavePortFace`, defaulting to zero;
  document the convention in the API doc.
- **Complex-symmetric vs Hermitian pivoting in `faer` sparse LU.** The
  driven matrix `K(ω) + j k₀ B` is complex-symmetric (not Hermitian) for
  real `ε_r`. `faer`'s complex sparse LU handles this — Phase 4.fem.eig.1's
  `ComplexInverseIterEigen` already exercises the same code path on the
  lossy-cavity gate. Pre-flight at impl time, but no new risk.
- **PEC corner case at port-face / sidewall edges.** Edges that lie on the
  intersection of `Γ_port` and `Γ_PEC` must be tangential-E-zero (the
  PEC sidewall) and NOT carry a modal-source contribution (the port term).
  The face-classification iterator must detect shared edges and apply PEC
  precedence. Test fixture in plan step E2.
- **Sweep frequency density.** A 50-point uniform sweep across 8–12 GHz
  (80 MHz spacing) is sufficient for fem-eig-003 because the ABC reflection
  spectrum is smooth. Resonant geometries (fem-eig-004 slot antenna) will
  need adaptive sweeping; v0 ships uniform sweeps only.

## 11. Dependencies

- **`yee-fem` extension** — new `open_boundary` module; new face-block
  helpers in `element.rs`; assembly layer gains face-iteration over
  `face_kinds`.
- **`yee-mom`** — re-export `NumericalCrossSection` via `yee_mom::ports`
  (already public; no API change). The FEM side consumes it through the
  existing `e_tangential_at` accessor.
- **`yee-core`** — `MaterialDatabase` is already there from Phase 4.fem.eig.1;
  no API change. Possibly one new error variant for `BadFaceKind`.
- **`faer`** — already in workspace; complex sparse LU already exercised by
  Phase 4.fem.eig.1.
- **No new external crate.** The Phase 1.3.1.1 modal source and Phase
  4.fem.eig.1 complex sparse LU together cover the entire v2 surface.

No strict ordering constraint relative to other Phase 4 sub-projects.
fem-eig-003 stands alone behind its own walking-skeleton-extend gate.

## 12. Phase numbering ladder

- **Phase 4.fem.eig.0** — closed-cavity walking skeleton (shipped):
  rectangular metallic cavity, first-order Nedelec, lossless real `ε_r`,
  fem-eig-001 passes at 0.09 % rel.err.
- **Phase 4.fem.eig.1** — lossy dispersive `ε(ω)` Newton tracker (shipped):
  single-pole Drude / Lorentz / Debye, fem-eig-002 passes at 1.3e-3 Re(f) /
  3e-3 Im(f).
- **Phase 4.fem.eig.2** — **this spec**: 1st-order Engquist–Majda ABC +
  modal wave-port driven analysis; fem-eig-003 WR-90 stub gate passes at
  ±0.5 dB across 8–12 GHz vs Pozar §3.3.
- **Phase 4.fem.eig.2.0.1** — cubic / barycentric modal-profile
  interpolation on port-FEM mesh mismatch, if fem-eig-003 hits the v0
  interpolation floor.
- **Phase 4.fem.eig.2.1** — driven sweep over Phase 4.fem.eig.1 dispersive
  Newton tracker (combined open-boundary + lossy dispersive).
- **Phase 4.fem.eig.2.5** — 2nd-order Engquist–Majda ABC / Higdon ABC /
  UPML / CFS-PML, replacing the 1st-order termination when (and only when)
  a published case shows up that 1st-order Engquist–Majda cannot meet at
  the required tolerance.
- **Phase 4.fem.eig.3** — dielectric-resonator antenna (`fem-eig-005` or
  re-anchored `fem-eig-003`) with the puck modelled dispersively and the
  surrounding air halo ABC-terminated. End-to-end DRA validation against
  Petosa ch. 3.
- **Phase 4.fem.eig.4+** — periodic / Floquet BCs, GPU sparse solve,
  FEM-BEM hybrid for finite-aperture radiating problems. Open-ended.

## 13. Lane

Spec file:

```
docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md
```

Implementation lane (declared here for the follow-up plan, not edited by
this spec):

- `crates/yee-fem/src/open_boundary.rs` *(new)* — `OpenBoundarySolver`,
  `WavePortFace`, `FaceKind`, `SParameterRow`, `SParameters`.
- `crates/yee-fem/src/element.rs` — `assemble_abc_face_block`,
  `assemble_port_face_block`, `assemble_tet_element_complex` gains optional
  `abc_faces: &[FaceId]` arg.
- `crates/yee-fem/src/assembly.rs` — face-iteration over `face_kinds`,
  complex-symmetric driven-system assembly.
- `crates/yee-fem/src/solve.rs` — driven-solve helper (single complex LU
  back-substitution per swept frequency, reusing the Phase 4.fem.eig.1
  `faer::sparse::FaerLuSolver<Complex64>` surface).
- `crates/yee-fem/validation/README.md` — `fem-eig-003 (WR-90 stub + ABC)`
  row.
- `crates/yee-validation/{src,tests}/...` — fem-eig-003 driver.
- Out-of-lane (do not touch in the implementation PR): `yee-cli`,
  `yee-gui`, `yee-mom` (consume `NumericalCrossSection` via the existing
  public API only), `yee-mesh`, `yee-cuda`. The Python binding (`yee-py`)
  is plan step E6, optional.

## 14. References

- Engquist, B. and Majda, A., "Absorbing boundary conditions for the
  numerical simulation of waves", *Math. Comp.* 31 (1977), pp. 629–651 —
  the canonical 1st- and 2nd-order ABC derivation; v0 ships the 1st-order
  variant.
- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014. Ch. 10 (driven FEM analysis), §10.4 (ABC face contributions),
  §10.5 (wave-port modal decomposition), §10.7 (S-parameter extraction).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012. §3.3
  (waveguide TE/TM modes, propagation constants), §6.3 (rectangular
  cavity), §8.4 (cavity-coupled filters).
- Higdon, R. L., "Absorbing boundary conditions for difference
  approximations to the multi-dimensional wave equation", *Math. Comp.* 47
  (1986), pp. 437–459 — deferred to Phase 4.fem.eig.2.5.
- Berenger, J.-P., "A perfectly matched layer for the absorption of
  electromagnetic waves", *J. Comput. Phys.* 114 (1994), pp. 185–200 — the
  PML reference; deferred to Phase 4.fem.eig.2.5.
- Sacks, Z. S. et al., "A perfectly matched anisotropic absorber for use as
  an absorbing boundary condition", *IEEE Trans. Antennas Propag.* 43
  (1995), pp. 1460–1463 — UPML reference; deferred to Phase 4.fem.eig.2.5.
- `docs/superpowers/specs/2026-05-18-phase-4-fem-eigenmode-design.md` —
  Phase 4.fem.eig.0 spec; this spec strictly extends it.
- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-1-dispersive-design.md`
  — Phase 4.fem.eig.1 spec; the complex sparse LU surface used here is the
  same one shipped there.
- `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-cross-section-eigensolver-design.md`
  — Phase 1.3.1.1 spec; `NumericalCrossSection::e_tangential_at` is the
  modal-profile source for the wave-port face.
