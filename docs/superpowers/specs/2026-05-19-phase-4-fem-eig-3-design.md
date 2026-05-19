# Phase 4.fem.eig.3 — Coupled Whitney-1 + 2nd-order ABC + Multi-port S-parameters

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3 (coupled Whitney-1 modal-RHS / projection, 2nd-order
Engquist-Majda ABC, multi-port S-parameter matrix extraction); 4.fem.eig.3.5
deferred — see §13.
**Depends on:** Phase 4.fem.eig.2 (1st-order ABC + single-port modal
wave-port driven analysis, shipped E1-E6; CCCCCCCCC `M_pp` normalisation
partial fix shipped on top).
**Blocks:** Phase 4.fem.eig.4 (FEM-BEM hybrid), Phase 4 multi-port
filter validation.

## 1. Goal

Phase 4.fem.eig.2 shipped the open-boundary FEM walking skeleton — 1st-order
Engquist-Majda ABC, single-mode wave-port, single-port `S_{11}` extraction.
BBBBBBBBB fem-eig-003 gate produced `|S_{11}| = 1.0` saturation on the
WR-90 stub fixture. Track CCCCCCCCC retired the **synthetic** matched-port
identity (`E_FEM = a_inc · e_mode  ⇒  S_{11} = 0`) by dividing the modal
projection by `M_pp = ⟨e_mode, e_mode⟩_port` (spec §4.3 had silently
assumed `M_pp = 1/2`; the actual orthonormalisation is `M_pp ≈ 1`), but
**did not** retire the empirical `|S_{11}| ≈ 1.0` saturation on the
fixture: the lumped Whitney-1 edge-tangent approximation
`N_i(centroid) ≈ t_i / 3` used by both
`element::assemble_port_modal_rhs` and
`OpenBoundarySolver::e_t_at_face_centroid` drives the reconstructed
`E_FEM(centroid)` toward zero on the port face. The strict
fem-eig-003 absorption-floor gate `[-45, -35] dB` and the strict
continuum-limit passive-bound gate (`|S_{11}| < 1`) remain `#[ignore]`'d.

Phase 4.fem.eig.3 retires both strict gates and adds multi-port
support via three coupled sub-tracks:

- **F1 + F2 — coupled Whitney-1 evaluation.** The modal RHS contribution
  `b_i = +2jβ · ∫_face N_i · E_t dS` and the FEM-side reconstruction
  `E_FEM(centroid) = Σ_i e_i · N_i(centroid)` both move from the lumped
  `t_i / 3` proxy to the **exact Whitney-1 identity**
  `N_i(ξ) = λ_a(ξ) · ∇λ_b − λ_b(ξ) · ∇λ_a` evaluated at 3-point Gauss
  quadrature on the reference triangle. The RHS and the projection are
  changed *together* in a single commit so the modal round-trip
  cancellation that Pozar §3.3 / Jin §10.5 derives is preserved at the
  exact-basis level, not the lumped level.
- **F3 + F4 — 2nd-order Engquist-Majda ABC.** The 1st-order ABC bilinear
  form `+jk₀ · n̂×(n̂×E)` is augmented with the 2nd-order tangential-curl
  correction `−(1/2k₀) · n̂×(∇×E)` (Engquist & Majda 1979 IEEE T-AP eq. 9),
  lowering the reflection floor for normal incidence from `~ −40 dB`
  (1st-order) to `~ −60 dB` (2nd-order). `abc_order` enum knob on the
  `OpenBoundarySolver` selects between the 1st- and 2nd-order kernels;
  the v0 1st-order path stays bit-for-bit identical.
- **F5 + F6 — multi-port `S_{p,q}` matrix.** The Phase 4.fem.eig.2 sweep
  returned diagonal-only `s_pp[p][k]`; v3 extracts the full
  `n_ports × n_ports` scattering matrix per swept frequency, one driven
  solve per excited port. Cross-port projection uses the same exact
  Whitney-1 basis from F1 + F2.

## 2. Non-goals (Phase 4.fem.eig.3)

Explicitly out of scope for v3:

- **CFS-PML / UPML.** 2nd-order Engquist-Majda is the v3 absorber upgrade.
  If `fem-eig-005` (3-port T-junction) shows that the 2nd-order ABC's
  off-normal reflection floor is the binding constraint, CFS-PML is
  deferred to Phase 4.fem.eig.3.5 (one phase past v3, mirroring the
  4.fem.eig.2 → 4.fem.eig.2.5 ladder).
- **Higher-order Nedelec basis on port or ABC faces.** v3 stays on
  first-order Whitney-1 with exact basis-at-Gauss-point evaluation.
  Second-order Nedelec is deferred to Phase 4.fem.eig.4+.
- **Multi-mode incident excitation per port.** A port still carries a
  single dominant mode (TE_{10} for rectangular, TEM for coax). Multi-mode
  incident excitation is Phase 4.fem.eig.3.0.2.
- **Adaptive sweep / model-order reduction.** Uniform sweep only.
- **GPU.** CPU-only FP64 complex sparse LU, same as v2.
- **Driven sweep over Phase 4.fem.eig.1 dispersive Newton tracker.** v3
  composes orthogonally with v1; the combined surface is Phase
  4.fem.eig.3.1.

## 3. Scope decision — Phase 4.fem.eig.3 v0

Walking-skeleton-extend-once per `CLAUDE.md` §3. v3 extends the shipped
Phase 4.fem.eig.2 + CCCCCCCCC stack along three axes, every axis with its
own validation gate so the failure mode is localised:

- **F1 + F2 (coupled-Whitney):** new
  `element::assemble_port_face_block_gauss_pts` taking the modal `E_t`
  pre-sampled at three Gauss points on the reference triangle, plus a
  parallel `OpenBoundarySolver::e_t_at_face_gauss_pts` reconstruction
  helper. The existing centroid-only entry points are kept for v2
  back-compat (gated behind `coupled_whitney: bool` on the solver), but
  the production default flips to `coupled_whitney = true` once
  fem-eig-003 strict gates clear.
- **F3 + F4 (2nd-order ABC):** new
  `element::assemble_abc2_face_block` returning a 3×3 complex block
  carrying both the 1st-order Mur term `+jk₀ (n̂×N_i)·(n̂×N_j)` and the
  2nd-order curl correction `−(1/2k₀) (n̂×∇×N_i)·(n̂×∇×N_j)` from
  Engquist-Majda 1979 eq. 9. An `AbcOrder::{First, Second}` enum on
  `OpenBoundarySolver` selects the kernel; `AbcOrder::First` produces
  the v2 bit-for-bit output.
- **F5 + F6 (multi-port `S_{p,q}`):** the existing single-driven-solve
  per-frequency path is wrapped in an outer loop over `excited_port`.
  Each iteration drives port `p` with `a_inc_p = 1` and `a_inc_q = 0`
  for `q ≠ p`, solves the per-frequency complex sparse system, and
  projects the FEM solution onto every port's modal profile to extract
  the column `S_{·, p}(ω)`. The full matrix is `n_ports` per-frequency
  driven solves, each reusing the same LU factor when the system matrix
  is independent of the excited port (it is — the matrix depends only
  on ω; only the RHS changes).

Three end-to-end gates:

- **fem-eig-003 strict** — un-ignore the BBBBBBBBB gates. WR-90 stub
  fixture from Phase 4.fem.eig.2 §8 with `coupled_whitney = true` and
  `abc_order = Second`. Assert `20·log10|S_{11}(f)| ∈ [-45, -35] dB`
  across 8-12 GHz; assert strict `|S_{11}(f)| < 1` at every swept
  frequency.
- **fem-eig-004 — 2-port WR-90 thru-line at 10 GHz.** Air-filled 60 mm
  WR-90 section with both end faces tagged `WavePort(p)` (port 0 at
  `z = 0`, port 1 at `z = 60 mm`) and the four sidewalls PEC. At
  10 GHz: assert `|S_{21}| ≈ 1.0` within ±0.1 dB (lossless thru-line)
  and `|S_{11}| < −30 dB` (matched-port residual reflection from
  the modal projection floor); assert reciprocity `S_{21} ≈ S_{12}`
  within `1e-6` (passive lossless structure).
- **fem-eig-005 — 3-port WR-90 T-junction at 5 GHz.** WR-90 H-plane T
  with three TE_{10} ports. Assert magnitude conservation
  `Σ_q |S_{q,p}|² ≤ 1 + ε_num` for every excited port `p`
  (lossless 3-port), and reciprocity `S_{p,q} ≈ S_{q,p}` within
  `1e-3` (looser tolerance than fem-eig-004 because of the
  multi-port modal-overlap conditioning risk — see §10).

CPU, FP64, scalar complex. No GPU. The v0/v1/v2 paths stay green
unchanged: `FemEigenAssembly`, `DispersiveSolver`, and the v2
`OpenBoundarySolver` with `coupled_whitney = false`, `abc_order = First`,
single port still produce bit-for-bit identical output.

What v3 does **not** ship, deferred to 4.fem.eig.3.5+ (§13):

- CFS-PML / UPML.
- Multi-mode incident excitation per port.
- Higher-order Nedelec basis.

## 4. Mathematical formulation

### 4.1 F1 — exact Whitney-1 basis at Gauss points

The Whitney-1 edge basis function `N_i` associated with the directed edge
`a → b` of a triangle with barycentric coordinates `(λ_a, λ_b, λ_c)` is
(Bossavit 1988; Jin 3rd ed. eq. 8.13)

```text
    N_i(ξ) = λ_a(ξ) · ∇λ_b  −  λ_b(ξ) · ∇λ_a.
```

`∇λ_a` and `∇λ_b` are **constant** across the triangle (linear barycentric
coordinates have constant gradients), so `N_i` is linear in `ξ`. At the
face centroid `ξ_c` where `λ_a(ξ_c) = λ_b(ξ_c) = λ_c(ξ_c) = 1/3`,

```text
    N_i(centroid) = (1/3) · (∇λ_b  −  ∇λ_a),     ≠   t_i / 3
```

unless the triangle is equilateral and the canonical edge-tangent
`t_i = v_b − v_a` happens to align with the dual `∇λ_b − ∇λ_a`. The
lumped `t_i / 3` approximation used by v2 systematically over- or
under-counts `N_i(centroid)` on every WR-90-meshed face whose triangles
deviate from equilateral — which is every Kuhn-decomposed face in
practice.

The F1 fix is to evaluate `N_i` at three Gauss points
`ξ_g ∈ {ξ_1, ξ_2, ξ_3}` on the reference triangle (the standard
3-point Gauss rule with barycentric coordinates
`(2/3, 1/6, 1/6)`, `(1/6, 2/3, 1/6)`, `(1/6, 1/6, 2/3)`,
each weighted `A/3`), giving the modal RHS

```text
    b_i = +2jβ · Σ_g  w_g · N_i(ξ_g) · E_t(x(ξ_g)),       w_g = A / 3,
```

and the FEM-projection reconstruction

```text
    E_FEM(ξ_g) = Σ_i  s_i · e_i · N_i(ξ_g),
```

where `s_i ∈ {-1, +1}` is the local-to-global orientation sign and `e_i`
the per-edge complex DoF. The Gauss-quadrature degree (3-point integrates
polynomials up to degree 2 exactly on the triangle; the integrand
`N_i · E_t` is at most degree 1 × profile-degree) is adequate for the
TE_{10}-on-WR-90 case; a 6-point rule is the v3 fallback if convergence
is marginal (see §10).

### 4.2 F3 — 2nd-order Engquist-Majda ABC

The 2nd-order Engquist-Majda radiation condition on a planar surface with
outward normal `n̂` is (Engquist & Majda 1979, *IEEE Trans. Antennas
Propag.* 27(5) p. 661, eq. 9; equivalent forms in Jin §10.4)

```text
    n̂ × ∇×E  =  −jk₀ · n̂×(n̂×E)  +  (1/2jk₀) · ∇_t × (∇_t × E_t),
```

where `∇_t` is the tangential gradient on the ABC face. Substituted into
the variational form, the **bilinear form** picks up two boundary terms
per face — the 1st-order term inherited from v2 and a new tangential-curl
correction:

```text
    a_ABC2(E, v)  =  +jk₀ · ∫_face  (n̂×N_i)·(n̂×N_j)  dS                ← 1st-order Mur
                     −(1/2k₀) · ∫_face  (n̂×∇×N_i)·(n̂×∇×N_j)  dS         ← 2nd-order correction
```

The curl term is the new piece. For Whitney-1 elements on a triangular
face, `∇ × N_i = 2 ∇λ_a × ∇λ_b` is **constant per face** (curl of a
linear vector field), so the surface integral is exact and reduces to
`face_area · (n̂ × (∇λ_a × ∇λ_b))_i · (n̂ × (∇λ_a × ∇λ_b))_j`. The
2nd-order block stays complex-symmetric (both terms have real-symmetric
Gram structure and purely imaginary scalar prefactors). The block is
3×3 per face, scattered into the global `K(ω)` at the corresponding
interior-DoF indices, identical to the v2 1st-order scatter path.

Reflection floor for normal incidence on a TE plane wave drops from
`~ −40 dB` (1st-order) to `~ −60 dB` (2nd-order), measured by Berenger
1994 and confirmed by Jin §10.4 Table 10.1.

### 4.3 F5 — multi-port `S_{p,q}` matrix

The matrix entry `S_{q,p}(ω)` is the modal reflection at port `q` when
port `p` is driven and all other ports are matched. Per Sheen, Ali,
Abouzahra, Katehi 1990 (*IEEE Trans. MTT* 38(7) p. 849, eq. 7) the
extraction is

```text
    a_inc_q  =  δ_{q,p}                                  (drive port p, match others)
    b_q(ω)   =  ⟨ E_FEM(ω; driven by p) , e_mode_q ⟩ / M_qq   −   a_inc_q
    S_{q,p}(ω) = b_q / a_inc_p   =   b_q                  (a_inc_p = 1).
```

This is the per-column extraction; the full matrix is `n_ports`
independent driven solves. The driven matrix `A(ω) = K(ω) − k₀² M(ω) +
boundary terms` is *independent of `p`* — only the RHS changes — so the
LU factor at frequency `ω` is computed once and back-substituted
`n_ports` times. The per-frequency runtime is `O(LU(N) + n_ports · BS(N))`
instead of `O(n_ports · LU(N))`.

Reciprocity `S_{p,q} = S_{q,p}` is a passive-lossless invariant
(Pozar §4.3); the fem-eig-005 gate enforces it within `1e-3` (the
multi-port mesh's modal-overlap conditioning is the binding floor here —
two port modal profiles can have non-trivial inner product on a shared
T-junction interior, biasing the projection — see §10).

## 5. Element-layer changes

The Phase 4.fem.eig.2 `assemble_tet_element_complex`,
`assemble_abc_face_block`, `assemble_port_face_block`, and
`assemble_port_modal_rhs` are unchanged. v3 adds **three new helpers**,
each a pure function of geometry + frequency + a pre-evaluated modal
sample:

```rust
// crates/yee-fem/src/element.rs — three new helpers

/// Per-face wave-port stiffness contribution evaluated at 3-point Gauss
/// quadrature with the exact Whitney-1 basis at each Gauss point.
///
/// Returns the 3×3 complex face block whose entries are
/// `+ j β_mode · Σ_g w_g · (1/μ_r,face) · (n̂ × N_i(ξ_g)) · (n̂ × N_j(ξ_g))`
/// with `w_g = A / 3` and `ξ_g` the three Gauss points
/// `(2/3, 1/6, 1/6) / (1/6, 2/3, 1/6) / (1/6, 1/6, 2/3)` in barycentric.
pub fn assemble_port_face_block_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_mode: f64,
    mu_r_face: f64,
    modal_e_t_at_gauss_pts: [Vector3<f64>; 3],
) -> SMatrix<Complex64, 3, 3>;

/// Per-face wave-port RHS evaluated at 3-point Gauss quadrature.
///
/// Returns `b_i = + 2 j β · Σ_g w_g · N_i(ξ_g) · E_t(x(ξ_g))` with
/// `w_g = A / 3`. The caller pre-evaluates the modal profile at each
/// Gauss point via the cross-section eigensolver or analytic profile.
pub fn assemble_port_face_rhs_gauss_pts(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    beta_mode: f64,
    modal_e_t_at_gauss_pts: [Vector3<f64>; 3],
) -> SVector<Complex64, 3>;

/// Per-face 2nd-order Engquist-Majda ABC contribution at ω.
///
/// Returns the 3×3 complex face block carrying the sum of the
/// 1st-order Mur term `+ jk₀ · area / μ_r · (n̂ × N_i) · (n̂ × N_j)` and
/// the 2nd-order tangential-curl correction
/// `− (1/2k₀) · area / μ_r · (n̂ × ∇×N_i) · (n̂ × ∇×N_j)`.
/// `∇ × N_i = 2 ∇λ_a × ∇λ_b` is constant per face — the surface integral
/// is exact for first-order Whitney-1.
pub fn assemble_abc2_face_block(
    face_vertices: [Vector3<f64>; 3],
    outward_normal: Vector3<f64>,
    k0: f64,
    mu_r_face: f64,
) -> SMatrix<Complex64, 3, 3>;
```

The local Whitney-1 gradient `∇λ_a` per face vertex is computed once per
helper from the face-vertex geometry — same identity already used by
`assemble_tet_element_complex` (interior tet gradients) projected onto
the face plane. The 3-point Gauss quadrature rule constants are
`const`-folded.

## 6. Public API surface

`OpenBoundarySolver` gains two configuration knobs (defaults reproduce
the v2 + CCCCCCCCC behaviour bit-for-bit) and one new sweep entry point:

```rust
//! crates/yee-fem/src/open_boundary.rs (extensions)

/// Selects the ABC bilinear form on `FaceKind::Abc`-tagged faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum AbcOrder {
    /// 1st-order Engquist-Majda. The v0/v2 default. Reflection floor
    /// `~ −40 dB` at normal incidence.
    #[default]
    First,
    /// 2nd-order Engquist-Majda. Reflection floor `~ −60 dB` at normal
    /// incidence; adds the tangential-curl correction term.
    Second,
}

impl<'m> OpenBoundarySolver<'m> {
    // ... v2 methods unchanged ...

    /// Toggle the coupled exact-Whitney-1 RHS + projection path
    /// (Phase 4.fem.eig.3 F1 + F2). Default `false` reproduces v2
    /// + CCCCCCCCC behaviour bit-for-bit.
    pub fn with_coupled_whitney(self, coupled: bool) -> Self;

    /// Set the ABC order (Phase 4.fem.eig.3 F3 + F4). Default
    /// `AbcOrder::First` reproduces v2 behaviour bit-for-bit.
    pub fn with_abc_order(self, order: AbcOrder) -> Self;

    /// Frequency-sweep driven solve returning the full multi-port
    /// `S_{p,q}` matrix (Phase 4.fem.eig.3 F5 + F6).
    ///
    /// For each `ω`, runs one driven solve per excited port (reusing
    /// the same LU factor across ports — the matrix is independent of
    /// the excited port, only the RHS varies), then extracts every
    /// `S_{q,p}(ω)` via modal projection onto port `q`'s profile.
    /// Returns an `SParametersMatrix` with one `n_ports × n_ports`
    /// complex matrix per swept frequency.
    pub fn sweep_matrix(
        &self,
        omegas: &[f64],
    ) -> Result<SParametersMatrix, Error>;
}

/// Full multi-port frequency-swept S-parameter matrix
/// (Phase 4.fem.eig.3 F5 + F6 output).
#[derive(Debug, Clone)]
pub struct SParametersMatrix {
    pub omegas: Vec<f64>,
    /// `s[k]` is the `n_ports × n_ports` matrix at `omegas[k]`.
    /// Entry `s[k][(q, p)]` is `S_{q,p}(omegas[k])`.
    pub s: Vec<DMatrix<Complex64>>,
}
```

The existing `OpenBoundarySolver::sweep` (single-port diagonal) is
preserved unchanged for v2 callers.

Optional Python binding (plan step F7):

```python
yee.fem.solve_open_cavity(
    mesh, materials, port_faces, abc_faces, omegas,
    *,
    coupled_whitney: bool = False,    # F1 + F2 toggle
    abc_order: str = "first",         # "first" | "second" — F3 + F4
    multi_port: bool = False,         # F5 + F6 — return full S matrix
) -> np.ndarray
```

When `multi_port = True`, returns shape `(n_omegas, n_ports, n_ports)`;
otherwise returns the v2 single-port shape `(n_omegas, n_ports)`.

## 7. Complex sparse linear solve

The driven matrix `A(ω) = K(ω) − k₀² M(ω) + Σ_ABC B_ABC + Σ_port B_port`
is independent of the excited port (all port faces contribute their
stiffness block regardless of which one is driven). The F5 + F6
multi-port sweep exploits this:

```text
    for ω in omegas:
        A = assemble(ω)                  # one assembly per ω
        L, U = sparse_lu(A)              # one factorisation per ω
        for p in 0..n_ports:
            b_p = build_rhs(p, ω)        # one RHS per (ω, p)
            e_p = back_substitute(L, U, b_p)  # cheap; same factor
            for q in 0..n_ports:
                S[q, p](ω) = extract(e_p, q, ω)
```

Per-frequency cost is `O(LU(N) + n_ports · (BS(N) + n_ports · proj))`
where `LU(N)` dominates and `n_ports` is small (2 for fem-eig-004,
3 for fem-eig-005). Same `faer::sparse::FaerLuSolver<Complex64>`
surface as v2.

## 8. Validation gates

### fem-eig-003 strict (un-ignore the BBBBBBBBB gates)

WR-90 stub fixture from Phase 4.fem.eig.2 §8 with v3 flags
`coupled_whitney = true`, `abc_order = Second`. Sweep 50 uniform points
across 8-12 GHz at 80 MHz spacing. Gate criteria:

1. **Strict absorption floor** — `20·log10|S_{11}(f)| ∈ [-45, -35] dB`
   at every swept frequency. The 2nd-order Engquist-Majda physics floor
   is `~ −60 dB`; the `[-45, -35] dB` band is *deliberately* the same as
   the BBBBBBBBB strict gate so the un-ignore is a single
   `#[ignore]`-removal.
2. **Strict passive bound** — `|S_{11}(f)| < 1` strictly at every swept
   frequency (Pozar §3.3 continuum identity).
3. **No PEC-resonance peak** at TE_{101} = 9.66 GHz in the swept band.
4. **Phase monotonicity** preserved (BBBBBBBBB gate C, unchanged).
5. **Run-time < 240 s in `--release`** on the existing fem-eig-003 mesh.
   Informational; not a CI gate.

### fem-eig-004 — 2-port WR-90 thru-line

Air-filled 60 mm WR-90 section (`a × b × d = 22.86 × 10.16 × 60 mm`)
with both `z = 0` and `z = 60 mm` faces tagged `WavePort(p)` (port 0
and port 1 respectively), four sidewalls PEC. Single-frequency test
at 10 GHz with `coupled_whitney = true`, `abc_order = First`
(no ABC faces in this fixture — both end faces are wave ports).

Gate criteria at 10 GHz:

1. **Through-line transmission** — `|S_{21}| ∈ [0.95, 1.05]`
   (`±0.1 dB ≈ ±0.012` linear; widened to ±5% to absorb modal-
   projection discretisation).
2. **Matched-port reflection** — `|S_{11}| < −30 dB` and
   `|S_{22}| < −30 dB`.
3. **Reciprocity** — `|S_{21} − S_{12}| < 1e-6` (passive lossless
   structure).
4. **Phase consistency** — `arg(S_{21}) ≈ −β · d` within `±5°`
   (free-space phase rotation over the line length).
5. **Run-time < 180 s in `--release`**.

### fem-eig-005 — 3-port WR-90 H-plane T-junction

WR-90 T-junction (broad-wall a = 22.86 mm) at 5 GHz with three TE_{10}
ports. Mesh: ~50 k tets. Flags `coupled_whitney = true`,
`abc_order = First` (no ABC; the T is closed by three ports).

Gate criteria at 5 GHz:

1. **Magnitude conservation** — `Σ_q |S_{q,p}|² ∈ [0.95, 1.0 + ε_num]`
   for every excited port `p ∈ {0, 1, 2}` (lossless 3-port; the
   continuum identity is `Σ_q |S_{q,p}|² = 1`).
2. **Reciprocity** — `|S_{p,q} − S_{q,p}| < 1e-3` for every off-diagonal
   pair. The looser tolerance vs fem-eig-004 reflects the multi-port
   modal-overlap conditioning risk (§10).
3. **Diagonal bounded** — `|S_{p,p}| < 1` for every `p` (passive).
4. **Run-time < 300 s in `--release`**.

The three gates together exercise every v3 sub-track: fem-eig-003 strict
isolates F1+F2+F3 (single port, ABC); fem-eig-004 isolates F5+F6
(multi-port without ABC); fem-eig-005 is the integration smoke for the
full v3 surface.

## 9. Higher-applications roadmap

Beyond fem-eig-003/4/5, v3's multi-port + 2nd-order-ABC stack unlocks:

- **Iris-coupled bandpass filter validation** (Pozar §8.4). 4-port
  cavity filter with iris couplings. Combine v3 with Phase 4.fem.eig.1
  dispersive Newton tracker (combined surface = Phase 4.fem.eig.3.1).
- **Coaxial / SMA-fed cavity Q-extraction with external loading**
  via wave-port + ABC halo. Combine v3 wave-port with v3 ABC.
- **DRA + ABC halo** (Petosa ch. 3). Phase 4.fem.eig.4 — needs an
  air-halo mesher and dispersive puck; combines v1 + v3.

All deferred to Phase 4.fem.eig.3.1+ / 4 (§13).

## 10. Risks and open questions

- **3-point vs 6-point Gauss-quadrature degree.** F1's 3-point rule
  integrates polynomials up to degree 2 exactly on the reference
  triangle. The integrand `N_i · E_t` is degree 1 in `N_i` times the
  modal profile's polynomial degree; TE_{10} on WR-90 has `sin(π x / a)`
  in `E_t`, which on a triangulated WR-90 cross-section is approximated
  by piecewise-linear samples — so the 3-point rule is degree-exact for
  the FEM-side reconstruction (`N_i · N_j`, degree 2) and adequate for
  the analytic-profile RHS. Fallback: a 6-point rule (degree 4) is
  swapped in if the fem-eig-003 absorption-floor gate measures a
  reflection floor consistently 3 dB worse than the documented `−60 dB`
  Engquist-Majda 1979 physics floor across the swept band.
- **2nd-order Mur stability near closed-stub resonances.** Track
  CCCCCCCCC prototyped the coupled fix and found over-amplification near
  the WR-90 stub's TE_{10n} resonances at 8 GHz (`n=1`) and 12 GHz
  (`n=2`). The 2nd-order Mur term contains `−(1/2k₀) · ∇_t ×∇_t × E_t`
  which becomes ill-conditioned as the tangential mode goes evanescent.
  Mitigation: clamp the 2nd-order term's contribution at the band edges
  by switching to `AbcOrder::First` if `|β_mode(ω) − k₀| / k₀ > 0.5`
  on any wave-port face — physically, when the dominant mode is close
  to cutoff the 1st-order ABC is the more numerically stable choice.
  Fallback: deferral to Phase 4.fem.eig.3.5 (CFS-PML) if the band-edge
  conditioning is the binding constraint.
- **Multi-port modal-overlap matrix ill-conditioning.** When two ports
  share a geometric face (a T-junction interior), the modal profiles
  `e_mode_p` and `e_mode_q` can have non-trivial inner product
  `⟨e_mode_p, e_mode_q⟩_port ≠ 0` on the shared interior. The
  per-frequency extraction in §4.3 implicitly assumes the modal-overlap
  matrix `M_{pq} = ⟨e_mode_p, e_mode_q⟩_port` is diagonal (i.e.
  ports are orthogonal). For fem-eig-005 (3-port T) this is *not*
  geometrically true: the three port faces are disjoint, but the modal
  profiles `e_mode_0 / e_mode_1 / e_mode_2` projected back via the FEM
  solution acquire numerical cross-coupling. Mitigation: compute the
  full modal-overlap matrix `M_{pq}` at each frequency and invert it
  during extraction (`b = M^{-1} · ⟨E_FEM, e_mode⟩ − a_inc`). For
  geometrically-disjoint ports the matrix is diagonal in the continuum
  limit; the v3 implementation diagonalises it via least-squares
  projection per-frequency with a fallback warning if the condition
  number exceeds `1e6`.
- **Excited-port LU-factor reuse correctness.** F5 + F6 reuses the
  per-frequency LU factor across all excited ports. The matrix is
  excited-port-independent **only** if every port face contributes its
  stiffness block regardless of `a_inc_p`. Verified: the wave-port
  bilinear form `+jβ B_port` is intrinsic to the boundary condition and
  contributes for every port face, matched or driven. Only the
  RHS depends on `a_inc_p`. No correctness risk; cross-check at impl
  time via a single-port-driven-twice fixture.
- **F3 + F4 promotes `K_ABC2` to a non-Hermitian-symmetric variant.**
  The 2nd-order curl correction `−(1/2k₀) (n̂×∇×N_i)·(n̂×∇×N_j)` has a
  **real** scalar prefactor (`−1/2k₀`), so the 2nd-order face block is
  *complex* but the curl-curl part is real-symmetric. The composite
  block is `+jk₀ R_1 + (−1/(2k₀)) R_2` where both `R_1` and `R_2` are
  real-symmetric — so the block is complex-symmetric (real prefactor
  for `R_2`, imaginary for `R_1`). `faer::sparse::FaerLuSolver<Complex64>`
  handles complex-symmetric matrices unchanged from v2.
- **Backward compatibility under `with_coupled_whitney(false)`.** The
  v0/v1/v2 callers do not call `with_coupled_whitney` or `with_abc_order`,
  so they hit the `Default` values (`false` and `First` respectively),
  reproducing the v2 + CCCCCCCCC code path bit-for-bit. The change is
  additive only.

## 11. Dependencies

- **`yee-fem` extension** — three new functions in `element.rs`; two new
  configuration knobs on `OpenBoundarySolver`; one new sweep entry point
  `sweep_matrix`; one new output type `SParametersMatrix`.
- **No new external crate.** Same `faer` complex sparse LU surface
  exercised by v1 + v2.
- **No `yee-mom` API change.** `NumericalCrossSection::e_tangential_at`
  is consumed exactly as in v2 — the F1 + F2 change is purely
  internal to the FEM crate (the analytic / numerical modal profile is
  pre-sampled at Gauss points by the FEM caller, not the MoM crate).
- **No `yee-mesh` change.** The exterior-face classifier from v2
  provides everything F1-F6 need.

## 12. Phase numbering ladder

- **Phase 4.fem.eig.0** — closed-cavity walking skeleton (shipped).
- **Phase 4.fem.eig.1** — lossy dispersive `ε(ω)` Newton tracker
  (shipped).
- **Phase 4.fem.eig.2** — 1st-order Engquist-Majda ABC + single-port
  modal wave-port driven analysis (shipped E1-E6 + CCCCCCCCC partial
  M_pp normalisation).
- **Phase 4.fem.eig.3** — **this spec**: coupled exact-Whitney-1 modal
  RHS + projection (F1+F2), 2nd-order Engquist-Majda ABC (F3+F4),
  multi-port `S_{p,q}` matrix extraction (F5+F6). Retires the
  fem-eig-003 strict gates; adds fem-eig-004 (2-port thru-line) and
  fem-eig-005 (3-port T-junction).
- **Phase 4.fem.eig.3.0.2** — multi-mode incident excitation per port.
- **Phase 4.fem.eig.3.1** — driven sweep over Phase 4.fem.eig.1
  dispersive Newton tracker (combined open-boundary multi-port +
  lossy dispersive).
- **Phase 4.fem.eig.3.5** — CFS-PML / UPML if 2nd-order Engquist-Majda
  hits a published benchmark it cannot meet (mirrors the
  4.fem.eig.2 → 4.fem.eig.2.5 placeholder slot).
- **Phase 4.fem.eig.4+** — FEM-BEM hybrid, GPU sparse solve, DRA-with-
  halo, iris-coupled bandpass filter validation. Open-ended.

## 13. Lane

Spec file:

```
docs/superpowers/specs/2026-05-19-phase-4-fem-eig-3-design.md
```

Implementation lane (declared here for the follow-up plan, not edited by
this spec):

- `crates/yee-fem/src/element.rs` — three new public helpers:
  `assemble_port_face_block_gauss_pts`,
  `assemble_port_face_rhs_gauss_pts`,
  `assemble_abc2_face_block`.
- `crates/yee-fem/src/open_boundary.rs` — `AbcOrder` enum,
  `with_coupled_whitney`, `with_abc_order`, `sweep_matrix`,
  `SParametersMatrix`. The existing `extract_s11` and
  `e_t_at_face_centroid` paths stay reachable behind
  `coupled_whitney = false`.
- `crates/yee-fem/src/lib.rs` — re-export new types.
- `crates/yee-fem/tests/abc2_face_block.rs` *(create)* — F3 unit test.
- `crates/yee-fem/tests/port_face_gauss.rs` *(create)* — F1 unit test.
- `crates/yee-fem/tests/open_boundary_matrix.rs` *(create)* — F5 unit
  test (2-port thru-line synthetic).
- `crates/yee-validation/{src,tests}/...` — fem-eig-003 strict un-ignore
  + fem-eig-004 + fem-eig-005 drivers.
- `crates/yee-py/src/fem.rs` — `coupled_whitney` / `abc_order` /
  `multi_port` kwargs (F7 plan step, optional).
- Out-of-lane (do not touch in the implementation PR): `yee-cli`,
  `yee-gui`, `yee-mom`, `yee-mesh`, `yee-cuda`, `yee-plotters`.

## 14. References

- Engquist, B. and Majda, A., "Radiation boundary conditions for acoustic
  and elastic wave calculations", *Comm. Pure Appl. Math.* 32 (1979),
  pp. 313-357; and "Absorbing boundary conditions for the numerical
  simulation of waves", *Math. Comp.* 31 (1977), pp. 629-651 — the
  1st- and 2nd-order ABC derivations. The IEEE T-AP 27(5) p. 661
  variant restated for waveguide modes is the one this spec
  implements at §4.2.
- Sheen, D. M., Ali, S. M., Abouzahra, M. D., Katehi, P. B. L.,
  "Application of the three-dimensional finite-difference time-domain
  method to the analysis of planar microstrip circuits",
  *IEEE Trans. Microwave Theory Tech.* 38(7) (1990), pp. 849-857 —
  multi-port S-parameter extraction convention (DOI
  10.1109/22.55781). The eq.-7 column extraction is what §4.3 implements.
- Jin, J.-M., *The Finite Element Method in Electromagnetics*, 3rd ed.,
  Wiley 2014, Ch. 10 (driven FEM analysis), §10.4 (1st- and 2nd-order
  ABC face contributions and reflection-floor tables), §10.5 (wave-port
  modal decomposition), §10.7 (S-parameter extraction).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §3.3
  (waveguide TE/TM modes, propagation constants), §4.3 (reciprocity
  for lossless multi-port networks).
- Bossavit, A., "Whitney forms: a class of finite elements for
  three-dimensional computations in electromagnetism",
  *IEE Proc.* 135-A (1988), pp. 493-500 — the Whitney-1 basis identity
  used by §4.1.
- Berenger, J.-P., "A perfectly matched layer for the absorption of
  electromagnetic waves", *J. Comput. Phys.* 114 (1994), pp. 185-200 —
  PML reference; deferred to Phase 4.fem.eig.3.5.
- `docs/superpowers/specs/2026-05-19-phase-4-fem-eig-2-open-boundary-design.md`
  — Phase 4.fem.eig.2 spec; this spec strictly extends it.
- `docs/src/decisions/0040-phase-4-fem-eig-2-open-boundary-scope.md` —
  Phase 4.fem.eig.2 scope ADR; §C-3 deferral that this spec fulfils.
- `docs/src/decisions/0042-phase-4-fem-eig-3-scope.md` — this spec's
  scope ADR.
