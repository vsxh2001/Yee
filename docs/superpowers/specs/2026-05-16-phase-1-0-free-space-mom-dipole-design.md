# Phase 1.0 — Free-Space MoM Dipole — Design

**Date:** 2026-05-16
**Status:** Approved; ready for writing-plans
**Repo base SHA at design time:** `ca1efee` (post Phase 0 walking-skeleton)
**Predecessor:** `docs/superpowers/specs/2026-05-16-phase-0-multi-agent-execution-design.md`

This spec describes the **first physics sub-project** under Phase 1 of the Yee electromagnetic-simulation workspace. It produces a small, self-contained MoM solver that validates the canonical mom-001 case (half-wave dipole impedance Z ≈ 73 + j42 Ω) and lays a reusable foundation for the rest of Phase 1.

---

## 1. Scope & Success Criteria

### Sub-project decomposition

Phase 1 in the workspace `ROADMAP.md` spans many independent subsystems (multilayer Green's functions, RWG basis, ports, de-embedding, surface roughness, GPU paths, Python bindings, GUI, Rerun streaming). This spec covers **one** sub-project: the free-space MoM kernel.

The other Phase 1 sub-projects each get their own spec → plan → implementation cycle.

### In scope

- RWG basis on a hand-meshed thin cylinder (mesh generator lives in test fixtures).
- Free-space dyadic Green's function.
- Duffy + Gauss quadrature for singular and near-singular triangle-pair integrals.
- Delta-gap port at the central edge of the cylinder.
- `faer` complex dense LU on CPU.
- `S₁₁(f)` and impedance `Z_in(f) = Z₀ (1 + S₁₁) / (1 − S₁₁)` over a 21-point sweep 130–170 MHz.
- Touchstone v1.1 `.s1p` export through existing `yee-io`.

### Out of scope

- Multilayer / dielectric Green's functions.
- RWG on non-cylinder geometry (microstrip, patch).
- GPU paths.
- Far-field / radiation pattern.
- Wave ports, TRL / SOLT de-embedding.
- Surface roughness models.
- Python bindings, GUI.
- Real Gmsh wiring of `yee-mesh::Session` — the cylinder mesh is a test fixture here.

### Done means

1. `cargo test -p yee-mom` includes a new `dipole_z_at_resonance` integration test. It constructs a dipole with `L = 1.0 m`, runs a sweep, and asserts at `f = 150 MHz`:

   ```
   |Z_in − (73 + j42)| / |73 + j42| ≤ 0.05
   ```

2. `cargo test -p yee-mom -- --include-ignored` adds a slower `dipole_full_sweep` test that walks all 21 frequencies and writes `tests/results/dipole.s1p` for human inspection. The file round-trips through `yee_io::touchstone::read` at `1 × 10⁻¹²` relative tolerance.
3. `PlanarMoM::run` no longer returns `Unimplemented` for a `TriMesh` representing a thin cylinder with a tagged central edge (port tag `1`).
4. All Phase 0 gates (1–9) stay green.
5. New validation case `mom-001` lands in `crates/yee-mom/validation/README.md` with reference value, tolerance, and the published source (Balanis, *Antenna Theory*, 4th ed., Ch. 8 §8.2).
6. `cargo doc --no-deps -p yee-mom` clean.

### Performance budget (informational, not a gate)

- Mesh: 24-around × 24-axial = 1152 triangles → ~1700 RWGs.
- System size: 1700 × 1700 complex-double = ~46 MB.
- Fill per frequency: < 60 s on a single CPU core with quadrature orders 5/7.
- LU per frequency: < 2 s.
- 21-frequency sweep budget: < 25 min single-threaded. Parallelize over frequencies with `rayon::par_iter` if budget is tight (mom-001 frequencies are independent).

---

## 2. Architecture Decisions Locked During Brainstorming

| # | Decision | Choice |
|---|----------|--------|
| D1 | First Phase-1 sub-project | Free-space MoM kernel (validates `mom-001`) |
| D2 | Formulation | RWG on triangle mesh of thin cylinder — production-shaped, reused by every downstream sub-project |
| D3 | Mesh source | Hand-coded cylinder mesher in `crates/yee-mom/tests/fixtures/cylinder.rs`; real Gmsh wiring deferred to a separate sub-project |
| D4 | Frequency scope | Narrow sweep 130–170 MHz, 21 points (exercises `SParameters` + Touchstone export end-to-end) |
| D5 | Extra capability | None — CPU only, no far-field, no GPU (walking-skeleton-of-physics) |
| D6 | Code organisation | Modular: `basis.rs` + `greens.rs` + `quadrature.rs` + `fill.rs` + `solve.rs` — every later sub-project replaces `greens.rs` only |

---

## 3. Module Architecture

```
crates/yee-mom/src/
├── lib.rs            # public surface unchanged: PlanarMoM, SParameters
├── basis.rs          # pub(crate) RwgEdge, RwgBasis
├── greens.rs         # pub(crate) FreeSpaceGreen
├── quadrature.rs     # pub(crate) GaussTriangle, DuffyTransform
├── fill.rs           # pub(crate) impedance_matrix
└── solve.rs          # pub(crate) delta_gap_rhs, s_parameters_sweep
crates/yee-mom/tests/
├── touchstone_roundtrip.rs   # existing
├── fixtures/
│   ├── mod.rs                # makes fixtures visible to integration tests
│   └── cylinder.rs           # thin-cylinder TriMesh generator + port-edge tagging
├── dipole.rs                 # mom-001 fast + ignored full-sweep tests
└── results/                  # gitignored; nightly outputs
```

`PlanarMoM`'s public surface does NOT grow in this sub-project. Module visibility is `pub(crate)` so we keep the API freeze for downstream consumers until at least one more sub-project ships.

---

## 4. Module APIs

```rust
// basis.rs
pub(crate) struct RwgEdge {
    pub v0: u32, pub v1: u32,                      // shared-edge vertex indices
    pub tri_plus: u32, pub tri_minus: u32,         // adjacent triangle indices
    pub free_plus: u32, pub free_minus: u32,       // opposite-vertex indices
    pub length: f64,
    pub port_tag: u32,                              // 0 = interior; non-zero = port id
}

pub(crate) struct RwgBasis {
    mesh: TriMesh,
    edges: Vec<RwgEdge>,
    centroids: Vec<Vector3<f64>>,
    normals: Vec<Vector3<f64>>,
    areas: Vec<f64>,
}

impl RwgBasis {
    pub fn from_mesh(mesh: TriMesh) -> Result<Self, yee_core::Error>;
    pub fn n_basis(&self) -> usize;
    pub fn n_tris(&self) -> usize;
    pub fn port_basis_indices(&self, port_tag: u32) -> impl Iterator<Item = usize> + '_;
    pub fn eval(&self, k: usize, tri: u32, bary: [f64; 3]) -> Vector3<f64>;
    pub fn div(&self, k: usize, tri: u32) -> f64;
}

// greens.rs
pub(crate) struct FreeSpaceGreen {
    pub k0: Complex64,
    pub eta0: f64,
}

impl FreeSpaceGreen {
    pub fn new(freq_hz: f64) -> Self;
    pub fn scalar(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
    pub fn scalar_smooth(&self, r1: Vector3<f64>, r2: Vector3<f64>) -> Complex64;
}

// quadrature.rs
pub(crate) struct GaussTriangle {
    pub points: Vec<[f64; 3]>,
    pub weights: Vec<f64>,
}

impl GaussTriangle {
    pub fn order_3() -> Self;
    pub fn order_5() -> Self;
    pub fn order_7() -> Self;
}

pub(crate) struct DuffyTransform { /* shared-vertex / shared-edge / same-triangle topology */ }

impl DuffyTransform {
    pub fn for_shared_vertex(/* topology */) -> Self;
    pub fn for_shared_edge(/* topology */) -> Self;
    pub fn for_same_triangle() -> Self;
    pub fn integrate<F>(&self, order: usize, f: F) -> Complex64
        where F: Fn(Vector3<f64>, Vector3<f64>) -> Complex64;
}

// fill.rs
pub(crate) fn impedance_matrix(
    basis: &RwgBasis,
    green: &FreeSpaceGreen,
) -> faer::Mat<Complex64>;

// solve.rs
pub(crate) fn delta_gap_rhs(basis: &RwgBasis, port_tag: u32) -> faer::Mat<Complex64>;
pub(crate) fn s_parameters_sweep(
    basis: &RwgBasis,
    port_tag: u32,
    freq_range: yee_core::FreqRange,
    z0_ref: f64,
) -> Result<yee_io::touchstone::File, yee_core::Error>;
```

---

## 5. Data Flow

```
TriMesh
  │
  ├──► RwgBasis::from_mesh
  │     · Enumerate interior + port edges
  │     · Cache centroids, normals, areas, edge length
  │     │
  │     ▼
  │   RwgBasis
  │
  ▼
For each f ∈ freq_range:
  ├── FreeSpaceGreen::new(f)
  ├── impedance_matrix(basis, green)  →  Z (n × n complex)
  ├── delta_gap_rhs(basis, port_tag) →  b (n × 1 complex)
  ├── faer LU on Z, solve Z * I = b →  I (n × 1 complex)
  ├── Compute port-edge current → V_port, I_port → Z_in = V_port / I_port
  └── S₁₁ = (Z_in − Z₀) / (Z_in + Z₀)
  ▼
yee_io::touchstone::File (n_ports = 1)
```

---

## 6. Test Strategy

### Per-module unit tests (`#[cfg(test)] mod tests`)

| Module | Test | What it asserts |
|--------|------|-----------------|
| `basis` | `edge_count_two_tri_mesh` | A 2-triangle mesh produces exactly 1 interior RWG edge. |
| `basis` | `divergence_sign_and_magnitude` | `div(k, tri_plus) = +1/A_plus`, `div(k, tri_minus) = −1/A_minus`. |
| `basis` | `port_lookup_returns_tagged_edges` | `port_basis_indices(1)` returns only edges with `port_tag = 1`. |
| `greens` | `scalar_at_quarter_wavelength` | `|G(R = λ/4)| matches analytical`, phase exact. |
| `greens` | `wave_number_relation` | `k0 = 2π / λ` with `λ = c/f`. |
| `quadrature` | `gauss_3_integrates_cubic_exact` | Integrate cubic polynomial over reference triangle exactly. |
| `quadrature` | `duffy_integrates_one_over_r` | Duffy self-triangle `∫1/R` matches Stratton-Chu analytical. |
| `fill` | `two_rwg_symmetric` | 2-RWG mesh → symmetric 2×2 Z; off-diagonal magnitudes match a closed-form near-field check. |
| `solve` | `s11_zero_when_z_equals_z0` | Synthetic Z_in = Z₀ yields S₁₁ = 0. |
| `solve` | `delta_gap_rhs_length_weighting` | RHS entry for port edge k equals `V × length_k`. |

### Integration tests (`tests/dipole.rs`)

| Test | Mark | What it does | Pass criterion |
|------|------|--------------|----------------|
| `dipole_z_at_resonance` | always-on | Build cylinder L = 1.0 m, r = 5 mm; run at single f = 150 MHz; compute Z_in | `|Z_in − (73 + j42)| ≤ 5%` rel |
| `dipole_full_sweep` | `#[ignore]` | 21 frequencies 130–170 MHz; export `.s1p`; round-trip through `yee_io::touchstone::read` | exit 0; `.s1p` produced; round-trip eq. at `1 × 10⁻¹²` rel |
| `condition_number_within_bound` | always-on | Build same cylinder; assert `cond(Z) ≤ 1 × 10⁶` at 150 MHz | guards mesh quality regression |

### Fixture

```rust
// crates/yee-mom/tests/fixtures/cylinder.rs
pub fn thin_cylinder(
    length_m: f64,
    radius_m: f64,
    n_axial: usize,
    n_around: usize,
) -> TriMesh;
// Triangulates the lateral surface of a cylinder (no end caps — mom-001
// reference is a thin-wire approximation; the ±5% tolerance accommodates
// the missing-end-cap error at L/r = 200). Edges on the two boundary
// circles are NOT interior and produce no RWG basis function.
// Tags the central axial ring of edges with port_tag = 1.
// Total: 2 * n_axial * n_around triangles.
```

---

## 7. Validation Gates

| Gate | Command | Pass |
|------|---------|------|
| Phase 0 gates 1–9 | unchanged | regression: all green |
| **`mom-001` fast** | `cargo test -p yee-mom dipole_z_at_resonance` | exit 0; assertion holds |
| `mom-001` sweep | `cargo test -p yee-mom -- --include-ignored dipole_full_sweep` | exit 0; `.s1p` produced |
| Per-module units | `cargo test -p yee-mom --lib` | all green |
| Condition number guard | `cargo test -p yee-mom condition_number_within_bound` | `cond(Z) ≤ 1 × 10⁶` |
| Documentation | `cargo doc --no-deps -p yee-mom` | warning-free |

The Phase 0 footnote G added to `ROADMAP.md` ("mom-001/002/003 reclassified as Phase 1") stays in place. After this sub-project lands, `mom-001` is a real Phase 1 gate; `mom-002` and `mom-003` remain footnoted until their own sub-projects.

---

## 8. Risk Register

| # | Risk | Likelihood | Impact | Mitigation |
|---|------|-----------|--------|------------|
| R1 | Singular integral handling buggy → wrong Z | Med | High | Duffy transform for same-tri, shared-edge, shared-vertex topologies handled separately; unit test against analytical `∫1/R` over reference triangle |
| R2 | Cylinder mesh aspect ratio extreme → ill-conditioned LU | Med | Med | 24 × 24 = 1152 tris, radius = L/200; record `cond(Z)` in `condition_number_within_bound` test; assert < `1 × 10⁶` |
| R3 | Delta-gap convention error (factor 2, sign) | Med | Med | Unit test on a synthetic 1-RWG "port" where `Z_in` is closed-form |
| R4 | Z₀ = 50 Ω conversion off by factor 2 | Low | Low | Unit test asserts `S(Z = Z₀) = 0` |
| R5 | `faer` complex LU performance for n ≈ 1700 | Low | Low | ~2 s/freq. If sweep too slow, parallelize over frequencies with `rayon` |
| R6 | mom-001 reference value (73 + j42 Ω) is geometry-specific | Confirmed | Med | Pin geometry exactly in fixture: L = 1.0 m, radius = L/200 = 5 mm, no end caps, delta-gap at central edge; cite Balanis Ch. 8 §8.2 in the validation README |
| R7 | Singularity subtraction misses near-singular outer-integration error | Med | Med | Quadrature order 7 outer integration when triangle centroids are within `3 × max(edge_len)`; cheap adaptivity |
| R8 | RWG sign convention (tri+ vs tri−) easy to flip during refactor | Med | High | Property test: assemble for a small mesh, assert Z is symmetric (reciprocal media). Skip the "div sums to zero on closed surface" check because the cylinder is open (no end caps). |

---

## 9. Next Step

After approval, invoke the `superpowers:writing-plans` skill to produce a detailed task-by-task implementation plan with TDD-shaped steps for each module, the cylinder fixture, the integration tests, and the validation README update.
