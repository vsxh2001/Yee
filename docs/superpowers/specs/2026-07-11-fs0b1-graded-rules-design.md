# FS.0b.1 — graded mesh rules + graded voxelization

**Phase:** FS.0b.1 (full-suite track). **Plan:**
`docs/superpowers/plans/2026-07-11-fs0b1-graded-rules.md`. **ADR:**
`docs/src/decisions/0210-fs0b1-graded-rules.md`.

## Why

FS.0b.0 (ADR-0208) shipped the graded CPU kernel (taper reflection measured
−52.7 dB), but nothing *generates* graded spacings from a `Layout` and the
voxelizer is uniform-only, so the payoff — refining only the
staircase-limited feature regions that hold the FS.0a residual
(max linear Δ|S| = 0.198, all in the stub's open-end skirt; the next uniform
pass costs ~2.4 h at 19 M cells, ADR-0204) — is unreachable. FS.0b.1 closes
that gap: a per-axis rule generator in `yee_engine::automesh` and a graded
rasterizer in `yee-voxel`, gated by a physics test that must reproduce the
uniform-converged stub-notch answer at a measured cell-count reduction.

## Design

### 1. Rule generator — `yee_engine::automesh::auto_spacings`

```rust
pub struct GradedMeshOptions {
    pub margin_m: f64,      // CPML margin each x/y side (metres)
    pub air_above_m: f64,   // air above the trace plane (metres)
    pub npml: usize,        // absorber depth in cells (uniform coarse inside)
    pub growth: f64,        // max cell-to-cell ratio (default 1.3)
    pub guard_m: f64,       // fine half-band around trace edges (default h)
}
impl GradedMeshOptions {
    pub fn for_board(layout: &Layout, f_max_hz: f64) -> Self; // FS.0a-shaped
}

pub struct AutoSpacings {
    pub dx: Vec<f64>, pub dy: Vec<f64>, pub dz: Vec<f64>, // primal widths
    pub x0_m: f64, pub y0_m: f64,   // domain origin (node 0)
    pub k_gnd: usize, pub k_top: usize, // ground / trace-plane layers
    pub coarse_m: f64, pub fine_m: f64, // the two rule outcomes, for reports
}
impl AutoSpacings {
    pub fn to_spacings(&self) -> crate::GradedSpacings; // JobSpec-ready
}

pub fn auto_spacings(
    layout: &Layout, f_max_hz: f64, opts: &GradedMeshOptions,
) -> Result<AutoSpacings, String>;
```

**Rules** (the FS.0a rulebook generalized per axis):

- **Coarse ceiling everywhere:** `coarse = auto_dx(layout, f_max_hz)` —
  λ/20-in-dielectric, h/3, feature/2, clamped to [1 µm, 1 mm] (unchanged).
- **Fine spacing:** `fine = min(min_feature/2, coarse/2)`. The
  `min_feature/2` term is the FS.0a feature rule where a small feature/gap
  binds below the ceiling; the `coarse/2` term guarantees at least one
  halving at every trace edge. Rationale is measured, not aesthetic:
  the brief-literal "min_feature/2 only" rule produces **no refinement at
  all** on the FS.0a stub board (min_feature/2 = 1.5 mm > coarse
  0.533 mm), i.e. exactly the uniform pass-0 mesh whose notch ADR-0204
  already measured at 5.100 GHz — 5.2 % from the converged 4.850 GHz,
  outside the ±2 % gate before a single new solve (recorded as iteration 0
  in ADR-0210). The FS.0a convergence trajectory showed one uniform
  halving (0.533 → 0.267 mm) converges the notch; `coarse/2` puts that
  halving only where the staircase error lives.
- **Fine bands (x and y):** an interval `edge ± guard_m` around every
  trace-AABB edge coordinate on that axis, plus the whole of every
  inter-trace axis-gap (the `min_feature_m` AABB idiom). Overlapping /
  near bands are merged; bands must clear the absorbers (error, not
  clamp, if geometry runs into the CPML margin).
- **Grading:** geometric ladder `fine·g^i` (g = `growth` ≤ 1.3, the
  compute-019-certified regime) between fine and coarse; consecutive-cell
  ratio ≤ g everywhere **including junctions** by construction
  (`fine·g^(m+1) ≥ coarse` ⇒ the ladder-top→coarse step < g).
- **z:** the substrate is `n_sub = ceil(h / (coarse/2))` cells of exactly
  `h / n_sub` (`k_gnd = 0`, `k_top = n_sub` — the ADR-0108 z-stack, no air
  gap at the ground); the air above grows geometrically from the substrate
  spacing to coarse and continues coarse until `air_above_m` is covered.
- **CPML uniformity:** the `npml` outermost x/y cells are exactly `coarse`
  (bit-equal — the same f64 constant), satisfying the FS.0b.0
  `validate_cpml_layers` scope rule; `JobSpec::dx_m` stays `coarse`, the
  nominal spacing the σ_max recipe assumes.

**Axis mesher** (private `mesh_axis`): march left→right — `npml` coarse
cells; per merged fine interval: up-ladder (if leaving a fine band), coarse
fill (`floor`), down-ladder, then fine cells `ceil((b − p)/fine)`; tail
up-ladder + coarse `ceil` through the far absorber. The fine band may start
up to one coarse cell early (the `floor` leftover) — coverage is only ever
extended, never clipped. Total length: starts at the domain min exactly,
overshoots the max by < one coarse cell (the same `ceil` behaviour as the
uniform voxelizer).

### 2. Graded voxelizer — `yee_voxel::voxelize_microstrip_graded`

```rust
pub struct GradedVoxelGrid {
    pub dx_m: Vec<f64>, pub dy_m: Vec<f64>, pub dz_m: Vec<f64>,
    pub x0_m: f64, pub y0_m: f64,
    pub k_gnd: usize, pub k_top: usize,
}
pub struct GradedMicrostripModel {
    pub dims: (usize, usize, usize),
    pub eps_r_cells: Array3<f64>,      // [nx+1, ny+1, nz+1]
    pub pec_mask_ex: Array3<bool>,     // [nx,   ny+1, nz+1]
    pub pec_mask_ey: Array3<bool>,     // [nx+1, ny,   nz+1]
    pub port_cells: Vec<(usize, usize, usize)>,
    pub k_gnd: usize, pub k_top: usize,
    pub x_nodes_m: Vec<f64>, pub y_nodes_m: Vec<f64>, pub z_nodes_m: Vec<f64>,
}
pub fn voxelize_microstrip_graded(
    layout: &Layout, grid: &GradedVoxelGrid,
) -> GradedMicrostripModel;
```

Rasterizes onto per-axis **coordinate arrays** (cumulative sums of the
spacings from `x0/y0/0`): cell-centre point-in-polygon tests against the
true centres `node[i] + d[i]/2`, PEC/ε assignment per true cell (identical
loop structure and z-stack to `voxelize_inner`), port cell lookup by
coordinate (`partition_point` over nodes — the graded generalization of
`floor((x − x0)/dx)`). Returns raw arrays plus the node coordinates (probe
and aperture placement needs them), **not** a `YeeGrid` (whose scalar `dx`
and dt are meaningless on a graded grid). The uniform entry points are
untouched — the voxel-001 z-stack pin and every downstream gate keep their
exact behaviour.

Precision note: for constant arrays the graded centres differ from the
uniform `x0 + (i + 0.5)·dx` by cumulative-sum rounding only (≲ 1 ulp of the
coordinate, ~1e-18 m at board scale), while geometry edges sit ≥ half a
cell from any centre in every generator this workspace emits — so the
classification, and hence the masks, are bit-identical (gate
`voxel-graded-001` verifies exactly this, it is not assumed).

### 3. The physics gate — `engine-graded-001`

`crates/yee-engine/tests/engine_graded_notch.rs`, `#[ignore]`, release.
The FS.0a S.6 stub-notch board (fixture copied from `board_automesh.rs`),
solved ONCE on the `auto_spacings` grid (DUT + through-line reference on
the **same** grid — the ADR-0204 same-physical-problem lesson), measured
with the launch-normalized double ratio (`sparams::forward_transfer`,
never the single ratio — ADR-0204). The JobSpec is built directly in the
test; the uniform `board.rs` fixture is not rewired (that integration is
FS.0b.2).

Asserts (targets from the ADR-0204 measured uniform-converged answer):

- notch frequency within ±2 % of 4.850 GHz (bins at 50 MHz: 4.80–4.90);
- notch depth ≤ −20 dB;
- `cells_graded / cells_uniform` < a pinned ceiling, where `cells_uniform`
  is the dx = 0.267 mm pass-2 grid built (not solved) via
  `two_port_board_job` with the convergence loop's exact rescaled options
  — **measured first, then pinned**; runtime reported.

Probe triples sit on uniform-coarse stretches of the graded x-axis
(scanned from the arrays; 12 coarse cells ≈ 6.4 mm spacing, the FS.0a
value) — `fit_standing_wave` assumes equally spaced probes, so a triple
must never straddle a taper.

### 4. Fast gates

- `voxel-graded-001` (`yee-voxel/tests/voxel_graded_001_uniform_bitexact.rs`,
  non-ignored): constant spacing arrays equal to dx ⇒ eps/PEC masks,
  port cells, and dims **bit-identical** to `voxelize_microstrip`'s
  (exact array comparison), plus graded-z and coordinate-lookup sanity.
- `automesh` unit tests (non-ignored): growth ratio ≤ 1.3 everywhere
  (programmatic scan of all three arrays, junctions included); fine bands
  cover every trace edge ± guard (every cell intersecting the band has
  width ≤ fine); inter-trace gaps are fine; a single-rect layout with no
  sub-coarse features degenerates to near-uniform (concrete assertion:
  every cell outside the four edge bands + tapers is bit-equal `coarse`,
  and coarse cells dominate the count); total-length bookkeeping
  (`x0 = domain min` exactly; `Σ spacings` covers the domain, overshoot
  < one coarse cell); absorber layers bit-equal coarse; fine bands inside
  the absorber margin are an error.

## Out of scope (FS.0b.2+)

Rewiring `board.rs` / `converge_two_port` to graded passes; GPU graded
kernels; graded NTFF; finite-board / lifted-ground graded variants
(`voxelize_finite_board` stays uniform); non-AABB polygon edge extraction;
per-band fine spacings (one global fine value per mesh).
