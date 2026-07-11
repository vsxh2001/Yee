# FS.0b.0 — Nonuniform (graded) grid support in the yee-compute CPU kernel

**Phase:** FS.0b.0 (full-suite track). **Plan:**
`docs/superpowers/plans/2026-07-08-fs0b-graded-grid.md`. **ADR:**
`docs/src/decisions/0208-fs0b-graded-grid.md`.

## Why

The FS.0a automesh convergence loop (ADR-0204) measured its residual error
(max linear Δ|S| = 0.198) concentrated entirely in staircase-limited regions
(a stub's open-end skirt), and the next *uniform* refinement pass costs
~2.4 h (19 M cells). Graded meshing refines only where the geometry needs it.
FS.0b.0 is the enabling physics step: per-axis nonuniform primal spacings in
the CPU kernel, gated bit-exact-on-uniform. Mesh *rules* and the automesh
integration are later increments (FS.0b.1+), out of scope here.

## Design (Taflove ch. 11 dual-step nonuniform Yee)

### Protocol

`yee_engine::JobSpec` gains

```rust
#[serde(default)]
pub spacings: Option<GradedSpacings>,
```

where `GradedSpacings { dx: Vec<f64>, dy: Vec<f64>, dz: Vec<f64> }` holds the
**primal cell widths** per axis (lengths `nx`, `ny`, `nz`, metres). `None` =
uniform `dx_m` everywhere — today's path, byte-identical wire format
(`serde(default)` keeps pre-FS.0b JSON deserializing).

`yee_compute` grows a mirror `GradedSpacings` type (no serde — the engine
translates, the `MaterialsSpec`/`Materials` precedent) with `validate`,
`validate_cpml_layers`, and `courant_limit` methods, plus
`CpuFdtd::set_spacings` (the `set_dispersive` idiom).

### Kernel

- **H updates** divide curl-E differences by the **primal** cell size at the
  H sample: adjacent E nodes sit on integer grid nodes, so their separation
  is one primal cell (`d_i`).
- **E updates** divide curl-H differences by the **dual** spacing at the E
  node: adjacent H samples sit at cell centres, so their separation is
  `dual_i = (d_{i−1} + d_i) / 2`; at a domain edge the single adjacent primal
  is used (`dual_0 = d_0`, `dual_n = d_{n−1}`). Interior bulk/CPML loops never
  touch the edge duals, but they are defined for the lumped-port arithmetic.
- Per-cell ε/μ/σ maps, PEC masks, and the update *order* are unchanged.
- **Lumped ports** use local spacings: a `ResistivePort` on `E_z` cell
  `(i, j, k)` takes transverse area `dual_x[i]·dual_y[j]` and length
  `primal_z[k]`; an `AperturePort`'s modal `V = Σ E_z·Δz` integral uses
  `primal_z[k]` per cell (its physical `height`/`area` stay caller-supplied,
  and the engine computes them as spacing sums when graded).

### One kernel, bit-exact on uniform (the implementation choice)

The existing kernel divides curl differences by the scalar `spec.dx/dy/dz`
**literally** (no precomputed inverses). We therefore store per-axis
**spacing** arrays (`primal`, `dual`) — *not* inverse-spacing arrays — and
keep the literal division: `Δ / dz_primal[k]` instead of `Δ / s.dz`. For the
uniform fill every array element is bit-equal to the scalar (`(d + d)/2 == d`
exactly in IEEE-754: `d + d` is an exact ×2, `/2` an exact halving), and
IEEE-754 division is a pure function of its operand bit patterns, so ONE
kernel serves both paths and bit-exactness on uniform holds **by
construction**. Precomputed inverses (multiply instead of divide) were
rejected: `a·(1/b)` differs from `a/b` in the last ulp. Gates `compute-018`
plus the pre-existing `compute-001/003/007` reference-parity gates verify the
claim.

### dt

Courant from the **minimum** spacing per axis, same formula and 0.9 factor
as the uniform path (`FdtdSpec::vacuum`):

```text
dt = 0.9 / (c₀ · sqrt(1/min_dx² + 1/min_dy² + 1/min_dz²))
```

`GradedSpacings::courant_limit` uses the identical expression shape as
`FdtdSpec::courant_limit`, so constant arrays produce a bit-identical dt.
An explicit `dt_s` is validated against the graded limit.

### CPML on a graded axis (scope decision)

For the walking skeleton, spacing must be **uniform within the CPML layers**
of any absorbing face (the `npml` outermost primal cells of that face's
axis); construction rejects violations with a clear error. Rationale: the
Roden–Gedney profile grading and the `sigma_max` recipe assume one cell size
per layer, and FS.0b.1+ mesh rules will never grade inside absorbers anyway.
`JobSpec::dx_m` remains the nominal spacing used by the
`CpmlConfig::for_spec` σ_max recipe — callers should keep the absorbing-layer
spacing equal to `dx_m` (documented, not enforced).

### Rejections (walking-skeleton scope)

- **GPU**: `spacings: Some(_)` + `backend: "gpu"` fails with the
  `ComputeError::Unsupported("graded grid (FS.0b) is not on the GPU yet")`
  message; `"auto"` falls back to CPU. The `GpuFdtd` API has no spacings
  parameter at this base, so the rejection lives at the engine translation
  layer (the ntff-rejection idiom in `run_job`); it moves into `gpu.rs` when
  the GPU kernel grows spacings (FS.0b.2+).
- **NTFF**: `ntff` + `spacings` is rejected (the accumulation samples a
  uniform scratch `YeeGrid`).
- **Dispersive ADE**: `CpuFdtd::set_spacings` and `set_dispersive` are
  mutually exclusive (the ADE fused E-step divides by the scalar spacings);
  asserted at attach time. Not reachable over the job protocol (no
  dispersive spec yet).

## Validation gates

1. **compute-018** (`tests/graded_uniform_bitexact.rs`, fast, non-ignored):
   a driven CPML scenario (soft source + resistive port + aperture port +
   probes) and a PEC-box scenario, each run BOTH ways — scalar `dx` vs
   constant `GradedSpacings` arrays — must produce **bit-identical** probe
   series and final fields (exact `f64` equality, the compute-007 idiom).
2. **compute-019** (`tests/graded_interface_reflection.rs`, `#[ignore]`,
   release gate): a pulse propagating along a graded x-axis (0.5 mm →
   geometric grading, ratio ≈ 1.12/cell ≤ 1.3, over 6 cells → 0.25 mm and
   back), free space, CPML at both x ends (uniform 0.5 mm inside the
   absorbers). Reference run: uniform 0.5 mm grid, identical upstream
   geometry, **same dt** (the graded Courant dt). The runs are bit-identical
   upstream until the pulse reaches the grading, so the probe **difference**
   isolates the grading-caused reflection; its peak within the window before
   any far-wall return, normalized by the incident peak, must sit below a
   floor **measured first, then pinned** (expected −40 dB or better;
   measured value recorded in ADR-0208).
3. Unit tests: dual-spacing arithmetic (edge cells, known graded axis,
   constant arrays reduce to bit-equal uniform values) and the validation
   errors (wrong array lengths, non-positive spacings, graded-inside-CPML
   rejection, protocol-level GPU/NTFF rejections, min-spacing dt).

## Out of scope (FS.0b.1+)

Mesh-grading rules and automesh integration; GPU graded kernels; graded NTFF;
graded dispersive ADE; σ_max per-layer recipes for absorbers with differing
face spacings; `yee-fdtd` reference-solver grading.
