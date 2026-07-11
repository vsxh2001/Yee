# ADR-0208: FS.0b.0 — graded (nonuniform) grid in the yee-compute CPU kernel

**Date:** 2026-07-11. **Status:** accepted. **Phase:** FS.0b.0
(`FULL-SUITE-ROADMAP.md`). **Spec:**
`docs/superpowers/specs/2026-07-08-fs0b-graded-grid-design.md`. **Plan:**
`docs/superpowers/plans/2026-07-08-fs0b-graded-grid.md`.

## Context

The FS.0a automesh convergence loop (ADR-0204) measured its residual error
(max linear Δ|S| = 0.198) concentrated entirely in staircase-limited regions,
and the next uniform refinement pass costs ~2.4 h (19 M cells). Graded
meshing refines only where needed; FS.0b.0 is the enabling physics step —
per-axis nonuniform primal spacings in the CPU kernel, Taflove ch. 11
dual-step style. Mesh rules and automesh integration are FS.0b.1+.

## Decision

- **Protocol:** `JobSpec.spacings: Option<GradedSpacings>`
  (`#[serde(default)]` — pre-FS.0b JSON still deserializes), with
  `GradedSpacings { dx, dy, dz: Vec<f64> }` = primal cell widths (lengths
  nx/ny/nz). `yee_compute::GradedSpacings` mirrors it (the
  `Materials`/`MaterialsSpec` precedent) and `CpuFdtd::set_spacings` attaches
  it (the `set_dispersive` idiom).
- **Kernel:** H updates divide curl-E differences by the **primal** cell
  width at the H sample; E updates divide curl-H differences by the **dual**
  spacing `(d_{i−1}+d_i)/2` (single adjacent primal at domain edges). Same
  mapping inside the CPML ψ corrections. Lumped ports use local spacings
  (resistive: dual transverse area × primal dz; aperture: per-cell primal dz
  in the modal V integral; the engine computes aperture height/area as
  spacing sums when graded).
- **One kernel, bit-exact-on-uniform by construction.** We chose ONE kernel
  that always divides by per-axis **spacing arrays** — *not* precomputed
  inverse arrays, and *not* a duplicated graded code path. The pre-existing
  kernel divides by the scalar `spec.dx/dy/dz` literally, so keeping literal
  division and filling the arrays with the scalar makes every divisor
  bit-equal (`(d+d)/2 == d` exactly in IEEE-754), and division is a pure
  function of operand bit patterns. Precomputed inverses were rejected
  because `a·(1/b)` differs from `a/b` in the last ulp; a separate graded
  path was rejected as ~500 lines of duplicated delicate kernel code. The
  claim is verified, not just argued: gate `compute-018` (probes + all six
  field components, exact `f64` equality, CPML and PEC scenarios) and the
  pre-existing reference-parity gates `compute-001/003/007` all pass
  unchanged.
- **dt** from the minimum spacing per axis, same formula and 0.9 factor as
  `FdtdSpec::vacuum` (`GradedSpacings::courant_limit` uses the identical
  expression shape, so constant arrays give a bit-identical dt); explicit
  `dt_s` is validated against the graded limit.
- **Scope rule:** spacing must be uniform within the CPML layers of every
  absorbing face (validated, clear error) — the Roden–Gedney grading and
  σ_max recipe assume one cell size per layer, and mesh rules will never
  grade inside absorbers. `dx_m` stays the nominal spacing feeding
  `CpmlConfig::for_spec`; keep absorbing layers at `dx_m`.
- **Rejections:** GPU + spacings → `ComputeError::Unsupported("graded grid
  (FS.0b) is not on the GPU yet")` (auto falls back to CPU). The check lives
  at the engine translation layer because `GpuFdtd` has no spacings input at
  this base — there is nothing to reject inside `gpu.rs` yet; move it there
  when the GPU kernel grows spacings. NTFF + spacings rejected (uniform
  scratch-grid sampling). Dispersive ADE + graded mutually excluded at
  attach time (the fused ADE E-step divides by the scalar spacings).

## Measured results (gates)

- **compute-018** (`graded_uniform_bitexact.rs`, fast, non-ignored): scalar
  vs constant-array runs — max |Δ| = 0 exactly on both probe series and all
  six field components, under CPML and PEC, with soft source + resistive
  port + aperture port. PASS (bit-identical). The protocol-level twin in
  `yee-engine` (`constant_spacings_match_scalar_path_bit_exactly`) also
  matches bit-exactly, including the default dt.
- **compute-019** (`graded_interface_reflection.rs`, `#[ignore]`, release):
  0.5 mm → 0.25 mm → 0.5 mm geometric taper (ratio 2^(1/6) ≈ 1.122/cell,
  ≤ 1.3) along x, free space, CPML all faces (uniform 0.5 mm inside the x
  absorbers), uniform-reference difference method at an upstream probe with
  a shared graded-Courant dt (upstream evolution bit-identical until the
  pulse reaches the grading; early-window Δ measured < 1e-12 of the
  incident). **Measured 2026-07-11: incident peak 2.721e-3, grading
  reflection peak 6.322e-6 → −52.68 dB** — under the −40 dB expectation.
  **Pinned at −48 dB** (~4.7 dB margin). Runtime ~5 s release (two 560-step
  ~260k-cell runs).
- Unit tests: dual-spacing arithmetic (primal `[1, 2, 4]` → dual
  `[1, 1.5, 3, 4]`; constant arrays reduce bit-exactly), validation errors
  (lengths, non-positive/NaN widths, graded-inside-CPML min/max faces,
  graded Courant `dt_s` rejection, NTFF + graded, GPU + graded message,
  auto → CPU fallback).

No negative results; nothing was weakened. One honest caveat: the far-wall
CPML return does **not** cancel between the compute-019 runs (different
downstream grids), so the gate's domain has a deliberately long downstream
tail (60 coarse cells) to keep that return outside the 560-step measurement
window — shortening the tail or lengthening the run invalidates the
measurement.

## Consequences

- The CPU kernel is graded-ready; FS.0b.1+ can add mesh-grading rules and
  wire `automesh` to emit `spacings` without touching kernel physics.
- GPU (`fdtd.wgsl`) still assumes uniform spacing; FS.0b.2+ must thread the
  primal/dual arrays through the WGSL bindings and move the Unsupported
  rejection into `gpu.rs` construction.
- `yee-fdtd` (the scalar reference) remains uniform-only; graded runs are
  gated against a uniform reference run plus the measured reflection floor,
  not against a graded reference solver.
