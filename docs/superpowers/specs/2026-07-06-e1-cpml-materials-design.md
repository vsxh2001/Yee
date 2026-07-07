# E.1 — CPML + per-cell materials on `yee-compute` (design)

**Date:** 2026-07-06
**Track:** ENGINE-STUDIO-ROADMAP phase E.1 (parent spec:
`2026-07-05-gpu-engine-web-studio-design.md`, ADR-0175)
**Plan:** `docs/superpowers/plans/2026-07-06-e1-cpml-materials.md`
**ADR:** ADR-0176 (written with the outcome)

## 1. Goal

Bring the E.0 walking skeleton up to "open-domain, real-materials" capability on both backends:

- Roden–Gedney 2000 CPML absorbing boundaries (the same formulation, profiles, and per-axis
  enable mask as `yee_fdtd::cpml`),
- per-cell ε_r / μ_r / σ maps (lossy CA/CB E-update, Taflove §3.7) and per-component interior
  PEC masks,
- the legacy outer-face PEC clamp (needed as the reflecting reference in the CPML gate),
- a Gaussian-in-time soft `E_z` source (CPU only; exactly `sources::gaussian_pulse_ez` — full
  source/port work remains E.2).

`yee-fdtd` stays the reference: the CPU backend must remain **bit-exact** against
`WalkingSkeletonSolver`'s step for every new arm; the GPU backend stays tolerance-gated against
the CPU backend.

## 2. Reference semantics being ported (verified against source)

- Full step (`WalkingSkeletonSolver::step_with_source`): `update_h` → CPML-H **or** legacy
  `apply_pec` → soft source → `update_e` → CPML-E **or** `apply_pec` → `apply_pec_mask` →
  clock++.
- CPML (`yee_fdtd::cpml`): ψ arrays sized per component (6 E-shaped + 6 H-shaped, order
  xy/xz/yx/yz/zx/zy); E profiles graded at `(d+1)/npml`, H profiles at `(d+0.5)/npml`;
  `b = exp(−(σ/κ+α)Δt/ε₀)`, `c = σ(b−1)/(σκ+κ²α)`; corrections applied as a second pass after
  the bulk update with coefficient `Δt/(ε₀ε_r)` (note: the CPML E pass ignores `sigma_cells` —
  mirrored as-is); `pml_depth` uses the component's own axis length (`ny+1` for E_x's y axis,
  `ny` for H_x's) and the per-axis enable mask.
- Per-cell arms (`yee_fdtd::update`): material arrays are `[nx+1, ny+1, nz+1]`, indexed by the
  component's own `(i,j,k)`; σ present ⇒ `E = CA·E + CB·curl` with
  `CA = (2ε₀ε_r − σΔt)/(2ε₀ε_r + σΔt)`, `CB = 2Δt/(2ε₀ε_r + σΔt)`.

## 3. Design

### 3.1 New public surface (`yee-compute`)

```rust
pub struct Materials {           // all optional; shapes validated at construction
    eps_r_cells, mu_r_cells, sigma_cells: Option<Vec<f64>>,   // [nx+1, ny+1, nz+1]
    pec_mask_ex, pec_mask_ey, pec_mask_ez: Option<Vec<bool>>, // per-component shapes
}
pub struct CpmlConfig { npml, m, sigma_max, kappa_max, alpha_max, axes }  // + for_spec()
pub enum Boundary { None, PecBox, Cpml(CpmlConfig) }
// Backends gain: with_config(spec, fields, Materials, Boundary),
// CPU additionally: step_with_gaussian_ez(source, t0, sigma) and probe access via fields().
```

`Boundary::None` preserves the E.0 raw-kernel semantics (compute-001 unchanged). `PecBox` is
the legacy reflecting clamp. Interior PEC masks apply whenever present (after the boundary
phase, as in the reference).

### 3.2 CPU backend

Flat-buffer, rayon-slab ports with per-cell arithmetic **identical** to the reference (same
match structure over `sigma_cells`/`eps_r_cells`, same op order). The CPML passes parallelize
over the outermost `i` index of the written component — ψ arrays share that component's shape,
so slabs of (field, ψ_a, ψ_b) zip together and stay disjoint; per-cell math is untouched, so
bit-exactness is preserved by the same argument as E.0.

### 3.3 GPU backend — arena-buffer refactor

Adding materials (3), masks (3), ψ (12), and profiles as separate bindings would blow WebGPU's
default limit of **8 storage buffers per stage**. The backend therefore moves to arena buffers:

| binding | content |
|---|---|
| 0 (uniform) | dims, npml, axis mask, flags, inv_d, field/ψ offsets |
| 1 | field arena: ex‥hz at fixed offsets |
| 2 | coefficient arena: `ca`, `cb`, `ce_cpml`, `ch` — four `[nx+1,ny+1,nz+1]` maps **always
      materialized host-side in f64** from scalars or per-cell arrays, then narrowed. This
      removes all material branching from WGSL (`ca = 1` reproduces the lossless add exactly). |
| 3 | ψ arena: 12 arrays (dummy 1-element when CPML off) |
| 4 | CPML profiles: b, c, κ, b_h, c_h, κ_h (6 × npml; dummy when off) |
| 5 | PEC mask arena: u32 per element, ex‥ez shapes (dummy when no masks) |

Kernels: the six update entry points fuse bulk + CPML correction (algebraically identical to
the reference's two passes: both read the same frozen opposite-family field and each ψ cell is
touched once per step); three mask-clamp entry points zero masked E cells after the E half.
`PecBox` on the GPU is implemented by zeroing the outer tangential E faces host-side at upload —
no kernel ever writes those faces, so the invariant holds for the whole run (documented on the
constructor).

### 3.4 Validation gates

- **compute-003 (CPU, bit-exact, heterogeneous)** — `tests/cpu_e1_reference_parity.rs`:
  24×20×22 grid, dielectric slab (per-cell ε_r), lossy block (σ), μ_r ≠ 1 region, an interior
  PEC sheet with a slot (masks), CPML npml=5 **and** a second PecBox scenario; driven Gaussian
  source, 30 steps, `WalkingSkeletonSolver` vs `CpuFdtd::with_config` — max |Δ| == 0.0 on all
  six components in both scenarios.
- **compute-004 (CPML reflection ≥ 30 dB)** — `tests/cpml_reflection.rs` in `yee-compute`:
  the `yee-fdtd` gate methodology reproduced on `CpuFdtd` (50³, npml=10, 300 steps, source at
  centre, probe (38,25,25), PEC-vs-CPML trace difference measurement) — reduction ≥ 30 dB.
- **compute-005 (GPU vs CPU, E.1 scenario)** — `tests/gpu_e1_parity.rs`: Gaussian-ball initial
  condition (no source needed on GPU in E.1), CPML + dielectric slab + lossy block + PEC-slot
  masks, 100 steps, family-relative tolerances as compute-002; plus an absorption sanity check
  (field energy decays vs the same run with PecBox). Self-skips without an adapter.

## 4. Out of scope (deferred)

Sources/ports as engine primitives, plane waves, lumped elements (E.2); dispersive ADE + NTFF
(E.5); GPU FP64 (E.3); performance tuning of the CPML region loops (E.4 — the reference walks
the full grid with an early-`continue`, and E.1 mirrors that shape).
