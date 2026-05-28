# Implementation Plan — Phase 1.validation.2: FDTD Aggregator Gate Integration

**Date:** 2026-05-27  
**Spec:** `docs/superpowers/specs/2026-05-27-phase-1-validation-2-fdtd-aggregator-gates-design.md`  
**ADR:** ADR-0074  
**Lane:** `crates/yee-validation/src/lib.rs`, `docs/src/decisions/0074-*.md`, `docs/src/SUMMARY.md`

---

## Steps

### Step 1 — Read the pattern + existing gate code

Read:
- Lines 1314–1496 of `crates/yee-validation/src/lib.rs` (the `fdtd-202` pattern)
- `crates/yee-fdtd/tests/cpml_reflection.rs` (for cpml-001 physics)
- `crates/yee-fdtd/tests/ntff_dipole.rs` (for ntff-001 physics)
- `crates/yee-fdtd/tests/dispersive.rs` (function `drude_slab_reflects_per_fresnel`)

Also check which functions are public in `crates/yee-fdtd/src/{update,sources,boundary}.rs`.

### Step 2 — Implement cpml-001 physics block

Immediately before `fn run_cpml_001()` (currently around line 1274 of lib.rs) add a
physics block `cpml001_run() -> f64` that:

1. Creates two `YeeGrid::vacuum(50, 50, 50, 1e-3)` grids.
2. Derives `dt` from the grid, sets `t0 = 20.0 * dt`, `sigma = 6.0 * dt`.
3. **PEC run**: `WalkingSkeletonSolver::new(grid)` — 300 steps of `step_with_source` at
   (25, 25, 25); records `E_z` probe trace at (38, 25, 25).
4. **CPML run**: `WalkingSkeletonSolver::with_cpml(grid, CpmlParams::for_grid(&grid, 10))`
   — same 300 steps; same probe.
5. Computes the late-time (steps 100–300) peak amplitude for each trace.
6. Returns `20 * log10(cpml_peak / pec_peak)` (negative dB = attenuation).

Update `run_cpml_001()` to call `cpml001_run()` and return
`CaseStatus::Passed` if result ≤ −30 dB (attenuation), else `CaseStatus::Failed`.

### Step 3 — Implement ntff-001 physics block

Add `ntff001_run() -> f64` (returns broadside/endfire dB ratio):

1. Creates `YeeGrid::vacuum(50, 50, 50, 1e-3)`.
2. `CpmlParams::for_grid(&grid, 10)`.
3. `WalkingSkeletonSolver::with_cpml(grid, params)`.
4. `NtffParams { f_probe: 15e9, box_margin_cells: 15, theta_rad: PI/2, phi_rad: 0.0 }`.
5. `NtffState::new(solver.grid(), ntff_params)`.
6. Loop 2000 steps calling `solver.step_with_source_and_ntff(25, 25, 25, t0, sigma, &mut ntff)`
   with `t0 = 12.0 * dt`, `sigma = 4.0 * dt`.
7. `ntff.far_field_at(PI/2, 0.0)` and `ntff.far_field_at(0.0, 0.0)`.
8. Returns `20 * log10(|broadside| / |endfire|)`.

Update `run_ntff_001()` to call `ntff001_run()` and return
`CaseStatus::Passed` if result ≥ 20 dB, else `CaseStatus::Failed`.

Special case: if `|endfire|` is 0 the ratio is infinite (perfect null) — return some
large dB value (e.g. 100 dB) and pass.

### Step 4 — Implement dispersive-001 physics block

Add `dispersive001_run() -> (f64, f64)` (returns `(gamma_measured, gamma_analytic)`):

This duplicates `drude_slab_reflects_per_fresnel` from `dispersive.rs`.

Constants:
```
N = 80, DX = 1e-3, NPML = 10, N_STEPS = 800
SRC = (20, 40, 40), PROBE = (30, 40, 40)
F_PROBE = 10e9   // 10 GHz
OMEGA_P = 2π × 20e9   // 20 GHz plasma frequency
SLAB_LO = 50, SLAB_HI = 70   // i range
```

Drude ε(ω): `ε(ω) = 1 − ω_p² / ω²`

Helper `run_trace_dispersive(n, dx, npml, n_steps, src, probe, t0, sigma, materials: Option<&MaterialMap>) -> Vec<f64>`:
- Creates `YeeGrid::vacuum(n, n, n, dx)` 
- Creates `CpmlState::new(&grid, CpmlParams::for_grid(&grid, npml))`
- If `materials` is Some, creates `DispersiveState::new(&materials)` (no, actually `DispersiveState::new(materials)`)
- Loop n_steps:
  - `update::update_h(&mut grid)`
  - `cpml.update_h(&mut grid)`
  - `sources::gaussian_pulse_ez(&mut grid, src.0, src.1, src.2, n as f64 * dt, t0, sigma)`
  - if dispersive: dispersive state update; else `update::update_e(&mut grid)`
  - `cpml.update_e(&mut grid)`
  - Record `grid.ez[probe]`

Main logic:
1. Build `MaterialMap` + `Material::Drude` for the slab.
2. Run vacuum reference (materials=None) → `trace_vac`.
3. Run Drude slab (materials=Some) → `trace_slab`.
4. `reflected = trace_slab - trace_vac` (element-wise subtraction).
5. DFT of `trace_vac` and `reflected` at f_probe.
6. Distances: `r_direct = (30-20)*DX = 10 mm`, `r_reflected = (30+20)*DX + 2*(50-30)*DX` 
   = 10 mm + 40 mm = 50 mm (method of images: source at x=20, reflected by wall at x=50,
   path length = x_probe + 2*(x_wall - x_probe) = 30 + 2*20 = 70 mm total? 
   Actually let me re-derive: source at x_s=20, probe at x_p=30, wall at x_w=50.
   Image source at x_i=80 (2*x_w - x_s = 100-20=80). 
   r_direct = |x_p - x_s| = 10 mm.
   r_reflected = |x_p - x_i| = 50 mm.
   The 1/r correction factor is r_reflected/r_direct = 5.)
7. `gamma_measured = |dft_reflected| / |dft_vac| * (r_reflected / r_direct)`.
8. Analytic: `ε = 1 - ω_p²/ω²` at ω = 2π * f_probe.
   At 10 GHz with ω_p = 2π * 20 GHz: ε = 1 - 4 = -3.
   n = sqrt(ε) = sqrt(-3) = j*sqrt(3) (evanescent in medium).
   R = (1-n)/(1+n). |R|² at normal incidence.
   Actually for complex n: |R| = |(1-n)/(1+n)| where n = sqrt(-3) = i*sqrt(3).
   n = i*1.732...; 1-n = 1 - i*1.732; 1+n = 1 + i*1.732.
   |1-n|² = 1 + 3 = 4; |1+n|² = 1 + 3 = 4. So |R| = 1 (total reflection for purely
   imaginary n — i.e., evanescent medium).
   
   Hmm, that gives |R| = 1, which doesn't match the test's "within 20%" claim.
   
   Wait, let me re-read the dispersive.rs test more carefully...

Actually I need to re-read the dispersive.rs test to get the exact parameters. Let me not guess and instead tell the agent to read the source test carefully and replicate it exactly.

### Step 4 — revised approach

Instead of trying to derive the exact physics from memory, just tell the agent to:
1. Read `crates/yee-fdtd/tests/dispersive.rs` carefully
2. Extract the exact constants and helper functions
3. Duplicate them in yee-validation with the same logic

The test file says tolerance 20% and ω_p = 2π × 20 GHz. I'll let the agent figure out the exact Fresnel formula.

### Step 5 — Add unit tests for the three gates

Mirror the pattern in `tests::fdtd_202_lossy_cavity_q_passes` and `run_all_includes_fdtd_202`:
- `fn cpml_001_passes()` — calls `run_cpml_001()`, asserts status == Passed
- `fn run_all_includes_cpml_001()` — calls `run_cpml_001()`, asserts id == "cpml-001"
- Same for ntff-001 and dispersive-001

### Step 6 — Write ADR-0074

Create `docs/src/decisions/0074-phase-1-validation-2-fdtd-aggregator-gates.md`.

### Step 7 — Update SUMMARY.md

Add ADR-0074 entry to the Decisions section.

### Step 8 — Verify

```bash
cargo test -p yee-validation --release -- --nocapture 2>&1 | grep -E "cpml-001|ntff-001|dispersive-001|PASS|FAIL"
cargo clippy -p yee-validation -- -D warnings
cargo fmt --check -p yee-validation
```

---

## WORKTREE

Branch: `feature/phase-1-validation-2-fdtd-gates`  
Base SHA: `a06fd9f11f821221afaea6ed63eb75496843932d`

## ESCAPE HATCH

If any gate takes > 30 s in debug mode → mark it Skipped in the aggregator with a
wall-time-gated note (same pattern as fdtd-201) and surface the finding.

If dispersive-001 Drude physics is significantly different from what this plan describes,
surface the discrepancy and stop — do not guess at physics.
