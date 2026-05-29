# fdtd-205 Ohmic Skin-Depth Penetration Gate — Implementation Plan

**Date:** 2026-05-28
**Spec:** `docs/superpowers/specs/2026-05-28-fdtd-205-ohmic-skin-depth-design.md`
**ADR:** 0078
**Phase:** 2.fdtd.9

---

## Overview

Implement the `fdtd-205` validation gate: verify that the FDTD CA/CB
Ohmic-loss E-update reproduces the exponential skin-depth penetration
profile inside a conducting half-space. Single-step increment: one new test
file + registration in `yee-validation`.

---

## Pre-flight

```bash
cargo test -p yee-fdtd --test cavity_q  # baseline must pass
cargo clippy --workspace --all-targets -- -D warnings  # lint baseline
```

---

## Step 1 — Write `crates/yee-fdtd/tests/ohmic_skin_depth.rs`

**Pattern:** `crates/yee-fdtd/tests/cavity_q.rs`

### Constants

```rust
const NX: usize = 5;
const NY: usize = 5;
const NZ: usize = 130;
const DX: f64 = 1.0e-3;       // 1 mm cells
const FREQ: f64 = 1.0e9;       // 1 GHz
const SIGMA: f64 = 2.5331;     // S/m  →  δ = 10 mm = 10 cells
const MU0: f64 = 1.2566370614_e-6;  // H/m
const Z_SURFACE: usize = 50;  // First conductor cell (z index)
const SRC_X: usize = 2;
const SRC_Y: usize = 2;
const SRC_Z: usize = 25;       // z = 25 mm, inside vacuum region
const N_TRANSIENT: usize = 6_000;
const N_MEASURE: usize = 2_000;
```

### Analytic helper

```rust
fn analytic_skin_depth(sigma: f64, freq: f64) -> f64 {
    let omega = 2.0 * std::f64::consts::PI * freq;
    (2.0 / (omega * MU0 * sigma)).sqrt()
}
```

Test: `analytic_skin_depth(SIGMA, FREQ) ≈ 0.01 m ±0.1%`.

### Simulation helper

```rust
fn run_skin_depth_sim() -> (f64, f64, f64) {
    // Returns (amp_surface, amp_1delta, amp_2delta)
    use std::f64::consts::PI;

    // --- Build grid ---
    let mut grid = yee_fdtd::YeeGrid::vacuum(NX, NY, NZ, DX);
    // set_sigma_box uses exclusive upper bounds:
    // i0..i1  →  0..NX+1  covers all NX+1 x-planes
    // j0..j1  →  0..NY+1  covers all NY+1 y-planes
    // k0..k1  →  Z_SURFACE..NZ+1  covers z=50..130 (i.e. cells 50..129)
    grid.set_sigma_box(0, NX + 1, 0, NY + 1, Z_SURFACE, NZ + 1, SIGMA);

    let dt = grid.dt;
    let mut solver = yee_fdtd::WalkingSkeletonSolver::new(grid);

    let mut amp_surface: f64 = 0.0;
    let mut amp_1delta: f64 = 0.0;
    let mut amp_2delta: f64 = 0.0;

    for n in 0..N_TRANSIENT + N_MEASURE {
        let t = n as f64 * dt;

        solver.update_h_only();
        // Inject sinusoidal E_z source (soft source)
        solver.grid_mut().ez[(SRC_X, SRC_Y, SRC_Z)] +=
            (2.0 * PI * FREQ * t).sin();
        solver.update_e_only();
        solver.apply_cpml_e();   // applies PEC when no CPML configured
        solver.advance_clock();

        // Record peak amplitude during measurement window
        if n >= N_TRANSIENT {
            let ez_s = solver.grid().ez[(SRC_X, SRC_Y, Z_SURFACE)].abs();
            let ez_1 = solver.grid().ez[(SRC_X, SRC_Y, Z_SURFACE + 10)].abs();
            let ez_2 = solver.grid().ez[(SRC_X, SRC_Y, Z_SURFACE + 20)].abs();
            amp_surface = amp_surface.max(ez_s);
            amp_1delta  = amp_1delta.max(ez_1);
            amp_2delta  = amp_2delta.max(ez_2);
        }
    }

    (amp_surface, amp_1delta, amp_2delta)
}
```

### Production test (NOT `#[ignore]`)

```rust
#[test]
fn skin_depth_ratios_match_analytic() {
    let delta = analytic_skin_depth(SIGMA, FREQ);
    // Guard: δ ≈ 10 mm = 10 cells
    assert!(
        (delta / DX - 10.0).abs() < 0.01,
        "δ should be 10 cells, got {:.4} cells",
        delta / DX
    );

    let (amp_s, amp_1, amp_2) = run_skin_depth_sim();

    assert!(
        amp_s > 1e-10,
        "amp_surface should be non-trivial (field reached surface)"
    );

    let ratio_1 = amp_1 / amp_s;
    let ratio_2 = amp_2 / amp_s;
    let target_1 = (-1.0_f64).exp();   // e^{-1}
    let target_2 = (-2.0_f64).exp();   // e^{-2}

    let err_1 = (ratio_1 - target_1).abs() / target_1;
    let err_2 = (ratio_2 - target_2).abs() / target_2;

    println!(
        "fdtd-205: δ_analytic={:.1} mm, amp_surface={:.4e},
  ratio_1δ = {ratio_1:.4}  (target e⁻¹ = {target_1:.4}, rel_err = {:.1}%)
  ratio_2δ = {ratio_2:.4}  (target e⁻² = {target_2:.4}, rel_err = {:.1}%)",
        delta * 1e3,
        amp_s,
        err_1 * 100.0,
        err_2 * 100.0,
    );

    assert!(
        err_1 < 0.10,
        "Gate A FAILED: ratio_1δ = {ratio_1:.4}, target = {target_1:.4}, \
         rel_err = {:.1}% (threshold 10%)",
        err_1 * 100.0,
    );
    assert!(
        err_2 < 0.15,
        "Gate B FAILED: ratio_2δ = {ratio_2:.4}, target = {target_2:.4}, \
         rel_err = {:.1}% (threshold 15%)",
        err_2 * 100.0,
    );
}
```

**Verification:** `cargo test -p yee-fdtd --test ohmic_skin_depth` → exit 0.

---

## Step 2 — Register fdtd-205 in `yee-validation/src/lib.rs`

### 2a. Add `SkinDepthResult` struct

Place near the other Result structs (e.g. near `CavityQResult`):

```rust
/// Result struct for [`fdtd205_run`].
#[derive(Debug, Clone)]
pub struct SkinDepthResult {
    pub id: &'static str,
    pub delta_analytic_m: f64,
    pub amp_surface: f64,
    pub amp_1delta: f64,
    pub amp_2delta: f64,
    pub ratio_1delta: f64,
    pub ratio_2delta: f64,
    pub rel_err_1delta: f64,
    pub rel_err_2delta: f64,
    pub passed: bool,
}
```

### 2b. Add `pub fn fdtd205_run() -> SkinDepthResult`

```rust
pub fn fdtd205_run() -> SkinDepthResult {
    use std::f64::consts::E;
    // --- constants (mirror ohmic_skin_depth.rs) ---
    const NX: usize = 5;
    const NY: usize = 5;
    const NZ: usize = 130;
    const DX: f64 = 1.0e-3;
    const FREQ: f64 = 1.0e9;
    const SIGMA: f64 = 2.5331;
    const MU0: f64 = 1.2566370614e-6;
    const Z_SURFACE: usize = 50;
    const N_TRANSIENT: usize = 6_000;
    const N_MEASURE: usize = 2_000;

    let omega = 2.0 * std::f64::consts::PI * FREQ;
    let delta = (2.0 / (omega * MU0 * SIGMA)).sqrt();

    let mut grid = yee_fdtd::YeeGrid::vacuum(NX, NY, NZ, DX);
    grid.set_sigma_box(0, NX + 1, 0, NY + 1, Z_SURFACE, NZ + 1, SIGMA);
    let dt = grid.dt;
    let mut solver = yee_fdtd::WalkingSkeletonSolver::new(grid);

    let mut amp_surface: f64 = 0.0;
    let mut amp_1delta: f64 = 0.0;
    let mut amp_2delta: f64 = 0.0;

    for n in 0..N_TRANSIENT + N_MEASURE {
        let t = n as f64 * dt;
        solver.update_h_only();
        solver.grid_mut().ez[(2, 2, 25)] +=
            (2.0 * std::f64::consts::PI * FREQ * t).sin();
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();

        if n >= N_TRANSIENT {
            amp_surface = amp_surface.max(solver.grid().ez[(2, 2, Z_SURFACE)].abs());
            amp_1delta  = amp_1delta.max(solver.grid().ez[(2, 2, Z_SURFACE + 10)].abs());
            amp_2delta  = amp_2delta.max(solver.grid().ez[(2, 2, Z_SURFACE + 20)].abs());
        }
    }

    let target_1 = 1.0 / E;
    let target_2 = 1.0 / (E * E);
    let ratio_1 = if amp_surface > 0.0 { amp_1delta / amp_surface } else { 0.0 };
    let ratio_2 = if amp_surface > 0.0 { amp_2delta / amp_surface } else { 0.0 };
    let rel_err_1 = (ratio_1 - target_1).abs() / target_1;
    let rel_err_2 = (ratio_2 - target_2).abs() / target_2;
    let passed = rel_err_1 < 0.10 && rel_err_2 < 0.15 && amp_surface > 1e-10;

    SkinDepthResult {
        id: "fdtd-205",
        delta_analytic_m: delta,
        amp_surface,
        amp_1delta,
        amp_2delta,
        ratio_1delta: ratio_1,
        ratio_2delta: ratio_2,
        rel_err_1delta: rel_err_1,
        rel_err_2delta: rel_err_2,
        passed,
    }
}
```

### 2c. Add private `fn run_fdtd_205() -> CaseResult`

```rust
fn run_fdtd_205() -> CaseResult {
    let r = fdtd205_run();
    let notes = format!(
        "fdtd-205: δ_analytic={:.1}mm ratio_1δ={:.4}(err{:.1}%) \
         ratio_2δ={:.4}(err{:.1}%)",
        r.delta_analytic_m * 1e3,
        r.ratio_1delta,
        r.rel_err_1delta * 100.0,
        r.ratio_2delta,
        r.rel_err_2delta * 100.0,
    );
    CaseResult {
        id: "fdtd-205".to_string(),
        status: if r.passed { CaseStatus::Passed } else { CaseStatus::Failed },
        notes,
    }
}
```

### 2d. Register in `run_all()`

Add `run_fdtd_205()` call inside `run_all()` alongside the other FDTD
cases (e.g. after the `run_fdtd_202_lossy_cavity_q()` call):

```rust
cases.push(run_fdtd_205());
```

---

## Step 3 — Lint and format

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
```

---

## Step 4 — Verification commands

```bash
# Must exit 0:
cargo test -p yee-fdtd --test ohmic_skin_depth
cargo test -p yee-validation
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check --all
```

---

## Step 5 — Commit

Commit message:

```
yee-fdtd,yee-validation: fdtd-205 Ohmic skin-depth penetration gate

Adds validation gate fdtd-205: CA/CB Ohmic-loss E-update reproduces the
exponential skin-depth penetration profile |E(z)| = |E₀| exp(-z/δ) inside a
conducting half-space. δ_analytic = 10 mm (σ = 2.533 S/m, f = 1 GHz). Gates:
|ratio_1δ - e⁻¹| < 10 %, |ratio_2δ - e⁻²| < 15 %. Runs in ~0.3 s, NOT
#[ignore]'d, registered in run_all(). Complements fdtd-202 (Q-factor temporal
decay) with the spatial dimension of the same Ohmic E-update. Reference:
Griffiths §9.4.1 / Taflove §3.7.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
```

---

## DoD checklist

- [ ] `cargo test -p yee-fdtd --test ohmic_skin_depth` exits 0
- [ ] `cargo test -p yee-validation` exits 0
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` exits 0
- [ ] `cargo fmt --check --all` exits 0
- [ ] `"fdtd-205"` appears in `run_all()` as `CaseStatus::Passed`
- [ ] `fdtd205_run().passed == true`
