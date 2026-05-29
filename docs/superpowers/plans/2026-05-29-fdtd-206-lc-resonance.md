# fdtd-206 series-LC resonant frequency gate — Implementation Plan

**Phase:** 2.fdtd.6.1  
**Spec:** `docs/superpowers/specs/2026-05-29-fdtd-206-lc-resonance-design.md`  
**ADR:** ADR-0080  
**Date:** 2026-05-29

---

## Step 1 — Constants + result struct in `yee-validation/src/lib.rs`

Add after the `fdtd-205` block (around line 1900):

```rust
// -----------------------------------------------------------------------
// fdtd-206 constants
// -----------------------------------------------------------------------
const FDTD206_NX: usize = 5;
const FDTD206_NY: usize = 5;
const FDTD206_NZ: usize = 40;
const FDTD206_DX: f64 = 1.0e-3;
const FDTD206_L_H: f64 = 1.0e-9;
const FDTD206_F0_HZ: f64 = 1.0e9;
// C = 1 / (4π² f₀² L)  →  25.330 pF
const FDTD206_C_F: f64 = 1.0 / (4.0 * PI * PI * FDTD206_F0_HZ * FDTD206_F0_HZ * FDTD206_L_H);
const FDTD206_R_OHM: f64 = 1.0;
const FDTD206_N_KICK: usize = 30;
const FDTD206_N_RING: usize = 5_000;
const FDTD206_DFT_N_BINS: usize = 1_000;
const FDTD206_DFT_F_LO_HZ: f64 = 0.5e9;
const FDTD206_DFT_F_HI_HZ: f64 = 1.5e9;
const FDTD206_TOL_F0_REL: f64 = 0.02;

/// Result from [`fdtd206_run`] (Phase 2.fdtd.6.1 series-LC resonance gate).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct LcResonanceResult {
    pub f_measured_hz: f64,
    pub f_analytic_hz: f64,
    pub rel_err: f64,
    pub passed: bool,
}
```

---

## Step 2 — Driver `fdtd206_run()` in `yee-validation/src/lib.rs`

```rust
/// Runs the fdtd-206 series-LC resonant frequency gate.
///
/// Sets up a 5×5×40 PEC box at dx=1 mm, places a series-LC port at the
/// centre cell (L=1 nH, C≈25.33 pF → f₀=1 GHz analytic, R=1 Ω → Q≈6.28),
/// kicks it with a broadband Gaussian, and extracts the ring-down frequency
/// via a DFT scan.  The gate asserts |f_measured − f₀| / f₀ < 2 %.
///
/// # Reference
/// f₀ = 1/(2π√LC) — Pozar §2.4 / Hayt & Kemmerly §14.1 (circuit theory).
pub fn fdtd206_run() -> LcResonanceResult {
    use std::f64::consts::PI;
    use yee_fdtd::{boundary, lumped::{LumpedRlcPort, SourceWaveform},
                   FdtdSolver, WalkingSkeletonSolver, YeeGrid};

    let mut grid = YeeGrid::vacuum(FDTD206_NX, FDTD206_NY, FDTD206_NZ, FDTD206_DX);
    let dt = grid.dt;

    let port_cell = (FDTD206_NX / 2, FDTD206_NY / 2, FDTD206_NZ / 2);
    let mut port = LumpedRlcPort::series_rlc(
        port_cell,
        FDTD206_R_OHM,
        FDTD206_L_H,
        FDTD206_C_F,
        SourceWaveform::None,
    );

    let mut solver = WalkingSkeletonSolver::new(grid);

    // Gaussian kick parameters: peak at step 10, sigma = 4 steps
    let t0_kick = 10.0 * dt;
    let sigma_kick = 4.0 * dt;
    let v_kick = 1.0_f64; // V/m

    // Kick phase
    for n in 0..FDTD206_N_KICK {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        let t = (n as f64) * dt;
        let kick = v_kick * (-(t - t0_kick).powi(2) / (2.0 * sigma_kick.powi(2))).exp();
        solver.grid_mut().ez[port_cell] += kick;
        solver.update_e_only();
        port.correct_e(solver.grid_mut(), n, dt);
        solver.advance_clock();
    }

    // Ring-down phase — record I_L
    let mut il_probe = Vec::with_capacity(FDTD206_N_RING);
    for n in FDTD206_N_KICK..(FDTD206_N_KICK + FDTD206_N_RING) {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.update_e_only();
        port.correct_e(solver.grid_mut(), n, dt);
        il_probe.push(port.inductor_current());
        solver.advance_clock();
    }

    // DFT scan: find peak frequency between F_LO and F_HI
    let df = (FDTD206_DFT_F_HI_HZ - FDTD206_DFT_F_LO_HZ) / (FDTD206_DFT_N_BINS as f64 - 1.0);
    let mut peak_amp = 0.0_f64;
    let mut f_peak = FDTD206_DFT_F_LO_HZ;
    for k in 0..FDTD206_DFT_N_BINS {
        let f = FDTD206_DFT_F_LO_HZ + (k as f64) * df;
        let omega = 2.0 * PI * f;
        let (mut re, mut im) = (0.0_f64, 0.0_f64);
        for (n, &il) in il_probe.iter().enumerate() {
            let phase = omega * (n as f64) * dt;
            re += il * phase.cos();
            im += il * phase.sin();
        }
        let amp = (re * re + im * im).sqrt();
        if amp > peak_amp {
            peak_amp = amp;
            f_peak = f;
        }
    }

    let rel_err = (f_peak - FDTD206_F0_HZ).abs() / FDTD206_F0_HZ;
    LcResonanceResult {
        f_measured_hz: f_peak,
        f_analytic_hz: FDTD206_F0_HZ,
        rel_err,
        passed: rel_err < FDTD206_TOL_F0_REL,
    }
}
```

---

## Step 3 — `run_fdtd_206_lumped_lc_resonance()` wrapper + `run_all` registration

Add the case runner (pattern: `run_fdtd_202_lossy_cavity_q`):

```rust
fn run_fdtd_206_lumped_lc_resonance() -> CaseResult {
    let t0 = Instant::now();
    let result = fdtd206_run();
    let wall_time_seconds = t0.elapsed().as_secs_f64();
    let status = if result.passed { CaseStatus::Passed } else { CaseStatus::Failed };
    let notes = format!(
        "f_measured={:.4e} Hz, f_analytic={:.4e} Hz, rel_err={:.4e} (gate < 2 %)",
        result.f_measured_hz, result.f_analytic_hz, result.rel_err
    );
    CaseResult { id: "fdtd-206".into(),
                 description: "Lumped series-LC resonant frequency (FDTD phase 2.fdtd.6.1)".into(),
                 status, notes, wall_time_seconds, plot_paths: Vec::new() }
}
```

Add to `run_all()` after `run_fdtd_205()`:
```rust
run_fdtd_206_lumped_lc_resonance(),
```

---

## Step 4 — Integration test `crates/yee-fdtd/tests/lumped_lc_resonance.rs`

Self-contained, pattern from `ohmic_skin_depth.rs`:
- No `yee-validation` dependency; implement the physics inline.
- `#[test] fn lumped_lc_resonance_f0_within_two_percent()`.
- Gate: same `rel_err < 0.02`.
- Detailed `eprintln!` diagnostic.

---

## Step 5 — Python wrapper

In `crates/yee-py/src/fdtd.rs`, add:
```rust
#[pyclass(name = "LcResonanceResult")]
#[derive(Clone)]
pub struct PyLcResonanceResult {
    inner: yee_validation::LcResonanceResult,
}

#[pymethods]
impl PyLcResonanceResult {
    #[getter] fn f_measured_hz(&self) -> f64 { self.inner.f_measured_hz }
    #[getter] fn f_analytic_hz(&self) -> f64 { self.inner.f_analytic_hz }
    #[getter] fn rel_err(&self) -> f64 { self.inner.rel_err }
    #[getter] fn passed(&self) -> bool { self.inner.passed }
    fn __repr__(&self) -> String {
        format!("LcResonanceResult(f_measured_hz={:.4e}, f_analytic_hz={:.4e}, \
                 rel_err={:.4e}, passed={})",
                self.inner.f_measured_hz, self.inner.f_analytic_hz,
                self.inner.rel_err, self.inner.passed)
    }
}

#[pyfunction]
pub fn run_lc_resonance(py: Python<'_>) -> PyResult<PyLcResonanceResult> {
    py.allow_threads(|| Ok(PyLcResonanceResult { inner: yee_validation::fdtd206_run() }))
}
```

In `crates/yee-py/src/lib.rs`, register:
```rust
m.add_class::<fdtd::PyLcResonanceResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_lc_resonance, m)?)?;
```

---

## Step 6 — Tutorial 16

File: `docs/src/tutorials/16-fdtd-lumped-lc-resonance-from-python.md`

Content:
- Brief physics (series LC resonance)
- Code snippet calling `run_lc_resonance()`
- Expected output with f_measured ≈ 1 GHz, rel_err ≈ 0.1%
- Reference to ADR-0080

---

## Step 7 — ADR-0080 + SUMMARY.md + ROADMAP.md

- `docs/src/decisions/0080-fdtd-206-lc-resonance.md`
- Add to `docs/src/SUMMARY.md` decisions list
- Add tutorial 16 entry to SUMMARY.md tutorials list
- Update ROADMAP.md Status snapshot

---

## Verification commands

```bash
# Core gate
cargo test -p yee-fdtd --test lumped_lc_resonance -- --nocapture

# Validation suite
cargo test -p yee-validation -- fdtd_206 --nocapture

# Full workspace
cargo test --workspace --no-default-features

# Lint
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check --all
```

All must exit 0.

---

## Escape hatch

If the DFT peak is more than 5 % from 1 GHz after a first implementation
attempt: inspect the `il_probe` time series (should oscillate at ~1 GHz
with decaying amplitude). If the signal is flat (no oscillation), the
kick is not coupling into the LC — check that `correct_e` is called AFTER
`update_e_only`, and that `SourceWaveform::None` is used in the port.
If the simulation diverges, check α_LC < 2 (verify DX and L values).
Surface finding and stop if blocked > 15 min.
