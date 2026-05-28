# Phase 2.fdtd.py.3 Implementation Plan  
## FDTD TF/SF Fresnel-Transmission Python Driver (fdtd-204 gate)

**Date:** 2026-05-27  
**Spec:** `docs/superpowers/specs/2026-05-27-phase-2-fdtd-py-3-fresnel-tfsf-python-design.md`  
**ADR:** 0076  
**Branch:** `feature/phase-2-fdtd-py-3-fresnel-tfsf`  
**Lane:** `crates/yee-py/src/**, crates/yee-validation/src/**, docs/**`

---

## Step 0 — Pre-flight

```bash
# Verify base
cargo check --workspace --no-default-features 2>&1 | tail -5
git log --oneline -3
```

## Step 1 — Validation gate in yee-validation

**File:** `crates/yee-validation/src/lib.rs`

1a. Add constants before the gate function:
```rust
/// fdtd-204 gate constants — TF/SF Fresnel transmission test.
const FDTD204_N: usize = 80;
const FDTD204_DX: f64 = 1.0e-3;   // 1 mm
const FDTD204_NPML: usize = 10;
const FDTD204_FREQ: f64 = 10.0e9;  // 10 GHz
const FDTD204_EPS_R: f64 = 2.2;    // slab ε_r
const FDTD204_SLAB_I0: usize = 50; // slab front face (inclusive)
const FDTD204_SLAB_I1: usize = 55; // slab back face (exclusive, 5-cell slab)
const FDTD204_N_STEPS: usize = 600;
const FDTD204_SETTLE: usize = 200;  // discard first 200 steps
```

1b. Add `fdtd204_run()` inline helper:
```rust
/// Inner driver for fdtd-204: TF/SF plane-wave Fresnel transmission.
/// Returns (t_measured, t_analytic).
pub(crate) fn fdtd204_run() -> (f64, f64) {
    use ndarray::Array3;
    use yee_fdtd::{CpmlParams, PlaneWaveDirection, PlaneWaveSource, WalkingSkeletonSolver, YeeGrid};

    let n = FDTD204_N;
    let dx = FDTD204_DX;
    let npml = FDTD204_NPML;
    let freq = FDTD204_FREQ;
    let eps_r = FDTD204_EPS_R;
    let slab_i0 = FDTD204_SLAB_I0;
    let slab_i1 = FDTD204_SLAB_I1;
    let n_steps = FDTD204_N_STEPS;
    let settle = FDTD204_SETTLE;

    // Geometry: TF box and probes
    let tf_i0 = 12_usize;
    let tf_i1 = n - npml - 1; // = 69
    let tf_j0 = 1_usize;
    let tf_j1 = n - 2;        // = 78
    let tf_k0 = 1_usize;
    let tf_k1 = n - 2;        // = 78
    let probe_inc  = (25_usize, n / 2, n / 2); // incident probe (vacuum, inside TF)
    let probe_trans = (62_usize, n / 2, n / 2); // transmitted probe (vacuum after slab)

    // Helper: run one simulation (vacuum or slab), return E_z trace at probe.
    let run_sim = |with_slab: bool, probe: (usize, usize, usize)| -> Vec<f64> {
        let grid = if with_slab {
            let mut eps_cells = Array3::<f64>::from_elem((n + 1, n + 1, n + 1), 1.0);
            for i in slab_i0..slab_i1 {
                for j in 0..=n {
                    for k in 0..=n {
                        eps_cells[(i, j, k)] = eps_r;
                    }
                }
            }
            YeeGrid::vacuum(n, n, n, dx).with_eps_r_cells(eps_cells)
        } else {
            YeeGrid::vacuum(n, n, n, dx)
        };
        let dt = grid.dt;
        let params = CpmlParams::for_grid(&grid, npml);
        let mut solver = WalkingSkeletonSolver::with_cpml(grid, params);
        let mut pw = PlaneWaveSource::new(
            tf_i0, tf_i1, tf_j0, tf_j1, tf_k0, tf_k1,
            PlaneWaveDirection::PlusX,
            freq, 50, dx, dt, 4,
        );
        let mut trace = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            solver.step_with_plane_wave(&mut pw);
            trace.push(solver.grid().ez[probe]);
        }
        trace
    };

    // Vacuum run: probe at probe_inc.
    let trace_vacuum = run_sim(false, probe_inc);
    // Slab run: probe at probe_trans.
    let trace_slab   = run_sim(true, probe_trans);

    // Measure peak amplitude after settling.
    let a_inc = trace_vacuum[settle..].iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
    let a_trans = trace_slab[settle..].iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
    let t_measured = if a_inc > 0.0 { a_trans / a_inc } else { 0.0 };

    // Analytic transfer-matrix prediction.
    let t_analytic = fdtd204_t_analytic(eps_r, slab_i1 - slab_i0, dx, freq);

    (t_measured, t_analytic)
}

/// Amplitude transmission coefficient for a lossless dielectric slab
/// (thickness `d_cells * dx`, ε_r = `eps_r`) at normal incidence,
/// frequency `freq`. Born & Wolf §1.6.2 / Pozar §2.3 transfer-matrix.
pub(crate) fn fdtd204_t_analytic(eps_r: f64, d_cells: usize, dx: f64, freq: f64) -> f64 {
    use std::f64::consts::PI;
    let c0 = yee_core::units::C0;
    let n2 = eps_r.sqrt();
    let d = d_cells as f64 * dx;
    let delta = 2.0 * PI * freq * n2 * d / c0; // phase depth
    let r12 = (1.0 - n2) / (1.0 + n2);
    let r23 = (n2 - 1.0) / (n2 + 1.0);
    let t12 = 2.0 / (1.0 + n2);
    let t23 = 2.0 * n2 / (n2 + 1.0);
    let exp_j2d = num_complex::Complex64::new(delta.cos() * 2.0 * delta.cos() - 1.0,
                                              2.0 * delta.cos() * delta.sin());
    // Use exact formula: e^{j2δ} = cos(2δ) + j sin(2δ)
    let cos2d = (2.0 * delta).cos();
    let sin2d = (2.0 * delta).sin();
    let exp2d = num_complex::Complex64::new(cos2d, sin2d);
    let _ = exp_j2d; // superseded above
    let ejd = num_complex::Complex64::new(delta.cos(), delta.sin());
    let denom = 1.0 + r12 * r23 * exp2d;
    let t = t12 * t23 * ejd / denom;
    t.norm()
}
```

1c. Add `run_fdtd_204()` function returning a `CaseResult`.

1d. Register in `run_all()`:
```rust
cases.push({
    id: "fdtd-204".into(),
    description: "TF/SF Fresnel transmission ε_r=2.2 slab (80³×600 steps, gate ≤5%)".into(),
    status: CaseStatus::Skipped,
    notes: "Wall-time ~5–15 min release; run via cargo test -p yee-validation \
            -- --ignored --release (fdtd_204_live_gate)".into(),
    wall_time_seconds: 0.0,
    plot_paths: Vec::new(),
});
```

Also add `#[ignore]` unit tests mirroring the cpml-001/ntff-001/dispersive-001 pattern.

**Pattern file:** `crates/yee-validation/src/lib.rs` — the `fdtd202_run()` function (lines ~1689+).

---

## Step 2 — Python binding in yee-py

**File:** `crates/yee-py/src/fdtd.rs`

2a. Add `PyFresnelTfsfResult` struct:
```rust
#[pyclass(name = "FresnelTfsfResult", module = "yee._yee")]
#[derive(Clone)]
pub struct PyFresnelTfsfResult {
    #[pyo3(get)]
    pub t_measured: f64,
    #[pyo3(get)]
    pub t_analytic: f64,
    #[pyo3(get)]
    pub rel_err: f64,
    #[pyo3(get)]
    pub passed: bool,
}

#[pymethods]
impl PyFresnelTfsfResult {
    fn __repr__(&self) -> String {
        format!("FresnelTfsfResult(t_measured={:.4}, t_analytic={:.4}, \
                 rel_err={:.4}, passed={})",
                self.t_measured, self.t_analytic, self.rel_err, self.passed)
    }
}
```

2b. Add `run_fresnel_tfsf()` function:
```rust
/// Run the fdtd-204 TF/SF Fresnel-transmission gate.
///
/// Parameters accept keyword-only overrides for testing; defaults match
/// the published fdtd-204 scenario (80³ grid, ε_r=2.2, d=5mm, 10 GHz, 600 steps).
#[pyfunction]
#[pyo3(signature=(*, nx=80, ny=80, nz=80, dx=1e-3, eps_r=2.2,
                  slab_i0=50, slab_i1=55, freq_hz=10e9, n_steps=600))]
pub fn run_fresnel_tfsf(
    nx: usize, ny: usize, nz: usize, dx: f64, eps_r: f64,
    slab_i0: usize, slab_i1: usize, freq_hz: f64, n_steps: usize,
) -> PyFresnelTfsfResult {
    // Inline the physics (mirrors yee-validation's fdtd204_run).
    // ... (implementation)
    let t_analytic = yee_validation::fdtd204_t_analytic(eps_r, slab_i1 - slab_i0, dx, freq_hz);
    let t_measured = /* ... run simulation ... */;
    let rel_err = (t_measured - t_analytic).abs() / t_analytic;
    PyFresnelTfsfResult { t_measured, t_analytic, rel_err, passed: rel_err < 0.05 }
}
```

**Pattern file:** `run_cavity_q` (lines ~697+ in `crates/yee-py/src/fdtd.rs`).

---

## Step 3 — Module registration

**File:** `crates/yee-py/src/lib.rs`

Add to module init:
```rust
m.add_class::<fdtd::PyFresnelTfsfResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_fresnel_tfsf, m)?)?;
```

**File:** `crates/yee-py/python/yee/__init__.py`

Add to imports and `__all__`:
```python
from ._yee import (
    ...
    FresnelTfsfResult,
    run_fresnel_tfsf,
)
__all__ = [
    ...
    "FresnelTfsfResult",
    "run_fresnel_tfsf",
]
```

---

## Step 4 — Pytest cases

**File:** `crates/yee-py/tests/test_fdtd.py`

```python
# ---------------------------------------------------------------------------
# Phase 2.fdtd.py.3 — fdtd-204 TF/SF Fresnel-transmission gate
# ---------------------------------------------------------------------------

def test_run_fresnel_tfsf_smoke():
    """Smoke test: 5-step run verifies API plumbing, not physics."""
    from yee import FresnelTfsfResult, run_fresnel_tfsf
    r = run_fresnel_tfsf(n_steps=5)
    assert isinstance(r, FresnelTfsfResult)
    assert hasattr(r, "t_measured")
    assert hasattr(r, "t_analytic")
    assert hasattr(r, "rel_err")
    assert hasattr(r, "passed")
    assert r.t_analytic > 0.0
    assert r.t_analytic < 1.0


def test_run_fresnel_tfsf_repr_smoke():
    """__repr__ contains FresnelTfsfResult and t_measured."""
    from yee import run_fresnel_tfsf
    r = repr(run_fresnel_tfsf(n_steps=5))
    assert "FresnelTfsfResult" in r
    assert "t_measured" in r
```

Note: the full physics gate (`n_steps=600`) is **not** run in the fast pytest
suite (wall-time ~5–15 min release); the `test_run_fresnel_tfsf_smoke` test
verifies the API is wired correctly. The physics gate is covered by the
`#[ignore]`-gated unit test in yee-validation.

---

## Step 5 — ADR + Tutorial + SUMMARY

**ADR:** `docs/src/decisions/0076-phase-2-fdtd-py-3-fresnel-tfsf-python.md`

**Tutorial:** `docs/src/tutorials/13-fdtd-fresnel-tfsf-from-python.md`

**SUMMARY.md changes:**
- Fix missing ADR-0075 entry (docs gap from py.2 session)
- Add tutorial 13 entry
- Add ADR-0076 entry

---

## Verification Commands

```bash
# Lint
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check --all

# Fast tests (< 60 s)
cargo test -p yee-validation --lib -- fdtd_204 --release
cargo test -p yee-py --test test_fdtd -- test_run_fresnel_tfsf --release

# Gate (slow, ~5-15 min)
cargo test -p yee-validation -- fdtd_204_live_gate --release --ignored

# Python smoke (requires maturin develop)
python -c "from yee import run_fresnel_tfsf; r = run_fresnel_tfsf(n_steps=5); print(r)"

# SUMMARY integrity
grep 'ADR-0075' docs/src/SUMMARY.md
grep 'ADR-0076' docs/src/SUMMARY.md
grep 'tutorial.*13' docs/src/SUMMARY.md
```

**Expected exit codes:** all 0.

---

## Escape Hatch

Blocked > 15 min on any single step → surface and stop. Do NOT:
- Touch `crates/yee-fdtd/src/**` (out of lane)
- Weaken the 5% gate without documenting the reason as a finding
- Modify other Python functions (only add new `FresnelTfsfResult` + `run_fresnel_tfsf`)
