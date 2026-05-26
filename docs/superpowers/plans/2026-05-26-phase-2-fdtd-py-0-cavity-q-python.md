# Implementation Plan — Phase 2.fdtd.py.0 Python Cavity-Q Driver

**Phase:** 2.fdtd.py.0  
**Date:** 2026-05-26  
**Spec:** [2026-05-26-phase-2-fdtd-py-0-cavity-q-python-design.md](../specs/2026-05-26-phase-2-fdtd-py-0-cavity-q-python-design.md)  
**ADR:** [ADR-0072](../../docs/src/decisions/0072-phase-2-fdtd-py-0-cavity-q-python.md)

---

## Worktree / Base

Branch: `feature/phase-2-fdtd-py-0-cavity-q` off current `main` HEAD  
Worktree: `worktrees/fdtd-py0/`  
Base SHA: (record the SHA at worktree creation time in the implementation commit)

---

## Lane

**ALLOWED paths (everything else is a finding, NOT a fix):**

```
crates/yee-py/src/fdtd.rs
crates/yee-py/src/lib.rs             (export PyCavityQResult + run_cavity_q)
crates/yee-py/tests/test_fdtd.py
crates/yee-validation/src/lib.rs
docs/src/tutorials/10-fdtd-lossy-cavity-from-python.md   (new file)
docs/src/SUMMARY.md                   (tutorial registration only)
```

---

## Pattern file

`crates/yee-py/src/fdtd.rs` — imitate `PyFdtdDriverConfig`, `PyFdtdDriver`,
`PyRadiationPattern` for the new `PyCavityQResult` and `run_cavity_q`.

`crates/yee-validation/src/lib.rs` lines 1270–1310 — `run_cpml_001()` and
`run_ntff_001()` show how a `CaseResult` is built.  Also look at
`run_fem_eig_001()` (lines 1315–1335) for a non-Skipped gate that calls
into a solver.

---

## Steps

### P1 — `PyCavityQResult` class in `crates/yee-py/src/fdtd.rs`

Add after `PyRadiationPattern`:

```rust
/// Result of a lossy-cavity Q-factor simulation.
///
/// Returned by [`run_cavity_q`].
#[pyclass(name = "CavityQResult", module = "yee._yee")]
pub struct PyCavityQResult {
    /// Q extracted from the log-linear ring-down fit.
    #[pyo3(get)]
    pub q_measured: f64,
    /// Analytic Q = 2π · f₁₀₁ · ε₀ / σ.
    #[pyo3(get)]
    pub q_analytic: f64,
    /// Analytic TE₁₀₁ resonant frequency (Hz).
    #[pyo3(get)]
    pub f101_hz: f64,
    /// |q_measured − q_analytic| / q_analytic.
    #[pyo3(get)]
    pub rel_err: f64,
    /// True iff rel_err < 0.05.
    #[pyo3(get)]
    pub passed: bool,
    /// Ring-down probe time series (n_ring samples).
    pub probe_vec: Vec<f64>,
}

#[pymethods]
impl PyCavityQResult {
    /// The ring-down probe time series as a numpy array.
    pub fn probe_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.probe_vec.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "CavityQResult(q_measured={:.4}, q_analytic={:.4}, rel_err={:.4e}, passed={})",
            self.q_measured, self.q_analytic, self.rel_err, self.passed
        )
    }
}
```

### P2 — Physics helpers and `run_cavity_q` free function

Add the physics helpers and the simulation runner immediately after P1 in
`fdtd.rs`.  The logic mirrors `crates/yee-fdtd/tests/cavity_q.rs` exactly —
no new physics, just a PyO3 wrapper:

```rust
use std::f64::consts::PI;
use yee_fdtd::boundary;
use yee_fdtd::{WalkingSkeletonSolver, YeeGrid};

const EPS0: f64 = 8.854_187_817e-12;
const C0: f64 = 299_792_458.0;

fn cavity_q_f101(nx: usize, nz: usize, dx: f64) -> f64 {
    let a = nx as f64 * dx;
    let d = nz as f64 * dx;
    0.5 * C0 * ((1.0 / (a * a)) + (1.0 / (d * d))).sqrt()
}

fn cavity_q_analytic(nx: usize, nz: usize, dx: f64, sigma: f64) -> f64 {
    2.0 * PI * cavity_q_f101(nx, nz, dx) * EPS0 / sigma
}

fn gaussian(t: f64, t0: f64, sigma_t: f64) -> f64 {
    let arg = (t - t0) / sigma_t;
    (-arg * arg).exp()
}

fn fit_log_decay(series: &[f64], dt: f64, t_start: f64) -> f64 {
    let n = series.len() as f64;
    let ts: Vec<f64> = (0..series.len()).map(|i| t_start + i as f64 * dt).collect();
    let ys: Vec<f64> = series.iter().map(|&v| v.abs().max(1e-30).ln()).collect();
    let t_mean = ts.iter().sum::<f64>() / n;
    let y_mean = ys.iter().sum::<f64>() / n;
    let num: f64 = ts.iter().zip(ys.iter()).map(|(&t, &y)| (t - t_mean) * (y - y_mean)).sum();
    let den: f64 = ts.iter().map(|&t| (t - t_mean).powi(2)).sum();
    -1.0 / (num / den)
}

/// Run a lossy rectangular PEC cavity simulation and return the Q-factor.
///
/// Builds a vacuum grid of `nx × ny × nz` cells at cell size `dx` metres,
/// fills it uniformly with electric conductivity `sigma` S/m, injects a
/// broadband Gaussian pulse for `n_src` steps, then records the ring-down
/// for `n_ring` steps and fits an exponential decay to extract Q.
///
/// Args:
///     nx: grid cells in x (default 20).
///     ny: grid cells in y (default 10).
///     nz: grid cells in z (default 20).
///     dx: cell size in metres (default 0.01 = 10 mm).
///     sigma: electric conductivity in S/m (default 2.96e-3 → Q ≈ 20).
///     n_src: number of source-injection steps (default 200).
///     n_ring: number of ring-down steps to record (default 6000).
///
/// Returns:
///     A [`CavityQResult`] containing q_measured, q_analytic, f101_hz,
///     rel_err, passed, and probe_array.
#[pyfunction]
#[pyo3(signature = (
    nx = 20,
    ny = 10,
    nz = 20,
    dx = 0.01,
    sigma = 2.96e-3,
    n_src = 200,
    n_ring = 6000,
))]
pub fn run_cavity_q(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    sigma: f64,
    n_src: usize,
    n_ring: usize,
) -> PyCavityQResult {
    let mut grid = YeeGrid::vacuum(nx, ny, nz, dx);
    grid.set_sigma_box(0, nx + 1, 0, ny + 1, 0, nz + 1, sigma);
    let dt = grid.dt;
    let mut solver = WalkingSkeletonSolver::new(grid);

    let src_i = nx / 4;
    let src_j = ny / 2;
    let src_k = nz / 4;
    let prb_i = nx * 3 / 4;
    let prb_j = ny / 2;
    let prb_k = nz * 3 / 4;

    let t0 = 12.0 * dt;
    let sigma_t = 4.0 * dt;

    for _ in 0..n_src {
        let t = solver.current_time();
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.grid_mut().ey[(src_i, src_j, src_k)] += gaussian(t, t0, sigma_t);
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
    }

    let mut probe = Vec::with_capacity(n_ring);
    for _ in 0..n_ring {
        solver.update_h_only();
        #[allow(deprecated)]
        boundary::apply_pec(solver.grid_mut());
        solver.update_e_only();
        solver.apply_cpml_e();
        solver.advance_clock();
        probe.push(solver.grid().ey[(prb_i, prb_j, prb_k)]);
    }

    let skip = n_ring / 3;
    let window = &probe[skip..];
    let t_start = (n_src + skip) as f64 * dt;
    let tau = fit_log_decay(window, dt, t_start);
    let f101 = cavity_q_f101(nx, nz, dx);
    let q_measured = PI * f101 * tau;
    let q_analytic = cavity_q_analytic(nx, nz, dx, sigma);
    let rel_err = (q_measured - q_analytic).abs() / q_analytic;

    PyCavityQResult {
        q_measured,
        q_analytic,
        f101_hz: f101,
        rel_err,
        passed: rel_err < 0.05,
        probe_vec: probe,
    }
}
```

### P3 — Register in `crates/yee-py/src/lib.rs`

In the `_yee` module builder, add:

```rust
m.add_class::<fdtd::PyCavityQResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_cavity_q, m)?)?;
```

Also add public re-exports in `yee/__init__.py` (if present) or trust the
`_yee` module auto-discovery.  Check how `FdtdRadiationPattern` is re-exported
(it appears at `yee.FdtdRadiationPattern` — mirror that pattern for
`CavityQResult` and `run_cavity_q`).

### P4 — Pytest cases in `crates/yee-py/tests/test_fdtd.py`

Append three new tests after the existing ones:

```python
def test_cavity_q_default_passes():
    """Default σ₀ = 2.96e-3 S/m → Q ≈ 20, rel_err < 5 %."""
    from yee import run_cavity_q, CavityQResult

    result = run_cavity_q()
    assert isinstance(result, CavityQResult)
    assert result.passed, (
        f"fdtd-202 gate failed: rel_err = {result.rel_err:.4f} "
        f"(q_measured={result.q_measured:.4f}, q_analytic={result.q_analytic:.4f})"
    )
    assert abs(result.q_analytic - 20.0) < 1.0, (
        f"analytic Q should be ~20, got {result.q_analytic}"
    )


def test_cavity_q_probe_array_shape():
    """probe_array should be a numpy array of length n_ring."""
    import numpy as np
    from yee import run_cavity_q

    result = run_cavity_q(n_ring=100)
    arr = result.probe_array()
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (100,)
    assert arr.dtype == np.float64


def test_cavity_q_repr_smoke():
    """__repr__ should include q_measured and passed."""
    from yee import run_cavity_q

    result = run_cavity_q()
    r = repr(result)
    assert "CavityQResult" in r
    assert "q_measured" in r
```

### P5 — `run_fdtd_202_lossy_cavity_q()` in `crates/yee-validation/src/lib.rs`

Add near the end of the existing FDTD case section (after `run_dispersive_001`).
Duplicate the needed physics helpers as private functions (`fdtd202_f101`,
`fdtd202_q_analytic`, `fdtd202_gaussian`, `fdtd202_fit_log_decay`,
`fdtd202_run`).  The constants mirror `cavity_q.rs`:

```rust
fn run_fdtd_202_lossy_cavity_q() -> CaseResult {
    // ... (implement using same logic as cavity_q.rs)
    // Returns CaseResult { id: "fdtd-202", status: Passed/Failed, ... }
}
```

Wire into `Report::run_all()`:

```rust
let cases = vec![
    run_mom_001(),
    run_mom_002(),
    run_mom_003(),
    run_cpml_001(),
    run_ntff_001(),
    run_dispersive_001(),
    run_fdtd_202_lossy_cavity_q(),   // ← new
    run_fem_eig_001(),
    // ...
];
```

Update the doc-comment on `run_all()` to mention fdtd-202 (~0.4 s).

### P6 — Tutorial `docs/src/tutorials/10-fdtd-lossy-cavity-from-python.md`

Write a self-contained tutorial covering:
1. Motivation: why Q-factor matters (resonator design, filter bandwidth).
2. Physics: TE₁₀₁ ring-down, Q = ε₀ω/σ, Taflove §3.7.
3. Code walkthrough: `run_cavity_q()` call, inspecting fields.
4. Plotting: matplotlib ring-down plot (code block).
5. Extending: higher-Q case (lower σ), or checking `run_cavity_q(sigma=sigma0/10)`.

Register in `docs/src/SUMMARY.md`:

```markdown
- [FDTD lossy cavity from Python](tutorials/10-fdtd-lossy-cavity-from-python.md)
```

---

## Verification commands

```bash
# V1 — build
cargo check --workspace

# V2 — clippy
cargo clippy --workspace --all-targets -- -D warnings

# V3 — format
cargo fmt --check --all

# V4 — validation unit test (fast, ~0.4 s)
cargo test -p yee-validation --lib -- run_fdtd_202 --nocapture

# V5 — maturin develop
maturin develop --manifest-path crates/yee-py/Cargo.toml

# V6 — pytest cavity_q cases
pytest crates/yee-py/tests/test_fdtd.py -k cavity_q -v

# V7 — run_all includes fdtd-202
cargo test -p yee-validation --lib -- run_all --nocapture 2>&1 | grep fdtd-202
```

---

## Escape hatch

Blocked > 15 min → surface and stop.  Do NOT spend more than 15 min on any
single step.  If `maturin develop` fails for a PyO3 reason (import path
mismatch, missing `__init__.py` re-export), note it and move on to the next
step — the Rust + validation changes are independently mergeable.
