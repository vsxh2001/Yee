# Phase 2.fdtd.py.2 — Implementation Plan

**Spec:** `docs/superpowers/specs/2026-05-27-phase-2-fdtd-py-2-dipole-pattern-python-design.md`  
**ADR:** `docs/src/decisions/0075-phase-2-fdtd-py-2-dipole-pattern-python.md`  
**Phase:** 2.fdtd.py.2  
**Date:** 2026-05-27  

---

## Steps

### S1 — `PyDipolePatternResult` struct + `run_dipole_pattern()` in yee-py

File: `crates/yee-py/src/fdtd.rs`

Add after the existing `PyCavityResonanceResult` section:

```rust
/// Result of a short-dipole NTFF radiation-pattern simulation (fdtd-203 gate).
///
/// Returned by [`run_dipole_pattern`]. The five scalar fields are exposed
/// as Python read-only properties via `#[pyo3(get)]`. The full θ-cut is
/// available through `.theta_deg_array()` and `.e_theta_array()`.
#[pyclass(name = "DipolePatternResult", module = "yee._yee")]
pub struct PyDipolePatternResult {
    /// `|E_θ|` at θ = 0° (endfire null; should be < 0.15).
    #[pyo3(get)]
    pub e_theta_0: f64,
    /// `|E_θ|` at θ = 45° (should be ≈ 0.707 = sin 45°, tolerance ±0.15).
    #[pyo3(get)]
    pub e_theta_45: f64,
    /// `|E_θ|` at θ = 90° (broadside peak; normalized to 1.0 by FdtdDriver).
    #[pyo3(get)]
    pub e_theta_90: f64,
    /// `|E_θ|` at θ = 135° (should be ≈ 0.707 = sin 135°, tolerance ±0.15).
    #[pyo3(get)]
    pub e_theta_135: f64,
    /// `|E_θ|` at θ = 180° (endfire null; should be < 0.15).
    #[pyo3(get)]
    pub e_theta_180: f64,
    /// `true` iff all five gate criteria pass (see `run_dipole_pattern` docs).
    #[pyo3(get)]
    pub passed: bool,
    /// Internal: full θ array.
    pub theta_deg_vec: Vec<f64>,
    /// Internal: full `|E_θ|` array.
    pub e_theta_vec: Vec<f64>,
}

#[pymethods]
impl PyDipolePatternResult {
    /// Return the full θ array as a 1-D numpy `float64` array (37 points,
    /// 0° to 180° in 5° steps).
    pub fn theta_deg_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.theta_deg_vec.clone().into_pyarray(py)
    }

    /// Return the full `|E_θ|` array as a 1-D numpy `float64` array.
    ///
    /// Normalized so the maximum equals `1.0`.
    pub fn e_theta_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.e_theta_vec.clone().into_pyarray(py)
    }

    fn __repr__(&self) -> String {
        format!(
            "DipolePatternResult(e_theta_0={:.4}, e_theta_90={:.4}, \
             e_theta_180={:.4}, passed={})",
            self.e_theta_0, self.e_theta_90, self.e_theta_180, self.passed,
        )
    }
}
```

Then add the `run_dipole_pattern()` function:

```rust
/// Run a short-dipole FDTD simulation and return the far-field radiation pattern.
///
/// Builds a vacuum grid of `nx × ny × nz` cells at cell size `dx` metres,
/// drives a z-polarised sinusoidal current source over a 5-cell dipole at the
/// grid centre, runs for `n_steps` time steps with CPML on all six faces, and
/// returns the NTFF θ-cut of `|E_θ|` at `φ = 0` normalized to its maximum.
///
/// The simulation parameters match the `dipole_pattern.rs` integration test
/// (Phase 2.fdtd.4 driver).
///
/// # Reference
///
/// C. A. Balanis, *Antenna Theory*, 4th ed., §4.2: short-dipole far-field
/// `E_θ ∝ sin θ` — maximum at θ = 90° (broadside), nulls at θ = 0°/180°
/// (endfire).
///
/// # Gate criteria (fdtd-203)
///
/// | θ    | Expected | Gate               |
/// |------|----------|--------------------|
/// |  0°  | 0        | e < 0.15           |
/// | 45°  | 0.707    | |e − 0.707| < 0.15 |
/// | 90°  | 1.0      | |e − 1.0|  < 0.05  |
/// | 135° | 0.707    | |e − 0.707| < 0.15 |
/// | 180° | 0        | e < 0.15           |
///
/// # Arguments
///
/// * `nx`, `ny`, `nz` — grid cell counts (default 60 each).
/// * `dx` — cell size in metres (default 0.005 = 5 mm).
/// * `n_steps` — FDTD time steps (default 800).
/// * `source_freq_hz` — dipole drive frequency in Hz (default 1e9 = 1 GHz).
///
/// # Returns
///
/// A [`DipolePatternResult`] with scalar fields for the five sampled θ
/// values, a `passed` flag, and numpy-array accessors for the full cut.
#[pyfunction]
#[pyo3(signature = (
    nx = 60usize,
    ny = 60usize,
    nz = 60usize,
    dx = 5.0e-3_f64,
    n_steps = 800usize,
    source_freq_hz = 1.0e9_f64,
))]
pub fn run_dipole_pattern(
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    n_steps: usize,
    source_freq_hz: f64,
    py: Python<'_>,
) -> PyDipolePatternResult {
    use yee_fdtd::{FdtdDriver as RustDriver, FdtdDriverConfig as RustCfg, YeeGrid};

    let (theta_deg, e_theta_phi0) = py.allow_threads(|| {
        let grid = YeeGrid::vacuum(nx, ny, nz, dx);
        let ci = nx / 2;
        let cj = ny / 2;
        let ck = nz / 2;
        let cfg = RustCfg {
            n_steps,
            dipole_center_cells: (ci, cj, ck),
            dipole_length_cells: 5,
            source_freq_hz,
            ntff_surface_pad_cells: 4,
            cpml_thickness_cells: 10,
        };
        let pat = RustDriver::new(grid, cfg).run();
        (pat.theta_deg, pat.e_theta_phi0)
    });

    // Helper: find the nearest θ in the discrete array to `target_deg`.
    let sample = |target_deg: f64| -> f64 {
        theta_deg
            .iter()
            .zip(e_theta_phi0.iter())
            .min_by(|a, b| {
                (a.0 - target_deg)
                    .abs()
                    .partial_cmp(&(b.0 - target_deg).abs())
                    .unwrap()
            })
            .map(|(_, &v)| v)
            .unwrap_or(0.0)
    };

    let e0   = sample(0.0);
    let e45  = sample(45.0);
    let e90  = sample(90.0);
    let e135 = sample(135.0);
    let e180 = sample(180.0);

    let passed = e0 < 0.15
        && (e45 - 0.707_f64).abs() < 0.15
        && (e90 - 1.0_f64).abs() < 0.05
        && (e135 - 0.707_f64).abs() < 0.15
        && e180 < 0.15;

    PyDipolePatternResult {
        e_theta_0: e0,
        e_theta_45: e45,
        e_theta_90: e90,
        e_theta_135: e135,
        e_theta_180: e180,
        passed,
        theta_deg_vec: theta_deg,
        e_theta_vec: e_theta_phi0,
    }
}
```

### S2 — Register in yee-py lib.rs

In `crates/yee-py/src/lib.rs`, add after the `run_cavity_resonance` lines:

```rust
m.add_class::<fdtd::PyDipolePatternResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_dipole_pattern, m)?)?;
```

### S3 — fdtd-203 gate in yee-validation

In `crates/yee-validation/src/lib.rs`:

1. Add `run_fdtd_203_dipole_pattern()` function (returning `CaseStatus::Skipped`).
2. Call it from `Report::run_all()` alongside fdtd-201 and fdtd-201-x.

```rust
/// fdtd-203: short-dipole sin θ NTFF radiation-pattern gate (Phase 2.fdtd.py.2).
///
/// Wall-time ~30 s in release mode; registered Skipped so the default
/// `yee validate all` stays fast. Run via:
///   cargo test -p yee-fdtd --test dipole_pattern --release -- --include-ignored
/// or from Python (maturin develop --release):
///   from yee import run_dipole_pattern; assert run_dipole_pattern().passed
fn run_fdtd_203_dipole_pattern() -> CaseResult {
    CaseResult {
        id: "fdtd-203".into(),
        description: "FDTD short-dipole sin θ NTFF radiation pattern (Balanis §4.2)".into(),
        status: CaseStatus::Skipped,
        notes: "wall-time ~30 s release; cargo test -p yee-fdtd \
                --test dipole_pattern --release -- --include-ignored".into(),
        ..Default::default()
    }
}
```

### S4 — Tutorial 12

Create `docs/src/tutorials/12-fdtd-dipole-pattern-from-python.md` following
tutorial 11 as the pattern file. Key sections:
- Background: short dipole, sin θ pattern, Balanis §4.2
- Quick start (smoke with n_steps=5)
- Full gate run (n_steps=800)
- Plotting θ-cut vs analytic sin θ
- Parameters section (grid, steps, frequency knobs)
- Connection to the Rust gate (`dipole_pattern.rs`)

### S5 — SUMMARY.md

Add to the Tutorials section in `docs/src/SUMMARY.md`:
```markdown
- [FDTD dipole radiation pattern from Python](tutorials/12-fdtd-dipole-pattern-from-python.md)
```

### S6 — pytest coverage

In `crates/yee-py/tests/test_fdtd.py` (or create if it doesn't exist), add:

```python
def test_run_dipole_pattern_smoke():
    """Smoke: quick n_steps=5 run — struct round-trips, not a physics gate."""
    import yee
    r = yee.run_dipole_pattern(n_steps=5)
    assert hasattr(r, 'passed')
    assert hasattr(r, 'e_theta_90')
    assert hasattr(r, 'e_theta_0')
    arr = r.e_theta_array()
    assert len(arr) == 37
    assert r.theta_deg_array()[0] == 0.0
    assert r.theta_deg_array()[-1] == 180.0
```

---

## Verification command

```bash
# After maturin develop in the yee-py venv:
cd /home/user/Yee
source .venv/bin/activate
maturin develop -m crates/yee-py/Cargo.toml
python -c "
from yee import run_dipole_pattern
r = run_dipole_pattern(n_steps=5)
assert hasattr(r, 'passed'), 'struct missing'
print('smoke OK:', r)
"
cargo test -p yee-validation
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check --all
```

Expected exit codes: all 0.

---

## Pattern file

Imitate: `crates/yee-py/src/fdtd.rs` (the `run_cavity_resonance` block).
