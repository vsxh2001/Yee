# Phase 2.fdtd.py.1 Implementation Plan

**Spec:** [2026-05-26-phase-2-fdtd-py-1-cavity-resonance-python-design.md](../specs/2026-05-26-phase-2-fdtd-py-1-cavity-resonance-python-design.md)
**ADR:** [ADR-0073](../../docs/src/decisions/0073-phase-2-fdtd-py-1-cavity-resonance-python.md)

---

## Step 1 — `PyCavityResonanceResult` struct + helper

In `crates/yee-py/src/fdtd.rs`, after the existing `PyCavityQResult` impl block:

1. Add helper function `cr_analytic_f101(nx, nz, dx) -> f64`:
   `c0 / 2 * sqrt((1/(nx*dx))^2 + (1/(nz*dx))^2)`.
2. Add `#[pyclass(name = "CavityResonanceResult", module = "yee._yee")]` struct with
   fields: `f_extracted_hz`, `f_analytic_hz`, `rel_err`, `passed: bool`,
   `probe_vec: Vec<f64>` (internal).
3. Add `#[pymethods] impl PyCavityResonanceResult`:
   - `probe_array<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>>` — returns
     the E_y probe time series as numpy array.
   - `__repr__` string.

## Step 2 — `run_cavity_resonance()` function

In `crates/yee-py/src/fdtd.rs`, after `run_cavity_q`:

1. Add `#[pyfunction]` with signature matching spec §3.2 defaults.
2. Implement the algorithm from spec §3.3:
   - Build `YeeGrid::vacuum(nx, ny, nz, dx)`.
   - Compute source/probe indices, Gaussian params.
   - Step loop with `update_h_only` / `apply_pec` / inject / `update_e_only` /
     `apply_cpml_e` / `advance_clock`.
   - Accumulate probe time series.
   - DFT scan over n_freq_bins candidates in `[0.65·f_ref, 1.50·f_ref]`.
   - Return `PyCavityResonanceResult`.
3. Include `#[allow(deprecated)]` on `apply_pec` call (same pattern as `run_cavity_q`).

## Step 3 — Wire into lib.rs

In `crates/yee-py/src/lib.rs`:
- `m.add_class::<fdtd::PyCavityResonanceResult>()?;`
- `m.add_function(wrap_pyfunction!(fdtd::run_cavity_resonance, m)?)?;`

## Step 4 — Wire into __init__.py

In `crates/yee-py/python/yee/__init__.py`:
- Add `CavityResonanceResult` to the import list from `yee._yee`.
- Add `run_cavity_resonance` to the import list.
- Add both to `__all__`.

## Step 5 — Pytest cases

In `crates/yee-py/tests/test_fdtd.py` (or `tests/test_cavity_resonance.py`), add:
1. `test_cavity_resonance_passes_gate()` — call `run_cavity_resonance(n_steps=15_000)`;
   assert `result.passed`, `result.rel_err < 0.025`.
2. `test_cavity_resonance_result_fields()` — verify field types and probe_array shape.
3. `test_cavity_resonance_repr_smoke()` — verify `__repr__` contains key fields.

## Step 6 — Aggregator registrations

In `crates/yee-validation/src/lib.rs`:
1. Add `fn run_fdtd_201_cavity_resonance() -> CaseResult` returning a Skipped result with
   descriptive notes.
2. Add `fn run_fdtd_201x_cavity_higher_mode() -> CaseResult` returning a Skipped result.
3. Add both to `Report::run_all()` between `run_fdtd_202_lossy_cavity_q()` and
   `run_fem_eig_001()`.

## Step 7 — Tutorial 11

Create `docs/src/tutorials/11-fdtd-cavity-resonance-from-python.md` covering:
- Physics: PEC cavity TE₁₀₁ analytic formula (Pozar §6.3), Gaussian broadband excitation,
  DFT scan, grid dispersion.
- Code walkthrough: `run_cavity_resonance()` call, printing results, plotting probe series
  and DFT scan with matplotlib.
- Gate: ±2.5 % interpretation; refinement path to ±0.5 %.
- Reference to the `#[ignore]`-gated Rust test as the full-resolution gate.

Register in `docs/src/SUMMARY.md` after the line for tutorial 10.

## Verification

```bash
cargo check --workspace                                         # exit 0
cargo clippy --workspace --all-targets -- -D warnings          # exit 0
cargo fmt --check --all                                         # exit 0
cargo test --workspace                                          # exit 0
# After maturin develop:
python -c "
from yee import run_cavity_resonance, CavityResonanceResult
r = run_cavity_resonance(n_steps=15000)
assert r.passed, f'FAIL: rel_err={r.rel_err:.4f}'
print(f'PASS: f_extracted={r.f_extracted_hz*1e-9:.4f} GHz, rel_err={r.rel_err:.4e}')
"
cargo run --bin yee -- validate all 2>&1 | grep -E 'fdtd-201'
# expected: two SKIP lines
```
