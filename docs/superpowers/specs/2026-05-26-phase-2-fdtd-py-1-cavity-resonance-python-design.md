# Phase 2.fdtd.py.1 — FDTD Rectangular-Cavity Resonance Python Driver

**Date:** 2026-05-26
**ADR:** [ADR-0073](../../../docs/src/decisions/0073-phase-2-fdtd-py-1-cavity-resonance-python.md)
**Phase:** 2.fdtd.py.1
**Status:** Accepted

---

## 1. Context

Phase 2.fdtd.py.0 (ADR-0072) shipped `run_cavity_q()` / `CavityQResult`, exposing the
fdtd-202 lossy-cavity Q-factor gate to Python. This let notebooks and tutorials verify the
lossy E-update CA/CB scheme interactively and registered fdtd-202 as the first Phase-2 FDTD
gate in `Report::run_all()`.

The **fdtd-201** gate (TE₁₀₁ rectangular cavity resonant frequency, Phase 2.fdtd,
ADR-0062) and its sibling **fdtd-201.x** (TE₂₀₁ higher-order mode, ADR-0066) are both
implemented as `#[ignore]`-gated Rust integration tests in `crates/yee-fdtd/tests/` but are
not yet:
- callable from Python,
- registered in `Report::run_all()` (so `yee validate all` does not surface them), or
- documented in the tutorial sequence.

This spec covers the work to close all three gaps: a Python driver for fdtd-201, Skipped
registrations in the aggregator for both fdtd-201 and fdtd-201.x, and tutorial 11.

---

## 2. Goal

After this phase:

1. `from yee import run_cavity_resonance` works and returns a `CavityResonanceResult` whose
   `.passed` attribute is `True` (rel_err < 2.5 %).
2. `yee validate all` lists both `fdtd-201` and `fdtd-201-x` as `SKIP` (wall-time-gated)
   with informative notes.
3. Tutorial 11 (`docs/src/tutorials/11-fdtd-cavity-resonance-from-python.md`) shows the
   complete Python workflow — physics background, code, and interpreting the output.
4. All existing gates continue to pass (`cargo check --workspace`, `cargo clippy`,
   `cargo fmt --check`, `cargo test --workspace`).

---

## 3. Design

### 3.1 Where the code lives

The implementation follows the fdtd.py.0 pattern exactly:

- **`run_cavity_resonance()`** is implemented entirely in
  `crates/yee-py/src/fdtd.rs`, consuming only the public `yee_fdtd` API
  (`YeeGrid`, `WalkingSkeletonSolver`, `boundary::apply_pec`).
- No `yee-fdtd/src/` changes are required — the simulation logic is
  self-contained in the Python layer.
- The existing `#[ignore]`-gated Rust test in
  `crates/yee-fdtd/tests/cavity_resonance.rs` is **untouched** (that is the
  run-able validation gate; the Python function is a re-expression of it).

### 3.2 Python API

```python
from yee import run_cavity_resonance

result = run_cavity_resonance()          # defaults: 20×10×20, dx=10 mm, 30 000 steps
result = run_cavity_resonance(n_steps=15_000)  # faster, still passes ±2.5 %
```

**`PyCavityResonanceResult`** fields (all exposed as `#[pyo3(get)]`):

| Field | Type | Meaning |
|-------|------|---------|
| `f_extracted_hz` | `float` | DFT-scan peak frequency (Hz) |
| `f_analytic_hz` | `float` | Analytic TE₁₀₁ `(c/2)·√((1/a)²+(1/d)²)` (Hz) |
| `rel_err` | `float` | `|f_extracted − f_analytic| / f_analytic` |
| `passed` | `bool` | `rel_err < 0.025` (fdtd-201 gate ±2.5 %) |

`probe_array() -> np.ndarray`: E_y probe time series (length = n_steps).

**`run_cavity_resonance()` keyword arguments** (all optional):

| Argument | Default | Meaning |
|----------|---------|---------|
| `nx` | 20 | cells in x |
| `ny` | 10 | cells in y |
| `nz` | 20 | cells in z |
| `dx` | 0.01 | cell size (m) |
| `n_steps` | 30 000 | FDTD time steps |
| `n_freq_bins` | 400 | DFT scan resolution |

### 3.3 Algorithm (mirrors `cavity_resonance.rs`)

1. Build `YeeGrid::vacuum(nx, ny, nz, dx)` with PEC walls.
2. Source cell: `ey[(nx/4, ny/2, nz/4)]`; probe cell: `ey[(3*nx/4, ny/2, 3*nz/4)]`.
3. Gaussian parameters: `t0 = 12·dt`, `sigma_t = 4·dt`.
4. Step loop: `update_h_only()` → `apply_pec(grid_mut())` → inject Gaussian into E_y →
   `update_e_only()` → `apply_cpml_e()` → `advance_clock()` → record probe.
5. DFT scan: `n_freq_bins` candidates in `[0.65·f_ref, 1.50·f_ref]`; peak = argmax power.
6. Gate: `|f_extracted − f_analytic| / f_analytic < 0.025`.

### 3.4 Aggregator registrations

In `crates/yee-validation/src/lib.rs`:
- Add `run_fdtd_201_cavity_resonance()` → `CaseResult` with `status: CaseStatus::Skipped`
  (wall-time ~5–15 s release, too slow for the default CI path).
- Add `run_fdtd_201x_cavity_higher_mode()` → `CaseResult` with `status: CaseStatus::Skipped`
  (wall-time similarly gated, sibling TE₂₀₁ test).
- Add both to `Report::run_all()`.

### 3.5 Tutorial

`docs/src/tutorials/11-fdtd-cavity-resonance-from-python.md`:
- Physics: PEC cavity modes, TE₁₀₁ analytic formula, grid dispersion.
- Code: `run_cavity_resonance()` call, interpreting the result, plotting the probe time series
  and DFT scan.
- Notes on the ±2.5 % gate and the ±0.5 % refinement path.

Register in `docs/src/SUMMARY.md`.

---

## 4. Validation / DoD

All of the following must pass before merge:

| Check | Command | Expected |
|-------|---------|----------|
| Build | `cargo check --workspace` | exit 0 |
| Lint | `cargo clippy --workspace --all-targets -- -D warnings` | exit 0 |
| Format | `cargo fmt --check --all` | exit 0 |
| Unit tests | `cargo test --workspace` | exit 0 |
| Python smoke | `python -c "from yee import run_cavity_resonance, CavityResonanceResult; r = run_cavity_resonance(n_steps=15000); assert r.passed, f'rel_err={r.rel_err:.4f}'"` | exit 0 |
| Aggregator | `cargo run --bin yee -- validate all 2>&1 \| grep -E 'fdtd-201'` | prints SKIP for both |

---

## 5. Out of scope

- `yee-fdtd/src/` changes (no core FDTD changes needed).
- Making fdtd-201 or fdtd-201.x a non-ignored gate (wall-time constraint; keep Skipped).
- `fdtd-201.x` Python driver (TE₂₀₁ variant; can be Phase 2.fdtd.py.2 if desired).
- GUI integration.
