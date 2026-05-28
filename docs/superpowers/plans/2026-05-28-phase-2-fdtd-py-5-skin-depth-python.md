# Plan: Phase 2.fdtd.py.5 — Ohmic Skin-Depth Python Driver

**Spec:** `docs/superpowers/specs/2026-05-28-phase-2-fdtd-py-5-skin-depth-python-design.md`
**ADR:** ADR-0079

---

## Step 1 — Add `PySkinDepthResult` and `run_skin_depth` to `fdtd.rs`

File: `crates/yee-py/src/fdtd.rs`

Append a `Phase 2.fdtd.py.5` section (after the Phase 2.fdtd.py.4 block):

- `PySkinDepthResult`: `#[pyclass]` wrapping all 9 public fields of
  `yee_validation::SkinDepthResult` via `#[pyo3(get)]`.
  Implement `__repr__` printing `delta_analytic_m`, `rel_err_1delta`,
  `rel_err_2delta`, `passed`.
- `run_skin_depth()`: `#[pyfunction]` that calls
  `yee_validation::fdtd205_run()` and maps each field onto `PySkinDepthResult`.

Pattern file: the Phase 2.fdtd.py.4 `run_cpml_reflection()` block (lines
993–999 in `fdtd.rs`).

## Step 2 — Register in `lib.rs`

File: `crates/yee-py/src/lib.rs`

Add two lines after the `PyDispersiveDrudeResult` registration block:

```rust
m.add_class::<fdtd::PySkinDepthResult>()?;
m.add_function(wrap_pyfunction!(fdtd::run_skin_depth, m)?)?;
```

## Step 3 — Export from `__init__.py`

File: `crates/yee-py/python/yee/__init__.py`

Add `SkinDepthResult` and `run_skin_depth` to the `from yee._yee import (…)`
block and to `__all__`.

## Step 4 — Add 3 pytest cases to `test_fdtd.py`

File: `crates/yee-py/tests/test_fdtd.py`

Add a `# Phase 2.fdtd.py.5` section with:

1. `test_run_skin_depth_returns_result()` — type-identity + all field names
   present + `delta_analytic_m > 0`.
2. `test_run_skin_depth_passes_gate()` — `r.passed`, `rel_err_1delta < 0.10`,
   `rel_err_2delta < 0.15`.
3. `test_run_skin_depth_repr_smoke()` — `"SkinDepthResult" in repr(r)` and
   `"delta_analytic_m" in repr(r)`.

## Step 5 — Tutorial

File: `docs/src/tutorials/15-fdtd-ohmic-skin-depth-from-python.md`

Sections:
1. Background (δ = √(2/(ωμ₀σ)), Griffiths §9.4.1)
2. Running the gate (`run_skin_depth()` code block + expected output)
3. Gate criteria (Gate A / Gate B)
4. Result field reference (table)
5. Connection to the Rust gate (`cargo test -p yee-validation …`)
6. See also (links to tutorials 13, 14 + ADR-0078 + ADR-0079)

## Step 6 — Update `docs/src/SUMMARY.md`

Add:

```md
- [FDTD Ohmic skin-depth from Python](tutorials/15-fdtd-ohmic-skin-depth-from-python.md)
```

after the tutorial-14 entry, and:

```md
- [ADR-0079: Phase 2.fdtd.py.5 skin-depth Python driver](decisions/0079-phase-2-fdtd-py-5-skin-depth-python.md)
```

after the ADR-0078 entry.

## Step 7 — Verify

```bash
cargo check -p yee-py
cargo clippy -p yee-py -- -D warnings
cargo fmt --check -p yee-py
```

Then (after `maturin develop`):

```bash
python -c "from yee import run_skin_depth; r = run_skin_depth(); assert r.passed"
pytest crates/yee-py/tests/test_fdtd.py -k skin_depth -v
```

Expected: 3/3 pass.
