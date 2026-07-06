# ADR-0178: `yee.compute` — Python bindings for the GPU/CPU engine

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0175..0177 (engine track), Phase 1.frontend.0 (yee-py shape).

## Context

The engine (`yee-compute`) shipped E.0–E.5a with Rust-side gates but had no scripting surface.
`yee-py` already exposes solver-specific drivers (`yee.fem`, FDTD gate runners, surrogate).

## Decision

New `yee.compute` submodule following the established register-and-sys.modules idiom. Shape: a
config-accumulating `FdtdSim` builder (`set_eps_r_cells` / `set_mu_r_cells` / `set_sigma_cells`
as `(nx+1,ny+1,nz+1)` float64 numpy maps, `set_pec_mask` per E component in its staggered
shape, `set_boundary("none"|"pec"|"cpml", npml, axes)`, `add_gaussian_source`,
`add_resistive_port`, `add_probe`) and an immutable `FdtdResult` (`probe(i)` 1-D arrays,
`field("ex".."hz")` 3-D staggered arrays, `backend`). `run(n_steps, backend="cpu"|"gpu"|"auto")`
builds a fresh engine per call (the `PyFdtdDriver` idiom) and releases the GIL for the whole
solve (`Python::detach`); `"auto"` falls back to CPU only on `NoAdapter`.

## Gate

`crates/yee-py/tests/test_compute.py` (6 tests, runs in the existing CI `python-bindings` job):
shapes/dtypes/staggering, propagation smoke, bit-identical repeatability, material + CPML +
mask path (masked cells verifiably clamped), shape/enum validation errors, GPU↔CPU probe
agreement within 1e-3 relative (skips without an adapter — green on llvmpipe here), and
auto-fallback. Physics itself is gated in Rust (compute-001..010); these gate the *binding*.

## Consequences

- `import yee; yee.compute.FdtdSim(...)` scripts the engine end-to-end, GPU included — wheels
  gain the wgpu dependency tree (size/build-time cost accepted).
- The studio's Python-side experimentation path and notebook tutorials can target one API.
