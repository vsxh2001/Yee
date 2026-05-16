# yee-core — Roadmap

Stable foundation. Move slowly here; every breakage cascades.

## Phase 0 (months 0–6)
- [ ] `units` constants + tests against CODATA 2018 values
- [ ] `FreqRange` with iteration, validation, log-spaced variant
- [ ] `Error` enum covering invalid input, numerical failure, IO bridge
- [ ] `Solver` trait skeleton
- [ ] Doc coverage: every public item documented; `cargo doc` clean

## Phase 1 (months 6–18)
- [ ] `Material` trait family (PEC, dielectric, lossy dielectric, anisotropic)
- [ ] `Port` trait (lumped, wave, modal) — abstract; impls in solver crates
- [ ] `SParameters<N>` container — n-port, complex, with reference impedances
- [ ] `FarFieldPattern` container (θ/φ grid, gain, directivity)
- [ ] `Mesh` trait (consumed by yee-mesh; abstract here)
- [ ] Convergence-estimator trait

## Phase 2 (months 18–30)
- [ ] `TimeSpan` and `Dt` types for FDTD
- [ ] Dispersive material traits (Drude / Lorentz / Debye / multi-pole)
- [ ] Source descriptors (Gaussian, modulated Gaussian, plane wave)

## Phase 3+
- [ ] Surrogate dataset trait (parameters → outputs)
- [ ] Optimization figure-of-merit trait

## Out of scope (forever)
- CUDA bindings (live in yee-cuda)
- I/O parsers (live in yee-io)
- GUI / plotting (live elsewhere)

## Validation gates per phase
- All public items documented, `cargo doc --no-deps` clean.
- `cargo test -p yee-core` passes including doc tests.
- Constants verified against CODATA reference.
- Breaking API change requires CHANGELOG entry and downstream-crate update PR in the same commit.
