# Phase 1.validation.2 — FDTD Aggregator Gate Integration

**Date:** 2026-05-27  
**Status:** Proposed  
**Relates to:** ADR-0074  
**Paired plan:** `docs/superpowers/plans/2026-05-27-phase-1-validation-2-fdtd-aggregator-gates.md`

---

## 1. Problem

Three validation cases — `cpml-001`, `ntff-001`, and `dispersive-001` — are registered in
`Report::run_all()` but hardcoded to `CaseStatus::Skipped` with notes reading
*"deferred to Phase 1.validation.2"*.  They were deferred from Phase 1.validation.1 because
the physics helpers were not yet available in a form the aggregator could call without
depending on test-internal helpers.

The underlying FDTD integration tests (`cpml_reflection.rs`, `ntff_dipole.rs`,
`dispersive.rs`) **already pass** in `cargo test`; the aggregator just doesn't run them.

This means `yee validate all` silently skips 3 passing physics checks, creating a false
impression of gap coverage.

---

## 2. Decision

Port the physics from the three existing yee-fdtd integration tests directly into
`crates/yee-validation/src/lib.rs`, following the same pattern as `fdtd-202`
(physics helpers duplicated inline; no cross-crate test-module deps).

All three gates are fast enough (< 10 s debug, < 2 s release) to be
non-`#[ignore]`-gated in the aggregator.

---

## 3. Gate specifications

### 3.1 cpml-001 — CPML attenuates ≥ 30 dB vs PEC

**Source:** `crates/yee-fdtd/tests/cpml_reflection.rs`  
**Grid:** 50³ cells, dx = 1 mm  
**Runs:** two 300-step simulations sharing the same Gaussian parameters
- PEC reference: `WalkingSkeletonSolver::new(grid)` (hard walls)
- CPML candidate: `WalkingSkeletonSolver::with_cpml(grid, CpmlParams::for_grid(&grid, 10))`  
**Source:** `E_z` at (25, 25, 25); Gaussian, t₀ = 20 dt, σ = 6 dt  
**Probe:** `E_z` at (38, 25, 25) — 2 cells inside inner PML edge on +x face  
**Measurement:** ratio of late-time peak (steps 100–300) vs early-time peak (step 5) in dB  
**Gate:** CPML late-time peak / PEC late-time peak ≤ −30 dB  
**Wall time (estimated):** < 0.3 s debug, < 0.1 s release  
**Reference:** Roden & Gedney 2000; Taflove & Hagness §7.9

### 3.2 ntff-001 — NTFF broadside/endfire null ≥ 20 dB

**Source:** `crates/yee-fdtd/tests/ntff_dipole.rs`  
**Grid:** 50³ cells, dx = 1 mm  
**Run:** 2000 steps, CPML (npml = 10), NTFF DFT at f = 15 GHz  
**Source:** `E_z` at (25, 25, 25); Gaussian, t₀ = 12 dt, σ = 4 dt  
**NTFF surface:** box_margin_cells = 15 (= npml + 5)  
**Measurement:** |E_far(θ=π/2, φ=0)| vs |E_far(θ=0, φ=0)|  
**Gate:** broadside / endfire ≥ 20 dB  
**Wall time (estimated):** < 2 s debug, < 0.5 s release  
**Reference:** Balanis §4.2 (infinitesimal dipole, null along axis)

### 3.3 dispersive-001 — Drude slab Fresnel reflection within 20 %

**Source:** `crates/yee-fdtd/tests/dispersive.rs::drude_slab_reflects_per_fresnel`  
**Grid:** 80³ cells, dx = 1 mm  
**Runs:** two 800-step simulations with CPML (npml = 10)
- Vacuum reference: no material map
- Drude slab: 20-cell slab at i ∈ [50, 70), ω_p = 2π × 20 GHz, γ = 0 (lossless)  
**Source:** `E_z` at (20, 40, 40); probe: `E_z` at (30, 40, 40)  
**Measurement:** DFT at 10 GHz of (slab − vacuum) trace; 1/r correction for point-source
geometry; |Γ_measured| vs |Γ_analytic| = |(1 − n) / (1 + n)| where n = ε(ω)^(1/2)  
**Gate:** |Γ_measured − Γ_analytic| / Γ_analytic ≤ 0.20  
**Wall time (estimated):** < 5 s debug, < 1 s release  
**Reference:** Drude permittivity ε(ω) = 1 − ω_p²/ω²; Fresnel normal incidence

---

## 4. Implementation scope

**Lane:** `crates/yee-validation/src/lib.rs`, `docs/src/decisions/0074-*.md`,
`docs/src/SUMMARY.md`  
**Out of lane:** yee-fdtd (read-only), yee-core (read-only); no changes to test files

**Pattern file:** the `fdtd202_run()` + `run_fdtd_202_lossy_cavity_q()` block in
`crates/yee-validation/src/lib.rs` (lines ~1314–1496)

---

## 5. Definition of Done

All three gates are non-Skipped and Passed in `yee validate all`:

```
cpml-001     PASS   ...
ntff-001     PASS   ...
dispersive-001 PASS ...
```

Verification command (non-ignored, fast path):
```bash
cargo test -p yee-validation --release -- --nocapture 2>&1 | grep -E "cpml-001|ntff-001|dispersive-001|PASS|FAIL"
```

Lint floor:
```bash
cargo clippy -p yee-validation -- -D warnings
cargo fmt --check -p yee-validation
```

---

## 6. Risks

- **Wall time too slow for CI**: all three gates are estimated < 5 s debug; if a gate
  exceeds 30 s debug it should be `#[ignore]`-gated and marked Skipped instead.
- **API surface**: the test files use some lower-level functions (`update::update_h`,
  `sources::gaussian_pulse_ez`) that are `pub` in yee-fdtd. Verify they remain accessible
  from yee-validation.
- **Dispersive CPML interaction**: the Drude slab must not extend into the CPML region
  (ε_r must stay at its CPML value in the PML cells) — this constraint is documented in
  `heterogeneous_substrate.rs` and must be replicated correctly.
