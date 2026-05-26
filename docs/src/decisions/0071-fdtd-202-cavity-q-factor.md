# ADR-0071 — fdtd-202 Lossy-Cavity Q-Factor Validation Gate

**Status:** Accepted (2026-05-26)  
**Phase:** 2.fdtd.8

---

## Context

The Phase 2 FDTD roadmap lists "Resonant cavity Q-factor: rectangular cavity
TE/TM modes match analytical to ±0.5%" as a validation milestone.  The
fdtd-201 (TE₁₀₁ frequency, ADR-0062) and fdtd-201.x (TE₂₀₁ frequency,
ADR-0066) gates verify resonant *frequency* only; both use PEC (lossless)
walls and do not test field *decay*.

No prior FDTD gate exercises the lossy-medium E-update path.  The standard
Taflove §3.7 CA/CB formulation (Yee 1966, Taflove & Hagness eq. 3.56a–b) adds
a per-cell electric conductivity σ to the E update:

```
CA = (2 ε₀ ε_r − σ Δt) / (2 ε₀ ε_r + σ Δt)
CB = Δt / (ε₀ ε_r + σ Δt / 2)
E^{n+1} = CA · E^n + CB · curl_H^{n+1/2}
```

This is distinct from the ADE dispersive path (Drude/Lorentz/Debye in
`crates/yee-fdtd/src/dispersive.rs`), which handles frequency-dependent ε_r.
For a simple DC conductivity the CA/CB update is exact and requires no
auxiliary variables.

The analytic Q-factor for a PEC-walled cavity uniformly filled with ε_r = 1
and conductivity σ is:

```
Q = ε₀ · ω₁₀₁ / σ
```

which is clean, parameter-free, and independent of geometry beyond the choice
of mode (TE₁₀₁).

---

## Decision

1. **Add `sigma_cells: Option<Array3<f64>>`** to `YeeGrid` (per-cell electric
   conductivity, S/m).  Add `with_sigma_cells` builder and `set_sigma_box`
   helper, matching the existing `eps_r_cells` / `mu_r_cells` pattern.

2. **Modify `update_e`** in `crates/yee-fdtd/src/update.rs` to apply the CA/CB
   formulation when `sigma_cells` is `Some`.  When σ = 0 or `sigma_cells` is
   `None`, the update is bit-identical to the prior path.  A unit test
   (`sigma_zero_matches_lossless_update`) enforces this.

3. **Implement `crates/yee-fdtd/tests/cavity_q.rs`** (gate `fdtd-202`):
   - Geometry: 20 × 10 × 20 cells, dx = 10 mm (a = d = 0.20 m, b = 0.10 m)
   - Fill: ε_r = 1, σ₀ = 2.96 × 10⁻³ S/m (target Q = 20 at f₁₀₁ = 1.0607 GHz)
   - Source: Gaussian pulse (N_src = 200 steps) at cavity centre; TE₂₀₁ NOT
     excited (node at source location)
   - Ring-down: N_ring = 3000 steps; probe E_y at centre
   - Fit: log-linear regression over last 2/3 of ring-down → τ → Q_meas
   - Gate: `|Q_meas / Q_analytic − 1| < 0.05` (±5 %)
   - Wall time: < 1 s release; **NOT** `#[ignore]`-gated
   - Companion: `fdtd_202_q_factor_hi_q_ignored` (Q = 200, 30 000 steps,
     `#[ignore]`) for regression coverage at high Q

4. **Write ADR-0071** (this document) and register it in `docs/src/SUMMARY.md`.

The ±5 % gate is deliberately loose relative to the roadmap's ±0.5 %.  The
Yee scheme has second-order phase error in the resonant frequency, which
introduces a small bias in the extracted decay rate.  The ±5 % target is
achievable without mesh refinement and matches the fdtd-201 frequency gate
philosophy.  Tightening to ±0.5 % on a refined mesh is a follow-on.

---

## Consequences

- **Positive:** The lossy E-update path is exercised and validated for the
  first time.  The `sigma_cells` API unlocks lossy-dielectric FDTD simulations
  (e.g., modelling conductor-backed absorbers, skin-depth approximations) for
  future users without requiring the ADE machinery.
- **Scope preserved:** The ADE dispersive path (`dispersive.rs`) is not
  modified.  Combining σ + ADE (e.g., a Drude metal with background
  conductivity) is a future extension.
- **CI impact:** `fdtd_202_q_factor_lossy_cavity` runs without `#[ignore]` in
  < 1 s (3200 steps of a 4000-cell grid); the hi-Q companion is `#[ignore]`-gated.
- **Aggregator:** The fast gate can be registered in `yee-validation::run_all`
  in a follow-on (after the aggregator pattern is confirmed to handle tests
  that take < 1 s).

---

## References

- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., §3.7 (lossy
  medium update), §3.1 (Q definition).
- Pozar, *Microwave Engineering*, 4th ed., §6.7 (cavity Q).
- ADR-0062: fdtd-201 cavity resonance gate (frequency, lossless).
- ADR-0066: fdtd-201.x higher-order cavity mode gate.
- Spec: `docs/superpowers/specs/2026-05-26-fdtd-202-cavity-q-factor-design.md`
- Plan: `docs/superpowers/plans/2026-05-26-fdtd-202-cavity-q-factor.md`
