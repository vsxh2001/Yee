# fdtd-202 — Lossy-Cavity Q-Factor FDTD Validation Gate

**Phase:** 2.fdtd.8  
**Date:** 2026-05-26  
**Status:** proposed  
**ADR:** [ADR-0071](../../../docs/src/decisions/0071-fdtd-202-cavity-q-factor.md)  
**Plan:** [2026-05-26-fdtd-202-cavity-q-factor.md](../plans/2026-05-26-fdtd-202-cavity-q-factor.md)

---

## 1  Motivation

The Phase 2 FDTD roadmap lists "Resonant cavity Q-factor: rectangular cavity
TE/TM modes match analytical to ±0.5%" as a validation milestone.  The
fdtd-201/201.x gates already verify resonant *frequency*; they use PEC
(lossless) walls and do not test field *decay*.  A Q-factor gate exercises a
complementary path: the Yee leapfrog must correctly model power dissipation in
a lossy dielectric, and the measured exponential ring-down must match the
analytic formula.

No prior FDTD gate has touched the lossy-medium E-update.  Adding per-cell
electric conductivity σ (Taflove §3.7 CA/CB formulation) is the minimal
infrastructure needed; it is a small, independent addition to the FDTD lane
that does not interact with the existing ADE dispersive path.

---

## 2  Physics

### 2.1  Geometry

Rectangular PEC cavity, `a × b × d = 0.20 × 0.10 × 0.20 m` (same as
fdtd-201).  All six outer walls are PEC.  The interior is filled uniformly with
ε_r = 1, σ = σ₀.

### 2.2  Resonant frequency

The dominant TE₁₀₁ mode (vacuum, same as fdtd-201):

```
f₁₀₁ = (c/2) · sqrt((1/a)² + (1/d)²)
      = (3e8/2) · sqrt((1/0.2)² + (1/0.2)²)
      = 1.0607 × 10⁹ Hz   (analytic)
```

With ε_r = 1 (vacuum permittivity), the resonant frequency is unchanged by σ
(conductivity does not shift resonance in the small-loss limit).

### 2.3  Q-factor

For a PEC-walled cavity uniformly filled with ε_r = 1 and conductivity σ
(Taflove §3.7 / Jackson §8.5):

```
Q = ε₀ · ω₁₀₁ / σ
```

where ω₁₀₁ = 2π f₁₀₁.

Design point: σ₀ = 2.96 × 10⁻³ S/m gives Q_analytic = 20 at 1.0607 GHz.

Verification:  
```
Q_analytic = ε₀ · ω₁₀₁ / σ₀
           = 8.854e-12 · 2π · 1.0607e9 / 2.96e-3
           ≈ 20.00
```

### 2.4  Ring-down measurement

After the source turns off at step N_src, the TE₁₀₁ amplitude decays as:

```
|E_y(t)| = A · exp(−ω₁₀₁ t / (2Q)) = A · exp(−t / τ)
```

where τ = 2Q / ω₁₀₁.  The FDTD measures τ from the field time series by
fitting a line to log|E_y(t)| vs t (least-squares linear regression over the
late ring-down).  Then:

```
Q_measured = ω₁₀₁ · τ / 2 = π · f₁₀₁ · τ
```

The source is placed at the TE₁₀₁ field maximum (centre of the cavity, E_y at
(nx/2, 1, nz/2)) so TE₂₀₁ is NOT excited (it has a node there).  Higher modes
(TE₁₀₃ etc.) decay faster; the fit uses only the late ring-down (t ≥ 2τ) where
they are negligible.

### 2.5  Lossy E-update (Taflove §3.7 eq. 3.56)

For a cell with ε_r and conductivity σ:

```
CA = (2 ε₀ ε_r − σ Δt) / (2 ε₀ ε_r + σ Δt)
CB = Δt / (ε₀ ε_r + σ Δt / 2)

E^{n+1} = CA · E^n + CB · curl_H^{n+1/2}
```

When σ = 0: CA = 1, CB = Δt/(ε₀ ε_r) → standard lossless update.

---

## 3  Definition of Done (DoD)

**G1 — Infrastructure:**  
`YeeGrid` gains `sigma_cells: Option<Array3<f64>>` (per-cell electric
conductivity, S/m).  A `with_sigma_cells(Array3<f64>) -> Self` builder and a
`set_sigma_box(i0,i1,j0,j1,k0,k1, sigma: f64)` helper are added.

**G2 — Lossy E-update:**  
`update_e` in `update.rs` uses the CA/CB formulation (Taflove §3.7 eq. 3.56)
when `sigma_cells` is `Some`.  When all sigma values are 0 or `sigma_cells` is
`None`, the update is bit-identical to the existing path.

**G3 — Gate test:**  
`crates/yee-fdtd/tests/cavity_q.rs` implements fdtd-202:

- Grid: 20 × 10 × 20 cells, dx = 10 mm (a = d = 0.20 m, b = 0.10 m)
- σ₀ = 2.96 × 10⁻³ S/m (target Q = 20)
- Source: Gaussian pulse centred at step 50, width 20 steps, injected into
  E_y at (10, 1, 10) (cavity centre)
- Source-on: N_src = 200 steps
- Ring-down: N_ring = 3000 steps
- Probe: E_y at (10, 1, 10) recorded for all ring-down steps
- Fitting: linear regression of log|E_y(t)| vs t over the last 2/3 of the
  ring-down (t ≥ N_ring/3 steps into ring-down → fast modes well decayed)
- Gate: |Q_measured / Q_analytic − 1| < 0.05 (±5%)

This test runs in release mode in < 1 second (3200 steps of a 20×10×20 grid)
and is **NOT** `#[ignore]`-gated.  A second `#[ignore]`-gated companion test
`fdtd_202_q_factor_hi_q_ignored` uses Q = 200 (σ = 1/10 × σ₀, 30 000 steps)
to verify the measurement at a higher Q where the ring-down is slow.

**G4 — Lossless regression:**  
A unit test `sigma_zero_matches_lossless_update` verifies that `update_e` with
uniform σ = 0 gives bit-identical results to a vanilla no-sigma grid.

**G5 — Lint and CI:**  
`cargo clippy --workspace --all-targets -- -D warnings` and
`cargo fmt --check --all` both exit 0.

---

## 4  Non-scope

- No interaction with the ADE dispersive path (Drude/Lorentz/Debye): sigma
  applies only to non-dispersive (Vacuum) cells.  Combining sigma + ADE is a
  future extension.
- No per-component sigma (isotropic scalar σ only, like the existing ε_r).
- No µ conductivity (magnetic loss).
- No yee-validation aggregator wiring in this increment (the gate is fast but
  the aggregator already Skips slow tests; a follow-on can register it).
- No CLI or Python binding changes.

---

## 5  References

- Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., §3.7 (lossy
  medium Yee update), §3.1 (Q-factor from ring-down).
- Pozar, *Microwave Engineering*, 4th ed., §6.7 (cavity Q-factor).
- Pattern file: `crates/yee-fdtd/tests/cavity_resonance.rs` (fdtd-201).
