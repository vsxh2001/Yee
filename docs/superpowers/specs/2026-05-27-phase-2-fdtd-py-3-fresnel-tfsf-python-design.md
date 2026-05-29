# Phase 2.fdtd.py.3 Design Spec  
## FDTD TF/SF Fresnel-Transmission Python Driver (fdtd-204 gate)

**Date:** 2026-05-27  
**Phase:** 2.fdtd.py.3  
**ADR:** 0076  
**Paired plan:** `2026-05-27-phase-2-fdtd-py-3-fresnel-tfsf-python.md`

---

## 1. Context

Phases 2.fdtd.py.0–2 exposed three FDTD capabilities to Python:

| Phase       | Gate     | Physics                              | Python fn               |
|-------------|----------|--------------------------------------|-------------------------|
| py.0        | fdtd-202 | Lossy-cavity Q-factor ring-down      | `run_cavity_q`          |
| py.1        | fdtd-201 | Rectangular-cavity TE₁₀₁ resonance   | `run_cavity_resonance`  |
| py.2        | fdtd-203 | Short-dipole NTFF sin-θ pattern      | `run_dipole_pattern`    |

**fdtd-204** closes the next gap: the TF/SF plane-wave source (Phase 2.fdtd.5.3.2) and the per-cell ε_r infrastructure (MMMMMMMM) have both shipped, but they have **never been tested together** in a quantitative validation case. A normal-incidence plane-wave transmission experiment through a lossless dielectric slab is the canonical test for this combination.

---

## 2. Goal

Ship `run_fresnel_tfsf()` / `FresnelTfsfResult` in `yee-py` and register `fdtd-204` in the `yee-validation` aggregator — following the pattern established by py.0/1/2.

The gate validates:
1. `PlaneWaveSource` (TF/SF) + `YeeGrid::with_eps_r_cells` (per-cell ε_r) work together correctly.
2. The measured amplitude transmission coefficient through a dielectric slab matches the analytic transfer-matrix prediction within 5%.

---

## 3. Physics Design

### 3.1 Geometry

| Parameter           | Value                                               |
|---------------------|-----------------------------------------------------|
| Grid                | 80 × 80 × 80 cells, dx = 1 mm                       |
| CPML                | npml = 10, all six faces                            |
| TF box              | i₀=12, i₁=69, j₀=1, j₁=78, k₀=1, k₁=78            |
| Dielectric slab     | i ∈ [50, 55), ε_r = 2.2, everywhere else ε_r = 1.0 |
| TF/SF frequency     | 10 GHz (CW sinusoid, ramp_steps = 50)               |
| Probe (incident)    | (25, 40, 40) — vacuum, inside TF box                |
| Probe (transmitted) | (62, 40, 40) — vacuum after slab, inside TF box     |
| n_steps             | 600                                                 |
| Settling window     | steps 0–199 discarded; measurement over [200, 600)  |

### 3.2 Measurement Protocol

Two identical runs with the same `PlaneWaveSource`:

1. **Vacuum run**: no `eps_r_cells` — records E_z at probe_inc (25, 40, 40).
2. **Slab run**: `with_eps_r_cells` applied (ε_r = 2.2 in cells i ∈ [50, 55), 1.0 elsewhere including all CPML cells) — records E_z at probe_trans (62, 40, 40).

```
A_inc   = max |E_z_vacuum(probe_inc)| over steps [200, 600)
A_trans = max |E_z_slab(probe_trans)| over steps [200, 600)
t_measured = A_trans / A_inc
```

### 3.3 Analytic Reference

Standard transfer-matrix formula (Born & Wolf §1.6.2, Pozar §2.3) for
normal incidence on a slab of thickness d = 5 mm, ε_r = 2.2, f = 10 GHz:

```
n₂ = √2.2 ≈ 1.4832
δ  = 2π f n₂ d / c ≈ 1.553 rad  (≈ π/2 quarter-wave)

r₁₂ = (n₁ − n₂)/(n₁ + n₂) = (1 − 1.4832)/(1 + 1.4832) ≈ −0.1946
r₂₃ = (n₂ − n₃)/(n₂ + n₃) = (1.4832 − 1)/(1.4832 + 1)   ≈ +0.1946
t₁₂ = 2n₁/(n₁ + n₂) ≈ 0.8054
t₂₃ = 2n₂/(n₂ + n₃) ≈ 1.1946

t_slab = t₁₂ · t₂₃ · e^{jδ} / (1 + r₁₂ · r₂₃ · e^{j2δ})
|t_analytic| ≈ 0.927  (power transmission T_power = |t|² ≈ 0.859)
```

The slab is near-quarter-wave at 10 GHz, giving a sub-unity transmission
(27% field reduction) that is well measurable while avoiding the near-zero
minimum of a half-wave transformer.

### 3.4 Gate Criterion

`|t_measured / t_analytic − 1| < 0.05`  (5%)

Rationale: the `heterogeneous_substrate.rs` point-source Fresnel test uses
±5% for the same ε_r=2.2 material; the TF/SF case is comparable (and in
some ways cleaner since the incident field is a true plane wave).

### 3.5 Why This Avoids Fabry-Perot / CPML-Mismatch Problems

- **5-cell slab** ends at i=55, well before the CPML boundary at i=70: no
  dielectric-CPML impedance mismatch.
- **Probe at i=62** is inside the TF box (i₁=69) and 7 cells after the slab
  back face: both the TF "total-field" condition and the vacuum measurement
  condition hold simultaneously.
- **Transmission measurement** (not reflection): the probe in the transmitted
  vacuum region after the slab directly measures A_trans. For a two-run
  comparison, A_trans / A_inc = |t_slab| without any subtraction artefact.

---

## 4. Implementation Scope

### 4.1 Lane

`crates/yee-py/src/**, crates/yee-validation/src/**, docs/**`

Out-of-lane (report as finding, do NOT touch):
- `crates/yee-fdtd/src/**` (no new FDTD source changes needed)
- Other crates

### 4.2 Files to create or modify

| File | Action |
|------|--------|
| `crates/yee-validation/src/lib.rs` | Add `fdtd204_run()` inline helper + `run_fdtd_204()` + register `fdtd-204` Skipped in `run_all()` |
| `crates/yee-py/src/fdtd.rs` | Add `PyFresnelTfsfResult` struct + `run_fresnel_tfsf()` PyO3 fn |
| `crates/yee-py/src/lib.rs` | Register new class + fn in Python module |
| `crates/yee-py/python/yee/__init__.py` | Add to imports + `__all__` |
| `crates/yee-py/tests/test_fdtd.py` | 2 pytest cases (smoke + gate) |
| `docs/src/decisions/0076-phase-2-fdtd-py-3-fresnel-tfsf-python.md` | New ADR |
| `docs/src/tutorials/13-fdtd-fresnel-tfsf-from-python.md` | New tutorial |
| `docs/src/SUMMARY.md` | Add tutorial 13 + ADR-0075 (missing) + ADR-0076 entries |

---

## 5. DoD (Definition of Done)

All items machine-checkable:

1. `cargo clippy --workspace --all-targets -- -D warnings` exits 0
2. `cargo fmt --check --all` exits 0
3. `cargo test -p yee-validation --lib -- fdtd_204 --release` exits 0 (unit smoke)
4. `cargo test -p yee-validation -- fdtd_204_live_gate --release --ignored` exits 0 (gate passes, ≤5% rel err)
5. `cargo test -p yee-py --test test_fdtd -- test_run_fresnel_tfsf` exits 0 (pytest via maturin)
6. `grep 'FresnelTfsfResult' crates/yee-py/python/yee/__init__.py` exits 0
7. `grep 'run_fresnel_tfsf' crates/yee-py/python/yee/__init__.py` exits 0
8. `grep 'fdtd-204' crates/yee-validation/src/lib.rs` exits 0
9. `grep 'ADR-0075' docs/src/SUMMARY.md` exits 0
10. `grep 'ADR-0076' docs/src/SUMMARY.md` exits 0

---

## 6. Risks

| Risk | Mitigation |
|------|-----------|
| CW amplitude not stable in settling window | Increase `n_steps` to 800 and settling to 300 if needed |
| Gate rel_err > 5% | Loosen to 10% (still physically meaningful; document in ADR) |
| TF/SF + per-cell eps_r interaction unexpected | Check against `plane_wave_propagation.rs` quiet-zone contrast as a sanity baseline |
