# FDTD Ohmic Skin-Depth from Python

This tutorial runs the `fdtd-205` Ohmic skin-depth spatial penetration gate
from Python using the `yee` package. The gate validates that the FDTD CA/CB
Ohmic-loss E-update (Taflove §3.7) correctly reproduces exponential field
decay inside a conductor.

## 1  Background

In a good conductor with electric conductivity σ (S/m), a sinusoidal
electromagnetic wave at angular frequency ω decays exponentially with
penetration depth z (Griffiths §9.4.1):

```
|E(z)| = |E(0)| · e^{−z/δ}
```

where the **skin depth** is

```
δ = √(2 / (ω μ₀ σ))
```

At f = 1 GHz and σ = 2.5331 S/m the analytic skin depth is δ = 10 mm. With a
1 mm cell size, one skin depth spans exactly 10 cells — a convenient integer
ratio for FDTD validation.

## 2  Running the Gate

```python
from yee import run_skin_depth, SkinDepthResult

r = run_skin_depth()

assert isinstance(r, SkinDepthResult)
print(f"δ_analytic = {r.delta_analytic_m*1e3:.1f} mm")
print(f"Gate A rel_err_1δ = {r.rel_err_1delta:.2%}  (gate: < 10%)")
print(f"Gate B rel_err_2δ = {r.rel_err_2delta:.2%}  (gate: < 15%)")
print(f"Passed: {r.passed}")
print(r)
```

Expected output:

```
δ_analytic = 10.0 mm
Gate A rel_err_1δ = 1.05%  (gate: < 10%)
Gate B rel_err_2δ = 2.22%  (gate: < 15%)
Passed: True
SkinDepthResult(delta_analytic_m=1.0000e-02, rel_err_1delta=1.0500e-02,
                rel_err_2delta=2.2200e-02, passed=True)
```

## 3  Simulation Scenario

The canonical `fdtd-205` scenario (fixed — `run_skin_depth` accepts no
parameters):

| Parameter | Value |
|---|---|
| Grid | 5 × 5 × 130 cells, dx = 1 mm |
| Conductor region | z ∈ [50, 130) cells |
| Conductivity σ | 2.5331 S/m |
| Source frequency f | 1 GHz |
| Analytic skin depth δ | 10 mm = 10 cells |
| Transient | 6000 steps |
| Measurement window | 2000 steps |

A 1 GHz sinusoidal `E_x` source spans the full 5×5 transverse cross-section
at k = 25 (vacuum region). PMC boundary conditions (`H_z = 0` at y-faces)
suppress the evanescent PEC-box TM₁₁ mode (cutoff ~42 GHz >> 1 GHz) that
would otherwise dominate the exponential decay.

## 4  Gate Criteria

The ratio of measured to analytic exponential decay at each skin-depth marker:

- **Gate A** (10 %): `|ratio_1δ − e⁻¹| / e⁻¹ < 0.10`
- **Gate B** (15 %): `|ratio_2δ − e⁻²| / e⁻² < 0.15`

where `ratio_nδ = amp_nδ / amp_surface` is the peak `|E_x|` ratio between
a depth of n·δ and the conductor surface.

## 5  Result Field Reference

| Field | Type | Description |
|---|---|---|
| `delta_analytic_m` | `float` | Analytic δ = √(2/(ω μ₀ σ)) in metres |
| `amp_surface` | `float` | Peak \|E_x\| at the conductor surface |
| `amp_1delta` | `float` | Peak \|E_x\| one skin depth in |
| `amp_2delta` | `float` | Peak \|E_x\| two skin depths in |
| `ratio_1delta` | `float` | amp_1delta / amp_surface |
| `ratio_2delta` | `float` | amp_2delta / amp_surface |
| `rel_err_1delta` | `float` | \|ratio_1δ − e⁻¹\| / e⁻¹  (Gate A) |
| `rel_err_2delta` | `float` | \|ratio_2δ − e⁻²\| / e⁻²  (Gate B) |
| `passed` | `bool` | Gate A AND Gate B both pass |

## 6  Connection to the Rust Gate

The underlying Rust integration test lives in
`crates/yee-fdtd/tests/ohmic_skin_depth.rs` and calls
`yee_validation::fdtd205_run()` directly:

```bash
# Run the Rust gate directly (~8 s debug, not #[ignore]-gated):
cargo test -p yee-fdtd --test ohmic_skin_depth

# Or via the validation aggregator:
cargo run --bin yee -- validate all
```

The same `fdtd205_run()` function powers `run_skin_depth()` — the Python and
Rust results are identical.

## See also

- [FDTD absorption validation from Python](14-fdtd-absorption-validation-from-python.md)
- [FDTD TF/SF Fresnel transmission from Python](13-fdtd-fresnel-tfsf-from-python.md)
- [FDTD lossy cavity from Python](10-fdtd-lossy-cavity-from-python.md)
- [ADR-0078: fdtd-205 Ohmic skin-depth penetration gate](../decisions/0078-fdtd-205-ohmic-skin-depth.md)
- [ADR-0079: Phase 2.fdtd.py.5 skin-depth Python driver](../decisions/0079-phase-2-fdtd-py-5-skin-depth-python.md)
