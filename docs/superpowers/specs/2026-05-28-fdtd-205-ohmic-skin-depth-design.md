# fdtd-205 Ohmic Skin-Depth Penetration Gate — Design Spec

**Date:** 2026-05-28
**Phase:** 2.fdtd.9
**ADR:** 0078
**Status:** Proposed

---

## 1. Context

`fdtd-202` (ADR-0071) validated the CA/CB Ohmic-loss E-update in the
*temporal* dimension: a lossy cavity ring-down gives a Q-factor matching
the analytic `Q = ωε₀/σ` to 0.04 % (gate ±5 %).

No existing gate validates the *spatial* dimension of the same update: the
field amplitude inside a conducting half-space should decay as
`|E(z)| = |E₀| exp(−z/δ)`, where the skin depth is

```
δ = √(2 / (ω μ₀ σ))
```

(Griffiths §9.4.1 / Jackson §5.18 / Taflove §3.7). This is the
publishable-benchmark companion to `fdtd-202` and completes the Ohmic-loss
validation axis.

---

## 2. Physics

A plane wave propagating in the +z direction enters a good-conductor
half-space (`z ≥ z_surface`) with conductivity σ (SI unit S/m). In the
good-conductor approximation (`σ >> ωε₀`), the complex wave-vector is

```
k = (1 + j) / δ
```

The time-harmonic field decays as

```
E_z(z) = E₀ exp(−(z − z_surface)/δ) exp(−j(z − z_surface)/δ)
```

The *magnitude* of the field decays exponentially with skin depth δ:

```
|E(z_surface + δ)|  /  |E(z_surface)|  =  e^{−1}  ≈  0.3679
|E(z_surface + 2δ)| /  |E(z_surface)|  =  e^{−2}  ≈  0.1353
```

These ratios are the gate quantities; they are independent of the absolute
source amplitude and of the standing-wave pattern in the vacuum region.

---

## 3. Target parameters

| Parameter                 | Value                               |
|---------------------------|-------------------------------------|
| Frequency                 | 1 GHz                               |
| Cell size dx              | 1 mm                                |
| Conductivity σ            | 2.533 S/m (→ δ = 10 mm = 10 cells) |
| δ analytic                | 10 mm = 10 cells                    |
| Grid NX × NY × NZ         | 5 × 5 × 130                         |
| Vacuum region             | z = 0..49 (50 cells)               |
| Conductor region          | z = 50..129 (80 cells > 7δ)        |
| Source cell               | E_z at (2, 2, 25), sinusoidal 1 GHz |
| N_transient               | 6 000 steps (≈ 11.6 periods)        |
| N_measure                 | 2 000 steps (≈ 3.9 periods)         |
| Total N_steps             | 8 000 steps                         |
| Execution time (release)  | ≈ 0.3 s                             |

The conductor spans `z = 50..129`: 80 cells = 8δ. The field at `z = 130`
is `|E₀| × e^{−8} ≈ 3.4 × 10^{−4} × |E₀|`, so the PEC reflection from
the far wall is negligible.

### σ derivation

```
δ = 10 mm = 0.01 m
σ = 2 / (ω μ₀ δ²)
  = 2 / (2π × 10⁹ × 4π × 10⁻⁷ × (0.01)²)
  = 2 / (8π² × 10⁻²)
  = 2 / 0.78957
  = 2.5331 S/m
```

CA coefficient check (stability):

```
dt  ≈ 1.9259 × 10⁻¹² s  (dx / (c√3))
σΔt = 2.533 × 1.926 × 10⁻¹² = 4.878 × 10⁻¹²
2ε₀ = 1.771 × 10⁻¹¹
CA  = (2ε₀ − σΔt) / (2ε₀ + σΔt)
    = (17.71 − 4.878) / (17.71 + 4.878)
    = 0.568   (positive, |CA| < 1 → stable)
```

---

## 4. Gate criteria (fdtd-205)

Measurement: over the last 2 000 steps (`n = 6000..7999`), record the
**peak absolute value** of `E_z` at cells `(2, 2, 50)`, `(2, 2, 60)`,
and `(2, 2, 70)`:

```
amp_surface = max |E_z(2, 2, 50, n)| for n ∈ [6000, 7999]
amp_1delta  = max |E_z(2, 2, 60, n)| for n ∈ [6000, 7999]
amp_2delta  = max |E_z(2, 2, 70, n)| for n ∈ [6000, 7999]
```

Gate A (1δ): `|amp_1delta / amp_surface − e^{−1}| / e^{−1} < 0.10`
Gate B (2δ): `|amp_2delta / amp_surface − e^{−2}| / e^{−2} < 0.15`

Both gates must pass. Gate B has a slightly looser tolerance because the
amplitude at 2δ is smaller and rounding / phase-measurement uncertainty
is proportionally larger.

---

## 5. Complementarity with fdtd-202

| Gate      | What it validates                 | Observable        |
|-----------|-----------------------------------|-------------------|
| fdtd-202  | CA/CB Ohmic loss temporal decay   | Q-factor (ring-down τ) |
| fdtd-205  | CA/CB Ohmic loss spatial profile  | Skin-depth ratio |

These two gates together give high confidence in the full Ohmic E-update.

---

## 6. Implementation scope (LANE)

**Allowed paths:**

```
crates/yee-fdtd/tests/ohmic_skin_depth.rs     (new file)
crates/yee-validation/src/lib.rs              (register fdtd-205)
docs/src/decisions/0078-fdtd-205-ohmic-skin-depth.md (already committed)
ROADMAP.md                                    (already committed)
```

No changes to `yee-fdtd/src/` are permitted.  The implementation must use
only the existing public API:
- `YeeGrid::vacuum(NX, NY, NZ, DX)`
- `YeeGrid::set_sigma_box(i0, i1, j0, j1, k0, k1, sigma)` — **exclusive**
  upper bounds (see grid.rs loop: `i0..i1.min(ni)`)
- `WalkingSkeletonSolver::new(grid)`
- `solver.update_h_only()` / `solver.update_e_only()`
- `solver.apply_cpml_e()` — no-op without CPML; internally calls
  `boundary::apply_pec` when no CPML is attached
- `solver.advance_clock()`
- `solver.grid_mut().ez[(i,j,k)] += ...` — direct field injection
- `solver.dt()` — time step in seconds

**Pattern file:** `crates/yee-fdtd/tests/cavity_q.rs` (source injection
style, sigma cells usage, step loop structure).

---

## 7. yee-validation integration

Add to `crates/yee-validation/src/lib.rs`:

```rust
pub fn fdtd205_run() -> SkinDepthResult { ... }
fn run_fdtd_205() -> CaseResult { ... }
```

Register in `run_all()` as a **Passed/Failed** case (NOT Skipped): it runs
in ~0.3 s and does not require `#[ignore]`.

Add a `SkinDepthResult` struct analogous to `CavityQResult`:

```rust
pub struct SkinDepthResult {
    pub id: &'static str,        // "fdtd-205"
    pub delta_analytic_m: f64,
    pub amp_surface: f64,
    pub amp_1delta: f64,
    pub amp_2delta: f64,
    pub ratio_1delta: f64,       // amp_1delta / amp_surface
    pub ratio_2delta: f64,       // amp_2delta / amp_surface
    pub ratio_1delta_target: f64, // e^{-1}
    pub ratio_2delta_target: f64, // e^{-2}
    pub rel_err_1delta: f64,
    pub rel_err_2delta: f64,
    pub passed: bool,
}
```

The `fdtd205_run()` function must be `pub` (follow the `fdtd204_t_analytic`
pattern, NOT the old private helper pattern).

---

## 8. Test file structure

The test file `crates/yee-fdtd/tests/ohmic_skin_depth.rs` contains:

1. **`fn analytic_skin_depth(sigma: f64, freq: f64) -> f64`** — analytic
   formula, testable in isolation.
2. **`fn run_skin_depth_sim(...) -> (f64, f64, f64)`** — shared helper
   that returns `(amp_surface, amp_1delta, amp_2delta)`.
3. **`#[test] fn skin_depth_ratios_match_analytic()`** — the production
   gate, NOT `#[ignore]`.

The Rust test file is independent of `yee-validation` (no cross-crate
import). The validation crate duplicates the driver code or calls an
internal helper from `yee-fdtd` if one is added. Given the existing
pattern (fdtd-202 has its driver in `cavity_q.rs` + copies the logic in
`yee-validation/src/lib.rs`), the simplest approach is to write the
driver inline in both places with shared constants.

---

## 9. DoD (Definition of Done)

Machine-checkable:

1. `cargo test -p yee-fdtd --test ohmic_skin_depth` exits 0
2. `cargo test -p yee-validation` exits 0 (fdtd-205 passes)
3. `cargo clippy --workspace --all-targets -- -D warnings` exits 0
4. `cargo fmt --check --all` exits 0
5. `crates/yee-validation/src/lib.rs` contains `"fdtd-205"` in `run_all()`
   with status `CaseStatus::Passed` (not `Skipped`)
6. `fdtd205_run().passed == true`

---

## 10. Risk and escape hatch

**Risk:** The 10% gate for Gate A might not hold due to FDTD numerical
dispersion with 10 cells/δ.

**Mitigation:** CA = 0.568 (positive, well-behaved). The FDTD numerical
dispersion error for the skin depth at 10 cells/δ is O((dx/δ)²) = 1% —
well within the 10% gate.

**Escape hatch (blocked > 15 min):** If Gate A fails, increase δ to
20 cells (σ = 0.633 S/m) and re-try; if still failing, gate at 20%
tolerance and document the finding in the ADR. Do NOT force the gate to pass
by weakening tolerance below the physics (the measurement should agree within
2× of the analytic).

---

## 11. References

- Griffiths, D. J. *Introduction to Electrodynamics*, 4th ed. §9.4.1 —
  "Electromagnetic Waves in Conductors", skin depth formula (eq. 9.25).
- Jackson, J. D. *Classical Electrodynamics*, 3rd ed. §5.18 — skin-depth
  formula and limiting forms.
- Taflove, A., & Hagness, S. C. *Computational Electrodynamics: FDTD*,
  3rd ed. §3.7 — CA/CB Ohmic-loss update (eqs. 3.36–3.37).
- Pozar, D. M. *Microwave Engineering*, 4th ed. §1.7 — complex propagation
  constant in a lossy medium (corroborates the good-conductor approximation).
