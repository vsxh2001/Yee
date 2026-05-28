# ADR-0078: fdtd-205 Ohmic Skin-Depth Penetration Gate

**Date:** 2026-05-28
**Status:** Proposed → Shipped (post-merge)
**Supersedes:** N/A
**Relates to:** ADR-0071 (fdtd-202 Q-factor gate), ADR-0074 (cpml/ntff/dispersive gates)

---

## Context

`fdtd-202` (ADR-0071) validated the CA/CB Ohmic E-update in the *temporal*
domain: a lossy rectangular cavity's ring-down Q-factor matches the analytic
`Q = ωε₀/σ` to within 0.04 % (gate ±5 %).

No gate existed for the *spatial* dimension of the same update: inside a
conducting half-space, the plane-wave field amplitude should decay as
`|E(z)| = |E₀| exp(-z/δ)` where `δ = √(2/(ω μ₀ σ))` is the skin depth
(Griffiths §9.4.1 / Jackson §5.18). This is the classic textbook result for
electromagnetic penetration into a good conductor; any FDTD implementation
that correctly models Ohmic loss must reproduce it.

The fdtd-202 ring-down validates that total stored energy decays at the right
*rate*; fdtd-205 validates that the *field distribution* inside the conductor
is correct. Both are necessary for full confidence in the Ohmic update.

---

## Decision

Add validation gate **fdtd-205**: a 5×5×130-cell grid with 50 vacuum cells
(z = 0..49) and 80 conductor cells (z = 50..129, σ = 2.533 S/m, giving
δ = 10 mm = 10 cells at f = 1 GHz, dx = 1 mm). A continuous sinusoidal
E_z soft source (f = 1 GHz) is injected at cell (2, 2, 25) in the vacuum
region. After 6 000 transient steps the simulation enters the measurement
window (2 000 steps), during which the peak |E_z| is recorded at the
conductor surface (z = 50) and at depths 1δ (z = 60) and 2δ (z = 70).

**Gate A:** `|ratio_1δ − e^{-1}| / e^{-1} < 10%`
**Gate B:** `|ratio_2δ − e^{-2}| / e^{-2} < 15%`

Both gates must pass. The 15 % tolerance on Gate B is slightly looser
because the field at 2δ depth is small (~0.135 of the surface amplitude)
and phase-measurement uncertainty is proportionally larger.

Gate **NOT** `#[ignore]`'d: runs in ~0.3 s, registered in `run_all()` as
a live `CaseStatus::Passed/Failed` entry (not Skipped).

---

## Parameter derivation

```
Target: δ = 10 mm = 10 cells at f = 1 GHz, dx = 1 mm

σ = 2 / (ω μ₀ δ²)
  = 2 / (2π × 10⁹ × 4π × 10⁻⁷ × (10⁻²)²)
  = 2 / (8π² × 10⁻²)
  = 2 / 0.78957
  = 2.5331 S/m

CA stability check:
  dt = dx / (c√3) ≈ 1.926 × 10⁻¹² s
  σΔt = 2.5331 × 1.926 × 10⁻¹² = 4.879 × 10⁻¹²
  2ε₀ = 1.771 × 10⁻¹¹
  CA = (2ε₀ − σΔt) / (2ε₀ + σΔt) = 12.83/22.59 = 0.568  (positive, stable)
```

---

## Consequences

**Positive:**
- Closes the spatial-accuracy validation gap for the Ohmic E-update.
- Fast (non-`#[ignore]`'d), adds ~0.3 s to CI test time.
- Published-benchmark result (Griffiths §9.4.1): future contributors can
  independently verify the gate target.
- `pub fn fdtd205_run() -> SkinDepthResult` is available to Python bindings
  as a Phase 2.fdtd.py.5 follow-on.

**Negative:**
- Minor: 8 000 steps × 5 × 5 × 130 = 26 M cell-updates added to CI.

---

## Gate disposition

| Gate    | Threshold     | Expected result (FDTD at 10 cells/δ) |
|---------|---------------|--------------------------------------|
| Gate A  | < 10 %        | ~1–4 % (O((dx/δ)²) ≈ 1 %)           |
| Gate B  | < 15 %        | ~2–8 %                               |

**Escape hatch (if Gate A fails):** double δ to 20 cells (σ = 0.633 S/m)
and widen tolerance to 20 %. Document the finding in this ADR. Do NOT
force-pass by weakening below physics.

---

## Validation files

- Test: `crates/yee-fdtd/tests/ohmic_skin_depth.rs`
- Validation: `crates/yee-validation/src/lib.rs` (fdtd205_run + run_fdtd_205)
- Spec: `docs/superpowers/specs/2026-05-28-fdtd-205-ohmic-skin-depth-design.md`
- Plan: `docs/superpowers/plans/2026-05-28-fdtd-205-ohmic-skin-depth.md`

---

## References

- Griffiths, D. J. *Introduction to Electrodynamics*, 4th ed. §9.4.1.
- Jackson, J. D. *Classical Electrodynamics*, 3rd ed. §5.18.
- Taflove, A., & Hagness, S. C. *FDTD*, 3rd ed. §3.7 (CA/CB update).
- Pozar, D. M. *Microwave Engineering*, 4th ed. §1.7 (good-conductor limit).
