# Filter F1.2.5 — Combline dimensional synthesis — Design Spec

**ADR:** ADR-0144 · **Date:** 2026-05-31 · **Status:** Accepted
**Maintainer pick:** combline (AskUserQuestion 2026-05-31), with the explicit caveat that
it ships **only** with a PROPER published-design gate (not a shallow λg/4 mirror). A
3-source research sweep (workflow `waws53n82`, confidence 0.9) confirmed the synthesis +
sourced a citeable Hong & Lancaster benchmark, so it is cleanly gateable.

## Problem

The studio recommends combline for narrow-band high-Q band-pass but has no engine; the
gallery card is "Soon". Combline = capacitively-loaded, short-circuited coupled
resonators (compact, high-Q, good spurious performance). This is the engine.

## Method (mirror `dimension_hairpin`; reuse the validated coupling/gap machinery)

A combline resonator is a short-circuited microstrip line (impedance `Z0`, electrical
length `θ0 < 90°` at `f0`) loaded by a shunt capacitor `C_L` at the open end. Adjacent
resonators couple via the line-to-line gap (coupled-line even/odd), exactly as
edge-coupled/hairpin. So the **coupling realization reuses** `solve_gap` /
`coupled_microstrip` / `coupling_coefficient` (validated `coupled_002`); the
combline-**distinct** pieces are the `θ0` short-circuited resonator + the loading cap.

### Synthesis (Hong & Lancaster §5.2.5)

- Coupling: `target_k[i] = FBW · m[i][i+1]` (the studio's existing normalized coupling
  matrix; equals H&L's `M_{i,i+1} = FBW/√(g_i·g_{i+1})`) → `solve_gap`. External Q
  `Qe = g0·g1/FBW` (already on the synthesized `CouplingMatrix`).
- Resonator: short-circuited line, electrical length `θ0` (a design choice, default
  **45° = λg/8** for compactness), physical length `L = θ0/β(f0)` with
  `β(f0) = 2π·f0·√ε_eff/c` (`ε_eff` at the synthesized width).
- Loading cap (H&L eq 5.43, per-resonator case): `C_L = cot(θ0)/(2π·f0·Z0)`.

### Types (mirror `HairpinDimensions`)

```rust
pub struct ComblineDimensions {
    pub line_width_m: f64,         // microstrip_width(Z0, εr, h)
    pub theta0_rad: f64,           // chosen resonator electrical length at f0 (< π/2)
    pub resonator_length_m: f64,   // L = θ0/β(f0)
    pub loading_cap_f: f64,        // C_L = cot(θ0)/(2π·f0·Z0) — the short-circuited end is grounded (via)
    pub gaps_m: Vec<f64>,          // N−1 inter-resonator gaps (solve_gap)
    pub target_k: Vec<f64>,        // N−1 target couplings (= FBW·m[i][i+1])
}
pub fn dimension_combline(
    project: &FilterProject, theta0_rad: f64, substrate: &Substrate,
) -> Result<ComblineDimensions, DimError>;
```

`θ0_rad` must be in `(0, π/2)` (else `cot ≤ 0` → non-physical cap → `DimError`).

## Changes

- `crates/yee-filter/src/dimension.rs` — `ComblineDimensions`, `dimension_combline`;
  re-export from the crate root. All public items documented.
- `crates/yee-filter/tests/` — the `dim_combline_001` gate (below).

## DoD (machine-checkable, NON-vacuous)

**Gate `dim_combline_001`:**

1. **Published-benchmark (H&L eq 5.46) — the non-tautological core.** Synthesize the
   5-pole 0.1 dB Chebyshev (the studio's `synthesize`/`prototype(Chebyshev{0.1}, 5)`)
   at FBW=0.1; assert the synthesized external Q and inter-resonator couplings match
   H&L's published combline design within tolerance: **Qe ≈ 11.468**, **M₁₂=M₄₅ ≈
   0.07975**, **M₂₃=M₃₄ ≈ 0.06077** (tol ±1% — tighten to ±1e-3 if the crate's g-values
   reproduce H&L's exactly). Second point (FBW=0.15 pseudocombline): **Qe ≈ 7.645**,
   **M₁₂ ≈ 0.11962**, **M₂₃ ≈ 0.09115**. These are specific published external numbers
   exercising the full g→Qe/M chain — a constant/wrong synthesis fails.
2. **Combline-distinct resonance consistency (first-principles, NOT the cap tautology).**
   For the dimensioned resonator (`θ0`, `Z0`, `L`, `C_L`), compute the loaded
   short-circuited-stub input susceptance **independently**
   `B(f) = −Y0·cot(β(f)·L) + 2π·f·C_L` with `β(f)·L = θ0·(f/f0)`, root-find `B(f)=0`
   over a band around f0, and assert the root equals **f0** within ±1% — re-deriving
   resonance from the admittance (catches a wrong cap / length / dispersion / sign), NOT
   asserting `C_L == cot/…` back.
3. Coupling gaps are bracketed/solved (no clamping); `θ0 ∈ (0, π/2)` enforced (a
   `θ0 ≥ π/2` input → `DimError`); all physical dims positive + finite; `C_L > 0`.
4. `cargo test -p yee-filter` green; `cargo clippy -p yee-filter --all-targets -- -D
   warnings` + `cargo fmt --check`; `cargo check --workspace`.

## Out of scope

The **studio lighting** of combline (a follow-on, like stepped-Z's engine→studio split);
discrete E-series selection of `C_L` (reuse the lumped BOM later); via/short-circuit 3-D
modelling; the rigorous Getsinger/Cristal self-/mutual-capacitance coupled-bar synthesis
(H&L eq 5.44 — this first-order engine reuses the proven `solve_gap` coupling realization
like hairpin, and documents that). EM verification (ADR-0133 wall) untouched.

## Why

Combline is the maintainer's chosen next technique. The research found a **proper,
non-tautological gate** (H&L eq 5.46 published Qe/M + a first-principles resonance check),
so it ships honestly — the compact, high-Q narrow-band band-pass realization, reusing the
validated coupling machinery + adding the combline-distinct short-circuited θ0 resonator
and loading cap.

## References
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications* (Wiley 2001),
  §5.2.5 Combline (eqs 5.42–5.46, design example 5.46) + §5.2.6 pseudocombline (eq 5.48).
- Matthaei, Young & Jones, *Microwave Filters…* Ch. 8 (combline; Getsinger/Cristal).
- Research workflow `waws53n82` (3-source, confidence 0.9).
