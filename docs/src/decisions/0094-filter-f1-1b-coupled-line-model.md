# ADR-0094: Filter Phase F1.1b.gate — coupled-microstrip even/odd model

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0086 (`yee-layout` + HJ single-line model), ADR-0093 (extraction),
`FILTER-DESIGN-ROADMAP.md` (F1.1b / F1.2)

---

## Context

F1.1b.1 (the FDTD coupled-resonator driver) needs a **validatable reference**
for the coupling coefficient `k` it extracts — there is no closed-form coupled-
microstrip electrical model in-tree (only the HJ *single*-line `microstrip_width`
/ `eps_eff`). The same model is also needed by F1.2 to pick **initial** gap/width
dimensions from a target `k` before the EM-in-the-loop refinement. So a static
coupled-microstrip even/odd model is a reusable prerequisite worth landing on its
own — and, crucially, it is **pure math, validatable against published data with
no FDTD run**.

## Decision

Add a `coupled` model to **`yee-layout`** (pure f64; WASM-safe; no new dep):

```rust
pub struct CoupledMicrostrip {
    pub z0e_ohm: f64, pub z0o_ohm: f64,   // even / odd characteristic impedance
    pub eps_eff_e: f64, pub eps_eff_o: f64, // even / odd effective permittivity
}
/// Static even/odd model for symmetric edge-coupled microstrip (a cited
/// closed-form model — Garg-Bahl / Hammerstad-Jensen coupled extension; the
/// implementation states which and its accuracy).
pub fn coupled_microstrip(w_m: f64, s_m: f64, h_m: f64, eps_r: f64) -> CoupledMicrostrip;

/// Coupler-style coupling coefficient `k = (Z0e − Z0o)/(Z0e + Z0o)`.
pub fn coupling_coefficient(m: &CoupledMicrostrip) -> f64;
```

Extends `yee-layout`'s existing HJ single-line model; stays WASM-safe (it is part
of the light flow, App.1). It does NOT touch FDTD.

> **Note on the two `k`s.** `coupling_coefficient` returns the *coupler* k
> `(Z0e−Z0o)/(Z0e+Z0o)`. The *resonator* coupling that F1.1b.1 extracts from the
> two split resonant frequencies relates to the even/odd phase velocities
> (`eps_eff_e`/`eps_eff_o`) — both are exposed so F1.1b.1 can derive the resonator
> reference; this ADR validates the underlying even/odd quantities.

## Consequences

**Ships:** the `coupled` model in `yee-layout`. Gates (crate tests, §4 — pure
math, NO FDTD): **(a)** `coupled-001` — `coupled_microstrip` for a *cited*
published worked example (a specific `εr, w/h, s/h`) reproduces the published
`Z0e`/`Z0o` within the model's stated accuracy (≈ 5–10 %; the test cites source +
numbers — do NOT invent a reference, surface if none is verifiable); **(b)**
`coupled-002` — monotonicity sanity: `Z0e > Z0o > 0`, `k ∈ (0,1)`, and `k`
decreases as the gap `s` increases (textbook coupling-vs-gap law) across ≥2 gaps.

**Not in scope:** any FDTD run; the F1.1b.1 driver; the resonator-k derivation's
own validation (F1.1b.1). New: none beyond `yee-layout`.

**Constraint (ADR-0089):** `yee-layout` stays WASM-safe — pure f64, no native dep.

---

## References
- Garg & Bahl, "Characteristics of coupled microstriplines," IEEE-MTT 1979;
  Hammerstad-Jensen 1980; Pozar §8.7 (coupled lines). (Implementation cites the
  exact model + the validation data point it uses.)
- `docs/superpowers/specs/2026-05-29-filter-f1-1b-coupled-line-model-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f1-1b-coupled-line-model.md`.
