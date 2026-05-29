# Filter Phase F1.1b.gate — coupled-microstrip even/odd model — Design Spec

**Phase:** F1.1b.gate · **ADR:** ADR-0094 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
A static even/odd electrical model for symmetric edge-coupled microstrip in
`yee-layout` — the validatable `k` reference F1.1b.1 (FDTD driver) needs and the
initial-dimensioning model F1.2 uses. Pure f64, WASM-safe, no FDTD, no new dep.

## API (`yee-layout`, new `coupled` module re-exported at crate root)
```rust
pub struct CoupledMicrostrip { pub z0e_ohm, pub z0o_ohm, pub eps_eff_e, pub eps_eff_o: f64 }
pub fn coupled_microstrip(w_m: f64, s_m: f64, h_m: f64, eps_r: f64) -> CoupledMicrostrip;
pub fn coupling_coefficient(m: &CoupledMicrostrip) -> f64;  // (z0e − z0o)/(z0e + z0o)
```

### Model
Implement a **cited** static coupled-microstrip even/odd model (e.g. Garg-Bahl
1979, or the Hammerstad-Jensen single-line `eps_eff`/`Z0` reused per mode with a
coupled even/odd correction; OR Akhtarzad). Document in the source which model +
its published accuracy. Reuse the existing `yee-layout` HJ single-line helpers
(`microstrip_width`/`eps_eff`) where the model builds on the single-line result.
Pure f64; doc every public item.

## DoD (machine-checkable; pure math, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-layout --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-layout` exit 0 (fast).
4. `coupled_001_vs_published` (`tests/`): for a **cited** published worked example
   (a specific `εr, w/h, s/h` with a literature `Z0e`/`Z0o`), `coupled_microstrip`
   reproduces `Z0e` and `Z0o` within the model's stated tolerance (state it, ≈
   5–10 %). The test MUST cite the source + the exact reference numbers in a
   comment. **Do NOT invent reference values** — if no verifiable published point
   is found (use WebSearch), STOP and surface (escape hatch), do not fabricate.
5. `coupled_002_monotonic` (`tests/`): `coupled_microstrip` over ≥2 gaps `s1<s2`
   (same w,h,εr): assert `z0e>z0o>0` and `eps_eff_e>0, eps_eff_o>0` for each, and
   `coupling_coefficient` is positive and **strictly decreases** as the gap grows
   (`k(s1) > k(s2) > 0`). Pure sanity from the textbook coupling-vs-gap law.

## Out of scope
Any FDTD run; the F1.1b.1 coupled-resonator driver; the resonator-k derivation +
its FDTD validation; asymmetric/offset coupled lines; full Kirschning-Jansen
dispersion (a static model suffices for the gate reference).
