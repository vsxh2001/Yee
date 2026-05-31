# Filter F1.2.3 — Stepped-Impedance Low-Pass Filter synthesis + dimensions — Design Spec

**ADR:** ADR-0137 · **Date:** 2026-05-31 · **Status:** Accepted
**Vision:** `2026-05-31-ideal-filter-design-app-vision.md` §5 (the recommended first
breadth topology — best value-per-effort: lights a gallery "Soon" card, adds the
**low-pass** response class, and makes the App.2.0 recommender's `SteppedImpedance`
recommendation backed by a real synthesizer). Mirrors the edge-coupled / hairpin
dimensional-synthesis pattern already in `yee-filter::dimension`.

## Problem

Yee's distributed synthesis is band-pass-only (edge-coupled, hairpin). The product
vision's highest value-per-effort breadth move is the **stepped-impedance low-pass
filter** — the simplest distributed topology (alternating high-Z / low-Z microstrip
line sections), a textbook design (Pozar §8.6), and the first **low-pass** capability.
The App.2.0 recommender already recommends `SteppedImpedance` for low-pass ≥ 500 MHz,
but there is no synthesizer behind it.

## Goal

Given a low-pass prototype (order + approximation) and the line-impedance choices
(Z₀, Z_high, Z_low) on a substrate, synthesize the **alternating line sections**
(electrical lengths + physical microstrip widths/lengths) and a placeable layout —
validated against the published Pozar §8.6 worked example.

## Method

Pure closed-form, mirroring `dimension_edge_coupled` in `crates/yee-filter/src/dimension.rs`.

### Synthesis (Pozar §8.6)

From the low-pass prototype g-values (`yee_synth::prototype(approx, order)` →
`[g0, g1, …, gN, g_{N+1}]`), each reactive element `g_k` (k = 1..N) becomes one short
transmission-line section, alternating shunt-capacitor (low-Z) / series-inductor
(high-Z), **starting with a shunt capacitor (low-Z)** (the standard prototype begins
with a shunt element):

- **Shunt capacitor → low-Z line** (`Z_low`): electrical length `βl = g_k · Z_low / Z₀`.
- **Series inductor → high-Z line** (`Z_high`): electrical length `βl = g_k · Z₀ / Z_high`.

(Derivation: a high-Z line of electrical length βl looks inductive with
`L = (Z_high/ω)·βl`; matching the prototype `L = g_k·Z₀/ωc` at ω=ωc gives
`βl = g_k·Z₀/Z_high`. Dually for the low-Z capacitive line.)

### Dimensions

For each section: physical width from `yee_layout::microstrip_width(Z, εr, h)`;
guided wavelength `λg = c / (f_c · √ε_eff)` with `ε_eff = yee_layout::eps_eff(width, h,
εr)` at that section's width; physical length `l = (βl / 2π) · λg`. Mirror the
`eps_eff` / `microstrip_width` usage in `dimension_edge_coupled`.

### Types (mirror `EdgeCoupledDimensions`)

```rust
pub struct SteppedSection {
    pub high_z: bool,            // true = series-inductor high-Z line; false = shunt-cap low-Z
    pub z_ohm: f64,              // Z_high or Z_low
    pub electrical_length_rad: f64,  // βl
    pub width_m: f64,
    pub length_m: f64,
}
pub struct SteppedImpedanceDimensions {
    pub sections: Vec<SteppedSection>,   // in order, source→load
    pub eps_r: f64,
    pub h_m: f64,
}
pub fn dimension_stepped_impedance(
    proto: &yee_synth::Prototype, f_c_hz: f64, z0: f64, z_high: f64, z_low: f64,
    sub: &yee_layout::Substrate,
) -> Result<SteppedImpedanceDimensions, DimError>;
pub fn dimension_stepped_impedance_layout(/* … */) -> Result<yee_layout::Layout, DimError>;
```

## Changes

- `crates/yee-filter/src/dimension.rs` — `SteppedSection`, `SteppedImpedanceDimensions`,
  `dimension_stepped_impedance`, `dimension_stepped_impedance_layout`; re-export from
  the crate root (mirror the edge-coupled re-exports). All public items documented.
- `crates/yee-filter/tests/` — the Pozar §8.6 gate (below).
- `crates/yee-layout/**` — ONLY if a missing helper is needed (prefer composing the
  existing `microstrip_width` / `eps_eff` / `Layout` / `Substrate`).

## DoD (machine-checkable)

1. **Gate `dim_stepped_001` (Pozar Example 8.6):** maximally-flat (Butterworth) N=6,
   f_c=2.5 GHz, Z₀=50 Ω, Z_high=120 Ω, Z_low=20 Ω → the six section **electrical
   lengths** (degrees), source→load, match Pozar's published table within **±1.0°**:
   `[11.85, 33.76, 44.28, 46.12, 32.41, 12.34]`. The test must derive βl from the
   formula (not hardcode the answer) and assert against these published values.
   Non-vacuous: six distinct values from real g-values; a constant fails. Also assert
   the alternation (section 0 is low-Z / `high_z == false`) and that physical lengths
   are positive and finite.
2. `cargo test -p yee-filter` green; `cargo clippy -p yee-filter --all-targets -- -D
   warnings` + `cargo fmt --check` clean; `cargo check --workspace` green.

## Out of scope

- **Studio low-pass wiring** (lighting the `SteppedImpedance` gallery card with a live
  flow) — the studio's Spec→Synthesis is band-pass-only; threading a low-pass response
  through it is a distinct follow-on increment. This increment is the synthesis +
  dimensions + gate (the validatable core the recommender now points at).
- Elliptic stepped-Z; stub low-pass; EM verification.

## Why

The simplest distributed topology, a textbook-validatable gate, the first low-pass
capability, and it makes the shipped recommender's `SteppedImpedance` recommendation
real — the vision's top breadth pick.
