# App.2.0 — Guided Technique-Recommender — Design Spec

**ADR:** ADR-0136 · **Date:** 2026-05-31 · **Status:** Accepted
**Vision:** `2026-05-31-ideal-filter-design-app-vision.md` §5 (the maintainer-chosen
next increment: the guided "recommend-a-technique" dual-UI entry — the Nuhertz
FilterQuick pattern, the most product-distinctive gap vs every free calculator).

## Problem

The studio's Technique stage is an **expert gallery** — the user must already know
they want edge-coupled vs lumped vs hairpin. Every free competitor is the same (pick a
calculator). The commercial tools (Nuhertz FilterQuick) add a **guided** entry: "tell
me your requirement, I recommend the topology." Yee has none. A novice with a spec but
no topology knowledge is stuck.

## Goal

Given a filter requirement (response, centre/cutoff frequency, fractional bandwidth),
**recommend a realization technique with a plain-language rationale and ranked
alternatives**, and surface it as a guided entry on the Technique stage that routes the
user into the recommended flow. Pure-domain decision logic (validatable), WASM-safe.

## Method

Two pieces, clean engine/UI split.

### 1. Engine (`yee-filter`, the validatable core)

A new pure-domain module:

```rust
/// A physical realization technique the studio can target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealizationTechnique {
    LumpedLc, EdgeCoupled, Hairpin, Combline, Interdigital, SteppedImpedance,
}

pub struct TechniqueRecommendation {
    pub primary: RealizationTechnique,
    pub rationale: String,                                   // why, citing the driver
    pub alternatives: Vec<(RealizationTechnique, String)>,   // ranked, one-line tradeoff
}

/// Recommend a realization technique from a spec. Pure; no panics on valid spec.
pub fn recommend_technique(spec: &FilterSpec) -> TechniqueRecommendation;
```

**Decision logic** (keyed on `response`, `f0_hz`, `fbw`; thresholds are documented
engineering judgment, pinned by the gate so they cannot silently drift):

- **Lowpass:** `f0_hz` (cutoff) ≥ 500 MHz → `SteppedImpedance` (distributed microstrip
  LP, alternating hi/lo-Z sections, Pozar §8.6); else `LumpedLc` (distributed sections
  impractically long below ~500 MHz). Alt: the other of the two.
- **Highpass:** `LumpedLc` (microstrip distributed high-pass is impractical; honest —
  Yee's distributed techniques are LP/BP-oriented). Alt: none distributed yet.
- **Bandpass / Bandstop:**
  - `f0_hz` < 500 MHz → `LumpedLc` (λ/4 ≈ 15 cm at 500 MHz — distributed resonators
    too large). Alt: EdgeCoupled (note: large board).
  - `f0_hz` ≥ 500 MHz:
    - `fbw` ≥ 0.20 → `EdgeCoupled` (parallel-coupled handles wide BP). Alt: LumpedLc.
    - 0.05 ≤ `fbw` < 0.20 → `EdgeCoupled` primary; `Hairpin` alt (folds the same
      resonators for a smaller board).
    - `fbw` < 0.05 → `Interdigital` primary (compact, high-Q, quarter-wave coupled —
      edge-coupled needs impractically tight gaps narrow-band); `Combline` alt
      (capacitively end-loaded, even more compact + tunable).
  - Bandstop additionally notes distributed band-stop (stub) is a future technique.

The **rationale** names the deciding factor ("At 5 GHz with 2% bandwidth, edge-coupled
gaps become impractically tight; an interdigital filter is compact and high-Q").

### 2. UI (`yee-studio-web`, a thin consumer)

A **Guided panel** at the top of the Technique stage (the expert gallery stays below =
dual-UI): a small form (response dropdown, f0/cutoff, fbw, optional stopband target) →
"Recommend" → shows the primary technique (highlighted), its rationale, and the ranked
alternatives. Each technique maps to the studio's gallery `Topology`; for **live**
techniques (EdgeCoupled, LumpedLc) a "Use this" routes into the flow (sets the topology
Signal + jumps to Spec/Synthesis); for **Soon** techniques it shows the recommendation
honestly + offers the nearest live alternative to proceed with. The recommendation also
seeds the editable `FilterSpec` from the form inputs.

Mapping `RealizationTechnique` → studio `Topology` is 1:1 for the six; live/soon status
is a UI concern (the engine stays pure-domain).

## Changes

- `crates/yee-filter/src/` — `RealizationTechnique`, `TechniqueRecommendation`,
  `recommend_technique`; re-export from the crate root; `#![warn(missing_docs)]` clean.
- `crates/yee-filter/` tests — the canonical spec→technique gate (table below).
- `crates/yee-studio-web/src/` — the Guided panel on the Technique stage + the
  `RealizationTechnique` → `Topology` map + routing for live techniques.

## DoD (machine-checkable)

1. `cargo test -p yee-filter` green, including a **non-vacuous** gate asserting each
   canonical case maps to its expected technique:
   - (Bandpass, 100 MHz, 0.05) → `LumpedLc`
   - (Bandpass, 2.4 GHz, 0.05) → `EdgeCoupled`
   - (Bandpass, 2.4 GHz, 0.25) → `EdgeCoupled`
   - (Bandpass, 5 GHz, 0.02) → `Interdigital`
   - (Lowpass, 1 GHz) → `SteppedImpedance`
   - (Lowpass, 50 MHz) → `LumpedLc`
   - (Highpass, 1 GHz) → `LumpedLc`
   Plus: every recommendation has a non-empty rationale; primary ∉ alternatives.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` + `cargo fmt --check` clean.
3. `cargo check --workspace` green.
4. `dx build --platform web --release` (working-dir `crates/yee-studio-web`,
   dioxus-cli 0.6.3 + wasm32) EXIT 0; the Guided panel is present in the built UI
   (the recommend form + a rendered recommendation block).

## Out of scope

Building the Soon techniques themselves (separate increments per the vision); auto-order
estimation tuning; ML/learned recommendation. This is a deterministic decision tree.

## Why

It turns the expert gallery into a **dual-UI** (guided + expert) — the single most
product-distinctive gap vs every free filter calculator, and pure validatable logic.
