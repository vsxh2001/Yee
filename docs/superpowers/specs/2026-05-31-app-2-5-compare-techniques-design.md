# App.2.5 — Compare techniques side-by-side — Design Spec

**ADR:** ADR-0142 · **Date:** 2026-05-31 · **Status:** Accepted
**Origin:** maintainer pick ("deepen the flows — optimize/compare") at the
complete-app boundary. The product vision (§1 dual entry, §5 P2 interactivity) wants
the user to *see their options*, not just pick one blind. Complements the App.2.0
recommender (recommend one) and the expert gallery (pick one) with **compare them all**.

## Problem

The Technique stage offers a guided recommendation + an expert gallery, but to compare
techniques the user must select each, walk its flow, and remember the board size /
verdict. There is no side-by-side view of "for *this* spec, here's how edge-coupled vs
hairpin vs lumped compare."

## Goal

For the current spec, synthesize **every live technique that realizes the spec's response
class** and show a side-by-side comparison — board size, PASS/FAIL, and the key graded
metrics — with the recommended technique marked and a "Use this" that routes into each.
Real engine output for every row; no fabricated numbers.

## Method

A pure, host-testable helper + a Compare panel on the Technique stage.

### Engine (`engine.rs`)

```rust
pub struct TechniqueComparison {
    pub technique: RealizationTechnique,
    pub realizable: bool,                 // false → the design failed to dimension
    pub board_w_mm: f64,
    pub board_h_mm: f64,
    pub pass: Option<bool>,               // None when not realizable
    pub order: usize,
    pub worst_passband_ripple_db: f64,
    pub worst_return_loss_db: f64,
    pub worst_stopband_rej_db: Option<f64>,
}
/// Synthesize every live technique that realizes `spec`'s response class and
/// collect a comparable metric row for each (real engine output).
pub fn compare_techniques(spec: &FilterSpec) -> Vec<TechniqueComparison>;
```

`compare_techniques` keys on `spec.response`:
- `Bandpass | Bandstop` → `[EdgeCoupled, Hairpin, LumpedLc]`: run `design_demo_from(spec,
  EdgeCoupled)` / `(.., Hairpin)` (metrics from `Designed.report` + `board_size_mm`,
  `realizable = layout.is_some()`), and `design_lumped_from(spec)` (metrics from
  `LumpedDesigned.verdict` + `board_size_mm`; `Err`/`None` → `realizable=false`).
- `Lowpass` → `[SteppedImpedance]`: `design_stepped_from(spec)` (metrics from the stepped
  fields + `board_size_mm`).
- `Highpass` → `[]` (no live technique yet).

Metrics are pulled directly from each design's existing graded structs (the same fields
`verify_view` uses) — `worst_passband_ripple_db`, `worst_return_loss_db`, stopband
rejection (min achieved over the stopband table, or `worst_stopband_rej_db` for lumped,
`None` when absent), `pass`, `order`. Pure (no signal reads).

### UI (`stages.rs`)

A **Compare** panel on the Technique stage (below the guided panel + the gallery): a
table — one row per `compare_techniques(spec)` entry — columns: technique name, board
size (mm × mm), PASS/FAIL chip (or "not realizable"), worst ripple, worst RL, stopband
rejection. The **recommended** technique (from `recommend_technique(spec).primary`,
mapped to live) is marked. Each realizable row has a "Use this" that sets the topology
signal + routes to Spec (reuse `route_into`). When only one technique realizes the
response (low-pass), show that single row with a note; when none (high-pass), an honest
"no live technique for this response yet" note. The spec is the live `spec` signal the
stage already holds.

## Changes

- `crates/yee-studio-web/src/engine.rs` — `TechniqueComparison`, `compare_techniques`
  (pure, documented) + a test.
- `crates/yee-studio-web/src/stages.rs` — the Compare panel in `technique_stage`.
- (No `main.rs` change expected — the Technique stage already has `topology`/`active`/
  `spec`.)

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): `compare_techniques`
   on a band-pass demo spec returns the three band-pass techniques (EdgeCoupled,
   Hairpin, LumpedLc) with REAL metrics pulled from each design (assert each row's
   metrics equal that technique's design's graded fields), and the rows are **not all
   identical** (e.g. the hairpin and edge-coupled board sizes differ — hairpin folds
   smaller — or the techniques differ). A low-pass spec returns exactly the
   `SteppedImpedance` row; a high-pass spec returns `[]`. A constant/empty
   `compare_techniques` fails.
2. `dx build --platform web --release` EXIT 0; the Compare panel renders the table.
3. Existing tests pass; `cargo clippy ... -D warnings` + `cargo fmt --check` clean;
   `cargo check --workspace` green.

## Out of scope

A live response-overlay plot of all techniques (a follow-on); tune-and-watch sliders
(the spec form already re-derives live); cross-response comparison (compare is within a
response class). No new physics.

## Why

Turns "pick a technique blind" into "see your options side-by-side for this spec" — the
product-vision dual-entry completion the maintainer chose, built entirely on the existing
per-flow engines, with a non-vacuous gate.
