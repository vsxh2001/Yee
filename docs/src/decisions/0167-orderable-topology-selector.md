# ADR-0167: Orderable-topology auto-selector (T3) — pick the topology that yields an orderable board

**Status:** Accepted (T3 — the ADR-0166 verdict'd follow-on)
**Date:** 2026-06-06
**Related:** ADR-0166 (T2 — `top_c_board` + `join_top_c_parts`), ADR-0165 (T1 — top-C synthesis + envelope),
ADR-0164 (JLCPCB track — `lumped_board` + `join_placed_parts` + the orderability capstone),
[[project-filter-design-final-goal]] (full pipeline spec→JLCPCB-orderable board + BOM).

---

## Context

Two lumped BPF topologies now reach the JLCPCB BOM/CPL path, covering **complementary** orderable regimes:

- **Alternating series/shunt `LumpedLadder`** (`synthesize_lumped` → `lumped_board` → `join_placed_parts`):
  orderable for **wideband** specs (ADR-0164 capstone: 1 GHz/70 % → zero blanks), blanks narrow-band (its
  *series* resonators want sub-pF/sub-nH).
- **Top-C-coupled `TopCNetwork`** (`synthesize_top_c_coupled` → `top_c_board` → `join_top_c_parts`):
  orderable for the **sub-GHz/moderate-band** regime (ADR-0166 gate: 0.5 GHz/20 % → zero blanks), blanks
  GHz-narrow (sub-pF coupling caps).

The user-facing goal is "give a spec, get an orderable board" — the user should NOT have to know which
topology their spec needs. Today they must call the right path by hand. T3 adds the **brain**: a selector
that, for a given spec, returns the topology that yields a fully-orderable board (or honestly reports that
neither lumped topology can — the distributed/planar track).

## Decision

Add an **orderable-topology selector** in `yee-filter` (pure-compute, gated; the CLI/studio wiring is the T4
follow-on so each brick stays one lane):

```rust
pub enum Topology { AlternatingLadder, TopCCoupled }

pub struct OrderableBoard {
    pub topology: Topology,
    pub board: LumpedBoard,
    pub parts: Vec<PlacedPart>,   // ready for jlcpcb_bom_csv / jlcpcb_cpl_csv
    pub fully_orderable: bool,    // true ⇔ every part resolved to a real LCSC #
}

pub fn synthesize_orderable(project: &FilterProject, footprint: Footprint)
    -> Result<OrderableBoard, LumpedError>;
```

**Policy (honest, deterministic):**
1. Try the **alternating ladder** (`synthesize_lumped` → `lumped_board` → `join_placed_parts`). If every part
   has an LCSC # → return `{ AlternatingLadder, …, fully_orderable: true }`.
2. Else try **top-C** (extract `(approx, n, f0, fbw, z0)` from the project → `synthesize_top_c_coupled` →
   `top_c_board` → `join_top_c_parts`). If every part resolves → return `{ TopCCoupled, …, true }`.
3. Else (neither fully orderable) → return the topology with the **fewer blanks** (ladder on a tie) with
   `fully_orderable: false` — an honest "no lumped topology is fully orderable for this spec; the distributed
   /planar track is the path." NEVER fake an orderable board.

The alternating ladder is tried first (the conventional/simplest topology) so wideband specs keep their
existing board; top-C is the fallback that rescues the narrow-band specs the ladder blanks on.

**Gate `topology-select-001`** (`crates/yee-filter/tests/`, pure-compute, non-`#[ignore]`'d), non-circular,
with EMPIRICALLY-CHOSEN discriminating specs:
- A **wideband** spec (e.g. 1 GHz/70 %/0402) → `AlternatingLadder` + `fully_orderable == true`.
- A spec in the **"ladder blanks but top-C is orderable"** window → `TopCCoupled` + `fully_orderable == true`.
  The implementer MUST find this spec empirically (probe the autopick on both topologies across a small
  (f0, FBW) grid) — it is the load-bearing proof that the selector's top-C fallback rescues a real spec the
  ladder can't make. If NO such spec exists in the realizable range (the two topologies never disagree),
  that is an honest surfaced finding (the selector is still correct, but the fallback is never exercised) —
  record it, do not fabricate a passing case.
- A **GHz-narrow** spec (e.g. 2 GHz/5 %/0402) → `fully_orderable == false` (honest: neither lumped topology
  resolves; distributed needed). The chosen topology = the fewer-blanks one; assert the blank set is real.

Every assertion runs the REAL `autopick`; the gate fails loudly on any mismatch between the returned
topology/flag and the actual orderability.

## Consequences

- The pipeline gains the "it just works" routing: one entry point (`synthesize_orderable`) returns an
  orderable board across the broadest spec range either lumped topology can cover, or an honest distributed
  pointer. This is the core of the spec→orderable-board goal.
- **T4 follow-on (noted, not in T3):** wire `synthesize_orderable` into `yee filter synth` (auto-route +
  report the chosen topology in the output / no-match note) and the studio export stage. Two lanes
  (`yee-cli`, `yee-studio-web`) — a separate brick.
- Scope T3: `crates/yee-filter/src/` (a `topology.rs` selector module + lib re-export) + `tests/`. Pure
  data/`f64`; WASM-safe. This ADR is the design record.
- **Not in scope:** the CLI/studio wiring (T4); a cost/size tie-breaker beyond fewer-blanks (future); the
  distributed/planar topology as a third selector arm (a separate track).

## Outcome (T3 — SHIPPED, merge `9da25b3`)

`yee_filter::{BoardTopology, OrderableBoard, synthesize_orderable}` shipped (+693, yee-filter only;
`topology.rs`). Policy as designed: ladder-first → top-C fallback (same resolved order
`project.prototype.order()`) → fewer-blanks (ladder on tie) with `fully_orderable=false`. Reviewer-confirmed
**honest — no path returns `fully_orderable=true` with a blank part** (both true-arms guarded by `blanks==0`
on the real `autopick`/`lcsc.is_some()` predicate); `BoardTopology` named to avoid the existing synthesis
`Topology{CoupledResonator}` enum.

**Gate `topology-select-001` (non-circular — independently recomputes each topology's orderability from the
lower-level synth/board/join paths, not the selector's own return):**
- wideband 1 GHz/70 %/0402 → `AlternatingLadder` + orderable (ladder 0/6).
- **The discriminating proof, 0.5 GHz/20 %/0402:** ladder blanks **4/6** (L 2.0 nH < floor, 91 nH > ceiling,
  C 1.2 pF < floor) but top-C resolves **0/10** → routes to `TopCCoupled` + orderable — the fallback
  genuinely rescues a spec the ladder can't make.
- GHz-narrow 2 GHz/5 %/0402 → ladder 4/6 & top-C 4/10 both blank → `fully_orderable=false`, real non-empty
  blank set (honest distributed pointer).

Reviewer APPROVE, no P0/P1/P2 (two P3s non-blocking: a pre-existing top-C `debug_assert` unreachable via the
public path; an intentional FR-4 substrate dup in the gate for test independence).

**T4 (the follow-on):** wire `synthesize_orderable` into `yee filter synth --jlcpcb` (auto-route the topology
+ report which topology was chosen and whether fully orderable in the output / no-match note), then the studio
export stage — so the user runs one command/click and gets an orderable board across the broadest spec range.

## References
- Topologies: `yee_filter::{synthesize_lumped, lumped_board, join_placed_parts}` (ADR-0164),
  `yee_filter::{synthesize_top_c_coupled, top_c_board, join_top_c_parts}` (ADR-0165/0166).
- Orderability check: `yee_filter::jlcpcb::autopick`; `PlacedPart.lcsc.is_some()`.
- Gate pattern: `jlcpcb-orderable-001` (ADR-0164), `top-c-board-001` (ADR-0166).
