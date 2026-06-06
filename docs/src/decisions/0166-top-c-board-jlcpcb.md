# ADR-0166: Top-C-coupled BPF → board + JLCPCB orderable upload set (T2)

**Status:** Accepted (T2 — the ADR-0165 verdict'd follow-on)
**Date:** 2026-06-06
**Related:** ADR-0165 (T1 — `yee_filter::top_c` synthesis + S21-validate + envelope; verdict GO for T2),
ADR-0164 (JLCPCB track — `jlcpcb_export` BOM/CPL, the orderability capstone + its narrow-band caveat),
F2.2 `board.rs` (`lumped_board` → `Layout` + `Placement`), [[project-filter-design-final-goal]].

---

## Context

ADR-0165 T1 shipped the top-C-coupled (capacitively-coupled) lumped BPF **synthesis** (`synthesize_top_c_coupled`
→ `TopCNetwork { shunt: Vec<(L,C)>, coupling_caps_farad: Vec<f64> }`) + a non-circular S21 validation, and
mapped the **realizability envelope**: top-C **extends** the JLCPCB-orderable regime to the sub-GHz / wider-band
corner ((0.2 GHz, 20 %) and (0.5 GHz, 20 %) synthesize to fully-orderable component sets) — the regime the
existing alternating-series/shunt `LumpedLadder` can NOT make orderable (its *series* resonators want
sub-pF/sub-nH there). T1's verdict was **GO for T2: wire the topology into the board + JLCPCB export**, scoped
to that sub-GHz/wider-band orderable regime.

The gap T2 closes: T1 proved the *component values* are orderable, but `TopCNetwork` is NOT yet placeable or
exportable. The JLCPCB path (`lumped_board` → `Placement` list → `join_placed_parts` → `jlcpcb_bom_csv` /
`jlcpcb_cpl_csv`) is built for the `LumpedLadder` (alternating series/shunt). Top-C has a **different
topology**: `N` shunt L–C resonators to ground + `N+1` **series** coupling caps in the through-arm between
(and at the ends of) the resonator nodes. So it needs its own board placement + value-join, then it reuses
the existing CSV emitters.

## Decision

Add the top-C → board → JLCPCB data path in `yee-filter` (data path only; CLI/studio topology-selection is
the T3 follow-on). Reuse the ADR-0164 `Placement` / `PlacedPart` / `jlcpcb_bom_csv` / `jlcpcb_cpl_csv`
machinery; add only the two top-C-specific pieces + a gate.

1. **`top_c_board(net: &TopCNetwork, substrate, footprint) -> LumpedBoard`** (`board.rs`): place the `N` shunt
   resonators along a line (each a `Lk` + `Ck` pad pair, `BranchKind::Shunt`) and the `N+1` series coupling
   caps in the through-arm between/at the ends of the nodes (ref-des `Cc1..Cc(N+1)`, `BranchKind::Series`),
   emitting the renderable `Layout` (pads) + the `Placement` list. Mirror `lumped_board`'s pad/placement
   construction (the PATTERN); the difference is the topology (shunt resonators + series coupling caps, vs
   alternating).
2. **`join_top_c_parts(placements: &[Placement], net: &TopCNetwork) -> Vec<PlacedPart>`** (`jlcpcb_export.rs`):
   join each placement's ref-des to its value — `Lk`/`Ck` → `net.shunt[k-1]`, `Cc j` → `net.coupling_caps[j-1]`
   — and `autopick` the LCSC part (reuse `autopick` + `PlacedPart` + `value_comment`). Unmatched values emit
   the same honest blank-LCSC `(NO BASIC PART)` row as `join_placed_parts` (never dropped/faked).
3. **Reuse** `jlcpcb_bom_csv` / `jlcpcb_cpl_csv` unchanged on the resulting `Vec<PlacedPart>`.

**Gate `top-c-board-001`** (`crates/yee-filter/tests/`): for a sub-GHz/moderate-band spec the T1 envelope
identified as orderable (e.g. Cheb 0.5 dB, N=3, **f0=0.5 GHz, FBW=20 %, 50 Ω, 0402**): `synthesize_top_c_coupled`
→ `top_c_board` → `join_top_c_parts` → the BOM is **FULLY orderable (zero blank LCSC #s)** — every shunt L/C
AND every coupling cap a real Basic part within tolerance; the CPL designators match the BOM; placements lie
within the board outline. Non-circular (the board/join consume the synthesized network, the gate asserts real
autopick coverage) + honest (if a cell the T1 probe called orderable does NOT fully resolve here, that is a
recorded discrepancy, not a weakened gate). This is the top-C analogue of `jlcpcb-orderable-001`.

## Consequences

- **Closes the narrow-band manufacturability gap** the ADR-0164 capstone flagged: the JLCPCB data path now
  has a topology (top-C) that yields a fully-orderable board in the sub-GHz/moderate-band regime where the
  alternating ladder blanks. The pipeline gains a second orderable topology covering a complementary regime.
- **T3 follow-on (noted, not in T2):** CLI `yee filter synth --topology top-c` + studio topology selection
  (ideally auto-pick the topology that yields an orderable board for the given spec) — the user-facing wiring.
- Scope T2: `crates/yee-filter/src/{board.rs, jlcpcb_export.rs, lib.rs}` + `crates/yee-filter/tests/`. Pure
  data / `f64` / `Complex64`; **WASM-safe** (the studio consumes it later). This ADR is the design record.
- **Not in scope:** the CLI/studio topology flag (T3); GHz-narrow top-C (T1 proved it still blanks on sub-pF
  coupling caps — distributed-only); J5 Gerber-completeness (ADR-0164, separate); distinct Q.

## Outcome (T2 — SHIPPED, merge `b764fc0`)

`yee_filter::top_c_board` + `join_top_c_parts` shipped (+725, yee-filter only). `top_c_board` realizes the
top-C schematic — `N+1` series coupling caps in-line on the through-arm interleaved with `N` shunt L‖C
resonators tapping to ground between them (`port1—Cc1—node1(L1‖C1)—Cc2—…—Cc{N+1}—port2`, reviewer-confirmed
electrically correct vs the `top_c_s21` ABCD cascade + the copper geometry) → renderable `Layout` +
`Placement` list. `join_top_c_parts` joins each ref-des (`Cc{j}` stripped before bare `C{k}`, panic-safe) →
value → `autopick`, reusing `jlcpcb_bom_csv`/`jlcpcb_cpl_csv` + the honest blank path.

**Gate `top-c-board-001` (non-circular, honest, non-vacuous):** a Cheb 0.5 dB/N=3/0.5 GHz/20 %/50 Ω/**0402**
top-C BPF → a **FULLY-ORDERABLE upload set, ZERO blank LCSC #s across all three arms** — coupling caps
0.96 pF→C1550 & 2.4 pF→C1559, shunt C 3.3 pF→C1565 & 4.3 pF→C1569, shunt L 16 nH→C27143 (every one a real
bundled Basic part within the 20 % band); CPL designators == BOM; placements within outline. **Closes the
narrow-band manufacturability gap** the ADR-0164 capstone flagged: top-C makes the sub-GHz/moderate-band
regime orderable where the alternating ladder blanks. Honest regime bound: **0402 required** (the 0603
RF-inductor grid jumps 12→22 nH, blanking the 16 nH shunt L — documented, not weakened). Reviewer APPROVE,
no P0/P1/P2 (two P3s non-blocking: E24 hardcoded — autopick re-snaps; bbox-outline check loose — ok for a
walking-skeleton board).

**T3 (the follow-on):** CLI/studio topology selection — ideally a `yee_filter` orderable-topology selector
(try the alternating ladder; if its BOM blanks, try top-C; return whichever is fully orderable + which
topology, honest "neither → distributed" otherwise) consumed by `yee filter synth` + the studio, so the
pipeline auto-emits an orderable board for the broadest spec range.

## References
- T1 synthesis: `yee_filter::top_c` (`synthesize_top_c_coupled`, `TopCNetwork`, `top_c_s21`) — ADR-0165.
- Reuse: `yee_filter::board` (`lumped_board`, `Layout`, `Placement`, `Footprint`, `BranchKind`, `PadSpec`),
  `yee_filter::jlcpcb_export` (`join_placed_parts`, `PlacedPart`, `value_comment`, `jlcpcb_bom_csv`,
  `jlcpcb_cpl_csv`), `yee_filter::jlcpcb::autopick`.
- Gate pattern: `jlcpcb-orderable-001` (ADR-0164 orderability capstone).
