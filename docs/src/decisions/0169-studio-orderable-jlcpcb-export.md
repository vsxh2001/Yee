# ADR-0169: Studio Export — auto-routed orderable JLCPCB assembly set (T5)

**Status:** Accepted (T5 — the ADR-0168 verdict'd follow-on)
**Date:** 2026-06-06
**Related:** ADR-0168 (T4 — CLI `--jlcpcb` auto-routes), ADR-0167 (T3 — `synthesize_orderable[_on]`),
ADR-0166 (T2 — `top_c_board`), ADR-0130 (Dioxus `yee-studio-web` is THE studio), ADR-0164 (J4 — studio
JLCPCB download buttons), [[project-filter-design-final-goal]] (the app is the deliverable).

---

## Context

The live Dioxus studio (`yee-studio-web`, `/Yee/studio/`) lumped-LC Export stage offers JLCPCB BOM/CPL
download buttons (J4, ADR-0164) — but they hardcode the **alternating ladder** at **0603**:
`join_placed_parts(&d.board.placements, &d.ladder, Footprint::Smd0603, ESeries::E24)`. So a narrow-band spec
the ladder can't make orderable (its series resonators want sub-pF/sub-nH) silently downloads a half-blank
BOM, and the T1-T4 top-C auto-route never reaches the app user — even though `synthesize_orderable_on`
(T3/T4) would route it to an orderable top-C realization.

Two coupled facts shape the fix: (1) the orderability-deciding topology is **spec-dependent** (ladder for
wideband, top-C for sub-GHz/moderate-band); (2) it is also **footprint-dependent** — the top-C narrow-band
envelope is orderable on **0402** (finer RF L/C value grid) but blanks on 0603 (the 0603 inductor grid jumps
12→22 nH). The studio's value-add over the CLI (which takes an explicit `--footprint`) is to **auto-find an
orderable realization** — search topology AND footprint — and surface it honestly.

## Decision

1. **Engine (`engine.rs`): `orderable_upload(project: &FilterProject, substrate: &Substrate) -> OrderableExport`** —
   searches footprints `[Smd0402, Smd0603, Smd0805]` (0402 first — finer grid, smaller board, the RF-lumped
   norm), calling `synthesize_orderable_on(project, substrate, fp)` for each; returns the FIRST
   `fully_orderable` result, else the one with the FEWEST total blanks (0402 on a tie). `OrderableExport`
   (pure data, WASM-safe):
   ```rust
   pub struct OrderableExport {
       pub topology_label: &'static str,   // "alternating ladder" | "top-C-coupled"
       pub footprint_label: &'static str,  // "0402" | "0603" | "0805"
       pub fully_orderable: bool,
       pub n_parts: usize,
       pub n_blank: usize,
       pub board: LumpedBoard,             // the auto-routed board (Gerber + CPL geometry)
       pub bom_csv: String,                // jlcpcb_bom_csv(&ob.parts)
       pub cpl_csv: String,                // jlcpcb_cpl_csv(&ob.board.placements)
   }
   ```
   Carry `orderable: OrderableExport` on `LumpedDesigned` (computed once in the lumped engine path).
2. **Export stage (`stages.rs`): a dedicated "JLCPCB orderable assembly set" subsection** replacing the two
   standalone hardcoded JLCPCB buttons. A self-consistent set — **all four files from the same auto-routed
   board**: JLCPCB BOM (`orderable.bom_csv`), JLCPCB CPL (`orderable.cpl_csv`), Gerber F.Cu + Edge.Cuts (from
   `orderable.board.layout`). A badge/fields: `"Orderable realization: {topology_label} · {footprint_label}"`
   and either `"✓ {n_parts}/{n_parts} parts on the JLCPCB Basic catalog"` or the honest
   `"⚠ {n_blank} of {n_parts} parts have no Basic match — narrow-band lumped is distributed-only (see the
   distributed techniques)"`. The existing "design board" Gerber/KiCad buttons (the displayed ladder board)
   stay, with a note clarifying the JLCPCB set is the auto-routed orderable realization (it may differ from
   the displayed ladder when the ladder isn't orderable; routing the whole lumped flow — incl. the finite-Q
   response — is a documented follow-on, blocked on a top-C lossy response).
3. **No fabricated orderability:** the badge reflects `fully_orderable` exactly; blanks are surfaced, never
   hidden.

## Validation

- **Engine unit test (`yee-studio-web`, fast, non-`#[ignore]`'d, NON-circular):** call `orderable_upload`
  (or the footprint-search helper) for three projects: a **wideband** spec → `topology_label == "alternating
  ladder"` + `fully_orderable`; a **0.5 GHz/20 %** spec → `topology_label == "top-C-coupled"` +
  `footprint_label == "0402"` + `fully_orderable` + `n_blank == 0` (the T3 discriminating cell, reached via
  the studio engine); a **2 GHz/5 %** spec → `!fully_orderable` + `n_blank > 0` (honest). Assert the
  `bom_csv`/`cpl_csv` are non-empty + the CPL/BOM designator sets are consistent. Mirrors `cli-jlcpcb-autoroute`
  at the studio-engine layer.
- The **`wasm-build` CI job** stays green (`synthesize_orderable_on` + the JLCPCB CSV emitters are pure-math,
  already WASM-safe; no new non-WASM dep). clippy/fmt clean.

## Consequences

- The deployed studio's lumped Export now delivers an **auto-routed, self-consistent, orderable** JLCPCB
  upload set — the app realizes the spec→production-ready-board goal for the broadest spec range, honestly
  bounded. This is the App-track surfacing of the T1-T4 engine work.
- Scope T5: `crates/yee-studio-web/src/{engine.rs, stages.rs}`. Pure data / WASM-safe. This ADR is the design
  record.
- **Not in scope:** routing the lumped DESIGN view + response/finite-Q to top-C (blocked on a top-C lossy
  `top_c_s21` finite-Q variant — a documented follow-on); a user-facing footprint picker (the search auto-
  picks); J5 Gerber-completeness (ADR-0164).

## Outcome (T5 — SHIPPED, merge `8167798`)

`yee-studio-web::{OrderableExport, orderable_upload}` shipped (+357/−43, 2 studio src files). `orderable_upload`
searches `[0402, 0603, 0805]` (0402-first) via `synthesize_orderable_on`, returns the first `fully_orderable`
result else the fewest-blanks one (0402 on tie); Err-skips a footprint; honest empty export if all error.
Carried on `LumpedDesigned.orderable`, computed in `design_lumped_from`. The Export stage's hardcoded-ladder/
0603 JLCPCB buttons are replaced by a **self-consistent "JLCPCB orderable assembly set"** — BOM + CPL + both
Gerbers ALL from the one auto-routed `orderable.board` (reviewer-confirmed: no display/export board mismatch)
+ a badge reflecting `fully_orderable` **exactly** (✓ N/N orderable, or ⚠ M-of-N blank → distributed-only;
no fabricated orderability). The displayed ladder board + finite-Q response stay ladder-based (documented
limit — routing the whole lumped flow is blocked on a top-C lossy response).

**Engine gate (non-circular, mirrors `cli-jlcpcb-autoroute`):** wideband 1 GHz/70 % → alternating ladder/
0402/orderable/0-blank; **0.5 GHz/20 % → top-C/0402/orderable/0-blank** (the rescue the old fixed-0603-ladder
path could NOT make orderable); 2 GHz/5 % → fewest-blanks/NOT-orderable/4-blank. 17 tests green; no new dep;
`cargo check --target wasm32-unknown-unknown` exit 0 (the `wasm-build` CI gate stays green). Reviewer APPROVE,
no P0/P1/P2 (two cosmetic P3s).

**⇒ The top-C arc is COMPLETE end-to-end (T1→T5): spec → orderable JLCPCB upload set now reaches the user via
BOTH the CLI (`yee filter synth --jlcpcb`) AND the deployed studio app** — auto-routing ladder↔top-C across
wideband + sub-GHz/moderate-band, honestly bounding GHz-narrow as distributed-only. **NEXT candidates:**
ADR-0164 **J5 Gerber-completeness** (F.Mask/F.Silk soldermask+silkscreen layers — the last real-fab-upload
artifact; lumped SMD has no vias → no drill needed); routing the whole lumped flow to top-C (needs a top-C
lossy/finite-Q response first); a `--topology` manual override.

## References
- Selector: `yee_filter::{synthesize_orderable_on, OrderableBoard, BoardTopology}` (ADR-0167/0168).
- CSVs/board: `yee_filter::{jlcpcb_bom_csv, jlcpcb_cpl_csv, lumped_board, top_c_board, LumpedBoard, Footprint}`.
- Studio: `crates/yee-studio-web/src/{engine.rs (LumpedDesigned, LUMPED_FOOTPRINT, design path), stages.rs (export_lumped)}`.
