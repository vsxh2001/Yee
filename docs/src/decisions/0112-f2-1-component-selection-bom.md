# ADR-0112: Filter Phase F2.1 — component selection + BOM

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0111 (F2.0 LC ladder), the lumped-LC → PCB goal,
[[project-lumped-lc-and-studio-redesign]]

---

## Context

F2.0 produces ideal L/C element values. The lumped-LC goal requires **component
choosing** + a **BOM** + **tolerance** — none of which exist. Real filters use
purchasable standard-value parts (IEC 60063 E-series), and the part tolerance is
the input to yield analysis (F2.4).

## Decision

Add `crates/yee-filter/src/parts.rs`: `ESeries{E24,E96}` (standard preferred
values, log-nearest selection, per-series tolerance), and
`select_components(&LumpedLadder, ESeries) -> Bom` mapping each ideal L/C to its
nearest E-series value with the recorded deviation %, grouping duplicates into
quantities. `BomLine` carries optional `esr_ohm`/`srf_hz` (defaulted `None`; a
real vendor-parts library with measured parasitics is a follow-on, F2.1b).
Pure-data/math, WASM-safe.

Gate `bom_001`: E-series selection is correct on textbook anchors (incl
log-nearest tie cases) and every result is an E-series member; for the F2.0
cheb N=5 ladder every chosen value is within the E24 quantization bound; the BOM
has `2·N` parts with correct duplicate-grouping. The gate validates *selection*,
not that the quantized response still passes the mask — that yield question is
F2.4.

## Consequences

**Ships:** ideal L/C → real standard parts + a BOM with per-part tolerance.
Feeds F2.4 (Monte-Carlo over the tolerances → yield) and the UI BOM panel.

**Gate:** `cargo test -p yee-filter` green incl. `bom_001`. Pure-math.

**Not in scope:** vendor DB + measured ESR/SRF (F2.1b); quantized-response/yield
(F2.4); footprints (F2.2); UI; series beyond E24/E96.

---

## References
- `docs/superpowers/specs/2026-05-30-f2-1-component-selection-bom-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-1-component-selection-bom.md`.
- IEC 60063 preferred number series (E24/E96).
