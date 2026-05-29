# ADR-0103: Filter Phase F1.4.1a — `yee-export` Gerber board outline (Edge.Cuts)

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0100 (F1.4.0 `layout_to_gerber` copper layer), ADR-0086
(`yee-layout` `BBox`), ADR-0089 (WASM-safe), `FILTER-DESIGN-ROADMAP.md` (F1.4)

---

## Context

F1.4.0 emits one copper layer. A fabricable board also needs a **board outline**
(the profile / KiCad `Edge.Cuts` layer) — the cut path the fab uses to route the
PCB shape. It is the next-smallest brick toward a real fab set and, like the
copper emitter, is pure-text + WASM-safe.

## Decision

Add `layout_to_gerber_outline` to `yee-export`:

```rust
/// Options for the board-outline Gerber.
pub struct OutlineOptions { pub layer_name: String, pub margin_mm: f64 }
impl Default for OutlineOptions { /* "Edge.Cuts", 1.0 mm */ }

/// Emit an RS-274X board-outline Gerber: a single closed rectangular contour
/// around the layout `bbox` expanded by `margin_mm` on each side, stroked (not
/// region-filled) with a thin aperture. mm / 4.6 fixed-point, same as the copper
/// emitter.
pub fn layout_to_gerber_outline(layout: &yee_layout::Layout, opts: &OutlineOptions) -> String;
```

Unlike the copper regions (`G36/G37` fill), the outline is a **stroked contour**:
`%FSLAX46Y46*% %MOMM*%`, a thin aperture (`%ADD10C,0.100*%` + `D10*`), a `D02`
move to the first corner, `D01` draws around the four `bbox±margin` corners and
back to the first, then `M02*`. Reuses the F1.4.0 metres→mm→4.6 conversion.

## Consequences

**Ships:** `layout_to_gerber_outline` + `OutlineOptions`. The copper (F1.4.0) +
outline (this) are the two layers a minimal single-sided board needs. The
CLI/studio wiring to write both as a bundle is a follow-on (kept out so the crate
stays pure/WASM-safe — file writing is the caller's job).

**Gates (`yee-export` tests):**
- **`gerber-003` (outline structure):** header + a `D02` move + four `D01` draws
  forming a closed contour + `M02*`; exactly one aperture; no `G36`/`G37` (it is
  stroked, not filled).
- **`gerber-004` (outline geometry):** the four emitted corner coordinates equal
  `bbox.min/max ± margin` (metres→mm→4.6 round-trip), confirming the profile
  encloses the layout with the requested margin.

**Constraint (ADR-0089):** stays WASM-safe (pure `String`, `yee-layout` only).

**Not in scope:** drill, multi-layer copper, soldermask/silkscreen, KiCad/STEP,
non-rectangular outlines, any file/bundle writer — F1.4.1b+.

---

## References
- ADR-0100 (copper emitter + coordinate model); Ucamco RS-274X (stroked contour
  vs `G36` region).
- `docs/superpowers/specs/2026-05-30-filter-f1-4-1-gerber-board-outline-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-f1-4-1-gerber-board-outline.md`.
