# ADR-0100: Filter Phase F1.4.0 — `yee-export` Gerber walking skeleton

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0086 (`yee-layout` geometry + `to_svg`), ADR-0097 (F1.2.0
dimensional synthesis → `Layout`), ADR-0089 (app architecture — WASM-safe light
flow), `FILTER-DESIGN-ROADMAP.md` (F1.4 export), the locked scope (outputs
KiCad/Gerber + STEP)

---

## Context

The pipeline now produces a physical `yee_layout::Layout` (F1.2.0). The user's
end goal is **spec → manufacturable**, and the locked scope names **Gerber +
KiCad + STEP** as outputs in a new `yee-export` crate (roadmap F1.4). Nothing
emits manufacturing files yet. Per the walking-skeleton convention, the first
brick is the minimal end-to-end pipe: a `Layout` → **single-copper-layer
RS-274X Gerber** emitter. Gerber is the de-facto PCB fab interchange format and
the simplest of the three outputs to land as a self-contained, validatable text
emission (KiCad footprints/PCB and STEP solids are later increments).

## Decision

Create a new **`yee-export`** crate (pure-Rust text emission; `serde` +
`yee-layout` deps only; **WASM-safe** — the app may export client-side per
ADR-0089, so no native dep):

```rust
/// RS-274X emission options (units, layer/aperture naming).
pub struct GerberOptions { pub layer_name: String, /* mm fixed-point assumed */ }
impl Default for GerberOptions { /* "F.Cu", … */ }

/// Emit a single-copper-layer RS-274X Gerber for a layout's polygons as filled
/// regions (G36/G37). Coordinates in millimetres, `%FSLAX46Y46*% %MOMM*%`.
pub fn layout_to_gerber(layout: &yee_layout::Layout, opts: &GerberOptions) -> String;
```

The emitter writes: the `%FSLAX46Y46*%` + `%MOMM*%` header, one aperture
definition, then each `Layout` polygon as a `G36* … D02/D01 … G37*` region
(metres → mm → 4.6 fixed-point integer), and `M02*`. Single layer (the trace
copper) for the walking skeleton.

## Consequences

**Ships:** `yee-export` + `layout_to_gerber`. First manufacturing-file output;
the brick F1.4.1+ (drill, board outline, KiCad, STEP) build on.

**Gates (`yee-export` tests — Gerber is I/O, so the gate is structural validity +
coordinate round-trip, analogous to the Touchstone round-trip gate, NOT a physics
benchmark):**
- **`gerber-001` (structure):** for a known small `Layout`, the output begins
  with the `%FS…%`/`%MO…%` header, contains exactly one `G36*`/`G37*` pair per
  polygon, defines ≥1 aperture, and ends with `M02*`.
- **`gerber-002` (coordinate round-trip):** parse the `X<int>Y<int>` coordinates
  back out of one region and assert they reproduce that polygon's vertices
  (metres→mm→4.6 and back) within the format's quantisation. This validates the
  coordinate emission — the part that actually matters for fab.

**Constraint (ADR-0089):** `yee-export` stays WASM-safe (pure string emission,
`serde` + `yee-layout` only — no native/FDTD dep).

**Cross-lane (sanctioned, new crate):** adding `crates/yee-export` to the root
`Cargo.toml` `[workspace] members` (and, if used, the workspace dep table) is the
one required edit outside the crate — standard for a new crate.

**Not in scope:** KiCad footprint/PCB, STEP/3-D, drill files, board outline,
multi-layer stack, soldermask/silkscreen — all F1.4.1+. No physics.

---

## References
- Ucamco, *The Gerber Layer Format Specification* (RS-274X / X2), region (G36/G37)
  + format (`%FS…%`) + units (`%MO…%`) statements.
- `docs/superpowers/specs/2026-05-30-filter-f1-4-0-yee-export-gerber-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-f1-4-0-yee-export-gerber.md`.
