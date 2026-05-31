# ADR-0145: Filter F1.2.6 — Combline board layout

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0144 (F1.2.5 — the combline synthesis engine; this is its layout
companion), ADR-0109/0138 (hairpin — `dimension_hairpin_layout` + the App.2.1 lighting
pattern this enables), ADR-0139 (stepped-Z — the primitive-composed layout pattern reused),
[[lumped-lc-and-studio-redesign]]

---

## Context

`dimension_combline` (F1.2.5) shipped, but combline has no board layout — so it cannot be
lit in the studio (which renders a `yee_layout::Layout`). Hairpin could be lit (App.2.1)
because `dimension_hairpin_layout` already existed; combline needs the equivalent.
Reusing `edge_coupled_bpf` would misrepresent combline: it draws staggered *open*
half-wave lines, whereas a combline is a **comb** — aligned resonator lines short-circuited
at a common spine and capacitively loaded at the open ends.

## Decision

Add `dimension_combline_layout(project, theta0_rad, substrate) -> Result<Layout, DimError>`
to `yee-filter::dimension`, composing the comb from `yee_layout` primitives (the same
approach `dimension_stepped_impedance_layout` used — no `yee-layout` edit):

- N parallel resonator-line rectangles (`line_width_m` × `resonator_length_m`), placed
  along +x at the **solved per-section gaps** (pitch = `line_width_m + gaps_m[i]`),
  short-circuit end at a common **ground spine**, **loading-cap pads** at the open ends,
  tapped input/output feeds + ports (ref impedance = spec `z0`).
- It calls `dimension_combline` for the dims (no physics recompute).

## Consequences

**Ships:** combline gets an honest board (the comb: aligned grounded resonators + cap
pads), the prerequisite for lighting combline in the studio end-to-end (App.2.7, next).
Composes proven primitives; no new physics, no `yee-layout` edit.

**Gate (`dim_combline_layout_001`, non-vacuous geometry):** the demo order-5 combline
layout has N resonator-line traces whose dimensions match `dimension_combline`
(`line_width × resonator_length`), placed at the solved per-section gaps (consecutive
x-pitch = `line_width + gaps_m[i]`, monotone — proving it consumes the real gaps, not a
uniform placeholder), exactly 2 ports (ref `z0`), positive/finite bbox + rects.

**Honest scope:** a first-class `yee-layout::combline_bpf` generator (a later refactor),
3-D via modelling, and the SMD-cap hybrid render are out of scope; the loading caps live
in the dimensions / the studio table, not the copper. The studio lighting (App.2.7) is
the next increment. EM-verify wall (ADR-0133) untouched.

---

## References
- `crates/yee-filter/src/dimension.rs` (`dimension_combline`,
  `dimension_stepped_impedance_layout`, `dimension_hairpin_layout`);
  `crates/yee-layout/src/lib.rs` (`Polygon`/`PortRef`/`BBox`/`Layout`).
- `docs/superpowers/specs/2026-05-31-f1-2-6-combline-layout-design.md`;
  `docs/superpowers/plans/2026-05-31-f1-2-6-combline-layout.md`.
