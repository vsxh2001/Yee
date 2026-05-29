# ADR-0086: Filter Phase F1.0 — `yee-layout` parametric planar-filter geometry

**Status:** Accepted
**Date:** 2026-05-29
**Related:** `FILTER-DESIGN-ROADMAP.md` (Phase F1); ADR-0084 (F0 synthesis core)

---

## Context

`FILTER-DESIGN-ROADMAP.md` Phase F1 takes the planar track to a first
manufacturable filter on the **FDTD** back-end (the MoM microstrip port is
ill-posed, ADR-0064). The first sub-step, **F1.0**, is the geometry walking
skeleton: turn explicit physical dimensions into a parametric microstrip
coupled-resonator **layout** — meshable and exportable later — with a
geometry-only gate (no EM yet, keeping walking-skeleton discipline). Dimensional
synthesis (coupling extraction + surrogate-BO with FDTD) and manufacturing
export are later F1 sub-steps that consume this crate's `Layout` type.

## Decision

New pure-geometry crate **`yee-layout`** (no EM, no new external dependency):

- **Types** (all `serde`): `Substrate { eps_r, height_m, loss_tangent,
  metal_thickness_m }`, `Point2`, `Polygon { verts }` (top-metal footprint),
  `PortRef { at: Point2, width_m, ref_impedance_ohm }`, `BBox`, and
  `Layout { substrate, traces: Vec<Polygon>, ports: Vec<PortRef>, bbox }`
  (top-metal-on-substrate; ground plane implied).
- **Microstrip line synthesis** (Hammerstad-Jensen closed form):
  `microstrip_width(z0, eps_r, h) -> w` and `eps_eff(w, h, eps_r)` — used to
  size feed lines and resonators.
- **Generators**: `edge_coupled_bpf(EdgeCoupledParams) -> Layout` (N parallel
  half-wave coupled resonators + feed lines) and `hairpin_bpf(HairpinParams) ->
  Layout` (N folded resonators + tapped feed). Parameters are explicit physical
  dimensions (lengths, widths, gaps, tap positions, substrate).
- **Preview**: `to_svg(&Layout) -> String` — a dependency-free top-view SVG for
  the interactive flow and CI artifacts.

The dims→geometry direction only; the coupling-matrix→dims mapping
(consuming `yee-filter`) is the later F1.2 dimensional-synthesis step.

## Consequences

**Ships:** `yee-layout` crate + the types/generators/HJ helpers/SVG, registered
in the workspace. Gates (crate tests, §4): `geo-001` (a generated edge-coupled
2-resonator layout has the expected trace+port counts, bbox within tolerance of
hand-computed extents, non-degenerate polygons, serde round-trip) and `geo-002`
(HJ `microstrip_width(50 Ω, ε_r=4.4, h=1.6 mm) ≈ 3.0 mm` and `eps_eff ≈ 3.3`,
both ±5% vs the published Hammerstad-Jensen reference).

**Not in scope (later F1 sub-steps):** any EM / FDTD meshing; coupling-matrix→
dimensions synthesis; manufacturing export (KiCad/Gerber — F1.4); the GUI.

**No new external dependency** (`serde` only; geometry in plain f64 metres).
Lane: `crates/yee-layout/**`, root `Cargo.toml`.

---

## References
- Hammerstad & Jensen, "Accurate Models for Microstrip Computer-Aided Design,"
  1980 (width/ε_eff synthesis); Pozar §3.8.
- `FILTER-DESIGN-ROADMAP.md` Phase F1; ADR-0064 (why FDTD not MoM for planar);
  `docs/superpowers/specs/2026-05-29-filter-f1-0-yee-layout-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f1-0-yee-layout.md`.
