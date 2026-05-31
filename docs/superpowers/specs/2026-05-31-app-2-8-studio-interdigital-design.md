# App.2.8 — Light the Interdigital technique in the studio — Design Spec

**ADR:** ADR-0150 · **Date:** 2026-05-31 · **Status:** Accepted
**Follows:** F1.2.7 (ADR-0148, `dimension_interdigital`) + F1.2.8 (ADR-0149,
`dimension_interdigital_layout`) — interdigital now has both engine + board, so lighting it is
a clean mirror of App.2.7 (combline lighting), which itself mirrored App.2.1 (hairpin).

## Problem

Interdigital is the **last** greyed studio gallery card. Its engine + layout shipped, but
`technique_status(Interdigital)` still returns `Soon(EdgeCoupled)` (a stand-in) and there is no
`Topology::Interdigital`. Light it as a live, routable technique — completing the
coupled-resonator gallery.

## Key insight (clean mirror of combline)

Interdigital is a coupled-resonator **band-pass** — same synthesis (coupling matrix / swept
response / PASS-FAIL) as edge-coupled / hairpin / combline; only the realization differs (full
λg/4 lines short-circuited at alternating ends, **no loading cap**). So everything
topology-independent in `Designed` is already correct; the new work is the geometry branch
(→ `dimension_interdigital` + `dimension_interdigital_layout`) + the routing/UI, with the
combline-distinct **loading cap** simply *absent* (`combline_loading_cap_f = None`) — the
interdigital resonator-table surfaces the λg/4 **resonator length** (already a column), not a
cap. The board renders via the existing distributed `layout_stage`. Compare + overlay
(App.2.5/2.6) **hardcode** the band-pass technique list, so interdigital must be added as a row
(as combline was in App.2.7).

## Method (`yee-studio-web`)

Add `Topology::Interdigital` and mirror every site combline occupies:

**`engine.rs`:**
1. `derive_geometry` — a `Topology::Interdigital` arm (mirror the `Combline` arm) calling
   `dimension_interdigital(project, &SUBSTRATE)` + `dimension_interdigital_layout(project,
   &SUBSTRATE)` (NO `theta0` param), `SolvedDistributed{ resonator_length_m, …,
   loading_cap_f: None }`.
2. `topology_name(Interdigital)` = `"interdigital (λ/4, alt. short)"`; `length_label` =
   `"resonator length (mm)"`.
3. The two distributed-flow match groups (`topbar_view` ~1147, `verify_view` ~1295):
   add `| Topology::Interdigital`.
4. `compare_techniques` (band-pass arm): add an interdigital row —
   `design_demo_from(spec, Topology::Interdigital)` + `from_distributed(
   RealizationTechnique::Interdigital, &interdigital)` — and update the index-based
   `compare_techniques_*` test (combline shifted indices the same way in App.2.7).
5. `overlay_curves`: interdigital shares the coupled-resonator **ideal** |S21| with
   edge-coupled / hairpin / combline (they differ only physically), so it joins the shared
   curve label — no new *distinct* curve (mirror how combline was handled; do not fabricate a
   separate ideal). Verify against the existing overlay test.

**`stages.rs`:**
6. `Topology::Interdigital` enum variant; `Stage::rail(Interdigital)` = `DISTRIBUTED`.
7. `technique_status(Interdigital)` → `Live(Topology::Interdigital)` (was `Soon`).
8. `topology_response(Interdigital)` → `Bandpass` (add to the band-pass group).
9. `technique_topology(Interdigital)` → `Topology::Interdigital`.
10. `technique_label(Interdigital)` → `"Interdigital"`.
11. The Interdigital gallery card → `selects: Some(Topology::Interdigital)`.
12. `layout_stage` distributed-render group: add `| Topology::Interdigital`.

**`main.rs`:** only if a `match Topology` site (StageCanvas) is non-exhaustive — interdigital
uses the DISTRIBUTED flow (like edge-coupled/hairpin/combline); confirm no `lumped`/`stepped`
branch mis-routes it.

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): `design_demo_from(demo_spec(),
   Topology::Interdigital)` returns a `Some(layout)` whose `layout_signature` **differs** from
   the edge-coupled, hairpin, AND combline layouts for the same spec (routes to the real
   `dimension_interdigital_layout`, not a clone), the coupling matrix / verdict are the **same**
   (shared synthesis), `combline_loading_cap_f` is **None** (no cap), and the surfaced
   resonator length > 0 (the λg/4 quantity). A stub/clone or a combline-with-cap fails.
2. `dx build --platform web --release` (crates/yee-studio-web) EXIT 0; the Interdigital card is
   live + routable (present in the built bundle).
3. Compare table includes an Interdigital row (5 band-pass rows: edge-coupled / hairpin /
   combline / interdigital / lumped); existing flows unregressed; `cargo clippy … -D warnings`
   + `cargo fmt --check`; `cargo check --workspace` green.

## Changes

- `crates/yee-studio-web/src/{engine.rs, stages.rs}` (+ `main.rs` only if a match needs the
  arm). NO `yee-filter` edits (engine + layout already exist + are gated).

## Out of scope

The interdigital via / 3-D short-circuit render (the alternating grounding is a layout/
fabrication detail surfaced in the comb geometry, not extra studio UI); precise tap/Qe→feed
(F1.2.1). EM-verify wall (ADR-0133/0147) untouched.

## Why

Lights the **last** gallery card → the coupled-resonator gallery is **complete**: edge-coupled
/ hairpin / lumped / combline / interdigital (band-pass) + stepped-impedance (low-pass), all
live + routable + in compare + overlay + export. Clean mirror of combline; no new physics (the
engine + layout shipped + are gated).
