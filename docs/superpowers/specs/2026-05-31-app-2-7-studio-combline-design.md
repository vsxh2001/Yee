# App.2.7 — Light the Combline technique in the studio — Design Spec

**ADR:** ADR-0146 · **Date:** 2026-05-31 · **Status:** Accepted
**Follows:** F1.2.5 (ADR-0144, `dimension_combline`) + F1.2.6 (ADR-0145,
`dimension_combline_layout`) — combline now has both engine + board, so lighting it is a
clean mirror of App.2.1 (hairpin lighting, which became trivial once its layout existed).

## Problem

Combline's synthesis + board shipped, but the studio gallery still greys **Combline** as
"Soon"; the App.2.0 recommender recommends it for narrow-band high-Q band-pass but can
only route the edge-coupled stand-in. Light it as a live technique.

## Key insight (clean mirror of hairpin)

Combline is a coupled-resonator **band-pass** technique — same synthesis (coupling
matrix / swept response / PASS-FAIL) as edge-coupled/hairpin; only the realization
differs (short-circuited θ0 resonators + loading caps). So everything topology-
independent in `Designed` is already correct; the only new work is the geometry branch
(→ `dimension_combline` + `dimension_combline_layout`) + the routing/UI, plus surfacing
the combline-distinct **loading cap C_L** (a single value — uniform θ0/Z0 gives the same
`C_L = cot(θ0)/(2π·f0·Z0)` for every resonator). The board renders generically via
`board_svg`. The compare table + response overlay (App.2.5/2.6) iterate the live
techniques, so combline joins both automatically.

## Method (`yee-studio-web`)

1. **`Topology::Combline`** added (stages.rs); `Stage::rail(Combline)` = `DISTRIBUTED`.
2. **Engine (engine.rs):** add a `Topology::Combline` arm to `derive_geometry`
   (mirror the `Hairpin` arm) → `dimension_combline(project, θ0, &SUBSTRATE)` +
   `dimension_combline_layout(project, θ0, &SUBSTRATE)` with **θ0 = π/4** (45° = λg/8, the
   compact default). Populate `layout`, `board_size_mm`, `line_eps_eff`, `dim_error`, and
   the resonator table (resonator length + gaps, like hairpin). Carry the **loading cap
   `C_L`** (a single f64 — surface it). `topology_name(Combline)` = "combline
   (capacitively-loaded)"; `length_label(Combline)` = "resonator length (mm)".
3. **Recommender / routing (stages.rs):** `technique_status(Combline)` →
   `Live(Topology::Combline)` (split it out of the current `Combline | Interdigital =>
   Soon`; Interdigital stays Soon); `topology_response(Combline)` = `Bandpass`;
   `technique_topology(Combline)` arm; the Combline gallery card → `selects:
   Some(Topology::Combline)`.
4. **Layout stage (stages.rs):** render the combline board (generic `board_svg`) + the
   resonator table + a **loading-cap line** ("loading cap C_L ≈ X pF per resonator, λg/8
   short-circuited resonators"). Export reuses the generic Gerber/KiCad from the `Layout`.

## Changes

- `crates/yee-studio-web/src/{stages.rs, engine.rs}` (+ `main.rs` only if a match needs
  the new arm). No yee-filter edits (the combline engine + layout already exist).

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): `design_demo_from(demo_spec(),
   Topology::Combline)` returns a `Some(layout)` whose `layout_signature` **differs** from
   both the edge-coupled and hairpin layouts for the same spec (the card routes to the
   real `dimension_combline_layout`, not a clone), the coupling matrix / verdict are the
   **same** (shared synthesis), and the surfaced loading cap `C_L > 0` and finite (the
   combline-distinct quantity is real). A stub/clone fails.
2. `dx build --platform web --release` (crates/yee-studio-web) EXIT 0; the Combline card
   is live + routable (present in the built bundle).
3. Existing band-pass / lumped / low-pass flows unregressed; `cargo clippy ... -D
   warnings` + `cargo fmt --check`; `cargo check --workspace` green.

## Out of scope

The combline SMD-cap hybrid render (the caps are surfaced as a value/table line, not
drawn as SMD footprints on the board — a polish follow-on); discrete E-series selection
of `C_L`; interdigital (a separate technique). EM-verify wall (ADR-0133) untouched.

## Why

Completes the maintainer's combline pick **end-to-end** — combline becomes a live,
routable studio technique (synthesis → board → Gerber/KiCad, in compare + overlay),
honestly surfacing the loading cap. Brings the gallery to four live band-pass techniques
(edge-coupled / hairpin / lumped / combline) + low-pass (stepped-Z). Clean mirror of
hairpin; built on the shipped combline engine + layout.
