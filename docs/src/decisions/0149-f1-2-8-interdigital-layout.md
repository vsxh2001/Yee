# ADR-0149: Filter F1.2.8 — interdigital board layout

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0148 (F1.2.7 `dimension_interdigital` engine — the physics this places),
ADR-0145 (F1.2.6 `dimension_combline_layout` — the pattern this mirrors), ADR-0146 (App.2.7
combline lighting — what the studio increment after this will mirror), [[lumped-lc-and-studio-redesign]]

---

## Context

`dimension_interdigital` (F1.2.7) ships the interdigital physics but no `Layout`; the studio
renders a `Layout`, so lighting interdigital needs a board generator. This is the exact analog
of `dimension_combline_layout` (F1.2.6).

## Decision

Ship `dimension_interdigital_layout` in `yee-filter`, composing the comb from `yee_layout`
primitives (no `yee-layout` edit), calling `dimension_interdigital` for the physics. The
interdigital comb differs from the combline comb in three concrete, drawable ways:

1. **Alternating-end grounding** — two ground rails (bottom + top); resonator `i` grounds to
   the bottom rail if even, the top rail if odd, with its open end gapped from the opposite
   rail (the resonator lines are alternately offset by `g_open`). Combline shorts all
   resonators at one common spine.
2. **No loading-cap pads** — the full λg/4 line resonates without a cap, so there are no SMD
   pads at the open ends (combline draws one per resonator).
3. **Full λg/4 lines** — `resonator_length_m` from the engine (combline's were λg/8).

Resonator x-positions are the combline positions (solved per-section gaps); `g_open = w` is a
neutral fixed default (precise end-coupling is an EM follow-on, like the hairpin fold spacing).
Tapped feeds + two `Z0` `PortRef`s mirror combline. Traces = N lines + 2 rails + 2 feeds (no
pads).

## Consequences

**Ships:** the interdigital board, gated by `dim_interdigital_layout_001` (non-vacuous): N
resonator lines matching the engine dims; **two** ground rails (the alternating-ground
structure); trace count `n+4` (no cap pads, vs combline's `2n+3`); even/odd `y`-origin
alternation with no resonator touching both rails (no accidental short); solved + symmetric
per-section pitch; two `Z0` ports. A combline-style single-spine/with-pads layout FAILS the
gate — so it is genuinely interdigital-specific, not a clone or a self-consistency tautology.
Unblocks the final increment: lighting interdigital in the studio (App.2.x).

**Gate honesty:** the geometry assertions distinguish interdigital from combline (two rails,
no pads, alternating offset); the engine-dim equality (parts 1, 5) pins the placement to the
real `dimension_interdigital` output (no recompute).

**Not in scope:** studio lighting (next, final). Precise `g_open` / tap / Qe→feed (F1.2.1).
Via / 3-D short modelling (fabrication annotation). The alternating-ground even/odd coupling
refinement (deferred EM follow-on, as combline's cap-coupling interaction).

---

## References
- `crates/yee-filter/src/dimension.rs` (`dimension_interdigital_layout`, mirroring
  `dimension_combline_layout`); `crates/yee-filter/tests/dim_interdigital_layout_001.rs`;
  `yee_layout::{Polygon, PortRef, Point2, BBox, Layout}`.
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, §5 (interdigital
  λg/4 short-circuited-at-alternating-ends resonators).
- `docs/superpowers/specs/2026-05-31-f1-2-8-interdigital-layout-design.md`;
  `docs/superpowers/plans/2026-05-31-f1-2-8-interdigital-layout.md`.
