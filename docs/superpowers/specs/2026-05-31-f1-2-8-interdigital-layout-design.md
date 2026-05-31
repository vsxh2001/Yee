# F1.2.8 — Interdigital board layout — Design Spec

**ADR:** ADR-0149 · **Date:** 2026-05-31 · **Status:** Accepted
**Follows:** F1.2.7 (ADR-0148, `dimension_interdigital` engine) + F1.2.6 (ADR-0145,
`dimension_combline_layout` — the pattern this mirrors). The board companion of the
interdigital engine + the prerequisite for studio lighting.

## Problem

`dimension_interdigital` ships the physics (λg/4 line widths, lengths, solved gaps) but no
`Layout`. The studio renders a `Layout`, so lighting interdigital needs a board generator —
`dimension_interdigital_layout`, the exact analog of `dimension_combline_layout`.

## Key insight: interdigital comb ≠ combline comb

Both are aligned coupled-line combs at the solved per-section gaps, but the interdigital
realization differs in three concrete, drawable ways (combline → interdigital):

1. **Alternating-end grounding (the "finger" structure).** Combline shorts *all* resonators
   at one common ground spine (y = 0). Interdigital shorts adjacent resonators at
   **alternating** ends — so there are **two** ground rails (bottom + top), and resonator `i`
   connects to the bottom rail if `i` is even, the top rail if `i` is odd; its opposite
   (open) end is gapped from the other rail.
2. **No loading-cap pads.** Combline draws a `w×w` cap pad at each open end (the SMD `C_L`
   mounts there). Interdigital has **no cap** (the full λg/4 line resonates on its own) →
   **no pads**.
3. **Full λg/4 lines** (`resonator_length_m` from the engine; combline's were λg/8).

## Method (`yee-filter`, `crates/yee-filter/src/dimension.rs`)

`pub fn dimension_interdigital_layout(project: &FilterProject, substrate: &Substrate) ->
Result<Layout, DimError>` — calls `dimension_interdigital` for the physics (no recompute),
composes the comb from `yee_layout` primitives directly (`Polygon::rect` / `PortRef` /
`BBox::from_polygons` / `Layout`), exactly as `dimension_combline_layout` does (no
`yee-layout` edit; there is no first-class `interdigital_bpf` generator). NO `theta0`
parameter (interdigital is θ = π/2 fixed).

Let `n = gaps_m.len() + 1`, `w = line_width_m`, `l = resonator_length_m`, and the open-end
coupling gap `g_open = w` (a neutral fixed default — like the hairpin fold spacing — since
mapping it to a precise end-coupling is an EM follow-on, not first-order). Resonator left
edges `x_i` are the combline positions: `x_0 = 0`, `x_i = x_{i-1} + w + gaps_m[i-1]` (pitch =
the **solved** per-section gap). `comb_right = x_{n-1} + w`.

**Traces** (count = `n + 2 + 2` = N lines + 2 rails + 2 feeds — note **no** cap pads, vs
combline's `2n+3`):

- **Bottom ground rail** — `Polygon::rect(0, −w, comb_right, w)` (y ∈ [−w, 0]), spanning the
  comb x-range. Grounds the **even** resonators (vias = a fabrication annotation, not copper,
  as in combline).
- **Top ground rail** — `Polygon::rect(0, l + g_open, comb_right, w)` (y ∈ [l+g_open,
  l+g_open+w]). Grounds the **odd** resonators.
- **N resonator lines**, alternately offset so the grounded end touches its rail and the open
  end is gapped `g_open` from the opposite rail:
  - **even `i`** (grounded bottom): `Polygon::rect(x_i, 0, w, l)` — shares the y = 0 edge with
    the bottom rail; open top at y = l, gapped `g_open` below the top rail (at y = l+g_open).
  - **odd `i`** (grounded top): `Polygon::rect(x_i, g_open, w, l)` — top at y = l+g_open shares
    the top rail's edge; open bottom at y = g_open, gapped `g_open` above the bottom rail (at
    y = 0).
- **Tapped feeds + ports** — neutral defaults mirroring combline/hairpin: feed width = `w`,
  feed length = `l`, tapped up the *first* (i=0, grounded bottom) and *last* (i=n−1)
  resonators from their grounded ends at a neutral `tap_y`; each ends in a `PortRef` at the
  spec `Z0`. (Tap height = a neutral fraction; qe→tap dimensioning is deferred, F1.2.1.)

`bbox = BBox::from_polygons(&traces)`. Propagates every `DimError` from
`dimension_interdigital`.

## DoD — gate `dim_interdigital_layout_001` (machine-checkable, non-vacuous)

Mirrors `dim_combline_layout_001`, asserting the interdigital-distinct geometry:

1. **N resonator lines, dims from the engine.** Exactly `n = gaps_m.len()+1` resonator-line
   rects, each `w` wide × `l` long, with `w == dims.line_width_m` and `l ==
   dims.resonator_length_m` (= the λg/4 engine value). No recompute drift.
2. **Two ground rails (the alternating-ground structure — interdigital-DISTINCT).** Exactly
   two horizontal rail bars spanning `[0, comb_right]` in x: one below (y < 0) and one above
   (y > l). Combline has exactly one. Assert both present + their span.
3. **No loading-cap pads (interdigital-DISTINCT).** Total trace count `== n + 4` (N lines + 2
   rails + 2 feeds), NOT combline's `2n+3` — i.e. no `w×w` pad at any open end. Assert the
   count and that no extra `w×w` square sits at an open end.
4. **Alternating connectivity.** Even-index resonators share the y = 0 edge (touch the bottom
   rail) and are open at the top (gapped `g_open` below the top rail); odd-index resonators are
   offset to share the top rail and open at the bottom. Assert the even/odd `y`-origin pattern
   (even y0 = 0, odd y0 = g_open) and that no resonator touches *both* rails (no accidental
   short → cavity).
5. **Solved per-section pitch + symmetry.** Centre-to-centre pitch `i→i+1` = `w + gaps_m[i]`
   (the real solved gaps); symmetric for the symmetric Chebyshev coupling (`gaps_m` palindrome,
   as combline's gate asserts).
6. **Two `Z0`-referenced ports.** `ports.len() == 2`, each `ref_impedance_ohm == spec.z0_ohm`.

## Changes

- `crates/yee-filter/src/dimension.rs` (`dimension_interdigital_layout`) +
  `crates/yee-filter/src/lib.rs` (re-export) +
  `crates/yee-filter/tests/dim_interdigital_layout_001.rs` (new gate). NO `yee-layout` edit;
  NO studio. Compose from `yee_layout` primitives like combline.

## Out of scope

Studio lighting (the next + final increment, App.2.x — mirror App.2.7). Precise open-end
coupling-gap (`g_open`) dimensioning + tap/Qe→feed (F1.2.1; EM follow-on). Via / 3-D
short-circuit modelling (a fabrication annotation, as in combline). The grounded-alternating
even/odd-mode coupling refinement (deferred, like combline's cap-coupling interaction).

## Why

Gives interdigital a board, the prerequisite for lighting the last greyed studio gallery card.
Clean mirror of the shipped + gated combline layout; the only new geometry is the two-rail
alternating-ground comb with no cap pads — gated by a non-vacuous geometry check that
distinguishes it from the combline comb (two rails, no pads, λg/4 lines).
