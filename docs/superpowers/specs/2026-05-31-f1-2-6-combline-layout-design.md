# Filter F1.2.6 — Combline board layout — Design Spec

**ADR:** ADR-0145 · **Date:** 2026-05-31 · **Status:** Accepted
**Follows:** F1.2.5 (ADR-0144 — the combline synthesis engine `dimension_combline`,
shipped). This is its **board-layout companion**, the prerequisite for lighting combline
in the studio (the studio renders a `yee_layout::Layout`; combline has none yet, unlike
hairpin's `dimension_hairpin_layout`).

## Problem

`dimension_combline` returns the combline dimensions (line width, resonator length θ0/β,
loading cap, gaps) but no placeable board. Lighting combline in the studio needs a
`yee_layout::Layout`. Reusing `edge_coupled_bpf` would be **wrong**: it draws staggered
*open* half-wave lines, whereas combline is a **comb** — aligned resonator lines all
short-circuited at a common spine, capacitively loaded at the open ends. An honest
render needs the comb geometry.

## Method (compose `yee_layout` primitives, like `dimension_stepped_impedance_layout`)

`dimension_combline_layout(project, theta0_rad, substrate) -> Result<Layout, DimError>`:
call `dimension_combline` for the dims, then compose the comb from `Polygon::rect` /
`PortRef` / `BBox::from_polygons` / `Layout` (no `yee-layout` edit — same approach the
stepped-Z layout used):

- **Resonator lines:** N parallel vertical rectangles, each `line_width_m` wide (x) ×
  `resonator_length_m` long (y), short-circuit end at `y = 0`, open end at
  `y = resonator_length_m`. Place along +x: resonator `i` left edge at
  `x_i = Σ_{j<i}(line_width_m + gaps_m[j])` (centre-to-centre pitch = width + solved gap).
- **Ground spine:** a horizontal rectangle at the short-circuit end (`y` just ≤ 0,
  a `line_width_m`-tall bar) spanning all N lines' x-range — the comb spine (grounded via
  vias; vias are a fabrication annotation, not separate copper here).
- **Loading-cap pads:** a small square pad (≈`line_width_m`) at each open end
  (`y = resonator_length_m`) where the SMD loading cap `C_L` mounts (the cap value lives
  in the dimensions / the studio table, not the copper).
- **Feeds + ports:** tapped feed lines to the first and last resonator (neutral defaults,
  mirroring hairpin/edge-coupled), with a `PortRef` (ref impedance = spec `z0`) at each.
- `bbox = BBox::from_polygons(&traces)`.

## Changes

- `crates/yee-filter/src/dimension.rs` — `dimension_combline_layout`; re-export from the
  crate root. Documented; reuses `dimension_combline` (no recompute of the physics).
- `crates/yee-filter/tests/` — the `dim_combline_layout_001` geometry gate.

## DoD (machine-checkable, NON-vacuous)

**Gate `dim_combline_layout_001`** (a geometry gate, like `geo-003` / the stepped-Z layout):
build the layout for the demo combline (order-5 band-pass, θ0=45°) and assert:
1. The layout contains the **N resonator-line traces** (N = order) plus the ground spine
   + the N cap pads + the 2 feeds — and the N resonator lines have the right dimensions:
   each is `line_width_m` × `resonator_length_m` (matching `dimension_combline`), to a
   tight tolerance.
2. The resonators are laid out with the **solved per-section gaps**: consecutive
   resonator-line x-positions differ by `line_width_m + gaps_m[i]` (centre-to-centre),
   asserted against `dimension_combline`'s own `gaps_m` — proving the layout consumes the
   real solved gaps, not a uniform placeholder. NB the pitches are **symmetric** about
   the centre (not monotone): the demo Chebyshev coupling is symmetric
   (M₁₂=M₄₅, M₂₃=M₃₄ → gaps `[g0, g1, g1, g0]`), so the gate asserts mirror symmetry +
   not-all-equal (a uniform placeholder fails).
3. Exactly 2 `ports` (in/out feeds), ref impedance = spec `z0`; `bbox` width + height
   positive and finite; all trace rects have positive extent.
4. `cargo test -p yee-filter` green; `cargo clippy -p yee-filter --all-targets -- -D
   warnings` + `cargo fmt --check`; `cargo check --workspace`.

A wrong/empty/uniform-placeholder layout fails (the per-section-gap + dimension-match
assertions are the non-vacuity).

## Out of scope

The **studio lighting** of combline (the next increment, App.2.7 — consumes this layout,
mirroring App.2.1 hairpin-lit); a `yee-layout::combline_bpf` first-class generator (this
composes primitives in yee-filter, like stepped-Z — a yee-layout generator is a later
refactor); 3-D via modelling; the SMD-cap hybrid render polish. EM (ADR-0133) untouched.

## Why

Gives combline a board, the prerequisite for lighting it in the studio end-to-end —
honestly drawing the comb (aligned grounded resonators + cap pads) rather than
misrepresenting it as an edge-coupled board. Composes proven primitives; geometry-gated.
