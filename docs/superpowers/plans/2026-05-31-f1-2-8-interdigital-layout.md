# F1.2.8 — Interdigital board layout — Plan

**Spec:** `2026-05-31-f1-2-8-interdigital-layout-design.md` · **ADR:** ADR-0149

## Lane
`crates/yee-filter/src/dimension.rs` (new `dimension_interdigital_layout`),
`crates/yee-filter/src/lib.rs` (re-export it),
`crates/yee-filter/tests/dim_interdigital_layout_001.rs` (new gate). NO `yee-layout` edit
(compose from primitives, as combline did); NO studio; NO other crate. Out of lane → finding.

## Base / worktree
New worktree off `main` (re-fetch first; main is post-41d7bae). Branch
`feature/f1-2-8-interdigital-layout`.

## Pattern files (READ FIRST — edit ONLY in the worktree, never the main checkout)
- `crates/yee-filter/src/dimension.rs` — MIRROR `dimension_combline_layout` (~line 936): the
  `resonator_x` position loop (`x_i = x_{i-1} + w + gaps_m[i-1]`), `comb_right`, the
  `Polygon::rect` resonator lines, the ground bar, the tapped feeds + `PortRef`s, `bbox`,
  `Layout`. Reuse `dimension_interdigital` for the physics (no recompute). Apply the THREE
  interdigital distinctions (see spec): TWO rails (bottom y∈[−w,0] + top y∈[l+g_open,
  l+g_open+w]); resonator lines ALTERNATELY offset (even `rect(x_i,0,w,l)` grounded bottom,
  odd `rect(x_i,g_open,w,l)` grounded top); NO cap pads. `g_open = w` (neutral fixed default).
- `crates/yee-filter/tests/dim_combline_layout_001.rs` — MIRROR the gate structure; assert the
  interdigital-distinct geometry per the spec DoD (two rails, trace count = n+4 not 2n+3, no
  cap pads, even/odd y-origin alternation, no resonator touches both rails, solved+symmetric
  pitch, 2 Z0 ports).
- `crates/yee-filter/src/lib.rs` — copy the `dimension_combline_layout` re-export pattern.
- `yee_layout::{Polygon, PortRef, Point2, BBox, Layout}` — the primitives (already imported in
  dimension.rs for combline; no new yee-layout API).

## Steps
1. `dimension_interdigital_layout(project, substrate) -> Result<Layout, DimError>`: call
   `dimension_interdigital`; compute `resonator_x`/`comb_right`; push the bottom rail, top
   rail, N alternately-offset resonator lines (even y0=0, odd y0=g_open), 2 tapped feeds +
   2 `PortRef`s; `bbox`; return `Layout`. NO cap pads. Full doc comment (cite H&L §5; state the
   two-rail / no-cap / λg/4 distinctions vs combline).
2. Re-export from `lib.rs`.
3. `tests/dim_interdigital_layout_001.rs` — the gate (spec DoD parts 1–6). Make it
   non-vacuous: a combline-style single-spine/with-pads layout would FAIL parts 2/3/4.

## Verify (run FROM THE WORKTREE; expected EXIT 0; quote output)
- `cargo test -p yee-filter --test dim_interdigital_layout_001 -- --nocapture` → passes; quote
  the trace-count + two-rail assertions.
- `cargo test -p yee-filter` → no regression (combline-layout, interdigital engine, hairpin,
  edge-coupled, stepped, etc. all green).
- `cargo clippy -p yee-filter --all-targets -- -D warnings` ; `cargo fmt --check -p yee-filter`.
- `cargo check --workspace`.
- `git -C /home/hadassi/Code/Yee status --porcelain crates/` EMPTY (main untouched).

Commit (in the worktree): `yee-filter: interdigital board layout (F1.2.8, ADR-0149)` + the
Co-Authored-By trailer.

## Escape hatch
The gate MUST assert the interdigital-distinct geometry (two rails, no cap pads, alternating
offset) — NOT a self-consistency check, and NOT a combline clone (a combline-style layout must
FAIL it). Do NOT edit `yee-layout` (compose from primitives). Do NOT add cap pads. Do NOT
recompute physics (call `dimension_interdigital`). If a resonator would touch both rails
(accidental short), the offset is wrong — fix the geometry, don't loosen the gate. NEVER edit
the main checkout / another crate. Blocked > 30 min → stop + surface.

## Done when
`dimension_interdigital_layout` exists + re-exported; `dim_interdigital_layout_001` passes
(two rails, no pads, alternating offset, solved+symmetric pitch, 2 ports — non-vacuous);
existing gates unregressed; clippy/fmt/check clean; diff = `crates/yee-filter/**` only. Then I
(dispatcher) verify + adversarial code-review + merge. (Studio lighting = the final increment.)
