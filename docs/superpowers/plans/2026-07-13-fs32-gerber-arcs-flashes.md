# FS.3.2b implementation plan — Gerber import arcs + flashes

**Spec:** `docs/superpowers/specs/2026-07-13-fs32-gerber-arcs-flashes-design.md`
**Lane:** `crates/yee-export/**` (+ this spec/plan/ADR-0220).
**Base:** main @ `352eec5`.

## Steps

1. **State + error surface** (`crates/yee-export/src/import.rs`):
   - `pub const ARC_CHORD_TOL_M: f64 = 1.0e-6` with the sagitta formula
     documented.
   - New `GerberImportError` variants: `FlashInRegion`,
     `UnknownAperture(String)`, `UnsupportedAperture(String)`,
     `BadArc(String)`; extend `Display`.
   - Private `Aperture` enum: `Circle { d_m }`, `Rect { x_m, y_m }`,
     `Other(String)` (bookkept, rejected only if flashed). Parse
     `%ADD<code><template>,<params>*%` (mm → metres); anything
     unparseable stays `Other` with the raw text.
   - Parser state: `HashMap<u32, Aperture>`, `current_aperture:
     Option<u32>`, interpolation mode (`Linear`/`Cw`/`Ccw`), `g75_seen`.

2. **Word handling** (copper importer only):
   - `G75*` sets multi-quadrant; `G74*` → `UnsupportedCommand`.
   - `G01*`/`G02*`/`G03*` standalone set the mode; a `G01`/`G02`/`G03`
     prefix on a coordinate word sets the mode then parses the rest.
   - Coordinate parser gains non-modal `I`/`J` (default 0 per word).
   - `D<code>*` with code ≥ 10 selects the aperture (replacing the
     blanket "ignore numeric D words"); bare `D01*`/`D02*`/`D03*` use
     the modal point.
   - `D03`: reject in-region; look up the aperture; emit the flash
     polygon (circle tessellated / rect exact) into `polys` in file
     order.
   - `D01` with circular mode inside a region: tessellate
     (`arc_vertices`: centre = current + (I, J); sweep CW/CCW normalized
     to (0, 2π], start == end → 2π; `n = ceil(sweep/φ_max)`; interior
     vertices at uniform angle with radius linearly interpolated
     r0 → r1; final vertex pushed exactly). Outside a region a circular
     `D01` stays a stroked draw → `UnsupportedCommand`.

3. **Gate** `crates/yee-export/tests/gerber_arcs_flashes.rs`
   (`gerber-rt-003`, pattern: `gerber_roundtrip_import.rs`): the five
   numbered cases from the spec, with pinned segment counts (18 / 71 /
   50), exact-vertex asserts at ≤ 0.5 nm, on-circle + sagitta ≤ 1 µm
   checks, and the full rejection matrix.

4. **Verify:** `cargo test -p yee-export && cargo clippy -p yee-export
   --all-targets -- -D warnings && cargo fmt --check --all` → exit 0.
   Root-workspace clippy/fmt must stay green (lane touches one crate).

5. **Docs:** ADR-0220 (`docs/src/decisions/0220-fs32b-gerber-arcs-flashes.md`)
   with the measured/pinned numbers; SUMMARY.md line; FULL-SUITE-ROADMAP
   FS.3.2 cell is the dispatcher's merge-time edit (out of this lane's
   scope if contested — noted in the report).

## Escape hatch

Blocked > 20 min on one problem → reduce to circle-flash-only skeleton
(rect + arcs rejected by name) and surface the remainder.
