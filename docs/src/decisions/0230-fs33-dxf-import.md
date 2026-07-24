# ADR-0230: FS.3.3 — DXF import (closes FS.3)

**Date:** 2026-07-24 · **Status:** accepted · **Track:** FS.3 (`FULL-SUITE-ROADMAP.md`)
**Spec:** `docs/superpowers/specs/2026-07-24-fs3-3-dxf-import-design.md`
**Plan:** `docs/superpowers/plans/2026-07-24-fs3-3-dxf-import.md`
**Predecessors:** FS.3.0/3.1/3.2 Gerber import chain (ADR-0209/0217/0220/0229) —
subset parser with named rejections, `gerber_to_outline`/`gerber_to_layout`,
byte-identical round-trip, bit-identical measurement twin (`engine-import-twin-001`).

## Context

FS.3's last open item was a second layout-import path: DXF (ASCII R12+
group-code format), the other format real board tools export alongside
Gerber. The Gerber importer's subset-plus-named-rejections discipline —
support a reviewable geometry subset exactly, reject everything else with a
typed, tested error rather than guessing — is the house style; this ADR
applies it to DXF rather than inventing a new import philosophy.

## Decision

### 1. Subset boundary

`yee_export::dxf::dxf_to_outline(dxf: &str, opts: &DxfOptions) ->
Result<Vec<Polygon>, DxfImportError>` supports:

- Closed `LWPOLYLINE` entities (straight segments + bulge arcs).
- Closed R12-style `POLYLINE`/`VERTEX`/`SEQEND` chains (the pre-LWPOLYLINE
  fallback some CAD exporters still emit).
- Entities on any layer, with an optional `DxfOptions::layer` filter
  (non-matching layers are silently skipped, not rejected — a layer filter
  is a selection, not a subset-boundary violation).

Return type mirrors the Task-1 wording literally: `Vec<Polygon>` (one closed
outline polygon per entity), matching `gerber_to_polygons`'s shape and the
gate's comparison target (`Layout.traces: Vec<Polygon>`) — not
`gerber_to_outline`'s single-path `Vec<Point2>` board-profile return, which
is a different thing (a board profile, not a set of trace outlines). No
`dxf_to_layout` counterpart to `gerber_to_layout` was added: DXF, like
Gerber, carries geometry only, and Task 1's own DoD asked for
`dxf_to_outline`, not a second entry point; studio wiring is an explicit
spec non-goal, so an unrequested `dxf_to_layout` would be scope beyond what
was asked (a caller with a `Substrate` + `Vec<PortRef>` can wrap
`dxf_to_outline`'s output exactly as `gerber_to_layout` wraps
`gerber_to_polygons`'s, when that follow-on lands).

### 2. Named rejection matrix (every entry a typed, tested error)

`DxfImportError` (8 variants, `Display` + `Error`, mirroring
`GerberImportError`'s one-variant-many-offenders idiom for entity kinds):

| Rejection | Trigger |
|---|---|
| `UnsupportedUnits(String)` | `$INSUNITS` missing, or present but not `4` (mm) / `1` (inch) |
| `UnsupportedEntity(String)` | `CIRCLE`, `ARC`, `ELLIPSE`, `SPLINE`, `TEXT`, `MTEXT`, `DIMENSION`, `INSERT`, `LINE`, `3DFACE`, `POINT`, `SOLID`, `HATCH`, or any other kind |
| `OpenPolyline` | `LWPOLYLINE`/`POLYLINE` closed flag (group 70 bit 0) unset |
| `NonzeroElevation` | nonzero `Z` (group 30) on a vertex, or nonzero constant elevation (group 38) on an `LWPOLYLINE` |
| `BadValue(String)` | a coordinate/bulge/flag failed to parse, or a required group code was missing |
| `UnclosedPolyline` | a `POLYLINE`'s `VERTEX` chain never reaches `SEQEND` |
| `BadBulge(String)` | a degenerate (zero-length-chord) bulge segment |
| `NoOutline` | the file parsed but contained zero closed polylines |

Every variant is exercised in `dxf_rt.rs`'s
`out_of_subset_inputs_are_rejected_explicitly` (or the bulge/POLYLINE
positive-path tests for `BadBulge`/`UnclosedPolyline`'s siblings), so the
subset boundary stays machine-checkable, not just documented prose.

### 3. Units decision: strict, no lenient default

`$INSUNITS` must be present and exactly `4` (mm) or `1` (inch); a missing
header variable is `UnsupportedUnits("missing")`, not a silently-assumed
default. DXF's own default for an absent `$INSUNITS` is "unitless" —
defaulting to mm there would be exactly the kind of silent
misinterpretation this importer exists to prevent. This is also the
simpler implementation: one code path (`Some(4) | Some(1)` else reject),
no separate "missing" special case needing its own justification.

### 4. Tessellation contract: reused, not mirrored

The spec asked to reuse the Gerber arc-tessellation helper "if
shape-compatible; else mirror, and say which." It was reusable:
`crate::import::arc_vertices` (bumped `private` → `pub(crate)`, its only
change) is called directly from `dxf.rs` after converting a DXF bulge
(`tan(included_angle/4)`, signed for CCW/CW) into the `(center, ccw)` shape
`arc_vertices` already expects. Only the bulge→center conversion
(`bulge_vertices`, ~15 lines) is new DXF-specific code; the pinned
[`ARC_CHORD_TOL_M`] chord tolerance and the angle-stepping loop itself are
shared verbatim with the Gerber importer, so DXF and Gerber arcs cannot
drift apart into two implementations of the same geometry. Proven directly,
not just asserted: `bulge_ccw_quarter_matches_gerber_pinned_tessellation`
reproduces `gerber_arcs_flashes.rs`'s pinned `n = 18` segment count for an
identical r = 1 mm quarter arc, vertex-for-vertex, plus the CW mirror.

The bulge sign/center formula (`center = M + e·n`, `e = (s² − h²)/(2s)`,
`s = bulge·h`) was verified before being trusted, against an
independently-constructible case: a CCW quarter-turn of the unit circle
from `(1,0)` to `(0,1)` must have `bulge = tan(π/8)` and center `(0,0)` by
construction (not by DXF-spec interpretation) — the formula reproduces
`(0,0)` exactly.

## Gate `dxf-rt-001`

`crates/yee-export/tests/dxf_rt.rs`, six tests, all instant (no FDTD):

1. `dxf_rt_001_vertex_exact_vs_native_stub` — the S.6 stub-notch trace
   geometry (feed line + Hammerstad-corrected open stub), built via
   `yee_layout`'s own public `eps_eff`/`open_end_delta_l` helpers (the same
   formula `import_twin.rs::native_stub_layout` uses, without duplicating
   the algebra), hand-emitted as two `LWPOLYLINE` rectangles in a minimal
   DXF file, imported, and asserted vertex-exact (0.5 nm tolerance — the
   `gerber-rt-001` house tolerance) against the native `Polygon`s.
2. `bulge_ccw_quarter_matches_gerber_pinned_tessellation` /
   `bulge_cw_quarter_matches_gerber_pinned_tessellation` — the
   `gerber-rt-003` pinned r = 1 mm quarter-arc wedge (n = 18 segments),
   reproduced bit-for-bit from a DXF bulge instead of a Gerber G03 arc,
   both winding directions.
3. `polyline_vertex_chain_parses_closed_rectangle` — the R12
   `POLYLINE`/`VERTEX`/`SEQEND` fallback path, vertex-exact.
4. `layer_filter_skips_non_matching_layers` — `DxfOptions::layer` drops
   the non-matching entity, keeps the matching one.
5. `out_of_subset_inputs_are_rejected_explicitly` — the full rejection
   matrix (units missing/0/2, open polyline, nonzero elevation on both
   `LWPOLYLINE` and `VERTEX`, unclosed `POLYLINE`, no-outline, and the
   `CIRCLE`/`ARC`/`ELLIPSE`/`SPLINE`/`TEXT`/`INSERT` named-entity matrix).

Geometry-only, per the spec's non-goal: the FS.3.2c twin gate
(`engine-import-twin-001`) already proved outline→measurement bit-identity
for the Gerber path, and that measurement path is shared (both importers
terminate in the same `Vec<Polygon>`/`Layout` shape), so a DXF-sourced FDTD
twin would re-run the same 293 s solve to confirm what geometry equality
already implies — explicitly out of scope, covered transitively.

## Measured result

```
$ cargo test -p yee-export --release
... dxf_rt.rs: 6 passed; 0 failed
... gerber_001_structure / 002_roundtrip / 003_outline_structure /
    004_outline_geometry / arcs_flashes (7) / roundtrip_import (3):
    all pass, unmodified
... kicad_001_structure / kicad_002_geometry: pass, unmodified
... doc-tests yee_export: 1 passed
```

All pre-existing gerber-rt/kicad gates green and byte-for-byte unmodified —
the only change to `crates/yee-export/src/import.rs` was a visibility bump
(`private` → `pub(crate)`) on `arc_vertices`, no behavioural change.

## Tolerances pinned

Vertex-exact: `0.5e-9` m (the `gerber-rt-001`/`gerber-rt-003` house
tolerance, reused, not invented). Bulge tessellation: `ARC_CHORD_TOL_M`
(the pinned 1 µm chord tolerance already gated by `gerber-rt-003`),
inherited by construction since DXF bulges tessellate through the same
`arc_vertices` call.

## Verdict — FS.3 COMPLETE

**Import side**: Gerber (FS.3.0/3.1/3.2) + DXF (FS.3.3), both terminating
in the same `Vec<Polygon>` outline shape and both closed by a named
rejection matrix instead of silent mis-parse. **Export side**: Gerber
(FS.3.1) and KiCad (`kicad_001`/`kicad_002`) — **DXF export was never a
goal** (the spec's non-goals section names it explicitly; nothing in the
FS.3 gap analysis calls for round-tripping DXF back out, only for accepting
it as an input format alongside Gerber). With both import paths gated and
the FS.3.2c twin gate having already proven outline→measurement fidelity
transitively for any `Vec<Polygon>` source, **FS.3 (`FULL-SUITE-ROADMAP.md`
§3) is COMPLETE.**

## What remains (queued, not attempted here)

**Studio DXF wiring** — the studio `ImportPanel` (ADR-0209/FS.3.1c) accepts
Gerber today; wiring a DXF file picker to `dxf_to_outline` is a small
follow-on once a `dxf_to_layout` (or equivalent panel-side wrapper
supplying `Substrate`/`Vec<PortRef>`) is written — explicitly deferred by
the spec's non-goals, not attempted in this track. No other FS.3 follow-on
is queued.
