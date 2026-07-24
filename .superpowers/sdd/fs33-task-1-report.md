# FS.3.3 Task 1 report — `import_dxf` + gates (yee-export)

Plan: `docs/superpowers/plans/2026-07-24-fs3-3-dxf-import.md` (Task 1 +
Global constraints only). Spec:
`docs/superpowers/specs/2026-07-24-fs3-3-dxf-import-design.md`.

- Branch: `feature/fs3.3-dxf-import`
- head_before: `9b91d0f06553a689f7557239f401bdc885f53dfb`
- head_after: `8d18b097f78524e8f27596e653e9d56c6f35bbdc`
- Status: DONE

## What shipped

- `crates/yee-export/src/dxf.rs` (new module): a hand-rolled ASCII DXF
  group-code-pair scanner (no new dependency) that parses closed
  `LWPOLYLINE` (straight + bulge segments) and closed R12-style
  `POLYLINE`/`VERTEX`/`SEQEND` chains into `Vec<yee_layout::Polygon>` —
  `dxf_to_outline(dxf: &str, opts: &DxfOptions) -> Result<Vec<Polygon>,
  DxfImportError>`.
- `DxfOptions { layer: Option<String> }` — optional layer filter
  (unmatched layers are silently skipped, not rejected).
- `DxfImportError` (8 variants, `Display` + `Error`, one typed
  rejection per named-rejection-matrix item, mirroring
  `GerberImportError`'s style):
  `UnsupportedUnits`, `UnsupportedEntity`, `OpenPolyline`,
  `NonzeroElevation`, `BadValue`, `UnclosedPolyline`, `BadBulge`,
  `NoOutline`.
- `crates/yee-export/src/import.rs`: `arc_vertices` bumped from
  private to `pub(crate)` (doc comment updated to note the new
  caller) — the only change to the existing Gerber module. No
  behavioural change.
- `crates/yee-export/src/lib.rs`: wires `pub mod dxf;` +
  re-exports (`DxfImportError`, `DxfOptions`, `dxf_to_outline`),
  mirroring the existing `import` wiring.
- `crates/yee-export/tests/dxf_rt.rs` (new): gate `dxf-rt-001` +
  the full rejection matrix, 6 test functions, all green.

## Key decisions (the plan explicitly left these to the implementer)

**`$INSUNITS` default policy.** Chose the strict option: `$INSUNITS`
must be present and exactly `4` (mm) or `1` (inch); a missing header
variable is a **named rejection** (`UnsupportedUnits("missing")`), not
a silently-assumed default. Rationale: DXF's own default for an absent
`$INSUNITS` is "unitless" — guessing mm there would be exactly the
kind of silent misinterpretation this importer exists to prevent, and
the strict rule is also the simpler implementation (one code path,
`Some(4) | Some(1)` else reject — no separate "missing" special case
needing its own justification). Documented in the module doc comment.

**Arc-tessellation reuse ("reuse the Gerber helper if shape-compatible;
else mirror, and say which").** Reused, not mirrored: `arc_vertices`
(the angle-stepping tessellation loop + the `ARC_CHORD_TOL_M` pinned
chord tolerance) is called directly from `dxf.rs` after bumping its
visibility to `pub(crate)`. Only the bulge → `(center, ccw)` geometry
conversion is new DXF-specific code (`bulge_vertices` in `dxf.rs`,
~15 lines) — it computes the arc's center and direction from the DXF
bulge value (`tan(included_angle/4)`, signed for CCW/CW) and hands off
to the existing loop, so DXF and Gerber arcs are tessellated by the
exact same code, not two implementations that could drift. This is
proven directly: `bulge_ccw_quarter_matches_gerber_pinned_tessellation`
reproduces `gerber_arcs_flashes.rs`'s pinned `n = 18` segment count for
an identical r = 1 mm quarter arc, vertex-for-vertex.

**Entry-point return type.** The spec's prose said `dxf_to_outline`
should mirror `gerber_to_outline`'s output type "exactly"; literally
`gerber_to_outline` returns `Vec<Point2>` (a single board-profile
path), which is a different thing from what Task 1's own bullet asks
for ("closed chains → **outline polygons**", plural, compared against
"the native generator's **polygons**" in the gate). Went with the
literal Task-1 wording: `dxf_to_outline` returns `Vec<Polygon>` (one
closed outline per LWPOLYLINE/POLYLINE entity), matching
`gerber_to_polygons`'s shape and the gate's actual comparison target
(`Layout.traces: Vec<Polygon>`). No `dxf_to_layout` counterpart to
`gerber_to_layout` was added — Task 1's DoD only requires
`dxf_to_outline`, and studio wiring is an explicit Non-goal in the
spec; adding an unrequested second entry point would be scope beyond
what was asked.

**Bulge sign convention** — verified before trusting it, not assumed.
Derived `center = M + e·n` (`n` = the chord's right-hand normal,
`e = (s² − h²)/(2s)`, `s = bulge·h`) and checked it against an
unambiguous, independently-constructible case: a CCW quarter-turn of
the unit circle from `(1,0)` to `(0,1)` must have `bulge = tan(π/8)`
and center `(0,0)` by construction (not by DXF-spec interpretation) —
the formula reproduces `(0,0)` exactly. The same formula, run with the
sign flipped, is the CW gate test. Documented in `bulge_vertices`'s
doc comment so the derivation doesn't need re-deriving later.

**Section-scoping.** `scan_entities` tracks `SECTION 2 <name>` and only
emits entities found inside `ENTITIES`, so an `LWPOLYLINE` sitting
inside a `BLOCKS`-section block definition (never inserted) can't be
mistaken for model-space geometry — needed for the `INSERT` rejection
to be meaningful rather than accidentally bypassable.

## Rejection matrix (all typed, all tested in `dxf_rt.rs`)

Open polyline; `$INSUNITS` missing; `$INSUNITS` = 0 (explicit
unitless); `$INSUNITS` = 2 (unsupported, feet); nonzero elevation on
`LWPOLYLINE` (group 38) and on `POLYLINE`/`VERTEX` (group 30);
`POLYLINE` never reaching `SEQEND`; no closed polylines in the file;
`CIRCLE`, `ARC`, `ELLIPSE`, `SPLINE`, `TEXT`, `INSERT` — each a named
`UnsupportedEntity` rejection (mirrors `GerberImportError::
UnsupportedCommand`'s one-variant-many-offenders idiom).

## Gate `dxf-rt-001`

`dxf_rt_001_vertex_exact_vs_native_stub`: builds the native S.6
stub-notch trace geometry (feed line + Hammerstad-corrected open stub)
using `yee_layout`'s own public `eps_eff`/`open_end_delta_l` helpers
(the same formula `yee-engine/tests/import_twin.rs::native_stub_layout`
uses, without duplicating the algebra — those two helpers are already
public API in a crate `yee-export` already depends on), hand-emits it
as two `LWPOLYLINE` rectangles in a minimal DXF file, imports it, and
asserts vertex-exact equality (0.5 nm tolerance) against the native
`Polygon`s. Plus two bulge fixtures (CCW/CW) reproducing the
`gerber-rt-003` pinned quarter-arc wedge (r = 1 mm, n = 18) bit-exact,
a `POLYLINE`/`VERTEX` positive-path test, and a layer-filter test.

## Verification (real output, this session)

```
$ cargo test -p yee-export --release
... dxf_rt.rs: 6 passed; 0 failed
... gerber_001_structure: 1 passed
... gerber_002_roundtrip: 1 passed
... gerber_003_outline_structure: 1 passed
... gerber_004_outline_geometry: 1 passed
... gerber_arcs_flashes: 7 passed
... gerber_roundtrip_import: 3 passed
... kicad_001_structure: 1 passed
... kicad_002_geometry: 1 passed
... doc-tests yee_export: 1 passed
(all pre-existing gerber-rt/kicad gates green and unmodified)

$ cargo clippy --workspace --all-targets -- -D warnings
Finished, 0 warnings/errors

$ cargo clippy -p yee-compute --all-targets --no-default-features -- -D warnings
Finished, 0 warnings/errors

$ cargo fmt --check --all
(clean after one `cargo fmt --all` pass fixing 2 formatting nits pre-commit)

$ cargo doc -p yee-export --no-deps
12 pre-existing warnings, all in lib.rs (unresolved-link / bare-URL
nits predating this change); zero warnings from dxf.rs.
```

`yee-engine --test import_twin` was **not** run: this task did not
touch `yee_layout::Polygon`/`Point2` (the shared outline types) or any
existing outline/Layout construction code — only added a new module
using those types as-is, plus a visibility bump on a private helper
function. The global-constraints trigger ("if you touch shared outline
types") does not apply.

## Concerns / follow-ups (none blocking)

- `dxf_to_outline`'s malformed-`$INSUNITS`-value case (e.g. a
  non-numeric group-70 value under `$INSUNITS`) collapses to the same
  `"missing"` error text as a genuinely absent variable. Low value to
  distinguish further (hand-authored fixtures won't hit it; a
  malformed numeric header is already a rejection either way) — left
  as is rather than adding a third units-error shape.
- Task 2 (ADR-0230 + roadmap row) is out of scope for this task and
  was not started.
