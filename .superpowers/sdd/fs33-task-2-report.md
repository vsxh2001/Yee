# FS.3.3 Task 2 report — ADR-0230 + roadmap row (docs)

Plan: `docs/superpowers/plans/2026-07-24-fs3-3-dxf-import.md` (Task 2 +
Global constraints only). Spec:
`docs/superpowers/specs/2026-07-24-fs3-3-dxf-import-design.md`.

- Branch: `feature/fs3.3-dxf-import`
- head_before: `f583e6f75e34f6ab5a0626bf11caf2325a1a2bd5`
- head_after: (recorded after commit below)
- Status: DONE

## What shipped

Docs-only, per Task 2's lane exception (`docs/src/decisions/0230-*.md`,
`docs/src/SUMMARY.md`, `FULL-SUITE-ROADMAP.md`). No source changes — Task 1
already landed `crates/yee-export/src/dxf.rs` +
`crates/yee-export/tests/dxf_rt.rs` (gate `dxf-rt-001`, 6 tests, all green).

- **`docs/src/decisions/0230-fs33-dxf-import.md`** — new ADR, structured to
  mirror ADR-0229 (Context / Decision subsections / Gate section / Measured
  result / Tolerances pinned / Verdict / What remains). Content pulled from
  the Task 1 report and a direct read of the shipped `dxf.rs` +
  `dxf_rt.rs` (subset boundary, the 8-variant rejection matrix as a table,
  the strict-`$INSUNITS` decision and its rationale, the tessellation-reuse
  contract with the `arc_vertices` visibility bump, the bulge sign-formula
  verification method, the `Vec<Polygon>` return-type decision vs the
  Task-1-wording-vs-literal-`gerber_to_outline`-type ambiguity Task 1 already
  resolved). States **FS.3 COMPLETE** with the import-side (Gerber + DXF)
  / export-side (Gerber + KiCad only, DXF export explicitly not a goal)
  split, and queues studio DXF wiring as the only named follow-on (not
  started, not blocking).
- **`docs/src/SUMMARY.md`** — one new line, `ADR-0230` immediately after the
  existing `ADR-0229` entry (same list format, same relative path pattern).
- **`FULL-SUITE-ROADMAP.md`** — the `FS.3` row's cell text appended with an
  **`FS.3.3 SHIPPED` (ADR-0230)** paragraph (same house style/idiom as the
  FS.3.2c paragraph immediately preceding it in the same cell) covering the
  subset, the rejection matrix, the tessellation-reuse proof, and gate
  `dxf-rt-001`'s six tests; the row's third column (validation-line-item
  checklist) gained a DXF vertex-exact bullet; the row's **Status** column
  changed from `**FS.3.0 + 3.1 + 3.2 SHIPPED** (FS.3 remainder: DXF
  import)` to `**FS.3 COMPLETE** (FS.3.0 + 3.1 + 3.2 + 3.3 SHIPPED; DXF
  export and studio DXF wiring explicitly out of scope)`. The unrelated
  gap-table row at line 153 (`Layout import (Gerber/DXF in, not just out) |
  ... | FS.3`) references the phase, not its status, and needed no edit.

No new dependency. No `crates/**` edits. No FDTD/measured content invented —
every number in the ADR (gate name, tolerance, test count, entity list)
was read directly from the shipped `dxf.rs`/`dxf_rt.rs` source, not
guessed or copied loosely from the Task 1 report's prose.

## Verification (real output, this session)

```
$ cargo fmt --check --all
(clean, no output — docs-only change, but rerun per global constraints)

$ cargo test -p yee-export --release
dxf_rt.rs: 6 passed (dxf_rt_001_vertex_exact_vs_native_stub,
  bulge_ccw_quarter_matches_gerber_pinned_tessellation,
  bulge_cw_quarter_matches_gerber_pinned_tessellation,
  polyline_vertex_chain_parses_closed_rectangle,
  layer_filter_skips_non_matching_layers,
  out_of_subset_inputs_are_rejected_explicitly)
gerber_001_structure / 002_roundtrip / 003_outline_structure /
  004_outline_geometry / arcs_flashes (7) / roundtrip_import (3):
  all pass, unmodified
kicad_001_structure / kicad_002_geometry: pass, unmodified
doc-tests yee_export: 1 passed
(all green, unmodified — expected: Task 2 touched no yee-export source)

$ cargo clippy --workspace --all-targets -- -D warnings
Finished, 0 warnings/errors

$ cargo clippy -p yee-compute --all-targets --no-default-features -- -D warnings
Finished, 0 warnings/errors
```

`yee-engine --test import_twin -- --include-ignored` was **not** run: Task 2
is docs-only and touches no shared outline types (`yee_layout::Polygon`/
`Point2`) — the global-constraints trigger for that 5-min gate does not
apply, consistent with Task 1's own reasoning for the same skip.

`missing_docs`: no Rust source touched, nothing to check beyond the
workspace clippy run above (which would surface `missing_docs` violations
as warnings under `-D warnings` if any crate's lint level caught this
change — it didn't, since nothing in `crates/**` changed).

## Concerns / follow-ups (none blocking)

- None. FS.3 is now closed per the roadmap; the only named follow-on
  (studio DXF panel wiring) is explicitly deferred by the spec's own
  non-goals section, not silently dropped.

## Fix round

Reviewer-found critical issue (1):

1. **False completeness claim.** ADR-0230:62-65 (mirrored in the
   `FULL-SUITE-ROADMAP.md` FS.3 row and this report's own line 23) asserted
   "every variant is exercised in `dxf_rt.rs`'s
   `out_of_subset_inputs_are_rejected_explicitly`". Grepping the test file
   showed only 6 of `DxfImportError`'s 8 variants were actually triggered
   by any test — `BadValue` (unparseable coordinate/bulge text, or a
   `VERTEX` missing group 10/20) and `BadBulge` (a degenerate zero-chord
   bulge) were declared but never exercised. This contradicted the plan's
   binding invariant "every rejection is a typed error with a test" for
   two reachable-on-malformed-input paths, not dead code.

   **Fix: added the two missing rejection tests** (the reviewer's
   preferred option over softening the prose, and the only option
   consistent with the plan's binding invariant) — both cheap fixtures,
   added directly to the existing `out_of_subset_inputs_are_rejected_explicitly`
   test function (test count stays 6, matching the ADR's/roadmap's "6
   tests" gate description; only the assertion count inside that one test
   grew):
   - `BadValue`: a `POLYLINE`/`VERTEX` chain where a `VERTEX` omits its
     required group `10` (X) — hits `dxf.rs`'s
     `.ok_or_else(|| DxfImportError::BadValue("VERTEX missing 10".into()))`
     path exactly (the line the reviewer cited, `dxf.rs:371`).
   - `BadBulge`: an `LWPOLYLINE` whose first two vertices coincide (a
     bulge on a zero-length chord) — hits `bulge_vertices`'s
     `d <= import::ARC_CHORD_TOL_M` guard (`dxf.rs:283`, also cited by the
     reviewer).

   Then corrected the docs that had over-claimed coverage now that the
   claim is literally true:
   - `docs/src/decisions/0230-fs33-dxf-import.md` — reworded the false
     "every variant is exercised … (or the bulge/POLYLINE positive-path
     tests for `BadBulge`/`UnclosedPolyline`'s siblings)" sentence (which
     was both inaccurate and internally confusing — `BadBulge` was never
     a "sibling" of a positive-path test) to name the two added fixtures
     directly; expanded gate item 5's rejection-matrix description to list
     `BadValue`/`BadBulge` explicitly instead of only 6 of 8.
   - `FULL-SUITE-ROADMAP.md` FS.3 row — the gate `dxf-rt-001` description
     said "the full 8-variant named-rejection matrix" while its own
     parenthetical enumerated only 6; added the missing two so the count
     and the list agree.
   - `crates/yee-export/tests/dxf_rt.rs` module doc — updated to name all
     8 variants covered, not 6.

### Verification (real output, this session)

```
$ cargo test -p yee-export --release --test dxf_rt
running 6 tests
test bulge_ccw_quarter_matches_gerber_pinned_tessellation ... ok
test dxf_rt_001_vertex_exact_vs_native_stub ... ok
test layer_filter_skips_non_matching_layers ... ok
test bulge_cw_quarter_matches_gerber_pinned_tessellation ... ok
test polyline_vertex_chain_parses_closed_rectangle ... ok
test out_of_subset_inputs_are_rejected_explicitly ... ok
test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out

$ cargo test -p yee-export --release
(full crate) all green, unmodified elsewhere: gerber_001/002/003/004,
arcs_flashes (7), roundtrip_import (3), kicad_001/002, doc-tests (1) —
same pass counts as the original Task 2 report, plus the 2 new BadValue/
BadBulge branches inside dxf_rt's existing 6-test count.

$ cargo fmt --check --all
(clean, no output)

$ cargo clippy --workspace --all-targets -- -D warnings
Finished, 0 warnings/errors

$ cargo clippy -p yee-compute --all-targets --no-default-features -- -D warnings
Finished, 0 warnings/errors
```

`yee-engine --test import_twin -- --include-ignored` not run: this fix
round touched only `crates/yee-export/tests/dxf_rt.rs` (test-only) and
docs — no `crates/**` source and no shared outline types
(`yee_layout::Polygon`/`Point2`) changed, so the global-constraints
trigger for that gate still does not apply.

All 8 `DxfImportError` variants are now machine-checked: the ADR/roadmap
claim is no longer aspirational prose, it is what the test file does.

- Status: DONE (fix round applied)
- head_before (this fix round): `fb7281e` (Task 2's original docs commit)
- head_before (task-level, per dispatch brief): `f583e6f75e34f6ab5a0626bf11caf2325a1a2bd5`
- head_after: recorded after the commit below
