# Filter Phase F1.4.0 — `yee-export` Gerber walking skeleton — Plan

**Spec:** `2026-05-30-filter-f1-4-0-yee-export-gerber-design.md` · **ADR:** ADR-0100

## Lane
`crates/yee-export/**` (new crate) + the root `Cargo.toml` `[workspace] members`
line for it (SANCTIONED cross-lane — required for a new crate; see spec). Do NOT
edit `yee-layout` or any other crate's source — consume `yee-layout`'s public
`Layout`/`Polygon`/`Point2` API. Out-of-lane beyond the member line → finding.
Keep `yee-export` WASM-safe (pure text; deps = `yee-layout` + maybe `serde`; no
native/FDTD dep).

## Base
New worktree off current `main` (base SHA in the brief). Branch
`feature/filter-f1-4-0-yee-export-gerber`.

## Pattern files
- `crates/yee-layout/Cargo.toml` — a small WASM-safe crate's manifest shape
  (workspace-inherited fields + `[lints.rust] unsafe_code="forbid"` /
  `missing_docs="warn"`) to mirror for `yee-export/Cargo.toml`.
- `crates/yee-layout/src/lib.rs` — `Layout`, `Polygon`, `Point2` (the vertex
  accessors + that coords are metres) + `to_svg` for how it walks polygons.
- `crates/yee-io/` — house style for an I/O-emitting crate + its round-trip test
  shape (Touchstone round-trip is the precedent for "I/O gate = round-trip").

## Steps
1. Create `crates/yee-export/` (`Cargo.toml` + `src/lib.rs`); add it to root
   `Cargo.toml` `[workspace] members`.
2. `src/lib.rs`: `GerberOptions` (+ `Default`), `layout_to_gerber` per the spec's
   RS-274X structure. A small private helper for the metres→4.6-mm integer
   formatting (documented). Doc every public item.
3. `tests/gerber_001_structure.rs` + `tests/gerber_002_roundtrip.rs` per DoD 4–5.

## Verify (exit 0; nice -n 19, --jobs 2)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-export --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-export --jobs 2
```
Pure text — sub-second. Do NOT run `cargo test --workspace`, FDTD, mom-001.

## Escape hatch
Blocked > 15 min — the RS-274X format details fight (region syntax / coordinate
format ambiguity) such that you can't produce a structurally-valid file you're
confident in, OR `yee-layout`'s `Layout`/`Polygon` API doesn't expose the vertex
list needed → STOP and surface: the exact format question + a minimal sample of
what you emit, OR the missing `yee-layout` accessor (as a finding — do NOT edit
yee-layout). Do NOT fabricate a passing gate. Keep the skeleton single-layer; do
NOT scope-creep into drill/outline/KiCad.

## Done when
DoD 1–5 pass; `git diff --stat <base>..HEAD` shows only `crates/yee-export/**`,
the root `Cargo.toml` member line, `Cargo.lock`, + the 3 committed docs;
`yee-export` has no native/FDTD dep (WASM-safe).
