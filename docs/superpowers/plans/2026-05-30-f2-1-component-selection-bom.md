# Filter Phase F2.1 — component selection + BOM — Plan

**Spec:** `2026-05-30-f2-1-component-selection-bom-design.md` · **ADR:** ADR-0112

## Lane
`crates/yee-filter/**` ONLY (new `src/parts.rs`, `lib.rs` re-export, `tests/`).
Consume F2.0's `LumpedLadder`/`LcResonator`. Out of lane → finding. WASM-safe:
pure `f64`/data + serde, no native dep.

## Base
New worktree off current `main` (re-fetch first). Branch
`feature/filter-f2-1-bom`.

## Pattern files (MIRROR)
- `crates/yee-filter/src/lumped.rs` (just-shipped F2.0) — module-doc style, serde
  structs, `lib.rs` re-export, `LumpedLadder`/`LcResonator`/`LcBranch` shape.
- `crates/yee-filter/tests/lumped_001.rs` — gate idiom + the cheb_bpf fixture
  (f0=2e9, fbw=0.10, z0=50, cheb 0.5 dB, N=5) → `synthesize` → `synthesize_lumped`.
  Clone the setup for `bom_001`.

## Steps
1. `src/parts.rs`: `ESeries{E24,E96}` (the 24/96 standard mantissas as consts;
   `values_decade`, `nearest` log-nearest across decades via
   `10^(round? no — search the tiled values)`, `tolerance_pct`), `CompKind`,
   `BomLine`, `Bom` (+serde, `total_parts`), `select_components(&LumpedLadder,
   ESeries)->Bom` (L + C line per resonator, nearest-selected, duplicates grouped
   into qty). Document all public items.
2. `lib.rs`: re-export `ESeries, CompKind, BomLine, Bom, select_components`.
3. `tests/bom_001.rs` per DoD 3 (E-series textbook anchors incl log-nearest tie
   cases; selection within bound for the cheb ladder; total_parts==2N; qty
   grouping for the symmetric ladder).

## Verify (exit 0; pure-math, fast — NO FDTD)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-filter --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-filter --jobs 2
```
Do NOT run workspace/FDTD/mom-001.

## Escape hatch
Blocked > 20 min (E-series tiling/rounding edge cases fight the tests; the
log-nearest tie-break is ambiguous on a boundary value) → STOP + surface the
exact value + the two candidate E-series members. Do NOT touch
`synthesize_lumped`/`dimension_*`; do NOT add deps.

## Done when
DoD 1–3 pass; `bom_001` green; `git diff --stat <base>..HEAD` = only
`crates/yee-filter/**`; F2.0 `synthesize_lumped` + the distributed `dimension_*`
untouched; WASM-safe preserved.
