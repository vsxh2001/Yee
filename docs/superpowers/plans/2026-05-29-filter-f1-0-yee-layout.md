# Filter Phase F1.0 — `yee-layout` — Implementation Plan

**Spec:** `2026-05-29-filter-f1-0-yee-layout-design.md` · **ADR:** ADR-0086

## Lane
`crates/yee-layout/**` (new), root `Cargo.toml` (add member + workspace dep
entry if needed). Out of lane (any other crate, docs already committed) →
finding, not fix.

## Base
Worktree `worktrees/filter-f10`, branch `feature/filter-f1-0-yee-layout`,
base `53df105`.

## Pattern files
- `crates/yee-synth/` — exact shape for a new pure-math crate: `Cargo.toml`
  (`[lints.rust]` form for unsafe/missing_docs, `serde = { workspace = true }`),
  `src/lib.rs` doc style, `tests/` layout.
- Root `Cargo.toml` `members` list (append `"crates/yee-layout"`).

## Steps
1. `crates/yee-layout/Cargo.toml` — pure-math crate, dep `serde` only (workspace
   version, derive). `[lints.rust] unsafe_code="forbid"`, `missing_docs="warn"`.
2. `src/lib.rs` — the types (Substrate, Point2, Polygon, PortRef, BBox, Layout —
   all `#[derive(..., Serialize, Deserialize)]`), `microstrip_width`/`eps_eff`
   (HJ formulas from spec), `edge_coupled_bpf`/`hairpin_bpf` generators,
   `Layout::to_svg`. Doc every public item.
3. Root `Cargo.toml` — add `"crates/yee-layout"` to `members`.
4. Tests: `geo_001_edge_coupled.rs`, `geo_002_hammerstad.rs` per spec §DoD 4–5.
   Use the published HJ numbers verbatim (W≈3.0 mm, eps_eff≈3.3, ±5%).

## Verify (exit 0; nice -n 19, --jobs 2; NO --workspace)
```
nice -n 19 cargo fmt --check --all
nice -n 19 cargo clippy -p yee-layout --all-targets --jobs 2 -- -D warnings
nice -n 19 cargo test -p yee-layout --jobs 2
```
Pure-geometry crate — builds/tests in seconds.

## Math notes
- HJ `W/h` has TWO branches; compute the `<2` branch first, and if the result is
  ≥2 recompute with the `>2` branch (or pick by `B`/`A` test) — see spec. For
  FR-4 50 Ω the ratio is ~1.9 (the `<2` branch).
- `eps_eff = (εr+1)/2 + (εr−1)/2·(1+12 h/W)^(−1/2)` (no thickness correction
  needed at this fidelity).
- Signed polygon area via the shoelace formula; assert > 0 (CCW) — pick a
  consistent winding when emitting verts.

## Escape hatch
Blocked >15 min — HJ width misses 3.0 mm beyond ±5% (re-check the A/B branch
selection), or generator geometry overlaps degenerate → STOP, surface the
computed vs expected. Do NOT loosen the ±5% gate.

## Done when
DoD 1–6 pass; `git diff --stat 53df105..HEAD` shows only `crates/yee-layout/**`
+ root `Cargo.toml`/`Cargo.lock` + the 3 committed docs.
