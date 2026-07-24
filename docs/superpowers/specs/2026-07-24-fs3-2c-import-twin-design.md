# FS.3.2c — imported-board-vs-native-twin measurement gate

**Date:** 2026-07-24 · **Track:** FS.3 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-export/**`, `crates/yee-engine/**` (+ docs)
**Predecessors:** FS.3.0 Gerber import subset (ADR-0209, `gerber-rt-001` byte-identical
round-trip), FS.3.1 `gerber_to_outline` + studio import (`gerber-rt-002` corner-exact),
FS.3.2b arcs/flashes (ADR-0220). The FS.3 row's remaining validation:
**"an imported reference board measures within tolerance of its native-built twin."**

## Point

Round-trip byte-identity proves the file layer. This gate proves the whole chain —
Gerber file → import → outline → `Layout` → voxelize → full-wave measurement —
lands on the same physics as the natively-constructed board. That is what "layout
import works" means to a user.

## Deliverables

1. **Outline→Layout twin path** (`yee-export` or `yee-engine`, wherever
   `gerber_to_outline`'s output type lives): a documented helper that rebuilds a
   `Layout` from an imported outline + user-supplied stackup/port metadata (Gerber
   carries no stackup — that asymmetry is the documented API contract; mirror how
   the studio ImportPanel already frames it). If FS.3.1 already left such a helper,
   reuse it and say so — do not duplicate.
2. **Gate `engine-import-twin-001`** (yee-engine): take the S.6 stub-notch board
   generator (the `sparams_stub_notch` fixture) as the native twin; export its
   layout to Gerber bytes (`yee_export`), import back, rebuild the Layout via (1),
   assert the rebuilt Layout is **geometrically identical** (same trace polygons —
   cheap structural assert first), then run the SAME measurement pipeline on both
   (same grid derivation from each Layout) and compare notch frequency + depth.
   Expectation: with vertex-exact import the twins should measure **identically or
   near-identically** (same voxelization inputs ⇒ same grid ⇒ bit-identical fields);
   assert bit-identical S-curves if that is what measurement shows, else pin the
   measured delta with an honest margin and explain the divergence source in the
   test docs (a nonzero delta needs a root-cause comment, not a shrug).
   Budget: one fixture, uniform or graded per whichever existing stub gate is
   cheapest (~≤ 5 min release target); `#[ignore]` + blanket CI pickup.
3. **ADR-0229** + FS.3 roadmap row (FS.3 remainder after this: DXF import only).

## Constraints

- Existing gates unmodified/green (gerber-rt-001/002/003, studio-import-e2e-001,
  sparams_stub_notch, bit-exact suite). Import subset untouched unless a genuine
  defect surfaces (that would be a finding first).
- Honest: if the twins do NOT measure bit-identically, the delta's cause must be
  identified (rasterization? float path? polygon ordering?) before pinning.

## Non-goals

DXF import; arcs/flash boards as the twin fixture (rectangles prove the chain;
arc-twins can extend later); studio wiring (the studio already has import + echo).
