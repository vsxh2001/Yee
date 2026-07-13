# Plan: FS.6.2 — single-stub matching

**Spec:** `docs/superpowers/specs/2026-07-12-fs62-stub-match-design.md`

1. (a, DONE) `single_stub_match` + `StubMatch` in yee-layout; gates
   `geo_005_stub_match` (Pozar position + null contract).
2. (b, next) long-feed edge-fed patch fixture with probe triple;
   measure Γ + β at f₀; synthesize; regenerate layout with the stub;
   re-measure; gate `match-em-001` per spec. CI: antenna job step.
