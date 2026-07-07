# Plan — A.1–A.3 antenna path

**Spec:** `docs/superpowers/specs/2026-07-06-a1-a3-antenna-path-design.md`

1. **A.1**: `yee-layout::inset_fed_patch` + unit tests (inset depth vs hand-computed
   R_edge); `sparams::directional_reflection_db` + unit test; gate
   `engine-antenna-002` (one release solve); ADR-0191; commit + push.
2. **A.2**: `CpmlConfig::with_faces` (CPU + GPU honor faces; fast bit-exact/behavior
   tests); `NtffSpec`/far-field on the protocol (CPU host-adapter path, mirrors
   compute-010); gate `engine-antenna-003` (patch upper-hemisphere pattern);
   ADR-0192; commit + push.
3. **A.3**: gate `engine-antenna-004` (secant on synthesis f₀ from the no-ΔL crude
   seed, directional-S11 dip observable); ADR-0193; roadmap rows + footer;
   commit + push.

Each phase ships behind its gate; measured numbers recorded in the ADRs.
