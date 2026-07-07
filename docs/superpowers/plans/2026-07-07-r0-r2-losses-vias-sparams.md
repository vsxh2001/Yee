# Plan — R.0–R.2 losses, vias, complex S-parameters

**Spec:** `docs/superpowers/specs/2026-07-07-r0-r2-losses-vias-sparams-design.md`

1. **R.0**: `yee_voxel::substrate_sigma_cells` (+ unit test: σ value at a known cell,
   air cells zero); gate `yee-engine/tests/board_loss.rs` (`engine-loss-001`,
   `#[ignore]`, one release solve): 6 λ_g lossy line, α from two directional triples vs
   Pozar ±20 %. CI: append to the antenna-gates job? No — new `rf-tool-gates` steps in
   `compute-engine-gates` (board-flow, not antenna). ADR-0194. Commit + push.
2. **R.1**: `yee_voxel::with_via` (+ unit test: mask cells set, bounds); gate
   `yee-engine/tests/board_via.rs` (`engine-via-001`, 3 runs): notch vanishes with the
   via, control keeps it. ADR-0195 (may share 0194). Commit + push.
3. **R.2**: `sparams::complex_s21/complex_s11` (de-embedded via fitted β) + phase gate
   + `.s2p` export of an engine measurement through `yee-io`. ADR-0196. Commit + push.
4. Roadmap statuses + footers after each; continue to R.3 (GPU parity) next.
