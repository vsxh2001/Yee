# Plan — R.6 hairpin seed fidelity

**Spec:** `docs/superpowers/specs/2026-07-07-r6-hairpin-seed-fidelity-design.md`

1. `yee-filter::dimension`: `HairpinOptions` (fold_widths 2.0, resonator_z_ohm
   None, corner_widths 0.85) + `dimension_hairpin_opts`; `_with_fold`/`plain`
   delegate. `HairpinDimensions.feed_width_m`. Unit tests: corrected arm
   formula; Zr = 70 on h = 1.6 mm dimensions OK; defaults consistency.
2. Consumers: `dimension_hairpin_layout`, studio `design.rs`/`verify.rs`, BO
   gate `candidate_layout` use `feed_width_m`; `hairpin_dim_001` arm assert
   evolves (+κ·w term); studio `design_e2e` unrealizable case re-triggered via
   `fold_widths = 3.5` (corner correction makes the old case realizable).
3. Full-wave: re-run `engine-bpf-bo-001` (13 solves) — the seed re-measures
   with corrected arms; adapt asserts to the honest outcome. ADR-0201,
   roadmap R.6 row, commit + push.
4. Then R.4c prep: env-parameterized BO gate (dx / backend) + a gpu-nightly
   job wired but gated on `YEE_GPU_RUNNER_ENABLED`.
