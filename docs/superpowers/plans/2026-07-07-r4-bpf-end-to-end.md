# Plan — R.4 BPF end-to-end + surrogate BO

**Spec:** `docs/superpowers/specs/2026-07-07-r4-bpf-end-to-end-design.md`

1. **R.4a geometry**: `yee-layout::HairpinSectionParams` + `hairpin_bpf_sections`
   (per-section `gaps_m`; `hairpin_bpf`/geo-003 untouched) + unit tests (x-pitch
   accumulates each gap; ports at tap height).
2. **R.4a qe→tap**: `tap_offset_from_qe(qe, z0, zr, halfwave_m)` in
   `yee-filter::dimension` (formula in the spec; errors on unreachable qe / tap
   off the arm); `HairpinDimensions.tap_offset_m`; `dimension_hairpin_layout`
   switches to per-section + real tap. Unit tests vs hand-computed values +
   monotonicity (qe ↑ → t ↑ toward the fold).
3. **R.4a gate** `engine-bpf-verify-001` (2 release solves): synthesized hairpin
   vs `coupling_matrix_s_params` — measured f0 + BW(−3 dB) + rejection,
   tolerances set from the first honest run (S.8 pattern). CI: append to the
   engine-filter gates. ADR-0197 (may fold R.4b), roadmap rows, commit + push.
4. **R.4b BO**: objective closure (voxelize → submit → directional S21 → scalar
   error on f0/BW) over ≥2 knobs via `yee_surrogate::bo::minimize`
   (n_initial + n_iters ≲ 10 solves); gate `engine-bpf-bo-001` asserts the
   optimum beats the seed and lands within spec tolerance. ADR, roadmap,
   commit + push; continue to R.5 (studio flow).
