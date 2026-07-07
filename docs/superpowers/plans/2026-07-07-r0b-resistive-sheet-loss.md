# Plan — R.0b resistive-sheet conductor loss

**Spec:** `docs/superpowers/specs/2026-07-07-r0b-resistive-sheet-loss-design.md`

1. yee-compute: `Materials.sheet_r_ohm: Option<f64>` (+validate); CPU
   `apply_pec_mask` applies the sheet relation to masked ex/ey edges (ez
   stays PEC); GPU rejects with `Unsupported`. Gate compute-017 (energy
   decay + R=0 bit-exact PEC).
2. yee-engine: `MaterialsSpec.sheet_r_ohm` protocol field (serde default),
   validation, pass-through. yee-voxel: `surface_resistance_ohm(f, sigma)`.
3. Gate `engine-closs-001` (one release solve): R.0 fixture, tan δ = 0,
   engineered σ = 5.8e4; measured α vs Pozar `α_c = R_s/(Z0·W)`; tolerance
   from the first honest run. Runs under the blanket engine gates CI step.
4. ADR-0202, roadmap row, commit + push. Follow-ons queued: ground-plane
   sheet, WGSL kernel, broadband SIBC.
