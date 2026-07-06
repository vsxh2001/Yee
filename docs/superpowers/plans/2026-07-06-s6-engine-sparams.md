# Plan — S.6 S-parameters on the engine (walking skeleton)

**Spec:** `docs/superpowers/specs/2026-07-06-s6-engine-sparams-design.md`

1. **`yee_engine::sparams` module**: `single_bin_dft` + `transmission_db`, pure f64,
   fully documented; unit tests (known sinusoid; −6.02 dB scaled copy).
2. **Gate** `crates/yee-engine/tests/sparams_stub_notch.rs` — `engine-sparams-001`,
   `#[ignore]`'d, release-only:
   - Geometry: 3 λ_g feed line + λ/4 open stub (`L_s = λ_g/4 − ΔL_Hammerstad`) at
     mid-line, two 50 Ω `PortRef`s (drive + passive `v0 = 0` load), voxelized with the
     S.5-certified options (dx 0.3 mm, margins 34/34).
   - Two `JobSpec` runs (reference line, DUT) over `submit()`; probe E_z under the
     trace ~3 mm before the load; ~9000 steps for ring-down.
   - `transmission_db` over 3–7 GHz; assert notch within ±15 % of the closed-form
     5 GHz prediction, ≥ 8 dB deep, passband edges shallow.
3. **CI**: no new step — the S.5 `yee-engine gates` step (`--include-ignored`) picks
   the new gate up automatically.
4. **Verify**: fmt, clippy `-D warnings`, `cargo test -p yee-engine`, release gate run
   locally (~2 × ~90 s solves), record measured numbers.
5. **Ship**: ADR-0183, roadmap S.6 row + footer, SUMMARY.md, commit + push.
