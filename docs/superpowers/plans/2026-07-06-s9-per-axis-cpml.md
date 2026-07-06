# Plan — S.9 per-axis CPML on the job protocol

**Spec:** `docs/superpowers/specs/2026-07-06-s9-per-axis-cpml-design.md`

1. `yee-engine`: `BoundarySpec::Cpml { npml, #[serde(default…)] axes: [bool; 3] }`;
   `run_job` applies `.with_axes(axes)`. Fast serde tests (legacy JSON, round-trip).
   Update the one existing `Cpml` construction site (engine doctest none; grep).
2. Experiment: LPF gate with `axes: [true, true, false]`, release run; record the
   table. Decide PEC box vs CPML-xy on the measured ripple/absolute levels.
3. Adopt the winner in `engine_lpf_verify.rs` with justified asserts; consider the
   stub gate too if the improvement is large (re-run it if touched).
4. ADR-0186 (root cause + measurements + decision), roadmap S.9 row + footer,
   SUMMARY.md, ADR-0185 stays as-is (historical record).
5. fmt/clippy/fast tests, release gates, commit + push.
