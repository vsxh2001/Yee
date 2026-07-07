# Plan — R.2b measured Γ of a via-shorted line

**Spec:** `docs/superpowers/specs/2026-07-07-r2b-shorted-line-gamma-design.md`

1. Gate `crates/yee-engine/tests/board_short_gamma.rs` (`engine-sparams-003`,
   `#[ignore]`, one release solve): via-shorted line, complex Γ at a plane
   `d ≈ 12 mm` before the short; asserts mean |Γ| ±15 % of unity and the
   round-trip phase slope `−4π d √ε_eff/c` ±5 % (d cell-snapped).
2. Runs under the blanket `yee-engine gates` CI step automatically
   (include-ignored, non-antenna) — no new CI step.
3. ADR-0200, SUMMARY, RF-TOOL-ROADMAP R.2b row with measured numbers.
   Commit + push.
