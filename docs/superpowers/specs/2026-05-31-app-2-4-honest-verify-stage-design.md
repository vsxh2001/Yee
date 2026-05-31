# App.2.4 — Honest Verify stage — Design Spec

**ADR:** ADR-0141 · **Date:** 2026-05-31 · **Status:** Accepted
**Origin:** the studio's last "Soon" stage. The current `verify_stage()` shows **fake
"—" placeholder stats** under a "FDTD realized response" header + a stale "App.D.5"
roadmap note — misleading (it implies an EM result that does not exist). §10 honest-UI.

## Problem

`verify_stage()` (no-arg) renders three `"—"` stat cards labelled "insertion loss @
f0" / "in-band return loss" / "rejection @ 2.4 GHz" under "FDTD realized response", plus
a note promising an FDTD overlay "Soon". It is the only studio stage that shows
fabricated-looking placeholders, and it hard-codes a band-pass-specific "2.4 GHz". Full
full-wave EM verification of the physical board is the deferred research frontier (the
ADR-0133 cavity wall) — not something the WASM studio does.

## Goal

Replace the fake placeholders with an **honest, topology-aware** Verify stage that shows
the **real** synthesized/realized metrics the engine already computes for the active
flow, states clearly **what level** was verified (realized LC ladder vs synthesized
ideal response), and honestly frames full-wave EM as a separate native step — no
fabricated numbers, no stale roadmap claim.

## Key insight (the metrics already exist)

Every flow already carries its graded metrics:
- band-pass distributed (`Designed.report: MaskReport`): `pass`, `worst_passband_ripple_db`,
  `worst_return_loss_db`, `stopband: Vec<(f, achieved, required, met)>`.
- lumped (`LumpedDesigned.verdict: MaskVerdict`): `pass`, `worst_passband_ripple_db`,
  `worst_return_loss_db`, `worst_stopband_rej_db` — graded on the **realized** ladder
  (`ladder_s21`), a genuine circuit-level verification.
- low-pass stepped (`SteppedLowpassDesigned`): `pass`, `worst_passband_ripple_db`,
  `worst_return_loss_db`, `stopband: Vec<(f, achieved, required, met)>`.

The Verify stage only needs to **surface** these honestly per flow.

## Method (mirror the App.2.3 `topbar_view` pure-helper pattern)

A pure, host-testable helper + a topology-aware stage:

```rust
/// What level of verification the metrics represent.
pub enum VerifyLevel {
    /// The realized LC ladder graded vs the mask (lumped — a real circuit-level check).
    RealizedLadder,
    /// The synthesized ideal / coupled-resonator response graded vs the mask
    /// (distributed + low-pass — the physical board's EM response is a native step).
    SynthesizedIdeal,
}
pub struct VerifyView {
    pub level: VerifyLevel,
    pub pass: Option<bool>,            // None = the active flow's design is not realizable
    pub worst_passband_ripple_db: f64,
    pub worst_return_loss_db: f64,
    pub worst_stopband_rej_db: Option<f64>,  // None when the mask has no stopband points
}
pub fn verify_view(
    topology: Topology, designed: &Designed,
    lumped: Option<&LumpedDesigned>, stepped: &SteppedLowpassDesigned,
) -> VerifyView;
```

Branches: `LumpedLc` → `RealizedLadder` from `lumped.verdict` (`None` when `lumped` is
`None`); `SteppedImpedance` → `SynthesizedIdeal` from the stepped fields (stopband rej =
min achieved over `stopband`); `EdgeCoupled | Hairpin` → `SynthesizedIdeal` from
`designed.report` (stopband rej = min achieved over `report.stopband`).

`verify_stage(topology, designed, lumped, stepped)` renders the three **real** metrics
(ripple, return loss, stopband rejection — `"—"` only when genuinely absent, e.g. no
stopband points), the PASS / FAIL / "not realizable" chip, a clear **level** label
("Realized LC ladder" vs "Synthesized ideal response"), and an honest note: the studio
verifies at the circuit / synthesis level; **full-wave EM verification of the physical
board is a separate native step** (the deferred ADR-0133 research frontier), not run in
the browser. No fabricated stats; the band-pass-specific "2.4 GHz" hard-code is gone.
The `StageCanvas` `Stage::Verify` arm passes the active flow's signals.

## Changes

- `crates/yee-studio-web/src/engine.rs` — `VerifyLevel`, `VerifyView`, `verify_view`
  (pure, documented) + a test.
- `crates/yee-studio-web/src/stages.rs` — `verify_stage` takes the four args, renders
  `verify_view` honestly.
- `crates/yee-studio-web/src/main.rs` — the `Stage::Verify` call site passes the signals.

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): `verify_view` returns the
   active flow's REAL metrics — `LumpedLc` → `RealizedLadder` + `lumped.verdict`'s
   ripple/RL/rejection (and `None` pass when `lumped` is `None`); `SteppedImpedance` →
   `SynthesizedIdeal` + the stepped metrics; `EdgeCoupled` → `SynthesizedIdeal` +
   `designed.report` metrics. Assert the metrics equal the source structs' fields (not
   `"—"`, not a constant) and the level differs (lumped vs distributed). A fake/constant
   `verify_view` fails.
2. `dx build --platform web --release` EXIT 0; the built UI no longer contains the fake
   "FDTD realized response" `"—"` placeholder cards (grep the source / built bundle for
   the removed strings).
3. Existing tests pass; `cargo clippy ... -D warnings` + `cargo fmt --check` clean;
   `cargo check --workspace` green.

## Out of scope

Any actual EM run (the deferred wall); a quantized-component (E-series) re-simulation
distinct from the existing realized-ladder verdict; tuning / auto-optimize.

## Why

Removes the studio's only fabricated-looking stat block, surfaces the real per-flow
verification metrics that already exist, and frames full-wave EM honestly — completing
the last "Soon" stage without faking or reopening the EM wall.
