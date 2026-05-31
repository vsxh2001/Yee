# App.2.3 — TopBar shows the active flow's verdict — Design Spec

**ADR:** ADR-0140 · **Date:** 2026-05-31 · **Status:** Accepted
**Origin:** the App.2.2 code review (P2, pre-existing) — the studio TopBar summary +
PASS/FAIL chip always read the band-pass `designed` signal, so for the **lumped** and
**low-pass stepped-Z** flows it shows a stale band-pass `f0`/FBW/verdict while the user
designs a different filter. Violates the §10 honest-UI principle.

## Problem

`TopBar(designed)` (main.rs) reads only the band-pass `Designed` — its summary
(`f0`/FBW%) and `report.pass` chip. When the active topology is `LumpedLc` (its own
ladder verdict, possibly different from the distributed verdict) or `SteppedImpedance`
(a **low-pass** design — the chip should show a *cutoff*, no FBW, and the low-pass
verdict), the TopBar is stale/wrong. Pre-existing for lumped; App.2.2 added a second
affected flow (stepped).

## Goal

The TopBar summary + PASS/FAIL chip reflect the **active flow's** real design: band-pass
(edge-coupled / hairpin) → the distributed verdict; lumped → the lumped ladder verdict;
stepped-impedance → the low-pass verdict with a cutoff (no FBW). No faked or stale state.

## Method

A pure, **host-testable** helper computes the view; `TopBar` only renders it.

```rust
/// The TopBar summary line + PASS/FAIL for the active flow. `None` verdict =
/// the flow's design is not realizable (e.g. an unrealizable lumped ladder).
pub fn topbar_view(
    topology: Topology,
    designed: &Designed,                 // band-pass (edge-coupled / hairpin)
    lumped: Option<&LumpedDesigned>,     // None when the ladder is unrealizable
    stepped: &SteppedLowpassDesigned,    // low-pass
) -> (String, Option<bool>);
```

Branches:
- `EdgeCoupled | Hairpin` → band-pass summary `· {approx} · N={order} · {f0} GHz ·
  {fbw}%`, verdict `Some(designed.report.pass)`.
- `LumpedLc` → the same band-pass summary (lumped shares the band-pass spec), verdict
  `lumped.map(|l| l.verdict.pass)` (`None` → "not realizable").
- `SteppedImpedance` → low-pass summary `· {approx} · N={order} · cutoff {f_c} GHz` (no
  FBW), verdict `Some(stepped.pass)`.

`TopBar` gains `(topology, designed, lumped, stepped)`, calls `topbar_view`, renders the
summary chip + a PASS / FAIL chip — or a muted "geometry not realizable" chip when the
verdict is `None`. The `App` call site passes the three flow signals + topology (all
already in scope).

## Changes

- `crates/yee-studio-web/src/engine.rs` — `topbar_view` (pure, documented) + a test.
- `crates/yee-studio-web/src/main.rs` — `TopBar` signature + body (call `topbar_view`) +
  the `App` call site.

## DoD (machine-checkable)

1. **Non-vacuous host test** (`cargo test -p yee-studio-web`): build a spec, then assert
   `topbar_view` is topology-aware — the `SteppedImpedance` summary contains "cutoff"
   and NOT "%" (low-pass), while the `EdgeCoupled` summary contains "%" (band-pass), and
   the verdicts come from the respective flows (`stepped.pass` vs `designed.report.pass`;
   `LumpedLc` from `lumped.verdict.pass`). A constant/band-pass-only `topbar_view` fails.
2. `dx build --platform web --release` (crates/yee-studio-web) EXIT 0.
3. Existing tests pass (band-pass + lumped + stepped flows unregressed);
   `cargo clippy -p yee-studio-web --all-targets -- -D warnings` + `cargo fmt --check`
   clean; `cargo check --workspace` green.

## Out of scope

Restyling the TopBar; other stages; new engine physics.

## Why

A real honesty fix — the top-bar verdict must reflect the filter the user is actually
designing, across all three flows. Bounded, testable, no new physics.
