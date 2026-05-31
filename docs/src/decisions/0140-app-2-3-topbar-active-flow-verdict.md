# ADR-0140: App.2.3 — TopBar shows the active flow's verdict

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0139 (App.2.2 low-pass flow — the review that surfaced this P2),
ADR-0120 (the lumped flow — the other affected flow), [[lumped-lc-and-studio-redesign]]

---

## Context

The studio `TopBar` reads only the band-pass `Designed` signal for its summary line
(`f0`/FBW) and PASS/FAIL chip. The studio now has three flows — band-pass distributed
(edge-coupled / hairpin), lumped LC, and low-pass stepped-impedance (ADR-0139). For the
lumped flow (its own ladder verdict) and the stepped-impedance flow (a **low-pass**
design — cutoff, no FBW, its own verdict), the TopBar shows a **stale band-pass**
summary + verdict. The App.2.2 review flagged this as a pre-existing P2; it violates the
§10 honest-UI principle (the top-bar verdict must reflect the filter the user is
actually designing).

## Decision

Make the TopBar topology-aware via a pure, host-testable helper:

- `topbar_view(topology, &Designed, Option<&LumpedDesigned>, &SteppedLowpassDesigned) ->
  (String summary, Option<bool> verdict)` — edge-coupled / hairpin → the band-pass
  summary + `designed.report.pass`; lumped → the band-pass summary + `lumped.verdict.pass`
  (`None` → not realizable); stepped-impedance → a **low-pass** summary (`cutoff f_c`, no
  FBW) + `stepped.pass`.
- `TopBar` takes the three flow signals + topology, calls `topbar_view`, and renders the
  summary + a PASS / FAIL chip (or a muted "not realizable" chip when the verdict is
  `None`).

## Consequences

**Ships:** the top-bar summary + verdict reflect the active flow — the lumped ladder
verdict, and the stepped-impedance low-pass cutoff + verdict (no more stale band-pass
state). An honesty fix across all three flows.

**Gate:** a non-vacuous host test on `topbar_view` (the `SteppedImpedance` summary shows
a cutoff and no FBW `%`, the `EdgeCoupled` summary shows FBW `%`, and verdicts come from
the respective flows — a band-pass-only/constant helper fails) + `dx build` EXIT 0 + no
regression.

**Not in scope:** TopBar restyling; other stages; new physics. The pure-helper split
keeps the rendering logic testable on host (the `#[component]` itself is not host-tested).

---

## References
- `crates/yee-studio-web/src/{main.rs (TopBar), engine.rs (Designed / LumpedDesigned /
  SteppedLowpassDesigned)}`.
- `docs/superpowers/specs/2026-05-31-app-2-3-topbar-active-flow-verdict-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-3-topbar-active-flow-verdict.md`.
