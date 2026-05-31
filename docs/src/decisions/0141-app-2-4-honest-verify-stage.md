# ADR-0141: App.2.4 — Honest Verify stage

**Status:** Accepted
**Date:** 2026-05-31
**Related:** ADR-0140 (App.2.3 — the `topbar_view` pure-helper pattern this mirrors),
ADR-0133 (the full-board EM-verify cavity wall — the deferred frontier this stage is
honest about), ADR-0120 / ADR-0139 (the lumped + low-pass flows whose metrics it
surfaces), [[lumped-lc-and-studio-redesign]]

---

## Context

The studio's `verify_stage()` was the last "Soon" placeholder — and the only stage
showing **fabricated-looking** content: three `"—"` stat cards under a "FDTD realized
response" header (with a hard-coded "rejection @ 2.4 GHz") plus a stale "App.D.5" roadmap
note. Full-wave EM verification of the physical board is the deferred research frontier
(the ADR-0133 cavity wall) — not something the WASM studio runs. Meanwhile every flow
already computes its graded verification metrics (`MaskReport` / `MaskVerdict` / the
stepped low-pass fields), including the lumped flow's **realized-ladder** verdict (a
genuine circuit-level check via `ladder_s21`).

## Decision

Replace the fake placeholders with an honest, topology-aware Verify stage, mirroring the
App.2.3 pure-helper pattern:

- `verify_view(topology, &Designed, Option<&LumpedDesigned>, &SteppedLowpassDesigned) ->
  VerifyView` — pulls the active flow's real metrics (worst in-band ripple, worst return
  loss, worst stopband rejection, pass) and tags a `VerifyLevel`: `RealizedLadder`
  (lumped — graded on the realized ladder) vs `SynthesizedIdeal` (distributed + low-pass
  — the ideal/synthesized response; the board's full-wave EM response is a native step).
- `verify_stage` renders those real metrics + a PASS / FAIL / "not realizable" chip + the
  level label + an honest note: the studio verifies at circuit/synthesis level; full-wave
  EM verification of the physical board is a separate native step (the deferred ADR-0133
  frontier), not run in the browser.

## Consequences

**Ships:** the studio's last placeholder stage becomes honest + useful — real per-flow
verification metrics (no fabricated "—", no hard-coded "2.4 GHz", no stale roadmap
claim), a clear statement of *what* was verified (realized ladder vs synthesized ideal),
and an honest framing of full-wave EM as a native step. Completes the six-stage flow's
honesty.

**Gate:** a non-vacuous host test on `verify_view` (per-flow metrics equal the source
structs' fields; level differs lumped-vs-distributed; the lumped `None`/unrealizable
case) + `dx build` EXIT 0 (the fake "FDTD realized response" `"—"` cards removed) + no
regression.

**Not in scope:** any actual EM run (the deferred wall is untouched); a quantized
E-series re-simulation distinct from the existing realized-ladder verdict; tuning /
auto-optimize. The stage is honest that full-wave EM is elsewhere.

---

## References
- `crates/yee-studio-web/src/{engine.rs (topbar_view pattern, the flow metric structs),
  stages.rs (verify_stage), main.rs (StageCanvas)}`.
- `docs/superpowers/specs/2026-05-31-app-2-4-honest-verify-stage-design.md`;
  `docs/superpowers/plans/2026-05-31-app-2-4-honest-verify-stage.md`.
