# ADR-0062 — fdtd-201 rectangular-cavity resonance gate (breadth rotation)

**Status:** Accepted
**Date:** 2026-05-24
**Context Phase:** 2.fdtd validation milestone `fdtd-201`

## Context

The recent autonomous cycles have run deep in one subsystem — the MoM
cross-section eigensolver / microstrip wave-port chain (Phase 1.3.1.1
steps 4→5.8, then 1.3.1.2 quasi-TEM, then 1.3.1.2-B numerical port,
ADRs 0050–0061). That chain is at a clean stopping point: the mom-002
numerical port shipped and localized its residual to a coordinate-frame
mapping. The two in-theme continuations (the frame-mapping port; a
mom-003 edge-feed port) are both *more* microstrip-port depth and both
grind-risky — which, per the project's load-bearing lesson, means
bounded-experiment-only, not a default.

A read-only rotation survey (FDTD / GUI / surrogate, verifying actual
source vs the repeatedly-stale ROADMAP/CLAUDE.md wording) found: the
surrogate layer is already shipped end-to-end (GP/BO/AL/NSGA-II + Python
bindings; the only next item is large NN-backend work needing a new
toolchain); the GUI plotting gaps are real but lack a published-benchmark
contract (rendering features, golden-image tests only); and FDTD's
headline open items are the *deferred quagmires* (Q6 energy-balance at
75–79% drift; fdtd-007 wrong-reference). But FDTD also has a clean,
un-shipped, high-value milestone hiding under those: `fdtd-201`.

## Decision

Rotate to **`fdtd-201`** — the rectangular-cavity resonance validation
gate listed as the first Phase-2 FDTD milestone but never implemented
(no cavity-resonance test exists). Extract a PEC cavity's dominant TE₁₀₁
resonance from a time-domain FDTD run (off-centre Gaussian pulse →
interior probe → single-bin-DFT peak find, **no new dependency**) and
match the analytic Pozar §6.3 closed form. Ship it **loose-tolerance-first**
(≈±2–3%, the strict ±0.5% deferred because it is grid-dispersion-limited
on a coarse mesh), `#[ignore]`-gated for wall-time. Pure consumer of the
existing `yee-fdtd` public API: **tests + README only, no `src/` change.**

## Rationale

(1) **Genuine breadth rotation** — volumetric FDTD, a different subsystem
from the MoM-port work, reducing tunneling risk.

(2) **High value × dispatchability.** It satisfies CLAUDE.md §4 (a
published-benchmark validation case — the literal `fdtd-201` ROADMAP
milestone against the analytic reference the FEM side already proved at
fem-eig-001), closes the most conspicuous FDTD-validation gap, and
exercises real solver physics (eigenfrequency / grid dispersion) no
current test covers. Yet it is a single-crate, single-test, no-`src/`
pass — it *cannot regress the solver* and needs no new toolchain. It
beats the GUI increment (no benchmark contract) and the surrogate
increment (large, new-NN-toolchain) on value × dispatchability.

(3) **Sidesteps the known FDTD quagmires** (Q6 energy-balance, fdtd-007)
by depending on none of their deferred infrastructure.

(4) **Grind-risk bounded** by the loose-tol-first policy + the "no `src/`
edits; if the API can't extract a resonance, stop + surface" escape hatch
— consistent with how this repo handles every tolerance-sensitive gate.

## Consequences

* A new `#[ignore]`-gated `crates/yee-fdtd/tests/cavity_resonance.rs`
  matching analytic TE₁₀₁ within a documented loose tolerance, + the
  `fdtd-201` README row flipped to live.
* The strict ±0.5% (and the Q-factor extraction the README row also
  names) become documented follow-ons, not part of the first slice.
* No solver `src/` change → zero regression surface; the rest of the FDTD
  suite is untouched.
* If the public API cannot probe a clean resonance, the track stops with
  a surfaced API-gap finding (no solver code added under this ADR).

## References

* Pozar §6.3; `crates/yee-fdtd/validation/README.md:13` (the contract).
* fdtd-201 spec + plan (2026-05-24).
* The rotation survey + the bounded-experiment lesson (memory
  `step5-mixed-solver-dielectric-underweight`). ADR-0061 (the MoM-port
  track this rotates off).
