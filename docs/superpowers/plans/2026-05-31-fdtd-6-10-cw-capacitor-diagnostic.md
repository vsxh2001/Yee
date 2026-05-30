# Phase 2.fdtd.6.10 — CW capacitor steady-state diagnostic — Plan

**Spec:** `2026-05-31-fdtd-6-10-cw-capacitor-diagnostic-design.md` · **ADR:** ADR-0127

## Lane
`crates/yee-fdtd/**` ONLY (`tests/cap_cw_001.rs`). No `src` change unless the
verdict is a cap-bug. Out of lane → finding.

## Base
Worktree off `main`, branch `feature/fdtd-6-10-cap-cw`. (Shipped: merge `f053164`.)

## Steps (as executed)
1. Built `cap_cw_001.rs`: CW sinusoid @ 2 GHz, ~200 cycles, Hann ramp; Probe 1
   (isolated arm, asserted) + Probe 2 (on-guide, recorded). Sliding-window
   `Z = V_T/I` + `V_C` envelope.
2. Ran in the bounded container; read the steady-state cap `Z` + `V_C`.
3. Verdict MEASUREMENT-LIMIT (cap correct) → no `src` change; recorded in ADR-0127.

## Verify (as run, container)
- fmt + clippy -D warnings → exit 0.
- `cap_cw_001 --ignored` → ok; no-regression (`aperture_port_001`,
  `lumped_rlc_twoway_001`) → ok.

## Done when
Verdict recorded (ADR-0127), gate green, no regression. The verdict drives F2.3-d
(the CW per-frequency drive).
