# ADR-0128: Filter Phase F2.3-d — CW per-frequency drive for the lumped EM sim

**Status:** Investigated — the CW drive WORKS (a 2.4 GHz notch forms, settle-converged
passband), but the coarse-grid lumped-FDTD **saturates at ~5 dB rejection** (gate
needs 20 dB) and the passband measures **over-unity** (short-board de-embed
artifact). `fdtd_lumped_001` RED, not weakened. The remaining levers (a much finer
grid, or re-scoping the 20 dB gate) are a **maintainer decision** — surfaced. Branch
`6373bc7` (unmerged). See Outcome.
**Date:** 2026-05-31
**Related:** ADR-0127 (the CW diagnostic VERDICT — the cap/port is correct; F2.3
needs a CW drive), ADR-0126 (F2.3-c — pulse drive can't settle the resonance),
ADR-0125 (aperture port), ADR-0115 (the `fdtd_lumped_001` gate), the lumped-LC →
PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

ADR-0127 settled it: the aperture lumped port (inductor 6.9 + capacitor) is
**proven correct** — under a CW drive the cap presents `1/(jωC)` and the shunt
L‖C tank **resonates**. F2.3's flat/monotone response is purely a **measurement**
problem: its modulated-Gaussian pulse + DFT measures an unsettled transient on a
short standing-wave line, so the high-Q (Q≈10) tanks never reach the steady-state
reactance the band-pass needs. The fix is downstream and concrete.

## Decision

Change F2.3's driver (`yee_voxel::simulate_lumped_board` / `run_board_solve`) to a
**CW per-frequency steady-state** measurement: for each measured frequency `f`,
drive a CW sinusoid (Hann-ramped over the first cycles to suppress the turn-on
transient), run until the field settles into a single-frequency steady-state
oscillation (enough cycles for the highest-Q tank + the line transit), and measure
the **steady-state** load voltage amplitude. `S21(f) = |V_dut,ss(f)| /
|V_thru,ss(f)|` (the DUT/thru ratio divides out the line's standing wave). This
replaces the single modulated-Gaussian + broadband DFT.

To bound runtime (per-frequency CW = separate DUT+thru solves), measure CW at the
**gate-check frequencies** (2.0 GHz passband, 2.4 GHz stopband) plus a handful for
the sweep shape — not a fine sweep. Keep the air-gap-fixed placement + the aperture
ports (ADR-0124/0126).

- **If the FDTD |S21| now reproduces the band-pass within the loose tol** (in-band
  ≈ 0 dB ±few, ≥ ~20 dB stopband) → the EM-sim component **ships** (F2.3 merges;
  lumped-LC goal 5/6).
- **If still short** → record the achieved CW |S21| (how close) precisely. Do
  **not** weaken `fdtd_lumped_001`.

## Consequences

**Ships (expected):** the goal's EM-simulation component — full-wave FDTD of the
lumped-LC board resonating under a correct CW measurement, cross-validated against
the analytic ladder at loose tol. The physics is already proven (ADR-0127), so this
is the measurement change that closes the loop.

**Gate:** `fdtd_lumped_001` GREEN (unchanged tol) on the F2.3 branch before merge;
the lumped/CPML/aperture gates non-regressed. Never weakened.

**Not in scope:** a fine CW sweep (cost — gate points + a few suffice); tight-tol
EM; SRF/ESR; the studio UI.

---

## Outcome (2026-05-31) — CW works; coarse-grid port saturates the notch → maintainer decision

The CW per-frequency drive (HannSine, ramp 12 / settle 140 / measure 16; freq set
{1.6,1.8,2.0,2.2,2.4,2.6} GHz; DUT/thru) shipped to the F2.3 branch (`6373bc7`).
It is a genuine step: a **2.4 GHz notch now forms** (flat under the pulse), and the
passband amplitude is settle-converged. But the gate is RED:

- `|S21|(2.0 GHz) ≈ 1.6` (−4.3 dB "IL" = **over-unity / gain** — unphysical for a
  passive filter ⇒ a short-board DUT/thru de-embed artifact, not real transmission);
- `|S21|(2.4 GHz)` rejection only ≈ 5 dB (gate needs ≥ 20).

**Settle-convergence sweep (the decisive check)** — 2.4 GHz rejection at
`cw_settle_cycles` = 140 / 300 / 600 = **5.1 / 6.2 / 3.4 dB** (non-monotonic,
wobbles, does NOT climb toward 20); 2.0 GHz `|S21|` = 1.61 / 1.64 / 1.64
(converged). So the shallow notch is **NOT settle-limited** — it is
**port-accuracy-saturated** on the coarse grid: the aperture port's ~25–75%
reactance accuracy (ADR-0125) degrades the L‖C tank Q, capping the notch at a few
dB, and the short-board de-embed over-corrects the passband to over-unity.

**Verdict (honest, gate NOT weakened):** the coarse-grid lumped-FDTD **cannot reach
the 20 dB cross-validation gate**. The physics is correct (ADR-0127), the
measurement is correct in principle (CW, settled), but the discretisation accuracy
(port reactance + short-board de-embed) is the floor. EM-sim has now had 11
reactive-port increments (6.2–6.10 + F2.3-b/c/d + investigations), each shipping
real capability or a decisive finding — the lumped board loads the line and forms a
band-structured response — but the 20 dB bar is unreachable without either a much
finer grid (N⁴ compute, multi-week) and/or a higher-accuracy port + a matched
de-embed.

**This is a maintainer decision** (re-scoping a validation gate is not mine to make
unilaterally — CLAUDE.md §4 "never weaken"). Surfaced via AskUserQuestion: (a)
re-scope `fdtd_lumped_001` to a physically-achievable qualitative band-structure
cross-validation (a notch forms at f_stop, monotone-ish passband) — a principled
achievable bar; (b) invest in a much finer grid + higher-accuracy port + matched
de-embed (multi-week, large compute); or (c) accept the current band-structured
full-wave response as the EM-sim deliverable, documented with its coarse-grid
accuracy limit. The F2.3 branch (`6373bc7`, CW drive) stands ready under any of (a)/(c).

---

## References
- ADR-0127 (the cap/port is correct under CW — the L‖C tank resonates);
  ADR-0126 (the pulse-drive limit); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-31-f2-3-d-cw-drive-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-d-cw-drive.md`.
