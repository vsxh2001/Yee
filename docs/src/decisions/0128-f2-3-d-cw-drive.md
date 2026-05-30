# ADR-0128: Filter Phase F2.3-d — CW per-frequency drive for the lumped EM sim

**Status:** Accepted
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

## References
- ADR-0127 (the cap/port is correct under CW — the L‖C tank resonates);
  ADR-0126 (the pulse-drive limit); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-31-f2-3-d-cw-drive-design.md`;
  `docs/superpowers/plans/2026-05-31-f2-3-d-cw-drive.md`.
