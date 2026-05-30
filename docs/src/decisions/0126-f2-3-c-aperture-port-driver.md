# ADR-0126: Filter Phase F2.3-c — wire F2.3 onto the aperture lumped port

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0125 (the aperture port, shipped — kills the O(dx²) collapse),
ADR-0124 (F2.3-b sheet placement: necessary, found the air-gap bug, but the
single-cell port couldn't resonate), ADR-0115 (the F2.3 EM-sim gate
`fdtd_lumped_001`), the lumped-LC → PCB goal, [[project-lumped-lc-and-studio-redesign]]

---

## Context

The aperture lumped port (ADR-0125, `LumpedRlcPort::aperture`) is shipped: it kills
the O(dx²) inductor collapse (dx-stable reactance) that prevented F2.3's L‖C tanks
from resonating. F2.3's driver `simulate_lumped_board` currently places each
element as a single-edge full-width *sheet* (ADR-0124, `bbc7e26`) — which loads the
line but can't resonate. Wiring it onto the aperture port is the decisive
end-to-end EM-sim test.

Known caveat (ADR-0125 Outcome): the KVL-branch **capacitor** under a single
Gaussian pulse reads near-short (sign correct, magnitude unvalidated) — it needs a
**CW / long-window steady state** to present `1/(jωC)`. F2.3's drive is a
modulated-Gaussian + DFT; for a linear system DFT-of-pulse = transfer function
*if* the record captures the full (slow) response, so a **long-enough `n_steps`**
may suffice — but if the capacitor still can't charge, a CW (single-frequency
steady-state) drive is the fallback.

## Decision

Change `yee_voxel::simulate_lumped_board` to place each ladder element via
`LumpedRlcPort::aperture(...)` over the full `(y,z)` port-face aperture (trace
width × substrate height) with one aggregate `R/L/C` per branch — replacing the
ad-hoc `C/N`,`N·L` single-edge sheet of ADR-0124. Keep the air-gap-fixed line-band
detection. Ensure `n_steps` is long enough for the capacitor's steady state
(extend the record; the band-pass needs the slow tail). Re-run `fdtd_lumped_001`
(unchanged, loose tol) in the bounded container.

- **If the FDTD |S21| now reproduces the band-pass within the loose tol** (in-band
  ≈ 0 dB ±few, ≥ ~20 dB stopband) → the EM-sim component **ships** (F2.3 merges;
  lumped-LC goal 5/6).
- **If the capacitor-under-transient limit still blocks resonance** → record the
  achieved |S21| (how close) → the next increment is a **CW single-frequency
  drive/de-embed** in F2.3 (or a CW cap validation first). Do **not** weaken
  `fdtd_lumped_001`.

## Consequences

**Ships (if it passes):** the goal's EM-simulation component — full-wave FDTD of
the lumped-LC board resonating + cross-validated against the analytic ladder at
loose tol. With F2.0/F2.1/F2.2/F2.4 → lumped engine 5/6 (only the maintainer-gated
UI merge remains).

**Gate:** `fdtd_lumped_001` GREEN (unchanged) on the F2.3 branch before merge; the
existing lumped/CPML/aperture gates non-regressed. Never weakened.

**Not in scope:** a CW drive (only if pulse + long window is insufficient); tight-tol
EM; SRF/ESR; the studio UI.

---

## References
- ADR-0125 (aperture port + the capacitor CW caveat); ADR-0124 (the air-gap fix +
  why the single-edge sheet couldn't resonate); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-30-f2-3-c-aperture-port-driver-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-3-c-aperture-port-driver.md`.
