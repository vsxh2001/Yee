# ADR-0126: Filter Phase F2.3-c — wire F2.3 onto the aperture lumped port

**Status:** Investigated — the aperture port made F2.3's elements **load the line**
(no longer inert), but the **shunt-tank capacitor reads a deepening near-short**
(a longer window makes it WORSE) so the band-pass doesn't form; `fdtd_lumped_001`
honestly RED (not weakened). Next = a CW single-frequency drive + a capacitor-arm
steady-state check. Branch `ab0f5a6` (unmerged). See Outcome.
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

## Outcome (2026-05-30) — loads the line; capacitor-arm steady-state blocks resonance

Wired onto the aperture port (merge of main was clean; ApertureSpec = trace band ×
full substrate height; one aggregate-R/L/C aperture port per branch; `n_steps`
4k→24k, probed 60k). Branch `ab0f5a6` (unmerged).

- **BEFORE (ADR-0124 single-edge sheet):** flat `|S21|≈1.0`, inert.
- **AFTER (aperture port):** strong, structured **loading** — the elements now
  genuinely couple (the O(dx²) inertness is gone). But **no band-pass**:
  `|S21|(2.0 GHz)=0.236` (12.5 dB IL, should be ≈0), `|S21|(2.4 GHz)=0.325`
  (9.8 dB rej, should be ≥20 — and *less* attenuated than the passband, wrong
  contrast); transmission rises monotonically toward the high edge.
- **Window test (24k→60k):** in-band loss got **WORSE** (12.5→22.1 dB) — a longer
  record does **not** help. For a *linear* element the DFT-of-pulse transfer
  function shouldn't degrade with window length, so this is **not** merely "pulse
  hasn't settled": the **capacitor arm reads a *deepening* near-short over time** —
  a steady-state / DC-windup behaviour in the shunt-tank cap, which dominates the
  L‖C and prevents resonance.

**Verdict (honest, gate NOT weakened):** `fdtd_lumped_001` RED. The aperture port
is a genuine step (elements load the line — the inertness of ADR-0124 is gone), but
EM-sim still doesn't ship. **Next increment:** a **CW single-frequency
steady-state drive** in F2.3 (the pulse + long-window path is demonstrably
insufficient) AND a **capacitor-arm steady-state diagnostic** — does the aperture
cap present `1/(jωC)` under a clean CW excitation, and is `V_C` bounded? This
distinguishes a *measurement* limit (pulse → CW fixes it) from a *cap-update*
windup bug (the "longer→worse" lead points at the latter). 

**Strategic note:** EM-sim is now ~10 reactive-port increments deep (6.2–6.9 +
F2.3-b/c + the de-risk investigations). Each shipped real, validated capability or
a decisive finding (stable two-way resistor port, per-axis CPML, the aperture port
that killed the O(dx²) collapse, the air-gap placement fix) — the lumped board now
loads the line — but full convergence to a validated band-pass keeps revealing the
next layer. The maintainer green-lit "the research track"; this depth warrants a
check-in on whether to keep investing (next: the CW cap drive), accept the
current full-wave-loads-the-line state, or re-scope EM-sim.

---

## References
- ADR-0125 (aperture port + the capacitor CW caveat); ADR-0124 (the air-gap fix +
  why the single-edge sheet couldn't resonate); ADR-0115 (the gate).
- `docs/superpowers/specs/2026-05-30-f2-3-c-aperture-port-driver-design.md`;
  `docs/superpowers/plans/2026-05-30-f2-3-c-aperture-port-driver.md`.
