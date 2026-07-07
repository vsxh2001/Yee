# R.2b — measured complex Γ of a known one-port (via-shorted line)

**Date:** 2026-07-07
**Track:** RF-TOOL-ROADMAP R.2b
**Related:** ADR-0195 (R.2 — whose THRU-Γ negative result motivated this),
ADR-0194 (R.1 vias — the short), ADR-0189 (the 3-probe split).

## Problem

R.2's negative result: on a THRU line, the plane-A "Γ" is the far port's
small residual reflection plus window-truncated multi-bounce — raw |Γ|
reached 2.4 and the exported S11 had to fall back to the exact-by-definition
zero. The engine has never measured a **known** complex reflection against
theory.

## Approach — measure the cleanest one-port there is

A via-shorted microstrip line: Γ = −1 exactly at the via plane; at a
reference plane `d` before it, `Γ(f) = −e^{−2jβ(f)d}` with β from
Hammerstad–Jensen. Unlike the THRU case the backward wave is as large as the
forward one, so the 3-probe split works with full signal on both arms — this
isolates "can the engine measure a complex Γ" from the separate THRU
de-embedding problem (which stays open; see below).

One release solve (`crates/yee-engine/tests/board_short_gamma.rs`): the
R.0/R.1 stack (FR-4 1.6 mm, dx = 0.3 mm), drive aperture port, probe triple
one λ_g from the launch, via `d ≈ 12 mm` (~0.73 λ_g at 5 GHz — the phase
wraps once, exercising the unwrap) beyond the first probe.

## Gate `engine-sparams-003`

1. **Magnitude**: mean |Γ| over 4–6 GHz within ±15 % of unity (a short
   reflects everything; the tolerance covers substrate leakage past the
   single-column via and fit residuals).
2. **Phase**: unwrapped-phase slope vs the round-trip closed form
   `dφ/df = −4π d √ε_eff / c` (twice the R.2 through-line slope), ±5 %,
   with `d` cell-snapped to the actual via placement.

## Explicitly still open after R.2b

Fixture de-embedding for the THRU S11 (removing launch/port residuals from a
*small* measured Γ) — a calibration-standards problem (SOL-style), separate
from this known-one-port validation and deferred until a use case forces it.
