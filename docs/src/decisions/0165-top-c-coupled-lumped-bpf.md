# ADR-0165: Top-C-coupled (capacitively-coupled) lumped BPF topology — de-risk first

**Status:** Accepted (track kickoff; de-risk brick first)
**Date:** 2026-06-06
**Related:** ADR-0164 (JLCPCB track; the narrow-band lumped BPF caveat — series resonators want sub-pF/
sub-nH, unrealizable), `yee_filter::{lumped, parts, jlcpcb}`, [[project-filter-design-final-goal]].

---

## Context

ADR-0164 (orderability capstone) showed the standard lumped BPF topology (alternating series/shunt L–C
resonators) is JLCPCB-orderable only for **wideband** filters; a **narrow-band** GHz BPF's *series*
resonators want sub-pF caps / sub-nH inductors below the discrete-part floor (~half the BOM unfillable).
The textbook fix for a manufacturable narrow-band lumped BPF is the **capacitively-coupled / top-C-coupled**
topology (Hong-Lancaster §8 / Zverev / Matthaei): `N` identical **shunt** L–C resonators tuned near `f0`,
coupled by `N+1` **series coupling capacitors** (admittance-inverters `J`). The shunt resonators are freely
realizable (pick a sane node `C`, get a realizable `L = 1/(ω0²C)`); the question is the **coupling caps**.

**Realizability concern (the de-risk target):** coupling caps scale as `Cij ≈ Jij/ω0` with `Jij ∝ FBW`, so
they SHRINK with frequency and narrow bandwidth. At 2 GHz / 10 % they may STILL be sub-pF (unrealizable) —
in which case top-C-coupled does NOT solve narrow-band *GHz* orderability (and the honest answer there is
the distributed/planar track), but likely DOES extend orderability to **lower-frequency** narrow-band
filters (where `Cij = J/ω0` is larger). So: **measure before building the full topology** — synthesize the
component values, validate the response, and probe the autopick-realizable envelope across `(f0, FBW)`.

## Decision

De-risk-first. **Brick T1 (this ADR's brick): synthesis + S21-validation + realizability probe** — NO new
ladder model / board / BOM wiring yet.

1. **Research + implement the top-C-coupled synthesis** (`yee_filter` or `yee-synth`): from the LP
   prototype g-values + `f0`, `FBW`, `Z0`, compute the shunt-resonator `L`/`C` and the `N+1` series
   coupling-cap values (the standard J-inverter / capacitively-coupled formulas — cite the source).
2. **Validate the synthesis (non-circular):** analyze the synthesized network's `S21` via an ABCD cascade
   (shunt-resonator admittance · series-coupling-cap impedance · …; reuse the `ladder_s21` Complex64 ABCD
   pattern) and assert it meets the target Chebyshev mask (passband ripple/return-loss + a stopband point)
   — the response from independent ladder analysis, NOT from the synthesis inputs.
3. **Realizability probe (the de-risk):** for a sweep of `(f0, FBW)` (e.g. 2 GHz & 0.5 GHz × 5–20 % FBW),
   feed the synthesized component values through `jlcpcb::autopick` and report the orderable-coverage
   envelope — WHERE (which f0/FBW) does top-C-coupled give a fully-orderable BOM (all caps incl. coupling
   ≥ the 1 pF floor, all inductors in range)? Honest finding either way (it extends the envelope to
   lower-freq narrow-band; GHz-narrow may still need distributed).

**Gate `top-c-coupled-001`:** synthesis S21 meets the Cheb mask for a representative narrow-band spec
(non-circular), AND the probe reports the orderable `(f0, FBW)` envelope (a measured table, not a claim).

## Consequences

- If T1 shows a useful orderable envelope (lower-freq narrow-band BPFs fully orderable), **T2** wires the
  topology into the ladder model + `lumped_board` + the JLCPCB export (so `--jlcpcb` emits an orderable
  narrow-band board) — a follow-on brick.
- If T1 shows top-C-coupled is *still* sub-pF-coupling-limited even at lower freq, that's an honest
  documented NO-GO (narrow-band lumped is fundamentally hard; distributed/planar is the path) — no T2,
  surfaced.
- Scope T1: `crates/yee-filter/src/` (a synthesis + ABCD-analysis module) + `crates/yee-filter/tests/`.
  Pure-math / WASM-safe. Research-first (cite the synthesis source). This ADR is the design record.
- **Not in scope (T1):** the ladder-model/board/BOM/CPL wiring (T2, only if T1 GO); distinct Q;
  the distributed/planar narrow-band path (a different track).

## References
- Synthesis: Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications* §8 (capacitively-coupled
  / J-inverter); Zverev, *Handbook of Filter Synthesis*; Matthaei-Young-Jones. (Cite the exact formulas in
  the module doc once researched.)
- Reuse: `yee_filter::lumped::ladder_s21` (Complex64 ABCD pattern), `yee_synth::prototype` (g-values),
  `yee_filter::jlcpcb::autopick` (the realizability probe).
