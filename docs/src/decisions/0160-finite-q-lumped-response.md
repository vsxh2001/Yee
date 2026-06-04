# ADR-0160: Finite-Q (realistic) lumped filter response + Cohn dissipation-loss gate

**Status:** Accepted
**Date:** 2026-06-04
**Related:** ADR-0158 (lumped-board CLI export), the lumped track (F2.x), [[project-filter-design-final-goal]]
(the app shows the filter response), `yee_synth::prototype` (g-values), `yee_filter::ladder_s21`.

---

## Context

The lumped filter response `yee_filter::ladder_s21` is **lossless** (documented at `lumped.rs:38,290`):
`LcResonator` carries only `l_henry` / `c_farad`, and the ABCD/S21 uses pure-reactive branches with the
lossless identity `|S11|² = 1 − |S21|²`. So the app can only show the *ideal* Chebyshev response —
flat-topped, zero midband insertion loss, infinitely sharp corners.

Real filters degrade from this because the L/C components have **finite unloaded quality factor `Q_u`**
(surface-mount inductors `Q ≈ 30–100`, capacitors higher). Finite `Q_u` produces **midband insertion
loss, rounded passband corners, and a shallower rejection skirt** — exactly what a user measures. For the
filter-design app ("see the parameters at the end"), the realistic finite-Q response is materially more
useful than the ideal one, and it is the honest curve to show alongside the mask.

This is a pure-analytic, fast, autonomous addition with a **published closed-form benchmark** to validate
against — Cohn's narrowband dissipation-loss formula.

## Decision

**1. Add a finite-Q lumped response.** A new `ladder_s21_lossy(ladder, f_hz, q_unloaded)` (the existing
lossless `ladder_s21` stays, behavior bit-identical — it is `q_unloaded = ∞`). Each resonator's
dissipation is modelled by a single unloaded `Q_u`:

- **series-branch resonator** (series L–C in a series ABCD arm): add a series resistance
  `R = ω₀·L / Q_u` so the branch unloaded Q is `Q_u = ω₀L/R`;
- **shunt-branch resonator** (parallel L–C in a shunt ABCD arm): add a shunt conductance
  `G = ω₀·C / Q_u` so the branch unloaded Q is `Q_u = ω₀C/G`.

`ω₀` is the band-centre. The lossy branch impedance/admittance is then complex with a real part, and S21
follows from the same ABCD cascade as the lossless path (no new topology — only the element model gains a
loss term).

**2. Validate against Cohn's dissipation-loss formula (published benchmark).** For a narrow-band
band-pass filter the midband insertion loss due to finite `Q_u` is (Cohn 1959; Hong-Lancaster §3.2):

```
IL₀ (dB) ≈ 4.343 · ( Σ_{k=1}^{n} g_k ) / ( Q_u · FBW )
```

where `g_k` are the low-pass prototype values and `FBW` the fractional bandwidth. Gate `lumped-q-001`
(fast, pure-compute, non-`#[ignore]`'d): synthesize the 3-pole 0.5 dB Chebyshev (`yee_synth::prototype`),
build the lossy ladder at a known `Q_u = 100`, measure `IL₀ = −20·log10|S21(f₀)|`, and assert it matches
Cohn's closed form within a narrowband tolerance (≤ 15 %, since FBW = 10 % is at the edge of the
narrowband approximation). For the 3-pole 0.5 dB Cheb, `Σg = g₁+g₂+g₃ ≈ 1.5963+1.0967+1.5963 = 4.289`, so
`IL₀ ≈ 4.343·4.289/(100·0.10) ≈ 1.86 dB`. **Non-circular:** Cohn's formula is the independent published
reference; the lossy ABCD S21 is the independent computation.

Sanity tripwires in the same gate: lossless (`Q_u = ∞`) gives `IL₀ ≈ 0` dB at f₀; finite-Q gives
`IL₀ > 0`; halving the loss (`Q_u → 2·Q_u`) roughly halves `IL₀` (the `1/Q_u` scaling Cohn predicts).

## Consequences

- The app gains a realistic insertion-loss curve (a follow-on can plot ideal vs finite-Q together);
  surfacing it in the UI is out of scope here (UI is framework-verdict-gated).
- One new published-benchmark validation case (`lumped-q-001`) — fast, in the routine CI matrix (no heavy
  EM). Strengthens the lumped track, which is the app's mask-clearing (analytic) design path.
- Scope: `crates/yee-filter/src/lumped.rs` (the lossy response) + `crates/yee-filter/tests/` (the gate).
  Small + well-defined ⇒ no separate spec/plan; this ADR is the design record.
- **Not in scope:** per-component distinct `Q_L`/`Q_C` (one lumped `Q_u` per resonator suffices for the
  Cohn benchmark; distinct-Q is a trivial later extension), predistortion synthesis, the UI plot.

## References
- Lossless response corrected: `yee_filter::ladder_s21` (`crates/yee-filter/src/lumped.rs`).
- g-values: `yee_synth::prototype` (`crates/yee-synth/src/lib.rs`), mirrored by `synth_001_gvalues.rs`.
- Method: S. B. Cohn, "Dissipation Loss in Multiple-Coupled-Resonator Filters," Proc. IRE, 1959;
  Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, §3.2 (the `4.343·Σg/(Q_u·FBW)`
  midband-loss formula).
