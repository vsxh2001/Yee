# F1.2.7 — Interdigital dimensional synthesis — Design Spec

**ADR:** ADR-0148 · **Date:** 2026-05-31 · **Status:** Accepted
**Follows:** F1.2.5 (ADR-0144, `dimension_combline`) — interdigital is a clean mirror:
the same coupled-resonator synthesis + `solve_gap` coupling, with the **θ = π/2 (full
λg/4) short-circuited resonator and NO loading cap** as the distinguishing realization.

## Problem

Interdigital is the last greyed studio gallery technique. It is a coupled-resonator
**band-pass** like edge-coupled / hairpin / combline, so it shares the synthesis (coupling
matrix → Qe/Mij) and the inter-resonator gap solve. What is interdigital-*specific* is the
resonator: a **straight λg/4 line short-circuited at one end, with adjacent resonators
grounded at *alternating* ends** (the interdigital finger structure) — and crucially **no
loading capacitor** (combline needs C_L only because its line is shortened to θ0 < π/2; the
interdigital λg/4 line resonates on its own).

This increment ships the **engine only** (`dimension_interdigital`), gated. The board
layout (the grounded-alternating comb) is F1.2.8; studio lighting is a later App increment —
the same 3-step decomposition combline used (engine → layout → lighting).

## Key insight: interdigital is the θ = π/2 limit of combline

A short-circuited microstrip stub of electrical length θ at f0 has input susceptance
`B(f) = −(1/Z0)·cot(θ·f/f0)`. Combline shortens the line to θ0 < π/2 (so `cot > 0`, an
inductive stub) and **adds** `C_L = cot(θ0)/(2π·f0·Z0)` to bring `B(f0) = 0`. **Interdigital
takes θ = π/2 (the full λg/4):** then `cot(π/2) = 0`, so `B(f0) = 0` with **no cap at all**.
`dimension_combline` deliberately *errors* at θ0 = π/2 (it would compute C_L = 0, non-physical
for a combline), so `dimension_interdigital` is a genuinely distinct function — same
coupling/gap machinery, different resonator (full λg/4, no cap).

## Method (`yee-filter`, `crates/yee-filter/src/dimension.rs`)

`pub fn dimension_interdigital(project: &FilterProject, substrate: &Substrate) ->
Result<InterdigitalDimensions, DimError>` — closed-form, a direct mirror of
`dimension_hairpin` / `dimension_edge_coupled` (NO θ0 parameter — θ is fixed at π/2 by
definition, unlike `dimension_combline`):

- **Line width** — spec-`Z0` Hammerstad-Jensen width (`yee_layout::microstrip_width`).
- **Resonator length** — `resonator_length_m = λg/4 = c / (4·f0·√ε_eff) = (π/2)/β(f0)`
  (`ε_eff` at the synthesized width via `yee_layout::eps_eff`). The factor-4 (λg/4) versus
  edge-coupled's factor-2 (λg/2); a *straight* λg/4 line (not folded like the hairpin arm).
- **NO loading cap** — `InterdigitalDimensions` has no `loading_cap_f` field (the structural
  contrast with `ComblineDimensions`).
- **Inter-resonator gaps** — identical to combline/edge-coupled: `target_k[i] = FBW ·
  m_{i,i+1}` realized by the shared `solve_gap` bisection over the monotone coupled-line
  coupling coefficient (no optimizer, no FDTD, no clamping). First-order: this reuses the
  same coupled-microstrip coupling model as the other techniques; the alternating-ground
  even/odd-mode refinement specific to interdigital is a deferred EM follow-on (the same
  scope boundary combline drew around the cap's effect on coupling).

`InterdigitalDimensions { line_width_m, resonator_length_m, gaps_m, target_k }` (the combline
struct minus `loading_cap_f`/`theta0_rad`, since θ is fixed and there is no cap).

**Errors** (mirror combline/hairpin): `UnsupportedTopology` if not `CoupledResonator`;
`OrderTooSmall` if N < 2; `GapNotBracketed` if a `target_k` is unreachable in the gap bracket.

## DoD — the gate `dim_interdigital_001` (machine-checkable, non-tautological)

Three parts, mirroring `dim_combline_001`:

1. **Published benchmark (H&L Qe/M) — the non-tautological synthesis core.** Synthesize the
   5-pole 0.1 dB Chebyshev band-pass at FBW = 0.10 and FBW = 0.15 (the H&L §5 worked
   coupled-resonator example) and assert the synthesized Qe / M₁₂ / M₂₃ reproduce H&L's
   *published* numbers (Qe ≈ 11.468, M₁₂ ≈ 0.07975, M₂₃ ≈ 0.06077 at FBW 0.10; Qe ≈ 7.645,
   M₁₂ ≈ 0.11962, M₂₃ ≈ 0.09115 at FBW 0.15) to < 1% (and the tighter 1e-3 band). This
   validates the synthesis chain feeding the interdigital dimensioner against the *book*,
   not against the synthesizer's own output. (Shared with combline's coupled-resonator
   synthesis — legitimately so; the coupling matrix is prototype-derived.)
2. **Interdigital λg/4 resonance (interdigital-DISTINCT).** From a `dimension_interdigital`
   result build the short-circuited-stub susceptance `B(f) = −(1/Z0)·cot((π/2)·f/f0)` (θ =
   π/2, **no cap term**) and root-find `B(f) = 0` over `[0.5 f0, 1.5 f0]`; assert the root
   equals f0 within ±1%. This is the interdigital analog of combline's resonance check — but
   the resonance comes from the λg/4 length alone (`cot(π/2) = 0`), NOT a cap. Catches a
   wrong length / dispersion / sign bug. Explicitly contrasts combline (which needs C_L).
3. **Dims solved + positive + structural.** N = 5 → 4 gaps; every solved gap re-evaluates to
   its `target_k` (< 1% via `yee_layout::coupling_coefficient`, no clamping); `line_width_m`
   and `resonator_length_m` finite and > 0; `resonator_length_m == (π/2)/β(f0) = λg/4`
   (closed-form, < 1e-9 rel); the struct carries **no loading cap**. `UnsupportedTopology` /
   `OrderTooSmall` error paths exercised.

## Changes

- `crates/yee-filter/src/dimension.rs` (`dimension_interdigital` + `InterdigitalDimensions`)
  + `crates/yee-filter/src/lib.rs` (re-export) + `crates/yee-filter/tests/dim_interdigital_001.rs`.
  NO studio / layout edits this increment.

## Out of scope

The grounded-alternating comb board layout (`dimension_interdigital_layout`, F1.2.8) and
studio lighting (a later App increment) — the next two steps. The alternating-ground even/odd
coupling refinement (a deferred EM follow-on, like combline's cap-coupling interaction). Tap
position / Qe→feed dimensioning (F1.2.1, shared across all coupled-resonator techniques).

## Why

Completes the coupled-resonator family's engine layer (edge-coupled / hairpin / combline /
interdigital) and unblocks the last gallery technique. Clean mirror of the shipped + gated
combline engine; the only new physics is the θ = π/2 λg/4 no-cap resonator, gated against
the published H&L synthesis + the closed-form quarter-wave resonance.
