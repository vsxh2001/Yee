# ADR-0173: Group delay — a core filter metric, from the complex S21 phase (T10)

**Status:** Accepted
**Date:** 2026-06-06
**Related:** ADR-0172 (T9 — `coupling_matrix_s_params` complex S, real phase — the enabler), `yee_filter::
{ladder_s_params_lossy}`, ADR-0163 (studio Verify/response), [[project-filter-design-final-goal]].

---

## Context

T9 gave the distributed response real phase (`coupling_matrix_s_params`); the lumped finite-Q response
(`ladder_s_params_lossy`) is already complex. **Group delay** — `τ_g = −dφ/dω` of `S21` — is a core filter
performance metric (flat group delay = linear phase = low signal distortion; a real design spec for comms
filters). The app does not surface it. Now that both response paths carry physical phase, group delay is a
clean, valuable addition.

## Decision

1. **`yee_filter::group_delay(s21: &[Complex64], freqs_hz: &[f64]) -> Vec<f64>`** — topology-agnostic:
   unwrap the `S21` phase (`φ_k = arg(s21_k)`, adding ±2π to keep adjacent samples continuous), then
   `τ_k = −dφ/dω` by central difference (`−(φ_{k+1}−φ_{k-1})/(ω_{k+1}−ω_{k-1})`, `ω = 2πf`; one-sided at the
   ends). Returns seconds. Pure-math, WASM-safe. Re-export from `lib.rs`.
2. **Gate `group-delay-001`** (`crates/yee-filter/tests/`, pure-compute, non-`#[ignore]`'d, NON-circular):
   feed a synthesized filter's `coupling_matrix_s_params` `S21` (a dense sweep) to `group_delay` and assert:
   - **Closed-form midband anchor (the load-bearing non-circular check):** `τ_g(f0)` matches the
     prototype sum-rule `τ_g(ω0) = (Σ_{k=1}^{N} g_k)/(FBW·ω0)` (`ω0 = 2π·f0`) within tolerance. **[T10
     correction — see Outcome:** the **doubly-terminated** lowpass prototype DC group delay is
     `τ_LP(Ω=0) = Σg_k/2` (NOT `Σg_k`); with the bandpass map's `dΩ/dω|_{ω0} = 2/(FBW·ω0)` the `2` and the
     `½` cancel, so the correct anchor is `Σg_k/(FBW·ω0)`. This ADR's original `2·Σg/(FBW·ω0)` was off by
     2×, caught + corrected by the gate.**] The `g_k` come from `proj.prototype` — independent of the phase
     path, so agreement validates the group-delay computation non-circularly.
   - **Causality:** in-band `τ_g > 0` at every sample.
   - **Symmetry:** `τ_g(ω)` is symmetric about `f0` (within tol) for the symmetric synchronous BPF.
3. **Studio Verify readout** (`yee-studio-web`): compute `group_delay` over the in-band sweep (distributed via
   `coupling_matrix_s_params`; lumped via `ladder_s_params_lossy` if clean) and surface a **"midband / max
   in-band group delay: X ns"** number in the Verify stage (a readout, NOT a new plot — bounded + the
   computation is gated engine-side). Carry a `group_delay_ns` field on the engine view(s).

## Consequences

- The app surfaces group delay — a core filter metric, newly enabled by the T9 complex phase. Consumes T9.
- Scope T10: `crates/yee-filter/src/lib.rs` (the fn) + `crates/yee-filter/tests/` (the gate) + `crates/
  yee-studio-web/src/{engine.rs, stages.rs}` (the readout). Pure-math, WASM-safe; keep `wasm-build` green.
- **Not in scope:** a group-delay PLOT (a readout suffices for T10; a plot is a later add); group-delay
  EQUALIZATION/optimization; the CLI group-delay readout (the studio is the deliverable surface).

## Outcome (T10 — SHIPPED, merge `e82e1c1`)

`yee_filter::group_delay(s21, freqs) -> Vec<f64>` shipped (+110 lib, +319 gate; studio `group_delay_ns` on
`Designed`/`LumpedDesigned` + a Verify-stage "group delay (midband, τ@f0): X ns" readout, `None`→"—" for the
no-complex-S21 flows — no fabrication). τ = −dφ/dω via unwrapped-phase central difference (reviewer-verified
vs numpy to 4.4e-23 on a many-times-wrapping pure-delay signal).

**Physics correction (the gate caught this ADR's 2× error; the agent resolved it from first principles, NOT by
weakening the tol; reviewer independently re-derived it):** the **doubly-terminated** lowpass prototype DC group
delay is `τ_LP(0) = Σg/2`, not `Σg` (N=1: a shunt `C=g₁` between 1 Ω terminations → `S21 = 2/(2+jΩC)`,
`τ = C/2 = g₁/2`; confirmed N=1/3/5 via an independent ABCD ladder). So the correct midband anchor is
`τ_g(ω0) = (Σg/2)·(2/(FBW·ω0)) = Σg/(FBW·ω0)` — the `2` (Jacobian) and `½` (prototype) cancel.

**Gate `group-delay-001` (non-circular — the `Σg` anchor reads `proj.prototype.g`, distinct from the
coupling-matrix `m` the phase comes from):** midband τ(f0) matches `Σg/(FBW·ω0)` to **0.00 % rel_err** (N=3:
3.4133 vs 3.4134 ns; N=5: 6.6937 vs 6.6939 ns; ratio 1.0000 across FBW {0.10, 0.02, 0.005}; tol 5 % NOT
weakened — the wrong 2× form is 50–100 % off); causality (in-band τ>0); symmetry about f0 (1.9 % vs 2.5 %, the
real f↔Ω-map asymmetry that shrinks with the window). yee-filter + studio 18/18 green; wasm32 exit 0; no new
dep. Reviewer APPROVE, no P0/P1/P2 (anchor independently re-derived; gate non-circular + un-weakened;
unwrap/sign correct; readout honest). One P3 nit (while-vs-if in the unwrap, harmless).

## References
- Sum rule: Pozar, *Microwave Engineering* §8 (group delay / prototype); Hong & Lancaster §3 (the lowpass
  prototype `Σg_k` group-delay-at-DC relation + the bandpass transformation `dΩ/dω`).
- Phase source: `yee_filter::{coupling_matrix_s_params (ADR-0172), ladder_s_params_lossy}`; `proj.prototype`
  (the `g_k`).
