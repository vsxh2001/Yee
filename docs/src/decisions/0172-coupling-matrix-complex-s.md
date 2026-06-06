# ADR-0172: Complex coupling-matrix S-parameters — real phase for the distributed response/`.s2p` (T9)

**Status:** Accepted
**Date:** 2026-06-06
**Related:** ADR-0171 (T8 — studio `.s2p`; the distributed export is magnitude-only/flat-phase), `yee_filter::
{ideal_response, CouplingMatrix}`, `yee_synth::coupling_design`, [[project-filter-design-final-goal]].

---

## Context

The distributed/coupling-matrix response in `yee_filter::ideal_response` is **magnitude-only**:
`Complex64::new(s21_sq.sqrt(), 0.0)` — `|S21|` from the lowpass characteristic function, **zero phase**. So the
T8 distributed `.s2p` and the response plot carry flat phase (the ADR-0171 honest caveat). Phase / group delay
is a core filter quantity (a user comparing to a VNA or doing time-domain / cascade analysis needs it). The
principled fix is the textbook **coupling-matrix → S-parameters** synthesis, and the inputs already exist on
every `FilterProject`: `proj.coupling: CouplingMatrix { m: Vec<Vec<f64>> (N×N, synchronous→zero-diagonal),
qe_in, qe_out }`.

## Decision

Add a complex coupling-matrix S-parameter evaluator to `yee-filter`; do NOT change `ideal_response` (it is the
magnitude model used by `check_mask` and the response plot — leaving it avoids semantic churn at those call
sites; the new complex `|S21|` is validated to agree with it).

1. **`yee_filter::coupling_matrix_s_params(coupling: &CouplingMatrix, freqs_hz: &[f64], f0_hz: f64, fbw: f64) ->
   Vec<(Complex64, Complex64)>`** — the complex `(S11, S21)` per frequency via the Hong-Lancaster §8.1 N×N
   form (RESEARCH the exact signs/normalization — Hong & Lancaster *Microstrip Filters* 2nd ed eq (8.30)–(8.31)
   / Cameron *Microwave Filters* / the standard `[A]⁻¹` formulation):
   - Normalized lowpass frequency `Ω = (1/FBW)·(ω/ω0 − ω0/ω)`, `ω = 2πf`.
   - `[A] = [q] + jΩ·[U] − j·[m]`, where `[U]` = N×N identity, `[m]` = the normalized coupling matrix, and
     `[q]` = diagonal with `q₁₁ = 1/qe_in`, `q_NN = 1/qe_out`, else 0.
   - `S21 = (2/√(qe_in·qe_out))·[A]⁻¹_{N1}`, `S11 = 1 − (2/qe_in)·[A]⁻¹_{11}` (confirm the sign convention
     against the source; the magnitude-agreement gate will catch a wrong sign/scale).
   - The inverse columns `[A]⁻¹_{·1}` come from solving `[A]·x = e₁` with a **hand-rolled dense complex
     Gaussian elimination with partial pivoting** (N ≤ ~10; NO LAPACK / external solver — pure `num_complex`,
     WASM-safe). Re-export from `lib.rs`.
2. **Gate `coupling-matrix-s-001`** (`crates/yee-filter/tests/`, pure-compute, non-`#[ignore]`'d, NON-circular):
   for a representative synthesized filter (e.g. Cheb N=3 and N=5), over the in/near-band sweep:
   - **Magnitude agreement (the load-bearing non-circular check):** `|S21|` from `coupling_matrix_s_params`
     matches `|S21|` from `ideal_response` (the independent characteristic-function route) within a tolerance
     across the band. Two independent synthesis routes agreeing validates the complex model (a wrong sign/
     scale/Ω-mapping breaks this).
   - **Losslessness:** `|S11|² + |S21|² ≈ 1` (the synthesized `m` is lossless) within tolerance at every sweep
     point.
   - **Phase is non-trivial + continuous:** `S21`/`S11` carry varying (non-zero, non-constant) phase across the
     band, and the phase has no spurious ±2π jumps between adjacent close samples (sanity, not over-fit).
3. **Wire the studio distributed `.s2p`** (`yee-studio-web`): the `export_distributed` `.s2p` builder uses
   `coupling_matrix_s_params(&d.project.coupling, &freqs, f0, fbw)` for the `(S11, S21)` pair (real phase)
   instead of `ideal_response` + the lossless quadrature. The matrix is passive (lossless), so `yee_io`
   re-import accepts it. The distributed response PLOT keeps `ideal_response` (a magnitude plot; the complex
   `|S21|` agrees per the gate).

## Consequences

- The distributed `.s2p` carries physical phase — a complete, VNA/time-domain-comparable S-parameter export,
  removing the T8 magnitude-only caveat for the distributed path.
- Scope T9: `crates/yee-filter/src/lib.rs` (the evaluator + the complex solve helper) + `crates/yee-filter/
  tests/` (the gate) + `crates/yee-studio-web/src/{engine.rs, stages.rs}` (the `.s2p` re-wire). Pure-math,
  WASM-safe; keep the `wasm-build` job green.
- **Not in scope:** changing `ideal_response`'s magnitude semantics (kept for the mask/plot); a lossy
  coupling-matrix (finite-Q resonator loss in `[A]` — the lumped path already has finite-Q; a coupling-matrix
  finite-Q is a later add); the CLI `.s2p` distributed phase (a noted follow-on if it shares the gap).

## Outcome (T9 — SHIPPED, merge `2b11e13`)

`yee_filter::coupling_matrix_s_params` shipped (+190 lib, +237 gate; studio `distributed_s2p_sweep` re-wired,
−45). The complex `(S11,S21)` come from the Hong-Lancaster §8.1 `[A]=[q]+jΩ[U]−j[m]` form, solved per-frequency
by a hand-rolled complex Gaussian elimination (partial pivoting; pure `num_complex`, WASM-safe, no LAPACK).
`ideal_response` (magnitude, used by mask/plot) UNCHANGED — additive.

**Key normalization (reviewer-validated as a domain LAW, not a curve-fit):** the stored normalized `m`
(`1/√(g_i·g_{i+1})`) pairs with the **FBW-scaled** external Q (`qe·FBW`), not the stored physical `Qe` — the
`Ω`-map's `1/FBW` forces the diagonal loading into the same normalized domain. The reviewer independently
solved for the optimal continuous `qe` scale across 24 fixtures (N×FBW×ripple) → `opt_s/FBW = 1.000000`
universally (fixture-invariant ⇒ principled), and matched the solve vs numpy to 2.3e-13 over 1000 systems.

**Gate `coupling-matrix-s-001` (non-circular — `[A]⁻¹` matrix-solve `|S21|` vs the INDEPENDENT characteristic-
function `ideal_response`):** max `|S21|` dev **2.09e-5** (N=3/5, tol 2e-3, NOT weakened, ~400× below the
wrong-convention gap); losslessness `|S11|²+|S21|²−1` = **6.66e-16**; phase span ~6.2 rad (non-flat vs
`ideal_response`'s 0), continuous. Studio distributed `.s2p` now carries real phase; matrix passive →
re-importable; 18/18 studio tests green. wasm32 check exit 0; no new dep. Reviewer APPROVE, no P0/P1/P2
(independently re-derived the normalization law + the solve + gate non-circularity).

**Honest follow-on:** the CLI's distributed/coupling-matrix `.s2p` (`yee_cli::filter`) is still magnitude-only
(`ideal_response` + the CLI's own `lossless_s_pair`) — wiring it to `coupling_matrix_s_params` for CLI↔studio
parity is a small clean follow-on.

## References
- Synthesis: Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications* 2nd ed §8.1 eq (8.30)–(8.31);
  R. J. Cameron, *Microwave Filters for Communication Systems* (coupling-matrix → S).
- Inputs: `yee_filter::CouplingMatrix { m, qe_in, qe_out }` (from `yee_synth::coupling_design`),
  `yee_filter::ideal_response` (the magnitude reference for the gate).
- Consumer: `yee-studio-web` `export_distributed` `.s2p` (ADR-0171 T8).
