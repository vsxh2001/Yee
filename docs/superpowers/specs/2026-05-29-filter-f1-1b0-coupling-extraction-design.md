# Filter Phase F1.1b.0 — coupling/Qe extraction algorithm — Design Spec

**Phase:** F1.1b.0 · **ADR:** ADR-0093 · **Date:** 2026-05-29 · **Status:** Accepted

## Goal
The pure-DSP extraction that turns an EM/measured response into a coupling
coefficient `k` (two split peaks) and an external `Qe` (ring-down decay) —
validated against analytic signals, no FDTD. The foundation the F1.1b.1 FDTD
driver feeds. `yee-filter` only; WASM-safe; no new dependency.

## API (`yee-filter`, new `extract` module, re-exported at crate root)
```rust
pub struct CouplingExtraction { pub f_lo_hz: f64, pub f_hi_hz: f64, pub k: f64 }

pub fn extract_coupling(freqs_hz: &[f64], mag: &[f64]) -> Option<CouplingExtraction>;
pub fn extract_q_ringdown(samples: &[f64], dt_s: f64, f0_hz: f64) -> Option<f64>;
```

### `extract_coupling`
- Require `freqs_hz.len() == mag.len() >= 5`.
- Find interior local maxima: index `i` (1..n-1) with `mag[i] > mag[i-1]` and
  `mag[i] > mag[i+1]`. Take the two with the largest `mag`; if fewer than two
  distinct maxima exist → `None`.
- Order their frequencies `f_lo < f_hi`. `k = (f_hi² − f_lo²) / (f_hi² + f_lo²)`.
- Return `CouplingExtraction { f_lo_hz, f_hi_hz, k }`.

### `extract_q_ringdown`
- Build the upper envelope: indices `i` (1..n-1) that are local maxima of
  `|samples[i]|`; skip those in the first 1/3 of the record (initial transient,
  matching `cavity_q.rs`). Require ≥ 3 envelope points with strictly positive
  magnitude.
- Least-squares fit `ln|env_k| = a − t_k/τ` (`t_k = i_k·dt_s`); slope `m < 0` →
  `τ = −1/m`. `Q = π · f0_hz · τ`. Return `None` if `m ≥ 0` (no decay) or the fit
  is degenerate (all-equal magnitudes / <3 points).

## DoD (machine-checkable; pure math, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 (fast).
4. `extract-001` (coupling): pick `f0 = 2e9`, `k_true = 0.04` → split
   `f_lo = f0/√(1+k)`, `f_hi = f0/√(1−k)` (Hong-Lancaster). Build `mag[i]` as the
   sum of two Lorentzians centred at `f_lo,f_hi` (modest Q so peaks resolve) over
   a 0.7–1.3·f0 sweep (≥ 401 pts). `extract_coupling` returns `Some`, and the
   recovered `k` is within `≤ 1e-2` of `k_true` (peak-bin limited); `f_lo < f_hi`.
   Negative control: a single Lorentzian → `None`.
5. `extract-002` (Q): `f0 = 2e9`, `τ = 5e-9`; `samples[n] = exp(−n·dt/τ)·sin(2π
   f0 n·dt)` with `dt` ~ `1/(40 f0)` over ~`8τ`. `extract_q_ringdown(samples, dt,
   f0)` returns `Some(q)` with `q` within `≤ 3 %` of `π f0 τ`. Negative control:
   a constant (no decay) → `None`.

## Out of scope
Any FDTD run; the coupled-resonator voxel driver; the published coupled-microstrip
`k` gate (all F1.1b.1). Hilbert-transform envelopes (local-max envelope suffices).
