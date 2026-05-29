# ADR-0093: Filter Phase F1.1b.0 — coupling/Qe extraction algorithm

**Status:** Accepted
**Date:** 2026-05-29
**Related:** ADR-0091 (`yee-voxel`), `FILTER-DESIGN-ROADMAP.md` §5/§5a (F1.1b),
ADR-0084 (`yee-filter`)

---

## Context

F1.1 (FDTD coupling/Qe extraction) is the EM-in-the-loop building block. F1.1a
shipped `yee-voxel` (geometry → FDTD grid). F1.1b is "drive a coupled-resonator
pair through `yee-fdtd`, extract `k` and external `Qe`". That increment has two
separable parts:
- **the extraction *algorithm*** — turn a frequency response (two split peaks)
  into `k`, and a ring-down time series into `Qe`. Pure DSP, validatable against
  *analytic* signals with no FDTD run.
- **the FDTD *driver*** — voxelize a coupled pair, place `LumpedRlcPort`s, run,
  feed the result into the extractor. Heavy (multi-minute FDTD), validated
  against a published coupled-microstrip reference.

Validating the extractor against synthetic analytic signals (known `k`/`Qe`)
decouples "is the algorithm right?" from "is the FDTD coupling accurate?", and
the algorithm is the foundation the FDTD driver depends on. So **F1.1b.0 = the
extraction algorithm** (this ADR; light, no FDTD) → **F1.1b.1 = the FDTD driver**
(separate, heavy).

## Decision

Add a `extract` module to **`yee-filter`** (pure math; WASM-safe; no new dep):

```rust
pub struct CouplingExtraction { pub f_lo_hz: f64, pub f_hi_hz: f64, pub k: f64 }

/// Coupling coefficient from the two split resonance peaks of a synchronously-
/// tuned coupled-resonator pair: find the two most prominent local maxima in
/// `mag` (paired with `freqs_hz`), order them `f_lo < f_hi`, and return
/// `k = (f_hi² − f_lo²) / (f_hi² + f_lo²)`. `None` if fewer than two peaks.
pub fn extract_coupling(freqs_hz: &[f64], mag: &[f64]) -> Option<CouplingExtraction>;

/// Loaded/external quality factor from a resonator ring-down: least-squares
/// log-linear fit of the decaying upper envelope (local maxima of |samples|,
/// skipping the initial transient) gives the time constant τ; `Q = π·f0·τ`
/// (Pozar §6.1; mirrors yee-fdtd's cavity-Q decay fit). `None` if degenerate.
pub fn extract_q_ringdown(samples: &[f64], dt_s: f64, f0_hz: f64) -> Option<f64>;
```

This is the inverse of synthesis (measured response → coupling); it lives next
to the forward `synthesize`/`CouplingMatrix` in `yee-filter`, stays WASM-safe
(pure f64), and is the exact API the F1.1b.1 FDTD driver will call.

## Consequences

**Ships:** the `extract` module + the two functions. Gates (crate tests, §4 —
**analytic signals, no FDTD**): `extract-001` — a synthetic two-Lorentzian
`|response|` at `(f_lo, f_hi)` derived from a known `k` (centre 2 GHz) → recovered
`k` within ≤1% (peak-resolution limited); `extract-002` — a synthetic
`e^{−t/τ}·sin(2π f0 t)` with known `τ, f0` → recovered `Q = π f0 τ` within ≤3%.
Negative controls: a single-peak response → `extract_coupling` returns `None`;
flat/no-decay samples → `extract_q_ringdown` returns `None`.

**Not in scope (F1.1b.1):** any FDTD run; the coupled-resonator voxel driver
(`yee-voxel` + `LumpedRlcPort`); the published coupled-microstrip `k` gate.

**No new dependency** (pure f64). Lane: `crates/yee-filter/**`.

---

## References
- Pozar, *Microwave Engineering* 4e, §6.1 (loaded Q / ring-down), §8.8 (coupled
  resonators, split-frequency `k`); Hong & Lancaster ch. 8.
- `yee-fdtd` `cavity_q.rs` decay-fit (the Q precedent); the F1.1 recon.
- `docs/superpowers/specs/2026-05-29-filter-f1-1b0-coupling-extraction-design.md`;
  `docs/superpowers/plans/2026-05-29-filter-f1-1b0-coupling-extraction.md`.
