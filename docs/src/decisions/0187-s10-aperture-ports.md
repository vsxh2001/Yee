# ADR-0187: S.10 aperture ports on the engine — the port-mismatch ripple solved

**Status:** Accepted
**Date:** 2026-07-06
**Related:** ADR-0186 (named port mismatch as the next fidelity lever), ADR-0125/0126
(the original aperture-port formulation on `yee-fdtd`), ADR-0177 (E.2 drive layer).
**Spec:** `docs/superpowers/specs/2026-07-06-s10-aperture-ports-design.md`

## What was tried and what won

1. **Naive series stack** (the spec's original idea): N single-cell resistive ports down
   the substrate column, `R/N` + `V₀/N` each — expressible in existing protocol
   primitives. **Measured worse than the single-cell port** (mid-passband dip to
   −11.3 dB, S11 mean +7.1 dB — non-physical): each sub-port independently damps its own
   cell toward a uniform field, over-constraining the quasi-TEM profile. Rejected;
   the helper was removed.
2. **The validated aperture formulation** (`yee_fdtd::LumpedRlcPort::aperture`,
   Phase 2.fdtd.6.9): ONE aggregate branch against the modal voltage `V = ∫E_z·dz`
   (averaged over the width columns), semi-implicit two-way solve with
   `β = dt·h/(2·ε₀·A)`, back-action distributed as a sheet current over the **physical**
   area `A = w·h`. **Ported verbatim into `yee-compute`** (`Drive::aperture_ports`,
   pure-R arm) and gated **bit-exact** against the reference — `compute-014`
   (`cpu_aperture_parity.rs`), full final field state, max |Δ| = 0.0.

## Protocol

`JobSpec` gains `#[serde(default)] aperture_ports: Vec<AperturePortSpec>` — the client
names the modal face (`i`, `j_lo..j_hi`, `k_top`) and the aggregate R/EMF; the engine
derives `h = k_top·dx`, `A = (j_hi−j_lo)·dx·h` and validates bounds (Error events).
**CPU-only**: `backend: "gpu"` + aperture ports → `ComputeError::Unsupported` (a new
variant) surfaced as an error event; `"auto"` falls back to CPU. The GPU kernel needs a
per-port column reduction each step — queued for the GPU track.

## Measured effect (LPF gate, CPML-xy walls, all else identical)

| | single-cell ports (S.9) | **aperture ports (S.10)** | ideal |
|---|---|---|---|
| band edge 0.8 GHz | +12.39 dB | **+2.42 dB** | 0 |
| passband mean @1 GHz | +1.32 dB | **−0.49 dB** | 0 |
| transition 2.5 GHz | +8.53 dB (ripple) | **−13.18 dB** | −10.13 |
| stopband 3.0 GHz | −3.38 dB (ripple) | **−21.15 dB** | −17.68 |
| passband S11 | +7.1 dB (non-physical) | **−9.2 dB** | matched |
| cutoff | 1.900 GHz | 1.700 GHz | 2.0 GHz |

The standing-wave ripple that has dogged every board measurement since S.6 is gone;
the response now tracks the ideal Butterworth across the band. The gate's asserts
tightened accordingly: passband mean **±3 dB** absolute (was ±6 relative-only) and a
new **return-loss ≤ −6 dB** assert. The measured cutoff (1.700 GHz, 15 % low, still in
the ±20 % gate) is now a clean **design-side** signal — the staircased ~2-cell high-Z
sections — exactly what F1.2.1's EM-in-the-loop refinement exists to close.

## Consequences

Board-level verify now produces trustworthy absolute |S21|/|S11|. Follow-ons: aperture
ports on the GPU backend (column-reduction kernel, certified on the nightly); reusing
them in the stub/S-parameter scenario gates (kept as-certified for reproducibility);
F1.2.1 consuming the clean measurement.
