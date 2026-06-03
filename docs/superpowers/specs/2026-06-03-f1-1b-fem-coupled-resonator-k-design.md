# F1.1b coupling-coefficient extraction via FEM driven-sweep — design

**Date:** 2026-06-03
**ADR:** [ADR-0155](../../src/decisions/0155-f1-1b-fem-coupled-resonator-k.md)

## Problem

The filter app must EM-validate its design engine by extracting the inter-resonator coupling `k`
from a full-wave sim of a coupled-resonator pair. FDTD-resonant extraction was abandoned (ADR-0108
cavity wall). The FEM frequency-domain driven sweep (ADR-0153/0154) is wall-free; a de-risk probe
proved it resolves the even/odd split and recovers k within 17–26 % of the analytic reference (GO).

## Goal

Productionize the probe (`spike/fem-coupled-k-probe`, `933940f`,
`crates/yee-fem/tests/coupled_k_probe.rs`) into a `yee-fem` coupled-resonator-k extraction API + a
non-circular `fem-coupling-001` gate, then a small k-vs-gap monotonicity sweep (K2).

## Non-goals

- Reopening the FDTD resonant-split (ADR-0108) / cavity wall (ADR-0133) / MoM port (ADR-0064) /
  the FEM strict-in-mask filter follow-on.
- The Layout→FEM-mesh app integration (auto-meshing a filter `Layout`) — a later F1.2.x step. K1/K2
  validate the physics on a directly-built coupled-pair geometry.
- Qe extraction (K3, a later increment).

## Architecture

### K1 — coupled-resonator-k API + `fem-coupling-001` gate

- **New `yee-fem` API** (promote the probe helpers from the test into `src`, mirroring how ADR-0154
  N1 promoted `microstrip_port_numerical`): a `CoupledResonatorGeom` (W, S, h, ε_r, f0 → λ_g/2
  length, box) + `pub fn coupled_resonator_k(geom, sweep) -> Result<CoupledKResult, Error>` that
  builds the two-λ_g/2-resonator mesh (`layered_microstrip_filter_mesh` + `TraceRect`), attaches two
  weakly-gap-coupled `microstrip_port_numerical_at` feeds, runs `sweep_matrix` with
  `with_coupled_whitney(true)`, finds the two |S21| peaks (reuse `yee_filter::extract_coupling`),
  and returns `CoupledKResult { f_lo, f_hi, k_fem, peaks_resolvable, valley_db, ... }`. `yee-filter`
  becomes a `yee-fem` dev-dep (test-side) — or a real dep if the API lives in `src` and calls
  `extract_coupling` (acyclic: yee-filter does not depend on yee-fem; confirm).
- **Gate** `fem-coupling-001` in `crates/yee-fem/tests/` (`#[ignore]`'d + `--release`): asserts
  (a) **two resolvable peaks** — `peaks_resolvable` true AND the valley is a real margin (e.g. ≥6 dB)
  below the shallower peak (a re-smearing tripwire); (b) `|k_fem − k_imp| / k_imp ≤ 0.30` where
  `k_imp = coupling_coefficient(&coupled_microstrip(W,S,h,ε_r))` (the synthesis-side reference); and
  it **reports** `k_eps` (the ε_eff even/odd-split reference, ~17 % in the probe) for traceability.
  Non-circular (KJ closed-form vs full-wave FEM). Do NOT weaken the tolerance to force green; the
  probe measured 17–26 %, so 30 % vs `coupling_coefficient` has honest headroom + catches a real
  regression.

### K2 — k-vs-gap monotonicity

- A small sweep over ≥3 gaps S (e.g. 1.5 / 2 / 3 mm): assert `k_fem(S)` is monotonic-decreasing AND
  each point within the K1 tolerance of `coupling_coefficient(S)`. A *curve* is far harder to pass
  by coincidence than a single point — strengthens the validation.

### K3 (later) — Qe extraction

- Deferred to a later increment; out of scope here.

## Data flow

`CoupledResonatorGeom` → build coupled-pair mesh (`layered_microstrip_filter_mesh`/`TraceRect`) +
two weak `microstrip_port_numerical_at` feeds → `sweep_matrix` per-ω complex LU → |S21|(f) →
`extract_coupling` (two peaks) → `k_fem = (f_hi²−f_lo²)/(f_hi²+f_lo²)` → grade vs
`coupling_coefficient`.

## Testing

- K1/K2 gates `#[ignore]`'d + a dedicated `--release` CI job (mirror the `fem-eigen` gate job;
  `libfontconfig1-dev` already installed there). Heavy runs boxed (`scripts/yee-box.sh` ≤14 g/3 cpu).
- A fast non-ignored unit test (debug-safe): geometry well-formedness + the analytic-k references
  (no FEM solve).

## Risks

- **Tolerance honesty (low):** the probe's 26 % vs `coupling_coefficient` is within the 30 % gate;
  the known systematic (peaks pulled low by mesh dispersion + feed-gap load) is a tightening lever
  for K2+, documented, not hidden.
- **Weak-coupling setup (low–med):** the feeds must stay weakly coupled (over-coupling smears the
  split). The probe's gap-coupled feed (1 dy cell) worked; K1 keeps it. If a geometry refactor
  breaks weak coupling, the two-peaks tripwire catches it (fails honestly, not silently).
- **`yee-filter` dep direction:** confirm `yee-filter` does not depend on `yee-fem` before adding
  the dep (acyclic). If it would cycle, keep `extract_coupling`'s tiny peak-finder inline in the
  gate test instead.
