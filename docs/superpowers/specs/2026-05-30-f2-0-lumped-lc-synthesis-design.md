# Filter Phase F2.0 — lumped-element LC ladder synthesis — Design Spec

**ADR:** ADR-0111 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

First brick of the **lumped-LC track** (new product goal: full lumped-LC filter
→ PCB with component choosing, EM sim, BOM, tolerance). F2.0 = closed-form
**LC ladder synthesis**: turn an abstract synthesized prototype (g-values) into
**ideal L/C element values** for a lumped band-pass filter. Pure `f64`, WASM-safe,
NO FDTD/parts/PCB — the foundation F2.1 (component selection/BOM), F2.3 (lumped EM
sim), F2.4 (tolerance) build on.

## Method (Pozar §8.3 / Hong & Lancaster ch. 3)

From a low-pass prototype `g0..g_{N+1}` (already produced by `yee-synth`), centre
`ω0 = 2π f0`, fractional bandwidth `Δ = FBW`, system `Z0`, the standard
**low-pass → band-pass** ladder transform maps each prototype element to a
series or shunt **LC resonator** (alternating along the ladder):

- **Series-branch resonator** (series L–C): `L_k = g_k·Z0/(ω0·Δ)`,
  `C_k = Δ/(ω0·Z0·g_k)`.
- **Shunt-branch resonator** (parallel L–C): `L_k = Z0·Δ/(ω0·g_k)`,
  `C_k = g_k/(ω0·Z0·Δ)`.

Each resonator is tuned to `ω0` (`L_k·C_k = 1/ω0²`). The first element is series
or shunt per the prototype convention (configurable; default shunt-first).

## Changes (`crates/yee-filter/**` ONLY)

- New `crates/yee-filter/src/lumped.rs`:
  - `pub enum LcBranch { Series, Shunt }`
  - `pub struct LcResonator { branch: LcBranch, l_henry: f64, c_farad: f64 }`
    (each is an L–C resonator tuned to f0).
  - `pub struct LumpedLadder { f0_hz, fbw, z0_ohm, resonators: Vec<LcResonator> }`
    (+ serde, like `EdgeCoupledDimensions`).
  - `pub fn synthesize_lumped(project: &FilterProject) -> Result<LumpedLadder, LumpedError>`
    applying the transform to `project.prototype` g-values. `LumpedError` for
    unsupported response (band-pass only for the skeleton).
  - A private ABCD-cascade helper `ladder_s21(&LumpedLadder, f_hz) -> Complex`
    (cascade each resonator's ABCD between `Z0` source/load) used by the gate.
- Re-export the public items from `lib.rs`.

## DoD (machine-checkable; pure-math, NO FDTD)

1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 — incl. a new gate `lumped_001`:
   - Synthesize the committed Chebyshev 0.5 dB N=5 fixture (f0=2e9, fbw=0.10,
     z0=50) via `yee_filter::synthesize` → `synthesize_lumped`.
   - Assert `resonators.len() == N` (5) and every resonator is tuned:
     `|L_k·C_k·ω0² − 1| < 1e-6`.
   - Compute the ladder `|S21|` via the ABCD cascade across the band and assert it
     **meets the same spec mask** the prototype does (passband ripple ≤ ripple_db
     with margin, in-band return loss ≥ mask RL, stopband point rejection ≥ mask
     dB) — i.e. the LC realization reproduces the synthesized design. This is the
     published-benchmark gate (self-consistent vs the synthesized response +
     textbook transform).
   - Sanity: element values are physical (L in nH–µH, C in pF–nF range for these
     numbers).

## Out of scope

Component/part selection + parasitics (F2.1); PCB/footprints (F2.2); FDTD lumped
sim (F2.3); tolerance/Monte-Carlo (F2.4); the Dioxus lumped UI track; low-pass/
high-pass/band-stop transforms (band-pass only for the skeleton). CLI/studio
wiring is a follow-on.

## Why this first

It is the pure-math, WASM-safe, immediately-validatable foundation of the whole
lumped track — mirrors the existing `synthesize` → `dimension_edge_coupled`
pattern, needs no new deps, and gives F2.1/F2.3/F2.4 their input (ideal L/C
values + the reference response).
