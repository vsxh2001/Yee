# Filter Phase F1.2.2 — hairpin dimensional synthesis — Design Spec

**ADR:** ADR-0109 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal
Add a SECOND filter topology to the synthesis→dimensions→layout pipeline:
closed-form **hairpin** band-pass dimensional synthesis, mirroring the shipped
F1.2.0 `dimension_edge_coupled`. Maps an abstract synthesized `CouplingMatrix`
to physical hairpin (U-folded half-wave) microstrip dimensions. Pure `f64`,
WASM-safe, NO FDTD/surrogate (the initial dimensioning that F1.2.1 BO later
refines). Directly advances the product goal's "topologies" (plural): today only
edge-coupled is dimensioned.

## Why hairpin reuses the edge-coupled coupling
A hairpin resonator is a half-wave line folded into a U; adjacent hairpins couple
through the **edge gap between their adjacent arms** — the SAME edge-coupled
mechanism `dimension_edge_coupled` already inverts. So the gap→k bisection is
identical; only the geometry (folded arms + tapped feed) differs. The
`yee_layout::hairpin_bpf` geometry generator already exists (Hong & Lancaster
ch. 6).

## Changes (`crates/yee-filter/**`; minimal `crates/yee-layout/**` only if needed)
- `crates/yee-filter/src/dimension.rs`:
  - `HairpinDimensions { line_width_m, arm_length_m, fold_spacing_m, gaps_m: Vec<f64>, target_k: Vec<f64> }` (+ serde, like `EdgeCoupledDimensions`).
  - `pub fn dimension_hairpin(project: &FilterProject, substrate: &Substrate) -> Result<HairpinDimensions, DimError>`:
    - `line_width_m = microstrip_width(z0, substrate)` (reuse).
    - `arm_length_m = λ_g/4 = c / (4·f0·√ε_eff)` — the U-folded half-wave resonator
      is two ≈λ/4 arms (vs edge-coupled's λ/2 straight length). Document the
      factor-4 (vs edge-coupled's factor-2) with the Hong & Lancaster reference.
    - `fold_spacing_m` — a documented closed-form choice (e.g. a few line widths;
      the two arms of ONE hairpin are weakly self-coupled — not the inter-resonator
      coupling — so a fixed sensible spacing is fine for the walking skeleton).
    - `gaps_m` — per adjacent-resonator coupling, the SAME gap-bisection as
      `dimension_edge_coupled`: solve `coupling_coefficient(coupled_microstrip(
      line_width, s, h, εr)) == target_k[i]` (`target_k[i] = fbw·m[i][i+1]`),
      bisection over `[GAP_MIN_M, GAP_MAX_M]`, `GAP_REL_TOL`, `GAP_MAX_ITERS`
      (reuse the existing consts).
  - `pub fn dimension_hairpin_layout(project, substrate) -> Result<Layout, DimError>`
    building a hairpin `Layout`. Prefer reusing `yee_layout::hairpin_bpf`. If
    `HairpinParams` carries only a single `coupling_gap_m` (it does today),
    EITHER (a) extend `HairpinParams`/`hairpin_bpf` to per-section gaps
    (`Vec<f64>`) in a minimal, backward-compatible `yee-layout` change, OR (b)
    pass a representative gap and document the uniform-gap walking-skeleton
    limitation in the doc-comment + ADR. Pick (a) if it's a small clean change;
    else (b). Surface the choice in the report.
- Wire `Topology::Hairpin` (the enum already exists) so `dimension_hairpin*` is
  selected when the project's topology is hairpin (if a dispatch point exists);
  otherwise expose the fns directly. Do NOT change `dimension_edge_coupled`.

## DoD (machine-checkable; pure-math, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0 (+ `-p yee-layout` if touched).
3. `cargo test -p yee-filter` exit 0 — sub-second; includes a new gate
   `hairpin_dim_001` (mirror `dim_001_inversion_roundtrip`): synthesize the
   committed Chebyshev N=5 fixture, `dimension_hairpin` on FR-4, assert each
   realized gap re-evaluates (`coupling_coefficient(coupled_microstrip(width,
   gap_i, h, εr))`) to its `target_k[i]` within < 1 %, AND assert
   `arm_length_m ≈ λ_g/4` and `gaps_m.len() == N-1`.
4. If `yee-layout` touched: `cargo test -p yee-layout` exit 0; the existing
   `hairpin_bpf` geometry gate stays green (do not weaken).

## Out of scope
FDTD validation of the hairpin (a later F1.1b-style gate); tapped-feed Qe
synthesis (reuse/defer); combline/interdigital (separate topologies, need
shorted-resonator + via models); CLI/studio wiring of `--topology hairpin` (a
follow-on, crosses lanes). Keep it closed-form + the per-section-gap round-trip
gate, exactly the dim-001 discipline.
