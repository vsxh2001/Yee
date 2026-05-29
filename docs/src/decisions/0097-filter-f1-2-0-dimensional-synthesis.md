# ADR-0097: Filter Phase F1.2.0 — closed-form edge-coupled dimensional synthesis

**Status:** Accepted
**Date:** 2026-05-30
**Related:** ADR-0084 (`yee-synth`/`yee-filter` synthesis + coupling matrix),
ADR-0086 (`yee-layout` HJ single-line model + `edge_coupled_bpf`),
ADR-0094 (`yee-layout::coupled_microstrip` even/odd model + coupler k),
`FILTER-DESIGN-ROADMAP.md` (F1.2)

---

## Context

The pipeline today goes spec → synthesis → **abstract** coupling matrix
(`yee-filter::CouplingMatrix` — normalized `m`, `qe_in`, `qe_out`). Nothing yet
turns that abstract network into **physical microstrip dimensions**. That mapping
is F1.2 ("dimensional synthesis"). The roadmap's eventual F1.2 is surrogate-BO
with EM-in-the-loop — heavy (multi-min FDTD per evaluation). But the *initial*
dimensioning needed to seed that loop is **closed-form and exact**:

- `yee-layout::coupled_microstrip` (ADR-0094) gives even/odd `Z0e`/`Z0o` →
  coupler `k = (Z0e−Z0o)/(Z0e+Z0o)`, and `coupled_002` proved `k` **strictly
  decreases** with the gap `s`. A strictly-monotonic function is invertible by
  **bisection** — so "find the gap that realizes a target `k`" is exact and
  cheap, no optimizer, no FDTD.
- `yee-layout::microstrip_width` already inverts `Z0 → width`; `eps_eff` gives
  the guided wavelength for the `λ_g/2` resonator length.

So F1.2.0 is the **light, closed-form** half of F1.2: map a `CouplingMatrix` to
edge-coupled microstrip dimensions by inverting the *already-validated*
`coupled_microstrip` + HJ models. The EM-in-the-loop refinement (and the
`qe`→feed-coupling dimensioning) is F1.2.1, gated on the F1.1b.1 FDTD driver.

## Decision

Add a `dimension` module to **`yee-filter`** (add a `yee-layout` dependency;
both are light-flow WASM-safe — `yee-layout` is `serde`-only, no cycle):

```rust
/// First-order physical dimensions of an edge-coupled half-wave microstrip BPF.
pub struct EdgeCoupledDimensions {
    pub line_width_m: f64,        // resonator/feed line width for the spec Z0 (HJ)
    pub resonator_length_m: f64,  // ≈ λ_g/2 at f0 (via eps_eff)
    pub gaps_m: Vec<f64>,         // N−1 inter-resonator coupled-section gaps
    pub target_k: Vec<f64>,       // the FBW·m_{i,i+1} each gap was solved for
}

/// Invert the validated coupled-microstrip model to size an edge-coupled BPF
/// from a synthesized coupling matrix + a substrate. Closed-form: line width
/// from HJ, resonator length from eps_eff, each inter-resonator gap by
/// bisecting `coupling_coefficient` (monotonic in gap) onto FBW·m_{i,i+1}.
pub fn dimension_edge_coupled(
    project: &FilterProject,
    substrate: &yee_layout::Substrate,
) -> Result<EdgeCoupledDimensions, DimError>;
```

Plus a convenience that assembles a `yee_layout::Layout` via the existing
`edge_coupled_bpf` from the synthesized dimensions.

**Mapping (first-order, narrowband — the standard initial-dimensioning
approximation):** for an `N`-pole edge-coupled half-wave filter, the
inter-resonator coupling `k_{i,i+1} = FBW · m_{i,i+1}` (= `yee-synth`'s
`FBW/√(g_i g_{i+1})`) is realized by a coupled section whose voltage coupling
`(Z0e−Z0o)/(Z0e+Z0o)` equals `k_{i,i+1}`. The line width is the spec-`Z0` HJ
width; the resonator length is `λ_g/2` at `f0`. This is an *initial* estimate to
seed EM refinement, not a final geometry.

## Consequences

**Ships:** the `dimension` module in `yee-filter` + a `yee-layout` dep. Pure
math, WASM-safe, no FDTD, no surrogate. First time the abstract coupling matrix
becomes concrete geometry.

**Gates (crate tests — pure math, NO FDTD):**
- **`dim-001` (inversion round-trip):** synthesize the committed Chebyshev 0.5 dB
  N=5 BPF, run `dimension_edge_coupled`, then evaluate `coupled_microstrip` on
  each `(line_width, gap_i)` and assert the recovered `coupling_coefficient`
  reproduces `target_k[i]` to < 1 % (bisection inverse correctness).
- **`dim-002` (physical sanity):** every gap > 0 and **strictly decreasing** as
  its target `k` increases (tighter coupling → smaller gap); `line_width_m`
  equals `microstrip_width(z0, εr, h)`; `resonator_length_m` ≈ `λ_g/2` within the
  model; all dimensions in a physically sane µm–mm range for the fixture.
- **`dim-003`:** the assembled `Layout` is non-degenerate (polygons, ports) and
  `EdgeCoupledDimensions` `serde` round-trips.

Correctness of the *physics* is inherited from `coupled_microstrip`
(`coupled-001` vs Steer) + `microstrip_width` (existing HJ gates); F1.2.0 adds
the *inversion* logic, which `dim-001` validates exactly. The end-to-end
**EM/published** validation of synthesized filter response is F1.3 (FDTD verify,
Swanson hairpin) — out of scope here.

**Not in scope:** `qe`→I/O feed-coupling dimensioning; surrogate-BO / EM-in-loop
refinement; the `N+1`-coupled-section J-inverter exactness — all F1.2.1.

**Constraint (ADR-0089):** `yee-filter` stays WASM-safe — the new `yee-layout`
dep is `serde`-only.

---

## References
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, ch. 8
  (coupling coefficients, edge-coupled resonator filters); Pozar §8.7.
- `docs/superpowers/specs/2026-05-30-filter-f1-2-0-dimensional-synthesis-design.md`;
  `docs/superpowers/plans/2026-05-30-filter-f1-2-0-dimensional-synthesis.md`.
