# Filter Phase F1.2.0 — closed-form edge-coupled dimensional synthesis — Design Spec

**Phase:** F1.2.0 · **ADR:** ADR-0097 · **Date:** 2026-05-30 · **Status:** Accepted

## Goal

Turn a synthesized `yee-filter::CouplingMatrix` into **physical microstrip
dimensions** for an edge-coupled half-wave BPF by inverting the already-validated
`yee-layout` models — closed-form, pure f64, WASM-safe, NO FDTD, NO surrogate.
The first stage that makes the abstract network concrete geometry; seeds the
later EM-in-loop refinement (F1.2.1).

## API (`yee-filter`, new `dimension` module re-exported at crate root)
```rust
pub struct EdgeCoupledDimensions {
    pub line_width_m: f64,
    pub resonator_length_m: f64,
    pub gaps_m: Vec<f64>,      // length N-1 (inter-resonator)
    pub target_k: Vec<f64>,    // length N-1; the FBW·m_{i,i+1} each gap solves
}
pub enum DimError { /* e.g. UnsupportedTopology, OrderTooSmall, GapNotBracketed */ }

pub fn dimension_edge_coupled(
    project: &FilterProject,
    substrate: &yee_layout::Substrate,
) -> Result<EdgeCoupledDimensions, DimError>;

/// Convenience: assemble a `yee_layout::Layout` from the synthesized dims via
/// the existing `edge_coupled_bpf`.
pub fn dimension_edge_coupled_layout(
    project: &FilterProject, substrate: &yee_layout::Substrate,
) -> Result<yee_layout::Layout, DimError>;
```

## Algorithm (all closed-form / bisection)
1. **Line width:** `w = yee_layout::microstrip_width(spec.z0_ohm, εr, h)`. (`εr`,
   `h` from the `Substrate`.)
2. **Resonator length:** `ε_eff = yee_layout::eps_eff(w, h, εr)`; `λ_g = c /
   (f0·√ε_eff)`; `resonator_length_m = λ_g / 2`. (`c = 299_792_458.0`.)
3. **Inter-resonator gaps:** for `i = 0..N-2`, `target_k[i] = spec.fbw ·
   coupling.m[i][i+1]` (the off-diagonal; equals `yee-synth`'s
   `k_{i,i+1}=FBW/√(g_i g_{i+1})` — cross-check this equality in a comment).
   Solve for the gap `s` such that
   `coupling_coefficient(coupled_microstrip(w, s, h, εr)) == target_k[i]` by
   **bisection** over a bracket `[s_min, s_max]` (e.g. `[0.01·h, 50·h]` or in
   metres a sane `[5 µm, 5 mm]`): `coupling_coefficient` is **strictly decreasing
   in `s`** (proven by `coupled_002`), so bisect to a relative tolerance
   (≤ 1e-4). If `target_k` is outside the bracket's achievable range, return
   `DimError::GapNotBracketed` (do NOT silently clamp).
4. Assemble into `EdgeCoupledDimensions`. The `_layout` helper maps these into
   `yee_layout::EdgeCoupledParams` + calls `edge_coupled_bpf` (read that struct's
   fields; fill width/length/gaps; if it needs a feed/end gap and `qe` is not yet
   mapped, use a documented placeholder = the first inter-resonator gap, and note
   `qe`→feed dimensioning is F1.2.1 — do NOT invent a `qe`→gap formula).

Only `Topology::CoupledResonator` is supported; other topologies →
`DimError::UnsupportedTopology`. `N < 2` → `DimError::OrderTooSmall` (no
inter-resonator coupling to realize).

## DoD (machine-checkable; pure math, NO FDTD)
1. `cargo fmt --check --all` exit 0.
2. `cargo clippy -p yee-filter --all-targets -- -D warnings` exit 0.
3. `cargo test -p yee-filter` exit 0 (sub-second).
4. **`dim-001` (inversion round-trip, `tests/`):** build the committed Chebyshev
   0.5 dB N=5 BPF spec (f0 = 2 GHz, FBW = 0.10, Z0 = 50 Ω) + a concrete substrate
   (e.g. εr = 4.4, h = 1.6 mm — FR-4, matching the existing `yee-layout` tests),
   `synthesize` → `dimension_edge_coupled`, then for each `i` assert
   `coupling_coefficient(coupled_microstrip(line_width, gaps_m[i], h, εr))`
   reproduces `target_k[i]` within **< 1 %** relative.
5. **`dim-002` (physical sanity, `tests/`):** for the same design — every
   `gaps_m[i] > 0`; gaps **strictly decrease** as `target_k` increases (sort by
   target_k, assert monotone); `line_width_m == microstrip_width(z0, εr, h)`
   (exact); `resonator_length_m` within ±2 % of `c/(2 f0 √eps_eff)`. **Range
   bounds (note — these are split by dimension kind):** the in-plane *coupling
   features* (`line_width_m`, every `gaps_m[i]`) lie in `[1 µm, 20 mm]`; the
   *resonator length* is a half guided-wavelength (≈ 41 mm at 2 GHz on FR-4, so
   it CANNOT fit a 20 mm cap) and is bounded by its own physical window
   `[1 mm, 200 mm]`. Do NOT apply the 20 mm cap to the resonator length.
6. **`dim-003` (`tests/`):** `dimension_edge_coupled_layout` returns a `Layout`
   with ≥ 1 polygon and the expected port count; `EdgeCoupledDimensions` `serde`
   round-trips byte-identically.

## Out of scope
`qe`→I/O feed-coupling dimensioning; surrogate-BO / any EM-in-loop; FDTD; the
exact `N+1`-section J-inverter synthesis; asymmetric/stepped lines — all F1.2.1+.
No `yee-layout` change (consume its existing public API only).
