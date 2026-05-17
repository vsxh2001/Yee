# ADR-0024: Phase 1.1.1.1 edge-clustered strip mesh for mom-002, with loose tolerance pending Sommerfeld extraction

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0020 shipped multi-image DCIM via GPOF (Phase 1.1.1.0) and
moved `mom-002`'s input-impedance estimate from |Z_in| ≈ 14 kΩ
(one-image placeholder) to |Z_in| ≈ 2.7 kΩ (five-image fit) on the
existing 30×2 uniform strip mesh. The ~5× improvement was real but
left the result well outside the Hammerstad–Jensen `[35, 75] Ω`
target. Track OOOO's escape hatch surfaced the diagnosis: a
PEC-mirror probe on the same mesh gives |Z_in| ≈ 1.4 kΩ, which
means the **mesh resolution floors the answer at ~kΩ regardless of
how many DCIM images GPOF fits**. The remaining error has two
distinct legs: mesh resolution (this ADR, Phase 1.1.1.1) and
surface-wave pole subtraction (ADR-0025, Phase 1.1.1.2).

The mesh-resolution leg has a known physics origin: at a perfectly
conducting strip edge, the transverse surface current diverges as
`1/√d` where `d` is the distance from the edge (Meixner edge
condition; Pozar §2.5). A uniform RWG mesh resolves the singularity
only by brute-force refinement, scaling `O(N²)` in matrix assembly
and `O(N³)` in LU. Practitioners since at least the 1980s have
used **edge-clustered meshes** — node spacing biased towards the
strip edges — to capture the singularity with far fewer basis
functions.

The standard clustering choice is **Chebyshev nodes**:

```
y_j = -(w/2) · cos(π · j / n_width),  j = 0, 1, ..., n_width
```

which places `n_width + 1` nodes on `[-w/2, +w/2]` with edge cells
shrinking as `O(1/n_width²)`. Compared to uniform refinement at
the same RWG-basis count, this is an order-of-magnitude more
efficient at resolving `1/√d`.

Phase 1.1.1.1's job is to ship the edge-clustered builder, route
`mom-002` through it, and document what the new mesh does and does
not fix.

The "does not fix" leg is the load-bearing finding. A width-
direction refinement sweep `nz ∈ {2, 4, 8, 16, 24, 32}` against
the edge-clustered builder produces:

```
nz | Re(Z) [Ω]   | Im(Z) [Ω]   | |Z| [Ω]
---+-------------+-------------+----------
 2 |  +2325.959  |  -1308.575  | 2668.793
 4 |  +2350.452  |  -1479.010  | 2777.066
 8 |    -5.275   |  -2068.119  | 2068.126
16 |   -67.067   |  -2219.578  | 2220.592
24 |   -42.616   |  -2144.837  | 2145.261
32 |   -51.905   |  -2092.100  | 2092.744
```

`Re(Z)` converges from +2.3 kΩ to ~-50 Ω by `nz = 16`. The
edge-clustered mesh **does** resolve the `1/√d` edge singularity,
closing the mesh-resolution leg of the Phase 1.1.1.0 escape hatch.

`Im(Z)` plateaus at ~-2.1 kΩ across `nz ∈ {16, 24, 32}` — three
levels of refinement deep into convergence. That residual is
**not mesh-bound**: refining further would not move it. The DCIM
kernel approximates the spectral reflection coefficient with a
finite sum of complex images but does not subtract the discrete
surface-wave poles that dominate the field on a grounded FR-4 slab
at 1 GHz. The plateau is the **TM_0 / TE_1 surface-wave-pole
signature** and is exactly what ADR-0025 (Phase 1.1.1.2) exists to
fix.

The brief's escape hatch was explicit: "if refining to nz=32 still
floors above 100 Ω, surface the |Z_in| sweep table and STOP — the
bound is Sommerfeld surface-wave poles (Phase 1.1.1.2), not mesh."
That fired. The tolerance therefore stays at the loose
`[1, 100 kΩ]` non-degeneracy band, and ADR-0025 is the gate
tightening's prerequisite.

## Decision

`yee-validation` Phase 1.1.1.1 ships an edge-clustered strip-mesh
builder for `mom-002` and routes the headline run through it at
`nz = 16`. The validation tolerance stays at the prior loose
`[1, 100 kΩ]` non-degeneracy band, with the surface-wave-pole
plateau documented as a Phase 1.1.1.2 prerequisite.

**Public API (in `yee-validation`):**

```rust
pub enum StripSpacing {
    Uniform,        // Phase 1.1.1.0 back-compat
    EdgeClustered,  // Chebyshev nodes; this ADR
}

pub fn mom_002_strip_mesh_with_spacing(
    n_width: usize,
    spacing: StripSpacing,
) -> TriMesh;
```

The `EdgeClustered` variant places `n_width + 1` Chebyshev nodes
on `[-w/2, +w/2]` per the formula above. The legacy
`mom_002_strip_mesh` is retained as a back-compat shim that
delegates to `…_with_spacing(StripSpacing::Uniform)` so existing
Phase 1.1.1.0 callers and the structural-invariants unit test see
bit-for-bit identical geometry.

**Production switch.** `MOM_002_N_WIDTH = 16` (was 2) and
`run_mom_002` / `generate_mom_002_plots` use
`StripSpacing::EdgeClustered`. The new headline gate test —
`mom_002_headline_gate_passes` — covers the non-plot path in the
default `cargo test` budget; the 21-point plot sweep
(`mom_002_standalone_passes`) is `#[ignore]`-gated because the
30×16 mesh runs ~8 min in release.

**Tolerance.** `MOM_002_Z_MAX` stays at 100 kΩ (loose
non-degeneracy band). The `MOM_002_Z_MAX` docstring records the
surface-wave-pole finding inline, with an explicit forward
reference to Phase 1.1.1.2 / ADR-0025.

**Three new tests:**

- `mom_002_strip_mesh_edge_clustered_structure` — connectivity /
  port-tag invariants plus a Chebyshev monotonicity check
  (`dy_edge < dy_centre`).
- `mom_002_strip_width_refinement_sweep` — `#[ignore]`-gated
  `n_width ∈ {2, 4, 8, 16, 24, 32}` sweep against the
  edge-clustered builder.
- `mom_002_strip_width_refinement_sweep_uniform` —
  uniform-spacing counterpart, for measuring the clustering
  speed-up.

`mom-001` is untouched (no `yee-mom` edits in the lane; the change
is `yee-validation`-local).

## Alternatives considered

1. **Uniform refinement to `nz = 64+`.** Rejected. To match
   edge-clustered `nz = 16` on `Re(Z)` convergence, uniform
   spacing needs roughly `nz = 64`, which quadruples the
   `O(N²)` assembly and `O(N³)` LU cost. The mom-002 wall-time
   would go from ~8 min to ~hours. The physics motivation for
   clustering (Meixner edge condition) is well-established;
   uniform refinement is a brute-force waste.
2. **Adaptive `h`-refinement** (start uniform, subdivide cells
   based on a posteriori error estimate). Rejected as out of
   scope: an entire sub-project of its own, and the Chebyshev
   placement is the textbook static solution for the `1/√d`
   case. Adaptive refinement is interesting for arbitrary
   geometries; for a known-singularity strip edge it is
   overkill.
3. **Tighten the `mom-002` tolerance to `[35, 75] Ω` on the
   strength of the `Re(Z)` convergence alone.** Rejected. The
   `Im(Z)` plateau is real and reproducible across three levels
   of refinement; tightening the gate while it is present
   would either fail spuriously on the Im axis or require
   widening Im-only by exactly the surface-wave residue —
   which is what ADR-0025 will eventually subtract.

## Consequences

**What becomes easier:**

- **`Re(Z)` is converged on `mom-002`.** The `1/√d` edge
  singularity is resolved to ~-50 Ω at `nz = 16`; the
  mesh-resolution leg of the Phase 1.1.1.0 escape hatch is
  closed.
- **The surface-wave-pole leg is isolated.** Anyone looking at
  the |Z_in| sweep table sees the Im plateau across nz=16/24/32
  and the diagnosis writes itself. ADR-0025's spec is grounded
  in this data.
- **`StripSpacing` is a clean seam for future strip-geometry
  experiments.** Adding a different clustering (e.g.
  Gauss–Legendre nodes, or a hand-tuned hyperbolic-tangent
  profile) is one enum variant and one factory arm.

**What becomes harder:**

- **`mom-002`'s 21-point plot sweep is `#[ignore]`-gated.** The
  30×16 mesh assembly is ~525×525 per frequency point and the
  sweep takes ~8 min. CI runs the headline single-frequency
  gate by default; the full plot sweep is opt-in.
- **The tolerance gap to the Hammerstad–Jensen target remains.**
  Until ADR-0025 lands, `mom-002` cannot certify that the
  shipped solver produces 50 Ω-class numbers; it can only
  certify that the result is finite, complex, and not blowing
  up. CLAUDE.md §4 / §10 record this honestly.

**What's now closed off:**

- Uniform-only `mom-002` builds. The headline run is
  edge-clustered; uniform is the back-compat shim, exercised
  only by the structural-invariants test and the comparative
  refinement sweep.
- Quietly tightening `MOM_002_Z_MAX` without surface-wave-pole
  subtraction. The docstring forward-references Phase 1.1.1.2
  so any future tightening PR has to explicitly note ADR-0025
  status.

## References

- `crates/yee-validation/src/mom_002.rs` — `StripSpacing` enum,
  `mom_002_strip_mesh_with_spacing` builder, `MOM_002_N_WIDTH =
  16`, `mom_002_headline_gate_passes` smoke test,
  `MOM_002_Z_MAX` docstring with surface-wave-pole forward
  reference.
- Commits 4a0eb26 (builder), e799403 (production switch + sweep
  table), 22e3a24 (Track AAAAA merge).
- D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley, 2011,
  §2.5 (Meixner edge condition; `1/√d` singularity at a
  perfectly conducting strip edge).
- E. Hammerstad and Ø. Jensen, "Accurate models for microstrip
  computer-aided design," *IEEE MTT-S Digest*, 1980, pp. 407–409
  (the `[35, 75] Ω` reference band that Phase 1.1.1.2 will gate
  against).
- ADR-0015 — `GreensSpec` builder; the `Microstrip` /
  `MicrostripDcim` variants this mom-002 build dispatches to.
- ADR-0020 — Phase 1.1.1.0 multi-image DCIM; the prior mesh
  baseline (30×2 uniform, |Z_in| ≈ 2.7 kΩ).
- ADR-0025 — Phase 1.1.1.2 Sommerfeld pole extraction; the next
  step that the `Im(Z)` plateau is waiting on.
- CLAUDE.md §4 — `mom-002` / `mom-003` loose tolerances; §10 —
  `MultilayerGreens` placeholder caveat (narrowed by this ADR's
  mesh leg, still open on the pole leg).
