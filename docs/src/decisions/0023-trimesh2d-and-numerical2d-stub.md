# ADR-0023: `TriMesh2D` in `yee-mesh` (plain-Rust, Gmsh-decoupled) and `ModalDistribution::Numerical2D(Box<…>)` stub

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0022 specced Phase 1.3.1.1 (numerical 2-D cross-section
eigensolver) and surfaced two out-of-lane findings: the
`ModalDistribution::Numerical2D` variant did not exist in `yee-mom`,
and the supporting 2-D triangular mesh type did not exist in
`yee-mesh`. Both are needed by the upcoming Nedelec edge-element
assembly (Phase 1.3.1.1 step 2) and by the validation fixtures
(hand-rolled WR-90 and septum-loaded WR-90).

Track BBBBB closes those two findings as a "step 0–1 stub" — the
*shape* of the API lands; the eigensolve is still
`Error::NotImplemented`. This is deliberate: the goal is to give
the follow-up steps a fixed target and to let consumer-side code
(`WavePort::with_numerical_cross_section`, downstream
`yee-py` builders eventually) be written ahead of the actual
eigensolve.

Two design questions need resolution before the stub lands:

1. **Where does the cross-section mesh type live, and how is it
   structured?** It cannot live behind the `gmsh` feature gate: the
   eigensolver validation fixtures (WR-90, septum-loaded WR-90)
   are hand-rolled and must build without the Gmsh SDK installed.
   It also should not live in `yee-mom` directly: mesh data is a
   `yee-mesh` concern by the workspace layout (CLAUDE.md §2), and
   future FDTD use cases (TF/SF surface meshes, NTFF probe
   surfaces) might want the same type.
2. **How is `ModalDistribution::Numerical2D` wrapped?** A bare
   `NumericalCrossSection` carries a `TriMesh2D` (hundreds to
   thousands of triangles, plus material-tag maps and cached
   eigensolve results). On a 64-bit target that is several
   hundred bytes per `Numerical2D` instance — comfortably over
   the `clippy::large_enum_variant` default threshold. The
   alternative variants (`Uniform`, `Te10(RectangularWaveguideTe10)`)
   are tiny (zero / ~24 bytes), so the enum size would be pinned to
   the largest variant.

The `clippy::large_enum_variant` decision is structural: every
`ModalDistribution` ever stored — on `WavePort`, in collections of
ports, in cached configuration — would carry the
`NumericalCrossSection`'s footprint. `Box`-wrapping the variant
moves the payload to the heap and brings the enum size back down to
the small-variant footprint, with a single allocation per
`Numerical2D` port (irrelevant for real workloads, which have
handfuls of ports per simulation).

## Decision

Track BBBBB ships **two changes**:

**(a) `TriMesh2D` in `yee-mesh`.** A plain-Rust 2-D triangular mesh
type, with no dependency on the `gmsh` feature:

```rust
pub type MaterialTag = u32;  // matches existing TriMesh tags

pub struct TriMesh2D {
    pub vertices: Vec<[f64; 2]>,
    pub triangles: Vec<[usize; 3]>,
    pub vertex_material: Vec<MaterialTag>,
    pub triangle_material: Vec<MaterialTag>,
}

impl TriMesh2D {
    pub fn new(...) -> Result<Self, Error>;
    pub fn area(&self, tri_index: usize) -> f64;        // signed shoelace
    pub fn centroid(&self, tri_index: usize) -> [f64; 2];
}
```

Construction enforces, by validation in `new`:

- `vertices.len() >= 3`, `triangles.len() >= 1`;
- every triangle index in range;
- **every triangle CCW with strictly positive signed area** (collinear
  / clockwise triangles are rejected);
- `vertex_material.len() == vertices.len()` and
  `triangle_material.len() == triangles.len()`.

Eight unit tests cover valid construction, out-of-range index,
collinear / CW rejection, area / centroid hand-calculations, and
material-length mismatches.

**(b) `ModalDistribution::Numerical2D(Box<NumericalCrossSection>)`
stub in `yee-mom`.**

```rust
pub struct NumericalCrossSection {
    pub mesh: TriMesh2D,
    pub eps_r_per_tag: HashMap<MaterialTag, f64>,
    pub mu_r_per_tag: HashMap<MaterialTag, f64>,
    // Cached eigensolve results; None until Phase 1.3.1.1 step 2+.
    pub(crate) beta_cache: Option<f64>,
    pub(crate) z_w_cache: Option<f64>,
    pub(crate) profile_cache: Option<Vec<(f64, f64, f64)>>,
}

pub enum ModalDistribution {
    Uniform,
    Te10(RectangularWaveguideTe10),
    Numerical2D(Box<NumericalCrossSection>),
}

impl WavePort {
    pub fn with_numerical_cross_section(
        mesh: TriMesh2D,
        eps_r_per_tag: HashMap<MaterialTag, f64>,
        mu_r_per_tag: HashMap<MaterialTag, f64>,
    ) -> Self { ... }
}
```

The stub's `solve()` returns `Error::NotImplemented` until the
eigensolve lands in Phase 1.3.1.1 step 2+. `WavePort::rhs` falls
back to `Uniform`'s delta-gap-equivalent behaviour for
`Numerical2D` ports until the cached profile is populated. The
mom-001 dipole gate is unchanged because it does not use
`WavePort`.

## Alternatives considered

1. **Put `TriMesh2D` behind the `gmsh` feature gate.** Rejected.
   The Phase 1.3.1.1 validation fixtures (WR-90, septum-loaded
   WR-90) must build on CI by default, which means without the
   Gmsh SDK. Feature-gating the mesh type forces the validation
   to either gate too, or to duplicate the type behind a different
   name. Plain-Rust avoids both.
2. **Put `TriMesh2D` in `yee-mom`.** Rejected. CLAUDE.md §2 puts
   mesh data in `yee-mesh`. Future use cases — TF/SF surface
   meshes, NTFF probe surfaces, possibly the Phase 2.fdtd.6.1
   calibrated TEM stripline — would want the same type. Owning it
   in `yee-mom` would force a downstream re-export or a
   duplicate type.
3. **Inline `Numerical2D` (no `Box`).** Rejected. Triggers
   `clippy::large_enum_variant` against the configured
   `-D warnings` floor (CLAUDE.md §3). The lint is correct on the
   merits: every `ModalDistribution` value would carry the
   `NumericalCrossSection`'s footprint regardless of variant. A
   `#[allow]` would be a workaround; `Box` is the fix.
4. **Use `Vec<f64>` for `vertices` (one flat array, stride 2)
   instead of `Vec<[f64; 2]>`.** Rejected. The
   `Vec<[f64; 2]>` shape is harder to misuse (no chance of a
   half-finished point), matches what FEM literature presents, and
   the allocation cost is identical.

## Consequences

**What becomes easier:**

- **The Phase 1.3.1.1 eigensolver work can start.** Surfaced
  findings #1 and #2 from Track QQQQ / ADR-0022 are closed. The
  follow-up assembly step (`A`, `B` sparse matrices over a
  `TriMesh2D`) has its inputs in place.
- **Cross-section meshes can be hand-rolled in tests** without
  Gmsh. The WR-90 fixture is a handful of vertices and two
  triangles; the septum-loaded variant is a handful more.
- **Other future cross-section-mesh use cases reuse the type.**
  TF/SF Huygens-surface meshes and NTFF probe meshes are
  candidates; the type is intentionally domain-agnostic.

**What becomes harder:**

- **`Numerical2D` ports incur one heap allocation per
  construction.** Irrelevant for real workloads (handfuls of
  ports); the `Box` cost is a one-time bookkeeping item.
- **The `solve()` stub returns `NotImplemented`.** Until Phase
  1.3.1.1 step 2 lands, any caller that constructs a
  `Numerical2D` port and calls `solve()` will hit the error.
  `with_numerical_cross_section` is therefore a *forward-compatible*
  builder, not a *usable* one, until the eigensolve lands.
- **CCW-winding enforcement at construction means** callers that
  hand-roll meshes have to get the winding right. The `Result`
  return type from `TriMesh2D::new` surfaces the error, and the
  test fixtures show the convention.

**What's now closed off:**

- A non-`Box`'d `Numerical2D` variant. The `Box` is part of the
  public API; un-boxing it later is a semver-meaningful change.
- A Gmsh-feature-gated cross-section mesh type. The `TriMesh2D`
  is plain-Rust by definition.
- A `TriMesh2D` that accepts collinear or clockwise triangles.
  Construction-time validation is mandatory; downstream code can
  rely on the invariants.

## References

- `crates/yee-mesh/src/lib.rs` — `TriMesh2D`, `MaterialTag`,
  validation in `TriMesh2D::new`.
- `crates/yee-mesh/tests/` — eight unit tests covering valid
  construction, out-of-range index, collinear / CW rejection,
  area / centroid, material-length mismatch.
- `crates/yee-mom/src/ports.rs` — `ModalDistribution::Numerical2D
  (Box<…>)`, `NumericalCrossSection`,
  `WavePort::with_numerical_cross_section`.
- Commits a7d865f (`TriMesh2D`), 8b14e9b (`Numerical2D` stub),
  661eb36 (Track BBBBB merge).
- ADR-0019 — Phase 1.3.1.0 TE10 (the other `ModalDistribution`
  variant).
- ADR-0022 — Phase 1.3.1.1 spec; the eigensolve that this ADR's
  stub is a placeholder for.
- CLAUDE.md §2 — workspace layout; §3 — `-D warnings` floor and
  the clippy convention this ADR's `Box` wrap respects.
