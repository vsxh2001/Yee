# ADR-0022: Phase 1.3.1.1 numerical 2-D cross-section eigensolver — spec + plan, deferred implementation

**Status:** Accepted (spec/plan only; implementation deferred to follow-up tracks)
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

ADR-0019 shipped the closed-form rectangular-waveguide TE10 mode
under `ModalDistribution::Te10` and explicitly deferred arbitrary-
cross-section wave-port support to Phase 1.3.1.1. The deferred work
is a numerical 2-D FEM eigensolver: given a triangular mesh of an
arbitrary waveguide cross-section (microstrip, CPW, GCPW,
septum-loaded WR-90), extract the dominant quasi-TEM / quasi-TE
mode (`β`, `Z_w`, transverse mode profile) and feed it into
`WavePort::rhs`.

The sub-project is large enough to need a written design before
code lands:

- **Mesh layer.** A plain-Rust 2-D triangular mesh type
  (`TriMesh2D`) is needed in `yee-mesh`, decoupled from the Gmsh
  FFI so that hand-rolled fixtures (WR-90, septum-loaded WR-90)
  can be constructed in tests without dragging in the SDK.
- **Assembly.** Nedelec edge elements for the transverse electric
  field `E_t` (curl-conforming, the right space for the lossless
  vector wave equation), plus nodal Lagrange for the axial `E_z`
  component. The generalised eigenproblem is sparse and complex-
  symmetric.
- **Solve.** Two paths: a dense `SymmetricEigen` fallback (always
  available, exercised on small fixtures), and a sparse
  shift-and-invert Arnoldi via `arpack-rs` behind a feature gate
  (the production path for realistic cross-sections).
- **Validation.** Two gates: the analytic TE10 cross-check from
  ADR-0019 (numerical β within 0.1%, `Z_w` within 1%, mode profile
  `L²` error < 1%), and a slow-wave loaded-WR-90 sanity case
  (`β > k_0` when a septum is added).

The total sub-project is on the order of 2000 LOC across `yee-mesh`,
`yee-mom`, optional `arpack-rs` wiring, and two validation gates.
That is large enough that the CLAUDE.md §3 walking-skeleton-first
principle bites: ship the **spec and the plan** as a Phase-X.0 step,
let consumer code (`ModalDistribution::Numerical2D` variant,
`WavePort::with_numerical_cross_section` builder) be written ahead
of the eigensolve so callers have a stable target, and then ship
the matrix-assembly + eigensolve as follow-up steps against the
already-fixed API.

This is the same pattern as Phase 1.3.1.0 / 1.3.1.1 itself
(ADR-0019 ships the analytic mode first; numerical comes later) and
as Phase 1.1.1.0 / 1.1.1.1 / 1.1.1.2 (ADR-0020 / 0024 / 0025 land
the multi-image DCIM, mesh refinement, and pole subtraction in
three separate sub-projects against the same `GreensSpec`
dispatch).

## Decision

Phase 1.3.1.1 is split into **spec + plan now, implementation in
follow-up tracks**. Track QQQQ delivers:

- `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-numerical-2d-eigensolver-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-17-phase-1-3-1-1-numerical-2d-eigensolver.md`
  — task-by-task TDD plan (eight steps; per-step files-touched /
  LOC budget / verification command).

The spec defines the abstraction surface that follow-up tracks
plug into:

```rust
pub enum ModalDistribution {
    Uniform,
    Te10(RectangularWaveguideTe10),               // Phase 1.3.1.0
    Numerical2D(Box<NumericalCrossSection>),      // Phase 1.3.1.1
}

pub struct NumericalCrossSection {
    pub mesh: TriMesh2D,                          // yee-mesh
    pub eps_r_per_tag: HashMap<MaterialTag, f64>,
    pub mu_r_per_tag: HashMap<MaterialTag, f64>,
    // Cached eigensolve results; populated by solve(), None until then.
    beta_cache: Option<f64>,
    z_w_cache: Option<f64>,
    profile_cache: Option<Vec<(f64, f64, f64)>>,  // (E_x, E_y, E_z) per node
}
```

`Numerical2D` is `Box`-wrapped to keep `ModalDistribution`'s enum
size bounded (ADR-0023 covers the `clippy::large_enum_variant`
rationale).

`WavePort::with_numerical_cross_section(mesh, eps_per_tag, mu_per_tag)`
is the builder. Until the eigensolve lands, the stub variant's
`solve()` returns `Error::NotImplemented` and `rhs` falls back to
`Uniform`'s delta-gap-equivalent behaviour.

**Approach (deferred to follow-up tracks):**

1. **Step 0–1 (Track BBBBB, ADR-0023).** `TriMesh2D` in
   `yee-mesh`; `ModalDistribution::Numerical2D` + stub in
   `yee-mom`. Closes the API surface so consumer code can land.
2. **Step 2 — sparse assembly.** Nedelec edge basis (`E_t`) +
   nodal Lagrange (`E_z`); assemble sparse `A` and `B` matrices
   for the generalised eigenproblem. `nalgebra-sparse`.
3. **Step 3 — dense fallback.** `SymmetricEigen` on the dense
   reduction; exercised on small fixtures, always-on path.
4. **Step 4 — sparse shift-and-invert Arnoldi.** `arpack-rs`
   behind a feature gate; production path for realistic meshes.
   Escape hatch: if `arpack-rs` is build-broken, drop the feature
   and ship dense-only.
5. **Step 5 — TE10 validation gate.** Analytic TE10 cross-check
   (β within 0.1%, `Z_w` within 1%, profile `L²` < 1%) on a
   hand-rolled WR-90 `TriMesh2D` fixture.
6. **Step 6 — slow-wave sanity gate.** Septum-loaded WR-90 case
   (`β > k_0` when the septum is added).
7. **Step 7 — `WavePort::rhs` wiring.** Replace the `Numerical2D`
   stub's fallback with the cached profile.

The spec surfaces three **out-of-lane findings** that Track QQQQ
flags but does not fix:

1. `ModalDistribution::Numerical2D` variant not yet present in the
   `yee-mom` enum (closed by Track BBBBB / ADR-0023).
2. `TriMesh2D` not yet in `yee-mesh` (closed by Track BBBBB /
   ADR-0023).
3. `yee-py` builder for numerical wave-ports deferred (closed by a
   later Python-bindings track, not this batch).

## Alternatives considered

1. **Ship the full Phase 1.3.1.1 eigensolver in one track.**
   Rejected. ~2000 LOC across two crates, two solver paths, two
   validation gates, and an optional FFI is past the working-set
   size where a single agent brief stays focused. The
   spec-then-implement pattern keeps each step small and
   verifiable.
2. **Skip the spec, write the eigensolver against an ad-hoc
   API.** Rejected. The `ModalDistribution::Numerical2D` variant
   has to slot into the same dispatch surface as `Uniform` and
   `Te10`; nailing that surface down before the matrix code lands
   prevents a costly retrofit. CLAUDE.md §3 ("sub-projects are
   decomposed before agents are dispatched") is explicit on this.
3. **Use a different mesh type (e.g. Gmsh-backed) for the
   cross-section.** Rejected. Cross-section meshes are tiny
   (hundreds to low thousands of triangles); the Gmsh SDK
   dependency is heavy for that scale and forces a feature-gate
   that the rest of the eigensolver does not need. `TriMesh2D` is
   purpose-built and trivial to hand-roll for fixtures.

## Consequences

**What becomes easier:**

- **Consumer code can be written ahead of the eigensolve.** The
  Track BBBBB stub (ADR-0023) gives the public API its final
  shape; follow-up tracks change the implementation without
  changing the surface.
- **Each follow-up step is small and verifiable.** Per-step LOC
  budgets and per-step verification commands in the plan mean a
  step that overruns its budget is visible immediately.
- **The TE10 cross-check (ADR-0019) becomes a literal bisection
  reference** for the numerical eigensolve. Any deviation > 0.1%
  on β is a bug in the numerical path, not a question of which
  reference is "right."

**What becomes harder:**

- **Microstrip / CPW / GCPW wave-ports remain delta-gap-equivalent
  until the eigensolve lands.** The stub returns
  `Error::NotImplemented` from `solve()`; only the rectangular-
  waveguide path (ADR-0019) and the `Uniform` default produce
  non-trivial modal behaviour.
- **The `Box<NumericalCrossSection>` indirection has a small
  per-port allocation cost.** For real workloads (handfuls of
  ports per simulation) this is irrelevant; for synthetic
  stress-tests it is measurable. ADR-0023 covers the rationale
  (clippy lint vs ergonomic cost).
- **The `arpack-rs` feature gate is a future test-matrix
  obligation.** CI will need to exercise the dense fallback by
  default and the sparse path on a feature-gated job. The
  escape hatch (drop arpack, ship dense-only) is in place if the
  crate proves unmaintainable.

**What's now closed off:**

- A `ModalDistribution::Numerical2D` variant with a different
  shape (e.g. inlined, non-`Box`'d). The Box wrap is part of the
  spec; changing it requires another ADR.
- An eigensolver that uses node-only Lagrange elements for `E_t`
  (spurious modes; well-known failure mode of nodal elements on
  the vector wave equation). Nedelec edge elements are mandated
  by the spec.

## References

- `docs/superpowers/specs/2026-05-17-phase-1-3-1-1-numerical-2d-eigensolver-design.md`
  — design spec.
- `docs/superpowers/plans/2026-05-17-phase-1-3-1-1-numerical-2d-eigensolver.md`
  — implementation plan.
- Commits 7802b80 (spec), fc7c28d (plan), b723bca (Track QQQQ
  merge).
- J. P. Webb, "Edge elements and what they can do for you," *IEEE
  Trans. Magnetics*, vol. 29, no. 2, pp. 1460–1465, March 1993
  (Nedelec edge elements vs nodal Lagrange on the vector wave
  equation).
- ADR-0019 — Phase 1.3.1.0 analytic TE10; this ADR's bisection
  reference.
- ADR-0023 — Track BBBBB stub; closes findings #1 and #2 surfaced
  here.
- CLAUDE.md §3 — walking-skeleton-first; §5 — multi-track
  orchestration; §6 — lane concept.
