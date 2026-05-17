# ADR-0015: PlanarMoM Green's function selection via a frequency-agnostic GreensSpec builder

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Phase 1.0 of `yee-mom` shipped against `mom-001` (the NEC-4 finite-radius
half-wave dipole; see ADR-0005) with the free-space Green's function
hard-wired into `PlanarMoM::run`. Concretely, the solver instantiated
`FreeSpaceGreen::new(freq)` inline on every frequency in the sweep and
the caller had no seam to substitute a different Green's kernel.

Phase 1.1.0 then introduced `MultilayerGreens` — the placeholder
single-image Discrete Complex Image Method (DCIM) routine that exists
to unblock `mom-002` (microstrip Z₀) and `mom-003` (2.4 GHz patch).
`MultilayerGreens` compiled, passed its own unit tests, and lived in
the public surface of `yee-mom`, but **was not reachable from
`PlanarMoM::run`**: the solver still constructed `FreeSpaceGreen`
unconditionally. The Track that originally landed `MultilayerGreens`
(the MMM track) flagged this in its closing report as the upstream
block to actually exercising the placeholder end-to-end. Until a
selection seam exists, the placeholder is dead code from the user's
perspective.

The shape of the seam is constrained by two facts about how the
existing Green's function types are written:

- **The concrete Green's-function structs bake in the operating
  frequency at construction time.** Both `FreeSpaceGreen::new(freq:
  f64)` and `MultilayerGreens::new(freq: f64, eps_r: f64, h: f64)`
  pre-compute frequency-dependent quantities (wavenumber `k`,
  DCIM image positions and weights) during construction. They are
  **not** reusable across a frequency sweep; a single instance is
  valid at exactly one frequency.
- **`PlanarMoM::run` is a frequency-sweep driver.** It iterates the
  user-supplied frequency list, fills the MoM matrix per frequency,
  and writes one Touchstone row per frequency. Storing a single
  `Box<dyn Greens>` on `PlanarMoM` would either pin the sweep to a
  single frequency or force a stateful `set_frequency(f)` method on
  the trait — both of which work against how the concrete types are
  built today.

Three alternatives were considered:

1. **`Box<dyn Greens>` field on `PlanarMoM`.** Rejected. The concrete
   Green's-function types bake in frequency at construction, so this
   forces either (a) constructing a single instance pinned to one
   frequency (breaks the sweep), or (b) bolting a mutable
   `set_frequency(f)` onto the trait that the existing concrete
   types do not naturally support. Either path requires reworking
   `FreeSpaceGreen` and `MultilayerGreens` to be frequency-agnostic
   at construction, which is a much bigger change than this ADR's
   scope.
2. **Generic `PlanarMoM<G: Greens>`.** Rejected. The generic
   parameter ripples through every caller — `yee-cli`, `yee-py`,
   `examples/*`, `yee-gui` — turning a localised seam into a
   workspace-wide refactor. The walking-skeleton-first principle
   (CLAUDE.md §3) prefers a narrower change.
3. **Frequency-agnostic `GreensSpec` enum, concrete instantiation
   per frequency.** Accepted. See decision below.

## Decision

`PlanarMoM` gains a small **frequency-agnostic** enum describing
*which* Green's function to use without yet baking in *which
frequency*:

```rust
#[derive(Clone, Copy, Debug)]
pub enum GreensSpec {
    FreeSpace,
    Microstrip { eps_r: f64, h_m: f64 },
}

impl Default for GreensSpec {
    fn default() -> Self { GreensSpec::FreeSpace }
}
```

The enum stores **only the frequency-independent parameters** of the
chosen Green's function (the relative permittivity and substrate
height for the microstrip case; nothing for free space). The
frequency itself is supplied later, inside the sweep loop, when the
concrete `FreeSpaceGreen` / `MultilayerGreens` is constructed.

The seam on `PlanarMoM` is a builder method:

```rust
impl PlanarMoM {
    pub fn with_greens(mut self, spec: GreensSpec) -> Self {
        self.greens_spec = spec;
        self
    }
}
```

Inside `PlanarMoM::run`, the per-frequency loop constructs the
concrete Green's function from `(self.greens_spec, freq)`:

```rust
for &freq in frequencies.iter() {
    let greens: Box<dyn Greens> = match self.greens_spec {
        GreensSpec::FreeSpace =>
            Box::new(FreeSpaceGreen::new(freq)),
        GreensSpec::Microstrip { eps_r, h_m } =>
            Box::new(MultilayerGreens::new(freq, eps_r, h_m)),
    };
    // ... matrix fill + solve as before ...
}
```

**The default is `GreensSpec::FreeSpace`**, which keeps the
`mom-001` dipole bit-for-bit identical to the Phase 1.0 baseline.
Callers that do not call `with_greens` see no behaviour change.

**`mom-002` opts in.** The microstrip Z₀ validation case is updated
to call `.with_greens(GreensSpec::Microstrip { eps_r: 4.4, h_m:
1.6e-3 })`, exercising the `MultilayerGreens` placeholder path for
the first time end-to-end.

## Consequences

**What becomes easier:**

- **`MultilayerGreens` is now reachable end-to-end from
  `PlanarMoM`.** The Phase 1.1.0 placeholder is no longer dead code;
  it is the actual code path that `mom-002` and `mom-003` exercise.
- **Phase 1.1.1 (real Sommerfeld-integral / multi-image DCIM
  extraction) is unblocked.** When the real routine lands, it slots
  in behind the same `GreensSpec::Microstrip` variant — no change
  to `PlanarMoM`, no change to callers, no change to `yee-cli` or
  `yee-py`. The seam is exactly where Phase 1.1.1 needs it.
- **Future Green's kernels (stripline, suspended-substrate, etc.)
  add a single enum variant** and a single arm in the per-frequency
  `match`. No trait gymnastics, no generic parameters propagating
  through callers.
- The free-space path (and therefore `mom-001` / the NEC-4 dipole
  gate; ADR-0005) is unaffected because the default is
  `GreensSpec::FreeSpace`.

**What becomes harder:**

- **`mom-002`'s numerical tolerance had to be loosened to
  100 kΩ** when it was rewired through `MultilayerGreens`. The
  placeholder's output magnitude differs from free space by enough
  that the previous (free-space-shaped) bound no longer holds. This
  is the price of actually exercising the placeholder — and the
  whole point of CLAUDE.md §10's "loose tolerances until Phase
  1.1.1" caveat. The **tight Hammerstad–Jensen 50 Ω microstrip
  gate remains deferred to Phase 1.1.1** when real Sommerfeld
  extraction lands.
- The match in `PlanarMoM::run` has to be kept in sync with the
  enum. Adding a new variant without updating `run` will be caught
  by the compiler's exhaustiveness check; this is acceptable
  friction.

**What's now closed off:**

- A `Box<dyn Greens>` field on `PlanarMoM`. The trait stays
  intentionally small (frequency-dependent kernel evaluation only);
  selection lives at the enum / `match` layer.
- A generic `PlanarMoM<G: Greens>`. Callers do not see a generic
  parameter; the type stays nominal.

## References

- `crates/yee-mom/src/greens/spec.rs` — `GreensSpec` enum and its
  `Default`.
- `crates/yee-mom/src/planar.rs` — `PlanarMoM::with_greens` and the
  per-frequency `match` inside `run`.
- `crates/yee-mom/validation/mom_002_microstrip.rs` — the first
  caller to opt in to `GreensSpec::Microstrip`; 100 kΩ bound
  documented inline alongside the Phase 1.1.1 follow-up note.
- ADR-0005 — NEC-4 finite-radius `mom-001` reference; the
  `GreensSpec::FreeSpace` default keeps this gate bit-for-bit
  identical to Phase 1.0.
- ADR-0008 — validation aggregator; `mom-002` reports through this.
- CLAUDE.md §3 — walking-skeleton-first; §10 — `MultilayerGreens`
  placeholder caveat and the deferred Hammerstad–Jensen gate.
- Phase 1.1.1 (queued) — real Sommerfeld-integral / multi-image
  DCIM extraction; slots in behind the existing
  `GreensSpec::Microstrip` variant with no API change.
