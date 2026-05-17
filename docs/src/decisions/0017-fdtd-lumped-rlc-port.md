# ADR-0017: yee-fdtd Phase 2.fdtd.6 ships a lumped RLC port with an energy-dissipation validation gate

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

Through Phase 2.fdtd.5 (ADR-0014 TF/SF slab) `yee-fdtd` had **sources**
(soft, hard, TF/SF slab), **boundaries** (CPML), **far-field
extraction** (NTFF), and **dispersive bulk materials** (Drude /
Lorentz / Debye via ADE). What it did **not** have was any way to
terminate a transmission-line geometry into anything other than the
material itself — no resistive load, no series-RLC load, no way to
couple the FDTD region to a circuit.

This blocks the next class of practical problems: filter networks,
matched / mismatched terminations on a feed line, baluns, antenna
feed-port circuit models, anything where a discrete component sits
between two metal traces and the goal is to see how it shapes the
field. The Taflove canonical solution is **lumped-element
incorporation** (Taflove & Hagness 3rd ed., §15.4 "Linear Lumped
Circuit Elements"): inside the FDTD update, a single edge `(i, j,
k)` carries the constitutive relation of a discrete element instead
of the usual `ε` / `σ` Yee update. For a series R-L-C this looks
like an auxiliary 2-state ODE (inductor current, capacitor voltage)
evolved alongside the field update on that edge.

Phase 2.fdtd.6 ships this. The design questions that needed
resolution before the brief was dispatched:

1. **What is the public API shape?** A single `LumpedRlcPort`
   struct with constructors for the common cases (pure resistor,
   series RLC), placed on a single edge in a single axis direction,
   was chosen over a more general "arbitrary 2-terminal element on
   an arbitrary path" formulation. Walking-skeleton-first
   (CLAUDE.md §3): ship the case that 80% of tutorials need, defer
   the general case to a follow-up.
2. **What is the validation gate?** This is the load-bearing
   question. The textbook gate for a lumped element on a
   transmission line is the **analytic reflection coefficient
   Γ = (Z_L − Z_0) / (Z_L + Z_0)** — inject a pulse, FFT the
   reflected wave, compare `|Γ|` to the analytic prediction. The
   problem: that gate requires a **calibrated TEM stripline
   geometry** with a known `Z_0`, and **`yee-fdtd` does not yet
   ship a calibrated stripline**. The TF/SF slab (ADR-0014) is a
   plane-wave source, not a transmission-line feed; nothing in the
   current `yee-fdtd` surface stands up a 50 Ω TEM line with the
   modal field profile correct to better than a few percent. A Γ
   gate against an un-calibrated `Z_0` is worse than no gate — it
   would fail by an amount that depends on the (uncalibrated)
   line, not on the lumped element's correctness, and chasing
   that failure would burn time on the wrong sub-problem.

So Γ-against-analytic is **deferred to Phase 2.fdtd.6.1**, when the
calibrated TEM stripline lands. The Phase 2.fdtd.6 gate has to be
something that tests the lumped element *in isolation*, without
depending on an upstream `Z_0` calibration.

The chosen gate is **energy dissipation**: charge a region with
field energy, run the simulation, and assert that the **lumped
resistor dissipates the expected amount of total energy**. This
gate has the property that it depends only on the local
constitutive relation on the resistor edge — it does not care
whether the surrounding region is a clean TEM line or a sloppy one.
It catches the bugs that would matter (wrong sign on the V-I
relation, off-by-one in the time-stepping of the auxiliary ODE,
units error in the resistance) without depending on the missing
upstream calibration.

The selected thresholds are **≥ 0.3% of initial field energy
dissipated globally** (the resistor must do *something*) and
**> 5× more energy dissipated at the resistor edge than on any
adjacent free-space edge** (the dissipation has to localise on the
resistor, not be a global numerical artefact). The first threshold
catches "resistor coefficient is zero / wired wrong"; the second
catches "the dissipation is real but the wrong edge has it."

## Decision

`yee-fdtd` Phase 2.fdtd.6 ships a `LumpedRlcPort` struct with two
constructors and an energy-dissipation validation gate.

**Public API:**

```rust
pub struct LumpedRlcPort {
    pub edge: (usize, usize, usize),
    pub axis: Axis,                   // X, Y, or Z
    pub r_ohms: f64,
    pub l_henries: f64,               // 0.0 for pure resistor
    pub c_farads: f64,                // ∞ (or sentinel) for pure resistor
    // ... internal auxiliary state ...
}

impl LumpedRlcPort {
    pub fn pure_resistor(edge: (usize, usize, usize),
                         axis: Axis, r_ohms: f64) -> Self { ... }

    pub fn series_rlc(edge: (usize, usize, usize), axis: Axis,
                      r_ohms: f64, l_henries: f64,
                      c_farads: f64) -> Self { ... }
}
```

The port is attached to a `FdtdGrid` like any other source / sink
and participates in the per-step update.

**Discrete update.** On the resistor edge, the standard Yee `E`
update is replaced by Taflove §15.4's discrete update incorporating
the R-L-C constitutive relation. The series-RLC case carries an
auxiliary 2-state ODE — **inductor current** `i_L` and **capacitor
voltage** `v_C` — evolved with a stable centered-difference scheme
at the same time step as the field update. The auxiliary state
lives on the `LumpedRlcPort` instance, not on the grid; one port =
one ODE state pair regardless of grid size.

**Validation gates.**

- `crates/yee-fdtd/tests/lumped_resistor_energy.rs`:
  - **Pure-resistor case.** Initialise a localised Ez pulse,
    place a `pure_resistor(50.0)` on an edge, run for ~5000 steps.
    Assert (a) `dissipated_global / E0 ≥ 0.003` and (b)
    `dissipated_on_resistor_edge > 5.0 × dissipated_on_max_other_edge`.
  - **Series-RLC case.** Same protocol with `series_rlc(50.0,
    1e-9, 1e-12)`. Same two assertions. Compiles, runs, and
    self-tests that the ODE state evolves in the expected
    direction (current ramps from zero, capacitor voltage tracks
    integrated current).

**What `LumpedRlcPort` does NOT do in Phase 2.fdtd.6:**

- No quantitative Γ-against-analytic gate. **Deferred to Phase
  2.fdtd.6.1**, paired with the calibrated TEM stripline geometry
  it requires.
- No multi-axis lumped element (one with current components on
  more than one Yee axis simultaneously). Single-axis only.
- No oriented-arbitrary 2-terminal element (one whose
  endpoints are not aligned to a grid axis). Single-edge,
  axis-aligned only.
- No parallel-RLC topology. Series only; parallel is a one-line
  algebra change in the ODE but is out of scope until a user
  asks.

These restrictions are listed in the docstring of `LumpedRlcPort`
with explicit forward-references to Phase 2.fdtd.6.1.

## Consequences

**What becomes easier:**

- **Resistive loads on transmission-line stubs simulate
  correctly** for the cases that 80% of tutorials need: a 50 Ω
  termination, a matched / mismatched line, a simple R-C low-pass
  network terminated into 50 Ω. These are the cases the
  energy-dissipation gate certifies as correct, and they are the
  same cases the upcoming Phase 2.fdtd.6.1 Γ gate will tighten.
- **Series-RLC compiles, runs, and self-tests.** The struct, the
  auxiliary ODE state, and the per-step update are all in place;
  the Phase 2.fdtd.6.1 follow-up is purely a *quantitative*
  validation upgrade, not a code-shape change. When the calibrated
  stripline lands, the same `series_rlc(...)` construction will
  pass the tighter Γ gate without source changes.
- The walking-skeleton-first principle (CLAUDE.md §3) is upheld:
  the minimum end-to-end pipe ships first, with a gate that is
  *honest* about what it's testing (local constitutive relation,
  not transmission-line behaviour).

**What becomes harder:**

- **Users get qualitative, not quantitative, Γ today.** A
  `pure_resistor(50.0)` on the current geometry will dissipate
  energy and reflect a wave, but the reflection coefficient
  cannot be compared to `Γ_analytic = 0` (matched) within better
  than the line's uncalibrated `Z_0` error. Users wanting that
  comparison have to wait for Phase 2.fdtd.6.1.
- **Phase 2.fdtd.6.1 is non-optional.** The energy-dissipation
  gate is *sufficient* for a placeholder but *not* sufficient to
  call the lumped-port subsystem "validated" in the sense the rest
  of `yee-fdtd` is validated (CPML against 30 dB reduction; TF/SF
  against textbook 10× contrast). The follow-up is on the queue,
  not a "nice-to-have."
- **No multi-axis / oriented-arbitrary / parallel-RLC**. Anyone
  who wants those today has to wait for Phase 2.fdtd.6.1.

**What's now closed off:**

- A "minimal viable" gate weaker than energy dissipation
  (e.g. only checking that the field on the resistor edge changes
  sign). Energy dissipation is the floor; the Γ gate is the
  ceiling; nothing in between is acceptable.
- Quietly extending the API to multi-axis or oriented-arbitrary
  ports inside Phase 2.fdtd.6. That work has a follow-up phase
  with its own validation requirements.

## References

- `crates/yee-fdtd/src/lumped/rlc_port.rs` — `LumpedRlcPort`
  struct, `pure_resistor` / `series_rlc` constructors, and the
  per-step discrete update on the lumped edge.
- `crates/yee-fdtd/tests/lumped_resistor_energy.rs` —
  energy-dissipation gate (≥ 0.3% global, > 5× local).
- `docs/src/theory/fdtd-details.md` — lumped-element section,
  including the explicit deferral of the Γ gate to Phase
  2.fdtd.6.1 and the rationale for the energy-dissipation
  surrogate.
- A. Taflove and S. C. Hagness, *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method*, 3rd ed., Artech
  House, 2005, §15.4 "Linear Lumped Circuit Elements" (the
  derivation of the discrete R / L / C edge update used here).
- ADR-0008 — validation aggregator; the new gate reports through
  this.
- ADR-0014 — TF/SF slab; the soonest the calibrated TEM stripline
  that Phase 2.fdtd.6.1 needs can borrow plumbing from.
- CLAUDE.md §3 — walking-skeleton-first; CLAUDE.md §4 —
  no-solver-feature-without-a-published-benchmark-gate policy and
  how the energy-dissipation surrogate satisfies it for Phase
  2.fdtd.6 with an explicit follow-up commitment.
- Phase 2.fdtd.6.1 (queued) — calibrated TEM stripline, Γ
  gate, multi-axis / oriented-arbitrary lumped ports, parallel-RLC
  topology.
