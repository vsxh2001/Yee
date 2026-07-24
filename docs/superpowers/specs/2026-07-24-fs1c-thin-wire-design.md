# FS.1c — thin-wire subcell (Holland–Simpson) + dipole gate vs NEC-4

**Date:** 2026-07-24 · **Track:** FS.1 (FULL-SUITE-ROADMAP §3) · **Lane:** `crates/yee-compute/**`, `crates/yee-engine/**` (+ docs)
**Reference in-repo:** mom-001 — NEC-4 finite-radius half-wave dipole `Z ≈ 87 + j41 Ω`
(L = 1 m, a = 5 mm; CLAUDE.md §4 — quote NEC-4 only, never the Balanis 73 + j42 wire limit).

## Why

FS.1's antenna catalog covers planar topologies; wire antennas (dipoles, monopoles,
wire feeds) need a wire much thinner than a cell. Naive one-cell PEC wires give a
radius-of-half-a-cell artifact. The standard cure is the **Holland–Simpson thin-wire
subcell model** (Holland & Simpson, IEEE Trans. EMC 1981; also Taflove §10.3 /
gprMax's implementation): the wire's in-cell inductance per unit length
`L = (μ₀/2π)·ln(dx/a)`-class correction modifies the E-update on wire-axis edges and
zeroes the radial E at the wire. Implement the simplest published variant that
supports a z-directed straight wire with a delta-gap source cell.

## Deliverables

1. **`yee-compute` thin-wire support (CPU)**: a `ThinWire { axis: z-only for now,
   i, j, k_lo, k_hi, radius_m, feed_k: Option<usize> }` attached via `Drive` or
   `Materials` (implementer picks the least-churn seam and documents why; the mask
   plumbing (PEC cells) and drive plumbing are both candidates). Update rule per
   Holland–Simpson: on wire edges the Ez update uses the in-cell inductance
   correction; radial E components at the wire are shorted per the model. Cite the
   exact equation source (Taflove eq. numbers or gprMax docs) in the module docs —
   research-first, don't invent a scheme. **GPU: named `Unsupported` rejection with
   a test** (walking skeleton; GPU port is a follow-on).
2. **Unit gate**: a coarse-vs-fine consistency check — the same physical wire
   (fixed L, a) at two grid resolutions gives resonant frequencies within a few %
   (the subcell model's whole point is grid-independence; a naive one-cell wire
   fails this badly — demonstrate the naive failure as the negative control if
   cheap).
3. **Gate `engine-thinwire-dipole-001`** (yee-engine): the mom-001 dipole (L = 1 m,
   a = 5 mm) in free space (open boundaries, CPML all faces), delta-gap fed at
   centre; measure input impedance at resonance from the feed-cell V/I (port
   records idiom, FS.2a) over a band around 143 MHz; compare Re/Im to NEC-4
   87 + j41 Ω. Tolerances: this is FDTD-vs-MoM-vs-NEC-4 — target ±10 % on Re(Z),
   ±20 % on Im(Z) (looser than mom-001's own gate; a subcell wire in a coarse grid
   is not a 176-segment MoM cylinder). Pin measured + margin; > 25 % on Re → STOP
   and root-cause (feed model, gap capacitance, CPML proximity), never widen.
   ALSO gate the resonant frequency (|Z| minimum / Im(Z) zero crossing) within
   ±5 % of the mom-001-observed resonance — frequency is the robust observable.
4. **ADR-0228** + FS.1 roadmap row (FS.1 then COMPLETE).

## Constraints

- Bit-exact + parity suites unmodified/green; a job with NO thin wires must be a
  provable no-op (existing results bit-identical — assert cheaply).
- Research-first: the implementer reads the published formulation before coding
  (gprMax docs/Taflove §10.3 summaries are acceptable sources; cite what was used).
- Runtime: free-space dipole at λ/20 ≈ 143 MHz is a big box — keep the grid modest
  (coarse λ/20, wire subcell is exactly what makes coarse OK); budget ≤ ~3 min
  release; `#[ignore]` + blanket CI pickup.

## Non-goals

Arbitrary-orientation/bent wires, wire junctions, monopole ground planes, GPU
kernel, loaded/insulated wires, NTFF pattern gate (S11/impedance only this
increment).
