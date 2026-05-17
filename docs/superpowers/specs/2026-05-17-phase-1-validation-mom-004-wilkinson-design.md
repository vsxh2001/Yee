# Phase 1 Validation — mom-004 Wilkinson Power Divider — Design

**Status:** Draft
**Owner:** TBD
**Phase:** 1 (validation case; gates v1.0)
**Depends on:** Phase 1.1.1.0 (multi-image DCIM, shipped at `f9e63c7`), Phase 1.1.1.1 (mesh refinement, in flight), Phase 1.3.1.0 (TE10 wave port, shipped)
**Blocks:** ROADMAP Phase 1 sign-off; tight microstrip 3-port validation; Phase 1.3.2 lumped-port milestone (see escape hatch).

## Assumption being challenged

The shipped `MultilayerGreens` (Phase 1.1.1.0, N=5 multi-image DCIM) lets `mom-002` recover the Hammerstad-Jensen 50 Ω characteristic impedance on a 30 mm microstrip line to within roughly ±10%. That is **a single-port** validation: one feed, one strip, one Z₀. It does not exercise multi-port S-parameter extraction, port-to-port isolation, or the lumped-impedance machinery that any real PCB design needs.

`mom-004` asks for the next obvious step on the same substrate stack: **a 2 GHz Wilkinson power divider on FR-4** (εr = 4.4, h = 1.6 mm). Geometry is canonical (Pozar 4th ed. §7.3), the closed-form answer for a lossless, matched, equal-split design is `S₁₁ = 0`, `S₂₁ = S₃₁ = −3 dB ∠−90°`, `S₂₃ = 0` (isolation), and `S₂₂ = S₃₃ = 0`. This is a small enough mesh (a Y-junction of three 50 Ω lines plus two λ/4 70.7 Ω arms) to fit comfortably in MoM-sized DoF counts, and it forces every code path that we eventually need for production designs:

- **Three-port S-matrix extraction** (mom-001/002 are one-port; mom-003 is one-port).
- **A 100 Ω isolation resistor between output ports.** This is the load-bearing element of the Wilkinson topology and is **not** a transmission line — it is a lumped impedance bridging two nodes.
- **Quarter-wave transformer arms at 70.7 Ω.** Width synthesis (Hammerstad-Jensen inverse) on FR-4 at 2 GHz gives `w ≈ 1.8 mm` (50 Ω trace ≈ 3.0 mm; 70.7 Ω trace ≈ 1.5 mm — these are guides, the validation case re-runs the synthesis at test-fixture-build time).

Substrate is kept identical to `mom-002` (FR-4, εr = 4.4, h = 1.6 mm) on purpose: the DCIM coefficients are already validated there, so any failure on `mom-004` localises to the multi-port / lumped-Z plumbing, not the Green's function.

## Scope

In:

- `WilkinsonGeometry::new(f_hz, eps_r, h_m)` factory that synthesises trace widths from Hammerstad-Jensen and builds a port-tagged 2-D triangle mesh of the Y-junction + two arms + lumped-resistor footprint.
- Three-port S-matrix extraction at the centre frequency `f₀ = 2.0 GHz`.
- Sweep `[1.5, 2.5] GHz` at 21 points, exported to `.s3p` Touchstone.
- Validation gate: `|S| values at f₀ within ±0.5 dB` of Pozar closed form; full sweep within `±1 dB` band.
- New row in `crates/yee-mom/validation/README.md`.

Out:

- Tight tolerances (the Wilkinson is **dispersive** on FR-4 — bandwidth narrows on lossy substrate; ±0.1 dB needs lumped lossy-resistor + dispersive-Green's, both Phase 1.1.2+).
- Differential / mixed-mode S-parameters.
- Radiation / spurious-mode capture (the substrate is treated as a lossless dielectric; the εr = 4.4 of FR-4 is a real value here).
- Mechanical / thermal effects.
- N-way Wilkinson (only the 2-way is in this milestone; N-way is Phase 1.4+).

## Approach

The Wilkinson is built as three concatenated microstrip sections plus one lumped load:

```
   Port 1 (50 Ω)
       │
   ┌───┴───┐
   │       │     <- Y-junction, mesh fan-out
   ▼       ▼
 70.7 Ω  70.7 Ω    <- two λ/4 arms, length ≈ 19.2 mm on FR-4 at 2 GHz
 (arm A) (arm B)
   │       │
   ▼       ▼
 Port 2  Port 3   <- both 50 Ω
   │       │
   └──[R=100 Ω]──┘  <- lumped isolation resistor across nodes (2,3)
```

The wavelength on a 70.7 Ω microstrip on FR-4 at 2 GHz is computed from the effective dielectric constant `εr_eff(70.7, FR-4) ≈ 3.05` (Hammerstad-Jensen): `λ_g = c / (f √εr_eff) ≈ 85.7 mm`, so `λ_g/4 ≈ 21.4 mm`. The 19 mm figure in the brief is a rough guide; the fixture computes the exact value from synthesis.

The 100 Ω resistor is a **lumped two-terminal port** spanning the two output nodes. Phase 1.3 shipped delta-gap and TE10 wave-ports as full ports of the network; the resistor here is internal — a load, not a measurement port. The implementation routes it as a `LumpedZ` element in `solve.rs`, post-processed into the impedance matrix as a Schur complement before the three-port S-matrix is computed. See escape hatch below if this surfaces a Phase 1.3.2 prereq.

Once the impedance matrix `Z` is filled, the 3-port S-matrix follows from:

```
   Z_port = (selection of rows/cols of Z corresponding to port basis fns, after isolation-resistor Schur reduction)
   S = (Z_port − Z₀ I) (Z_port + Z₀ I)⁻¹    with Z₀ = 50 Ω
```

References:

- Pozar, *Microwave Engineering* 4th ed., §7.3 (Wilkinson power divider closed form).
- Hammerstad & Jensen, "Accurate models for microstrip computer-aided design", IEEE MTT-S 1980 (width / εr_eff synthesis).
- Wadell, *Transmission Line Design Handbook*, 1991, §3.5 (FR-4 dispersion notes).

## Public API

```rust
/// Wilkinson divider validation fixture.
///
/// Phase 1 — gates `mom-004` against Pozar §7.3 at f₀.
pub struct WilkinsonGeometry {
    pub f_hz: f64,
    pub eps_r: f64,
    pub h_m: f64,
    /// 50 Ω input trace width (m), from Hammerstad-Jensen inverse synthesis.
    pub w_50: f64,
    /// 70.7 Ω arm trace width (m), from inverse synthesis.
    pub w_707: f64,
    /// Arm length λ_g/4 at f_hz (m).
    pub arm_len: f64,
    /// Isolation-resistor value (Ω). Always 100 Ω for the canonical Wilkinson.
    pub r_iso: f64,
}

impl WilkinsonGeometry {
    /// Build with Hammerstad-Jensen synthesis at the given (f, εr, h).
    pub fn new(f_hz: f64, eps_r: f64, h_m: f64) -> Self;

    /// Triangle mesh of the entire structure (Y-junction + arms),
    /// with port tags 1, 2, 3 on the three feed edges and a lumped-Z
    /// element tagged across the output nodes.
    pub fn mesh(&self) -> (TriMesh, LumpedZSpec);
}

/// Lumped impedance bridging two node ids.
pub struct LumpedZSpec {
    pub node_a: u32,
    pub node_b: u32,
    pub z_ohms: Complex64,
}
```

`yee-validation::run_mom_004` wires this into the existing `PlanarMoM` driver:

```rust
let geom = WilkinsonGeometry::new(2.0e9, 4.4, 1.6e-3);
let (mesh, r_iso) = geom.mesh();
let greens = GreensSpec::MicrostripDcim { eps_r: 4.4, h_m: 1.6e-3, n_images: 5 };
let result = PlanarMoM::run(mesh, greens, &[port1, port2, port3], &[r_iso], freq_sweep)?;
let s_at_f0 = result.s_at(2.0e9);
// gate: |S11| < -20 dB, |S21| within [-3.5, -2.5] dB, etc.
```

## Definition of done

1. `WilkinsonGeometry::new` and `::mesh` exist; `cargo build -p yee-mom` clean.
2. `LumpedZSpec` plumbed through `PlanarMoM::run` via a Schur-reduction path in `solve.rs`. The free-space mom-001 path is bit-for-bit unchanged (no `LumpedZSpec` → no Schur reduction).
3. **Validation gate (centre frequency).** At `f₀ = 2.0 GHz`:
   - `|S₁₁| ≤ −20 dB` (i.e. within ±0.5 dB of the ideal `−∞ dB`; finite floor allowed).
   - `|S₂₁| ∈ [−3.5, −2.5] dB` (within ±0.5 dB of −3 dB).
   - `|S₃₁| ∈ [−3.5, −2.5] dB` (symmetry).
   - `|S₂₃| ≤ −20 dB` (within ±0.5 dB of the ideal `−∞ dB` isolation).
   - `||S₂₁| − |S₃₁|| ≤ 0.1 dB` (equal split).
4. **Validation gate (band).** Across `[1.5, 2.5] GHz` (21 points): every S-parameter within `±1 dB` of the Pozar closed-form curve. The closed-form curve is computed in the test via a transmission-line ABCD-cascade; the gate is `max_abs_err ≤ 1.0 dB`.
5. Touchstone export: a `.s3p` file written to `tests/results/wilkinson.s3p`; round-trips through `yee_io::touchstone::read` at `1e-12` relative.
6. New row in `crates/yee-mom/validation/README.md`: `mom-004 / Wilkinson 2 GHz / Pozar §7.3 / ±0.5 dB f₀, ±1 dB band / multi-port + lumped-Z`.
7. `cargo doc --no-deps -p yee-mom -p yee-validation` warning-free.
8. mom-001 and mom-002 regression: still green.

## Lane (when implemented)

`crates/yee-mom/**` (LumpedZ plumbing in `solve.rs`, Schur reduction)
+ `crates/yee-validation/**` (fixture + gate)
+ `examples/**` (a `wilkinson_divider/` example binary mirroring `examples/microstrip_line/`).

No edits to `yee-cli`, `yee-gui`, `yee-fdtd`, `yee-cuda`.

## Verification

```bash
cargo build  -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo test   -p yee-validation --release run_mom_004
cargo fmt    --check --all
```

`mom-001` (`dipole_z_at_resonance`) and `mom-002` must remain green; this validation case is additive.

## Escape hatch

The **lumped 100 Ω isolation resistor** is the single non-canonical element. Phase 1.3 shipped two `Port` variants — delta-gap and TE10 wave-port — both of which are *external* ports (excitation + measurement). The isolation resistor is *internal*: it is neither an excitation nor a measurement, it is a load between two nodes that participates in the impedance matrix as a lumped Z spanning two RWG basis indices.

If the Schur-reduction implementation lands cleanly inside `solve.rs` (an N×N system + a K×K lumped block, reducing to an `(N-K)×(N-K)` effective system), this spec stands. **If it surfaces a wider need** (e.g. lumped capacitors / inductors, frequency-dependent lumped Z, lumped ports as a first-class `Port` variant), that's a new milestone:

- **Phase 1.3.2 — Lumped ports / loads.** A first-class `LumpedPort { z: f64 + j × ω L − j / ω C, nodes: (u32, u32) }` variant on the `Port` trait, with proper Schur reduction or augmented-system handling.

In that case, surface as a finding, complete mom-004 with the inline ad-hoc Schur reduction on this milestone, and open Phase 1.3.2 as a follow-up. **Do not** widen the lane to add a `LumpedPort` to `yee-mom/src/ports.rs` inside the mom-004 PR — that crosses the design boundary and would block the validation milestone on a larger refactor.

Blocked > 15 min on lumped-Z plumbing → surface and stop.

## References

- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §7.3 (Wilkinson power divider).
- Hammerstad, E. & Jensen, Ø., "Accurate models for microstrip computer-aided design", IEEE MTT-S Int. Microwave Symp. Digest, 1980, pp. 407–409.
- Wadell, B. C., *Transmission Line Design Handbook*, Artech House, 1991, §3.5.
- Wilkinson, E. J., "An N-way hybrid power divider", IRE Trans. Microwave Theory Tech., 8(1), 1960, pp. 116–118 (the original).
