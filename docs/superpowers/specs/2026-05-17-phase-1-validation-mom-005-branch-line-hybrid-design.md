# Phase 1 Validation — mom-005 Branch-Line (90°) Hybrid — Design

**Status:** Draft
**Owner:** TBD
**Phase:** 1 (validation case; gates v1.0)
**Depends on:** Phase 1.1.1.0 (multi-image DCIM, shipped at `f9e63c7`), Phase 1.1.1.1 (mesh refinement, in flight), mom-004 (multi-port + lumped-Z plumbing — shared infrastructure)
**Blocks:** ROADMAP Phase 1 sign-off; tight microstrip 4-port validation; in-phase / quadrature combiner work in Phase 1.4.

## Assumption being challenged

`mom-004` exercises a three-port microstrip structure with a single lumped impedance and a single quarter-wave transformer pair on FR-4. The branch-line hybrid (`mom-005`) is the natural next step: **four ports, all-microstrip (no lumped elements), four quarter-wave arms**, and the validation gate now spans both **amplitude balance** and **phase balance** between two output ports.

The branch-line is also a 90° hybrid coupler — its output ports have a 90° phase relationship by design, which is what makes it useful as a quadrature splitter. This means `mom-005` is the first validation case that gates **phase** rather than just magnitude. Any sign-convention drift in `delta_gap_rhs`, `WavePort::rhs`, or the S-matrix extraction will show up as a `0° / 180°` flip on either S₂₁ or S₃₁ relative to S₄₁; an exclusively magnitude-based gate (mom-001 / mom-002 / mom-003 / mom-004) would miss this.

Substrate is again FR-4 (εr = 4.4, h = 1.6 mm) to keep the DCIM coefficients consistent with mom-002 / mom-004. The geometry is a 4-port square ring: two horizontal arms at Z₀ = 50 Ω, length λ_g/4; two vertical arms at Z₀ = 35.36 Ω = 50/√2, length λ_g/4. Trace widths come from Hammerstad-Jensen inverse synthesis (a 35.36 Ω trace on FR-4 ≈ 4.6 mm — wider than the 50 Ω trace, which is ≈ 3.0 mm).

## Scope

In:

- `BranchLineGeometry::new(f_hz, eps_r, h_m)` factory — analogous to `WilkinsonGeometry` from mom-004.
- 4-port S-matrix extraction at `f₀ = 2.0 GHz`.
- Sweep `[1.8, 2.2] GHz` at 21 points, exported to `.s4p` Touchstone.
- Phase-balance gate (in addition to amplitude balance).
- New row in `crates/yee-mom/validation/README.md`.

Out:

- N-section / broadband hybrid topologies (only the canonical 1-section is in this milestone).
- Lange coupler (interdigital coupling structure; Phase 1.4+).
- Mixed-mode S-parameters (Phase 1.4+).
- Radiation effects (treated as lossless dielectric microstrip on a closed ground plane).

## Approach

The branch-line is a 4-port ring built from four λ_g/4 arms:

```
      Port 1 (50 Ω) ───── [Z₀ = 50, λ/4] ───── Port 4 (50 Ω, isolated)
              │                                       │
   [Z₀ = 35.36, λ/4]                       [Z₀ = 35.36, λ/4]
              │                                       │
      Port 2 (50 Ω) ───── [Z₀ = 50, λ/4] ───── Port 3 (50 Ω)
```

Pozar §7.5 gives the closed-form lossless ideal S-matrix for excitation at port 1:

- `S₁₁ = 0` (matched).
- `S₂₁ = −j / √2 = −3 dB ∠−90°` (through port).
- `S₃₁ = −1 / √2 = −3 dB ∠−180°` (coupled port).
- `S₄₁ = 0` (isolated port).
- Plus reciprocity and the obvious symmetries (`S_ji = S_ij`).

The phase-balance gate is the **difference** `arg(S₂₁) − arg(S₃₁) = 90°` (ideal), with the spec tolerance set to `±5°`. This is the load-bearing assertion of the hybrid: a 90° quadrature relationship at the output ports.

Note that no lumped elements appear — the mom-004 lumped-Z plumbing is **not** exercised by mom-005. This separation is intentional: if mom-004 reveals a Schur-reduction bug, mom-005 still gates the multi-port S-matrix path independently.

Wavelength on a 50 Ω microstrip on FR-4 at 2 GHz: `εr_eff(50, FR-4) ≈ 3.27`, so `λ_g ≈ 82.9 mm`, `λ_g/4 ≈ 20.7 mm`. On the 35.36 Ω arm, `εr_eff ≈ 3.46`, `λ_g ≈ 80.6 mm`, `λ_g/4 ≈ 20.2 mm`. The mesh fixture computes both lengths from the actual synthesis, not from these guides.

References:

- Pozar, *Microwave Engineering* 4th ed., §7.5 (90° hybrid / branch-line coupler).
- Reed & Wheeler, "A method of analysis of symmetrical four-port networks", IRE Trans. MTT, 4(4), 1956 (the original even/odd-mode decomposition).
- Mongia, Bahl, Bhartia, *RF and Microwave Coupled-Line Circuits*, 1999, Ch. 8.

## Public API

```rust
/// Branch-line (90°) hybrid validation fixture.
///
/// Phase 1 — gates `mom-005` against Pozar §7.5 at f₀.
pub struct BranchLineGeometry {
    pub f_hz: f64,
    pub eps_r: f64,
    pub h_m: f64,
    /// 50 Ω horizontal-arm trace width (m), from Hammerstad-Jensen synthesis.
    pub w_50: f64,
    /// 35.36 Ω vertical-arm trace width (m), from synthesis.
    pub w_35: f64,
    /// Horizontal-arm length λ_g/4 at f_hz on the 50 Ω trace (m).
    pub arm_len_50: f64,
    /// Vertical-arm length λ_g/4 at f_hz on the 35.36 Ω trace (m).
    pub arm_len_35: f64,
}

impl BranchLineGeometry {
    /// Build with Hammerstad-Jensen synthesis at the given (f, εr, h).
    pub fn new(f_hz: f64, eps_r: f64, h_m: f64) -> Self;

    /// Triangle mesh of the 4-port ring with port tags 1..=4 on the
    /// four feed edges. No lumped elements.
    pub fn mesh(&self) -> TriMesh;
}
```

`yee-validation::run_mom_005`:

```rust
let geom = BranchLineGeometry::new(2.0e9, 4.4, 1.6e-3);
let mesh = geom.mesh();
let greens = GreensSpec::MicrostripDcim { eps_r: 4.4, h_m: 1.6e-3, n_images: 5 };
let result = PlanarMoM::run(mesh, greens, &[p1, p2, p3, p4], &[], freq_sweep)?;
let s = result.s_at(2.0e9);
let phase_diff_deg = (s[(1,0)].arg() - s[(2,0)].arg()).to_degrees();
assert!((phase_diff_deg - 90.0).abs() <= 5.0);
```

## Definition of done

1. `BranchLineGeometry::new` and `::mesh` exist; `cargo build -p yee-mom -p yee-validation` clean.
2. **Validation gate (centre frequency amplitude).** At `f₀ = 2.0 GHz`:
   - `|S₁₁| ≤ −20 dB`.
   - `|S₂₁| ∈ [−3.5, −2.5] dB` (within ±0.5 dB of −3 dB).
   - `|S₃₁| ∈ [−3.5, −2.5] dB` (within ±0.5 dB of −3 dB).
   - `|S₄₁| ≤ −20 dB` (isolation gate — explicit per brief).
   - `||S₂₁| − |S₃₁|| ≤ 0.5 dB` (amplitude balance).
3. **Validation gate (centre frequency phase).** At `f₀ = 2.0 GHz`:
   - `|arg(S₂₁) − arg(S₃₁) − (−90°)| ≤ 5°`, modulo 360°.
4. **Validation gate (band).** Across `[1.8, 2.2] GHz` (21 points):
   - amplitude balance `||S₂₁| − |S₃₁|| ≤ 1 dB` at every swept frequency.
   - phase balance `|Δarg(S₂₁, S₃₁) − (−90°)| ≤ 10°` at every swept frequency.
5. Touchstone export: a `.s4p` file written to `tests/results/branch_line.s4p`; round-trips through `yee_io::touchstone::read` at `1e-12` relative.
6. New row in `crates/yee-mom/validation/README.md`: `mom-005 / branch-line 2 GHz / Pozar §7.5 / ±0.5 dB / ±5° / amp + phase balance`.
7. `cargo doc --no-deps -p yee-mom -p yee-validation` warning-free.
8. mom-001 / mom-002 / mom-004 regression: still green.

## Lane (when implemented)

`crates/yee-mom/**` (geometry / synthesis helpers; no new core types — multi-port S-matrix already lands with mom-004)
+ `crates/yee-validation/**` (fixture + gate)
+ `examples/**` (a `branch_line_hybrid/` example binary).

No edits to `yee-cli`, `yee-gui`, `yee-fdtd`, `yee-cuda`. No new lumped-element plumbing (relies on mom-004's path being merged but does not extend it).

## Verification

```bash
cargo build  -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo test   -p yee-validation --release run_mom_005
cargo fmt    --check --all
```

mom-001 / mom-002 / mom-004 must remain green.

## Escape hatch

The **phase-balance gate is new** with this milestone — no prior validation case exercises an absolute phase angle. Two failure modes are likely:

1. **Global phase reference drift.** If the S-matrix extraction has a port-edge midpoint convention that shifts globally with the mesh resolution, the absolute phase of `S₂₁` can drift while the **difference** `arg(S₂₁) − arg(S₃₁)` stays correct. In that case, weaken the gate from "absolute `arg(S₂₁) = −90°`" to "**difference** `arg(S₂₁) − arg(S₃₁) = −90°`", which is what the brief actually asks for and what the physics actually means. This is **not** a relaxation — it's the right statement of the gate.
2. **Sign-convention bug in `delta_gap_rhs`.** If a sign flip puts the output ports 180° out from the reference, the gate will fail by exactly 180°. Surface immediately — this is a real bug, not a tolerance issue, and silencing it by widening the gate would mask a real defect on every future multi-port case.

Blocked > 15 min on phase plumbing → surface and stop with both candidate diagnoses above.

## References

- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §7.5 (90° hybrid / branch-line coupler).
- Reed, J. & Wheeler, G. J., "A method of analysis of symmetrical four-port networks", IRE Trans. Microwave Theory Tech., 4(4), 1956, pp. 246–252.
- Mongia, R. K., Bahl, I. J., Bhartia, P., *RF and Microwave Coupled-Line Circuits*, Artech House, 1999, Ch. 8.
- Wadell, B. C., *Transmission Line Design Handbook*, Artech House, 1991, §3.5 (FR-4 dispersion).
