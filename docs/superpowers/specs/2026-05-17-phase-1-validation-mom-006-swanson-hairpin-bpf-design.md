# Phase 1 Validation — mom-006 Swanson 5-Pole Hairpin BPF — Design

**Status:** Draft
**Owner:** TBD
**Phase:** 1 (validation case; gates v1.0)
**Depends on:** Phase 1.1.1.0 (multi-image DCIM, shipped at `f9e63c7`), **Phase 1.1.1 complete** (full Sommerfeld extraction + surface-wave poles — see escape hatch), Phase 1.1.1.1 (mesh refinement, in flight), mom-005 (multi-port S-matrix path).
**Blocks:** ROADMAP Phase 1 sign-off (v1.0 cannot ship without mom-006); EM-CAD-grade filter validation; foundation for Phase 2 distributed-element synthesis tooling.

## Assumption being challenged

`mom-004` and `mom-005` are two- and four-port microwave-junction validation cases: they exercise the multi-port S-matrix path on simple topologies with two or four quarter-wave sections and at most one lumped impedance. **The Swanson 5-pole hairpin band-pass filter is qualitatively harder by ~1.5 orders of magnitude on every axis:**

- **Mesh size.** Each hairpin is a folded half-wavelength resonator with ~10 mm of trace per side. Five hairpins, plus tap-feed lines and inter-resonator gaps, mesh to roughly **5,000–10,000 RWG basis functions** at a useful refinement (versus ~500–1500 for mom-004 / mom-005). The complex-double impedance matrix at 10k DoFs is ~1.6 GB — this is the **first validation case that exercises the cuSOLVER GPU path** for real (CPU LU on faer would take ~30 minutes per frequency).
- **Substrate.** RT/Duroid 6006: `εr = 6.15` (vs 4.4 for FR-4) and `h = 1.27 mm` (vs 1.6 mm). The DCIM N=5 coefficients fitted for FR-4 do **not** transfer — they must be re-fit at the new (εr, h). This is the first validation case that exercises DCIM with **non-FR-4 substrate**, which is a precondition for any production use of the Green's function.
- **Inter-resonator coupling.** Hairpin resonators are coupled through their narrow-gap fringing fields. The strength of that coupling sets the filter bandwidth (~2% in the brief: passband `[1.97, 2.03] GHz`). Getting the resonant frequencies within ±0.5% means the **inter-resonator gap capacitance** must be predicted within roughly ±2-3%, which is at the edge of what the N=5 DCIM can do.
- **Surface waves.** RT/Duroid 6006 with εr = 6.15 and 1.27 mm substrate at 2 GHz has its TM0 surface-wave cutoff well below 2 GHz; the dominant TE1 surface-wave cutoff sits roughly above 30 GHz, so we are not exciting higher-order surface modes, but the TM0 contribution to inter-resonator coupling is non-negligible and **is the dominant error source if DCIM alone is used**. This is the load-bearing reason this spec is gated on Phase 1.1.1 (full Sommerfeld) being complete, not just 1.1.1.1.

This is the most aggressive validation case in Phase 1 and the one that demonstrates that Yee can do **EM-CAD-grade filter validation**, not just textbook junction analysis.

## Scope

In:

- `SwansonHairpinGeometry::new()` factory that produces a 5-pole hairpin BPF on RT/Duroid 6006 with center frequency ≈ 2 GHz, ≈ 2% bandwidth. The geometry is taken from Swanson's published design (see references — citation status surfaced below).
- Two-port S-matrix extraction at every swept frequency over `[1.5, 4.0] GHz`.
- Comparison to Swanson's published S-parameters / Sonnet reference within the documented tolerances.
- New row in `crates/yee-mom/validation/README.md`.
- A new runnable example binary `examples/swanson_hairpin/`.

Out:

- N-pole generalisation (only the canonical 5-pole structure is in this milestone; N-pole tooling is Phase 2+).
- Synthesis of the filter from spec (start with the published dimensions; synthesis is a Phase 2.cad sub-project).
- Loss budgeting / Q-factor extraction (this milestone is a lossless-substrate gate; lossy / dispersive substrates are Phase 1.1.2+).
- Mixed-mode or differential S-parameters.

## Approach

Five hairpin resonators side-by-side, fed by tap-feeds on the first and last hairpins:

```
   Port 1 ───tap─┐                                              ┌─tap─── Port 2
                 │                                              │
                 ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐ ┌──────┐  │
                 │      │ │      │ │      │ │      │ │      │  │
                 │  R1  │ │  R2  │ │  R3  │ │  R4  │ │  R5  │  │
                 │      │ │      │ │      │ │      │ │      │  │
                 └──────┘ └──────┘ └──────┘ └──────┘ └──────┘  │
                       gap     gap     gap     gap             │
```

Each hairpin is a half-wavelength resonator folded into a "U" shape; on εr = 6.15, h = 1.27 mm at 2 GHz, `λ_g/2 ≈ 36 mm`. The folded length is ≈ 18 mm with arm width chosen for ≈ 60 Ω characteristic impedance on the substrate. Inter-resonator gaps are computed from the filter prototype (Chebyshev or 0.01 dB ripple, depending on the published design) to give the ~2% bandwidth.

The mesh refinement that mom-006 needs is the one being shipped by Phase 1.1.1.1 (in flight at brief time). Adaptive refinement near gaps and tap-feed junctions is the load-bearing capability — uniform meshing at the resolution needed in the gap regions would give roughly 50k DoFs, well beyond the cuSOLVER capacity on a 16 GB GPU.

The cuSOLVER GPU LU path is Phase 1.5 (shipped); mom-006 is its first validation customer outside the synthetic GPU tests.

References:

- **Swanson, "Narrow-Band Microwave Filter Design", IEEE Microwave Magazine, 8(5), October 2007, pp. 105–114.** This is the **most likely** primary source — Daniel Swanson published a series of articles in IEEE Microwave Magazine in the mid-2000s on practical narrow-band filter design with explicit hairpin BPF examples and Sonnet-validated S-parameters. **STATUS — TBD: this citation must be verified before implementation begins.** The brief mentions a possible IEEE MTT publication 1992–2000; that may instead be Swanson's earlier work on EM-CAD validation of microstrip filters (his early IEEE MTT papers focus on direct-coupled-cavity filters and EM simulator validation). The hairpin-specific example may live in the 2007 Microwave Magazine article, in Swanson & Hofer, *Microwave Circuit Modeling Using Electromagnetic Field Simulation*, Artech House 2003 (book, Ch. 6), or in a different paper entirely. **Before any code lands**, the implementor must pull the actual published dimensions and S-parameters and pin the reference in `crates/yee-mom/validation/README.md`. Do not fabricate dimensions; if no source is found, surface and stop per the escape hatch below.
- Hong & Lancaster, *Microstrip Filters for RF/Microwave Applications*, Wiley 2001, §5.3 (hairpin resonator filter design — generic published procedure that produces equivalent geometries; can be a fallback reference if the original Swanson source cannot be located but is **not** what the brief asks for).
- Pozar, *Microwave Engineering* 4th ed., §8.3–8.4 (filter prototype synthesis).
- Sonnet Software documentation, "Hairpin filter on RT/Duroid 6006 — sample validation" (the brief's reference for the ±1 dB up-to-4-GHz tolerance band).

## Public API

```rust
/// Swanson 5-pole hairpin band-pass filter validation fixture.
///
/// Phase 1 — gates `mom-006` against Swanson's published S-parameters.
pub struct SwansonHairpinGeometry {
    /// RT/Duroid 6006: 6.15.
    pub eps_r: f64,
    /// RT/Duroid 6006: 1.27 mm.
    pub h_m: f64,
    /// Per-hairpin folded length (m). 5 entries.
    pub hairpin_len_m: [f64; 5],
    /// Per-hairpin arm width (m). 5 entries.
    pub hairpin_w_m: [f64; 5],
    /// Inter-resonator gap (m). 4 entries (between adjacent hairpins).
    pub gap_m: [f64; 4],
    /// Tap-feed position on the first / last hairpin (m from the bend).
    pub tap_offset_m: f64,
}

impl SwansonHairpinGeometry {
    /// Build with the published Swanson dimensions for ~2 GHz, ~2% BW on RT/Duroid 6006.
    ///
    /// **TODO at implementation time: pin published-source dimensions here.**
    /// Until then, the constructor returns dimensions from Hong & Lancaster
    /// §5.3 worked example with a deliberate `eprintln!("warning: ...")`.
    pub fn published() -> Self;

    /// Triangle mesh of all 5 hairpins + 2 tap-feeds, with port tags 1, 2.
    pub fn mesh(&self) -> TriMesh;
}
```

## Definition of done

1. `SwansonHairpinGeometry::published` and `::mesh` exist; dimensions are cited inline against a verified published source (escape hatch below if not).
2. **Validation gate (passband).** Across `[1.97, 2.03] GHz` (61 points, 1 MHz spacing):
   - `|S₂₁|` within `±1 dB` of the reference at every point.
   - return loss `|S₁₁| ≤ −20 dB` at every point (≥ 20 dB RL).
3. **Validation gate (pole frequencies).** The five in-band poles of `|S₁₁|` (local minima of the return loss) appear at frequencies within `±0.5%` of the reference values.
4. **Validation gate (wideband).** Across `[1.5, 4.0] GHz` (251 points, 10 MHz spacing):
   - `|S₂₁|` within `±1 dB` of the reference at every point.
5. Touchstone export: a `.s2p` file written to `tests/results/swanson_hairpin.s2p`; round-trips through `yee_io::touchstone::read` at `1e-12` relative.
6. The runnable example `examples/swanson_hairpin/` produces a PNG plot of `|S₂₁|` and `|S₁₁|` over `[1.5, 4.0] GHz` via `yee-plotters`.
7. New row in `crates/yee-mom/validation/README.md`: `mom-006 / Swanson 5-pole hairpin BPF / Swanson published S-params (cite once pinned) / ±1 dB passband / ±0.5% poles / ≥20 dB RL`.
8. `cargo doc --no-deps -p yee-mom -p yee-validation` warning-free.
9. mom-001 / mom-002 / mom-004 / mom-005 regression: still green.

## Lane (when implemented)

`crates/yee-mom/**` (no new core types — multi-port and lumped-Z paths already landed for mom-004 / mom-005)
+ `crates/yee-validation/**` (fixture + gate)
+ `examples/swanson_hairpin/**` (new example binary + plot).

Likely also touches `crates/yee-mom/Cargo.toml` to enable the `cuda` feature for the validation case in release builds, but this is **optional** — a sufficiently patient CI run on faer is acceptable if cuSOLVER is unavailable. The brief's escape hatch addresses the GPU dependency.

No edits to `yee-cli`, `yee-gui`, `yee-fdtd`.

## Verification

```bash
cargo build  -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo test   -p yee-validation --release run_mom_006   # may take 15–60 min on CPU; <5 min on GPU
cargo run    -p yee-validation --release --example swanson_hairpin
cargo fmt    --check --all
```

mom-001 / mom-002 / mom-004 / mom-005 must remain green.

## Escape hatch

**This is the most aggressive case in Phase 1.** Three layered escape paths, in order of severity:

1. **Citation not found.** If the implementor cannot pin Swanson's original published hairpin dimensions and reference S-parameters from a primary source (IEEE Microwave Magazine, IEEE T-MTT, or the Swanson & Hofer 2003 book) within **25 minutes** of literature search, surface and stop. Do **not** fabricate dimensions. Fallback option: cite Hong & Lancaster §5.3 worked example, mark `mom-006` as "validated against textbook hairpin, not against Swanson published" in `validation/README.md`, and surface as a finding requiring follow-up. This weakens the validation claim but ships the case.

2. **DCIM N=5 is insufficient.** If the run produces resonant frequencies off by ≥ 1% (i.e. > 2× the tolerance) and the residuals correlate with TM0 surface-wave coupling between hairpins, the **full Sommerfeld-integral Green's function (Phase 1.1.1 complete, including surface-wave pole extraction in Phase 1.1.1.2) is a hard prerequisite**. This is explicit in this spec's `Depends on` line: **mom-006 is gated on Phase 1.1.1 being complete, not just 1.1.1.1**. Surface and block the milestone until 1.1.1.2 lands; do not weaken the tolerances.

3. **Mesh size exceeds GPU memory.** If the refined mesh produces > ~12k RWG DoFs (impedance matrix > 2.3 GB), a 16 GB GPU cannot LU-factorize it without out-of-core support. Two paths: (a) refine adaptively only near gaps and accept slightly looser pole tolerances on the wing resonators; (b) fall back to CPU faer LU and accept the multi-hour CI wall-time, gated behind a `--features slow-validation` cargo flag so CI doesn't run it by default. Option (a) is preferred — surface the choice and stop for direction if neither produces a green gate.

Blocked > 25 min on any of the three → surface and stop.

## References

- **TBD: Swanson, D. G., the canonical hairpin-on-Duroid-6006 publication. Most likely IEEE Microwave Magazine 8(5), 2007 (Narrow-Band Microwave Filter Design); possibly Swanson & Hofer, *Microwave Circuit Modeling Using Electromagnetic Field Simulation*, Artech House 2003, Ch. 6. The brief mentions 1992–2000 IEEE MTT; that may instead be Swanson's earlier EM-CAD validation papers — those are about the validation methodology rather than this specific filter. Pin before implementation begins.**
- Hong, J.-S. & Lancaster, M. J., *Microstrip Filters for RF/Microwave Applications*, Wiley 2001, §5.3 (hairpin filter design procedure — fallback reference).
- Pozar, D. M., *Microwave Engineering*, 4th ed., Wiley 2012, §8.3–8.4 (filter prototype tables).
- Rogers Corporation, "RT/duroid 6006 Laminate Data Sheet" (substrate properties: εr = 6.15 ± 0.15, tan δ ≈ 0.0027 at 10 GHz; thickness 1.27 mm).
- Sonnet Software User's Guide, hairpin filter validation example (reference for the ±1 dB up-to-4-GHz Sonnet comparison band).
