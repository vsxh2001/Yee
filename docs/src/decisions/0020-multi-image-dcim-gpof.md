# ADR-0020: Multi-image DCIM via GPOF (Hua-Sarkar matrix pencil) for `MultilayerGreens`

**Status:** Accepted
**Date:** 2026-05-17
**Deciders:** Yee maintainers

## Context

`MultilayerGreens` shipped in Phase 1.1.0 as a one-image Discrete
Complex Image Method (DCIM) placeholder (CLAUDE.md §10 documents the
status). ADR-0015 then added the `GreensSpec` builder so `mom-002`
(microstrip Z₀) could actually route through it, but routing alone
exposed the limitation: a single complex image is too coarse a model
of a grounded dielectric slab's spectral reflection coefficient to
produce input impedance within an order of magnitude of the
Hammerstad–Jensen 50 Ω target. Phase 1.1.1 is the multi-step plan to
close that gap.

Phase 1.1.1 decomposes along orthogonal axes:

- **Phase 1.1.1.0 — multi-image DCIM.** Replace the single complex
  image with `N` images fitted to the spectral reflection coefficient
  via a numerical method. This ADR.
- **Phase 1.1.1.1 — edge-clustered strip mesh.** Resolve the
  `1/√d` RWG-current singularity at strip edges (ADR-0024).
- **Phase 1.1.1.2 — Sommerfeld surface-wave pole extraction.**
  Subtract discrete poles before fitting the smooth remainder
  (ADR-0025).

The choice for **this** ADR is the numerical method for fitting
`N`-image complex exponentials to a tabulated spectral kernel. The
candidate methods are well-established in the EM literature:

1. **Prony's method.** The classical approach. Numerically
   ill-conditioned on noisy or near-degenerate spectra; requires
   careful regularisation. Less common in modern DCIM practice.
2. **Generalised Pencil-of-Function (GPOF) / Hua–Sarkar
   matrix-pencil.** A closed-form linear-algebra pipeline: SVD of
   the Hankel matrix, complex eigendecomposition of the reduced
   pencil, Vandermonde least-squares for the amplitudes. No
   iteration. Hua & Sarkar 1989; Aksun 1996 popularised it for
   DCIM. Numerically equivalent to Prony in exact arithmetic, much
   more stable in floating point.
3. **Rational-function fitting (vector fitting and friends).** A
   different problem class — fits rational functions to frequency-
   domain data — and a less natural match for the
   complex-exponential ansatz that DCIM requires. The image-domain
   interpretation of rational poles is awkward.

GPOF wins on three axes simultaneously: it is the de-facto standard
in modern DCIM literature (Aksun 1996, Michalski–Mosig 1997), it
is a pure linear-algebra pipeline (no iteration, no convergence
tolerance to tune), and it is numerically equivalent to Prony in
exact arithmetic so the algorithm choice does not change the
mathematical content of the fit.

The contour for sampling the spectral reflection coefficient is the
**Aksun 1996 deformed contour** `k_{z0}(t) = k_0 (1 − j t)`, which
moves the integration path off the Sommerfeld branch cuts and away
from the surface-wave poles enough that GPOF can fit the smooth
remainder. The poles themselves are *not* removed by 1.1.1.0 —
that is 1.1.1.2's job (ADR-0025). 1.1.1.0 ships the fitting
machinery; 1.1.1.2 ships the pole-subtracted variant.

## Decision

`yee-mom` Phase 1.1.1.0 extends `MultilayerGreens` from one-image
to `N`-image complex-image approximation, with `N ∈ [1, 10]`. The
fitting method is **GPOF (Hua–Sarkar matrix pencil)**.

**New module:** `crates/yee-mom/src/gpof.rs` implements GPOF as a
closed-form linear-algebra pipeline:

1. Build the Hankel matrix `H` from the sampled spectral kernel.
2. Take the SVD `H = U Σ V*`, truncate to rank `N`.
3. Form the reduced pencil from the truncated singular vectors;
   take its complex eigendecomposition. The eigenvalues map to the
   image exponents `(α_n)`.
4. Solve a Vandermonde least-squares system for the amplitudes
   `(β_n)`.
5. Map `(α_n, β_n)` to the image `(b_n, a_n)` pairs via the
   Sommerfeld identity.

No iteration loop, no convergence tolerance. The numerical floor is
set by the SVD truncation, not by iteration count.

**New `MultilayerGreens` constructor:**

```rust
impl MultilayerGreens {
    pub fn new_microstrip_with_n_images(
        freq: f64,
        eps_r: f64,
        h_m: f64,
        n_images: usize,   // [1, 10]
    ) -> Self;
}
```

The `n_images = 1` path **preserves the Phase 1.1.0 placeholder
bit-for-bit** (a back-compat tripwire test enforces this). `n_images
> 1` samples the grounded-slab TE/TM reflection coefficients on the
Aksun contour, fits with GPOF, and produces the multi-image kernel.

**New `GreensSpec` variant:**

```rust
pub enum GreensSpec {
    FreeSpace,
    Microstrip { eps_r: f64, h_m: f64 },                    // Phase 1.1.0
    MicrostripDcim { eps_r: f64, h_m: f64, n_images: usize }, // Phase 1.1.1.0
}
```

with a `GreensSpec::microstrip_dcim(eps_r, h_m, n_images)`
convenience constructor mirroring `microstrip()`. The
`PlanarMoM::run` dispatch in ADR-0015 gets one new arm.

**Validation.** A synthetic three-image GPOF recovery test in
`gpof.rs` recovers the planted exponents and amplitudes to within
`1e-6` relative error. The `n_images = 1` tripwire confirms
bit-for-bit identity with the Phase 1.1.0 `vector_images` /
`scalar_images` entries. `mom-002` is routed through `n_images = 5`
and improves from |Z_in| ≈ 14 kΩ (one image) to |Z_in| ≈ 2.7 kΩ
(five images) on the same mesh — a ~5× improvement but still well
outside the Hammerstad–Jensen `[35, 75] Ω` target. **The
validation tolerance is therefore held at the prior loose
`[1, 100 kΩ]` non-degeneracy band**; the kΩ-vs-kΩ improvement is
captured in the case notes, not in the gate.

## Alternatives considered

1. **Prony's method.** Rejected. Numerically equivalent to GPOF in
   exact arithmetic but materially less stable on noisy or
   near-degenerate spectra. GPOF's SVD truncation is a built-in
   regulariser; Prony needs a hand-rolled one. The modern DCIM
   literature has converged on GPOF for this reason.
2. **Rational-function fitting (vector fitting).** Rejected as a
   conceptual mismatch: DCIM's ansatz is complex exponentials in
   the image domain, not rational functions in the frequency
   domain. Forcing the data through a rational-function fit and
   then converting to exponentials adds a conversion step that
   GPOF avoids.
3. **Iterative non-linear least squares (Levenberg–Marquardt on
   `(b_n, a_n)`).** Rejected. The literature has multiple decades of
   evidence that GPOF's closed-form solution is both faster and
   more robust than iterative fitting on this class of problem.
   Iteration introduces a convergence-tolerance knob with no
   physical justification.

## Consequences

**What becomes easier:**

- **`mom-002` exercises a credible multi-image DCIM kernel.** The
  Phase 1.1.0 one-image placeholder produced |Z_in| ≈ 14 kΩ; the
  multi-image fit produces ~2.7 kΩ. The mesh-resolution and
  surface-wave-pole legs of the remaining gap are now isolable
  (ADR-0024 and ADR-0025 respectively).
- **GPOF is available as a building block** for any future spectral
  fitting task in `yee-mom` (multilayer with multiple dielectric
  slabs, anisotropic substrates, etc.). The `gpof.rs` module is
  self-contained pure-Rust linear algebra and has no
  microstrip-specific assumptions.
- **The `n_images = 1` tripwire prevents accidental
  regression** of the Phase 1.1.0 baseline. Anyone refactoring
  `MultilayerGreens` in the future has a bit-for-bit reference
  point.

**What becomes harder:**

- **`mom-002` is still well outside the Hammerstad–Jensen target.**
  The Phase 1.1.1.0 escape hatch fired: a PEC-mirror probe on the
  same mesh gives |Z_in| ≈ 1.4 kΩ, which floors the achievable
  result at the *current mesh resolution* regardless of how many
  images GPOF fits. The remaining error is split between mesh
  (1.1.1.1) and surface-wave poles (1.1.1.2); GPOF alone cannot
  close either.
- **The `n_images` parameter is a knob.** `n_images = 5` is the
  production choice for `mom-002` but is not a universal default;
  larger `n_images` may overfit or hit the SVD-truncation floor on
  cleaner spectra. The docstring records the per-case choice.
- **The Aksun contour is hard-coded.** A future multilayer geometry
  with poles outside the `k_{z0}(t) = k_0 (1 − j t)` valid region
  will need a contour-selection seam. None of the planned Phase 1
  validation cases (mom-002 through mom-006) need this.

**What's now closed off:**

- Prony / Levenberg–Marquardt fitting paths inside `gpof.rs`. GPOF
  is the only fitter; if a future case needs a different one, it
  goes in a sibling module, not a flag on `gpof.rs`.
- Quietly tightening the `mom-002` tolerance to Hammerstad–Jensen
  on the strength of the multi-image improvement alone. That
  tightening requires both Phase 1.1.1.1 mesh refinement and Phase
  1.1.1.2 pole subtraction.

## References

- `crates/yee-mom/src/gpof.rs` — GPOF (Hua–Sarkar matrix-pencil)
  implementation; synthetic three-image recovery test.
- `crates/yee-mom/src/greens/multilayer.rs` —
  `new_microstrip_with_n_images`; `n_images = 1` back-compat
  tripwire.
- `crates/yee-mom/src/greens/spec.rs` — `GreensSpec::MicrostripDcim`
  variant; `microstrip_dcim` convenience constructor.
- `crates/yee-validation/src/mom_002.rs` — routes through
  `n_images = 5`; loose `[1, 100 kΩ]` gate retained.
- Commits 23decd6 (GPOF + N-image DCIM), da5d6d9 (`GreensSpec`
  variant), a7dd313 (`mom-002` routing), f9e63c7 (Track OOOO merge).
- Y. Hua and T. K. Sarkar, "Matrix pencil method for estimating
  parameters of exponentially damped/undamped sinusoids in noise,"
  *IEEE Trans. ASSP*, vol. 38, no. 5, pp. 814–824, May 1989.
- M. I. Aksun, "A robust approach for the derivation of
  closed-form Green's functions," *IEEE TMTT*, vol. 44, no. 5,
  pp. 651–658, May 1996.
- ADR-0015 — `GreensSpec` builder; this ADR adds one variant.
- ADR-0024 — Phase 1.1.1.1 edge-clustered strip mesh.
- ADR-0025 — Phase 1.1.1.2 Sommerfeld pole extraction; ships the
  pole-subtracted GPOF that this ADR's machinery feeds into.
- CLAUDE.md §10 — `MultilayerGreens` placeholder caveat; this ADR
  narrows it from "one-image DCIM" to "multi-image DCIM, no pole
  subtraction yet".
