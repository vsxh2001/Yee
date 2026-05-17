# Phase 1.1.1.0 — Multi-image DCIM (Discrete Complex Image Method, N ≥ 1 images)

**Status:** Draft  
**Owner:** TBD  
**Phase:** 1.1.1.0  
**Depends on:** Phase 1.1.0 (one-image DCIM, shipped), Phase 1.0 (mom-001, shipped), Track RRR (GreensSpec builder, shipped)  
**Blocks:** Phase 1.1.1 full Sommerfeld extraction; tight mom-002 / mom-003 tolerances

## Assumption being challenged

CLAUDE.md §10 records `MultilayerGreens` as a Phase 1.1.0 placeholder implementing **one-image** DCIM. The one-image approximation is the simplest possible Sommerfeld-integral surrogate: it places a single complex image charge below the ground plane to approximate the spectral-domain Green's function. Track MMM measured that this placeholder gives mom-002 a 14 kΩ input impedance vs. the Hammerstad-Jensen 50 Ω target — a factor of ~300 off. Not because the math is wrong, but because **one image is too few**.

The full fix is Phase 1.1.1 with proper Sommerfeld-integral evaluation. That's months of work: contour deformation in the complex spectral plane, Sommerfeld tails, surface-wave pole extraction, Prony / GPOF fitting of the spectral Green's function.

This spec asks: **can a multi-image DCIM (N = 5 to 10 images) close most of the gap with a fraction of the work?** Literature (Aksun 1996, Yang 1997) says yes — for moderate-thickness substrates, N = 5 images gets you within ~2-5% of full Sommerfeld evaluation. That's enough to get mom-002 tolerances inside ±10% of Hammerstad-Jensen, which would be a major improvement over current ±factor-of-2.

## Scope

In:
- Extend `MultilayerGreens` in `crates/yee-mom/src/multilayer.rs` to N images (configurable).
- Implement GPOF (Generalized Pencil of Function) fitting to extract complex image coefficients from a sampled spectral Green's function.
- Loose ±10% gate on mom-002 against Hammerstad-Jensen 50 Ω. Tighter ±3% gate still gated on Phase 1.1.1.

Out:
- Sommerfeld-integral evaluation (contour deformation, Bessel-J tails). Phase 1.1.1.
- Anisotropic / lossy / dispersive substrates. Phase 1.1.2+.
- Multilayer stacks with more than substrate + ground plane. Phase 1.1.3.
- Surface wave pole extraction. Phase 1.1.1.

## Approach

For a microstrip on a grounded dielectric slab (εr, h), the spectral-domain Green's function for the vector potential has the form:

$$\tilde{G}^A(k_\rho) = \sum_{n=1}^N \frac{a_n}{j2k_{z0}} \cdot e^{-jk_{z0} z_n}$$

where $(a_n, z_n)$ are complex coefficients/locations of the $n$-th image. The N=1 placeholder uses an arbitrary $(a_1, z_1)$ choice; multi-image DCIM fits N coefficients to match the true spectral Green's function $\tilde{G}^A(k_\rho)$ at a set of sample points.

GPOF fitting (Hua & Sarkar 1989, Aksun 1996):
1. Sample $\tilde{G}^A(k_\rho)$ at $M \ge 2N$ uniform points $k_{\rho,m}$ along a deformed integration contour.
2. Stack into a Hankel matrix $H \in \mathbb{C}^{(M-N) \times N}$.
3. SVD the matrix; the top-N right-singular vectors yield the image locations $z_n$ via a small generalized eigenproblem.
4. Solve a linear least-squares for the amplitudes $a_n$.

The actual spectral $\tilde{G}^A$ at substrate $(\epsilon_r, h)$ is given in Pozar §3.7 (microstrip) and Aksun §III in closed form — sampling is just plugging into a formula.

## Public API

```rust
//! Multi-image DCIM Greens for grounded dielectric slab.

pub struct MultilayerGreens {
    /// Substrate relative permittivity.
    pub eps_r: f64,
    /// Substrate thickness (m).
    pub h: f64,
    /// Number of DCIM image terms (default 5; N=1 reproduces Phase 1.1.0 behavior).
    pub n_images: usize,
    /// Fitted complex coefficients (built once at construction).
    coeffs: Vec<(Complex64, Complex64)>, // (a_n, z_n) per image
}

impl MultilayerGreens {
    /// Build with N-image DCIM fit for the given (eps_r, h) at the given freq.
    pub fn new_microstrip_with_n_images(eps_r: f64, h: f64, freq_hz: f64, n_images: usize) -> Self;

    /// Existing single-image constructor; preserved.
    pub fn new_microstrip(freq_hz: f64, eps_r: f64, h: f64) -> Self {
        Self::new_microstrip_with_n_images(eps_r, h, freq_hz, 1)
    }
}
```

`GreensSpec` (Track RRR) adds a sibling variant:

```rust
pub enum GreensSpec {
    FreeSpace,
    Microstrip { eps_r: f64, h_m: f64 },
    MicrostripDcim { eps_r: f64, h_m: f64, n_images: usize },
}
```

mom-002's wire-up in `yee-validation` updates to use `MicrostripDcim { eps_r: 4.4, h_m: 1.6e-3, n_images: 5 }`.

## Definition of done

1. `MultilayerGreens::new_microstrip_with_n_images` exists and builds N coefficient pairs without panic for $N \in [1, 10]$.
2. GPOF fit converges within 100 iterations (SVD-based, deterministic — no iteration needed actually, just a closed-form linear-algebra path).
3. `GreensSpec::MicrostripDcim` variant added.
4. `yee-validation`'s `run_mom_002` uses `n_images = 5`. The resulting `|Z_in|` at 1 GHz on the 30 mm microstrip strip mesh must be within `[35, 75] Ω` (i.e. ±50% of 50 Ω) — looser than the eventual ±3% Hammerstad-Jensen gate but far tighter than the current 100 kΩ degenerate bound.
5. Unit tests on the GPOF fitting machinery: synthetic 3-image data must be recovered to within `1e-6` relative error.
6. Verification chain green on `cargo build / clippy / test --release / fmt --check` for `-p yee-mom -p yee-validation`.

## Lane

`crates/yee-mom/**` (primary; multilayer.rs + lib.rs Greens-spec variant) + `crates/yee-validation/**` (mom-002 wire-up + tolerance update). No other crates.

## Verification

```bash
cargo build -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test -p yee-mom --release
cargo test -p yee-validation --release
cargo fmt --check --all
```

mom-001 must still pass (don't regress).

## Escape hatch

If GPOF fitting is numerically unstable at N ≥ 5 (small singular-value ratios → noise dominates), drop the test gate to "within ±100% of 50 Ω" — still a 100× improvement over the current 14 kΩ — and document. The full Sommerfeld path (Phase 1.1.1) is the real fix; this spec is the bridge, not the destination.

## References

- Aksun, "A robust approach for the derivation of closed-form Green's functions", IEEE TMTT 1996, 44(5).
- Yang, Mittra, Itoh, "Discrete complex image method for analyzing planar microstrip structures", *Time-Domain Methods for the Maxwell's Equations*, 1997.
- Hua & Sarkar, "Generalized Pencil-of-Function method for extracting poles of an EM system from its transient response", IEEE TAP 1989, 37(2).
- Pozar, *Microwave Engineering* 4th ed., §3.7 (Microstrip).
- Mosig, "Integral Equation Technique for Planar Geometries", 1989.
