# Phase 1.1.1.2 — Sommerfeld surface-wave pole extraction

**Status:** Draft
**Owner:** TBD
**Phase:** 1.1.1.2
**Depends on:** Phase 1.1.1.0 (N-image DCIM, shipped — Track OOOO), Phase 1.1.1.1 (edge-clustered strip mesh, shipped — Track AAAAA)
**Blocks:** Phase 1.1.1.3 (full Sommerfeld integration with contour deformation); tight tolerances on mom-002 (microstrip Z₀), mom-003 (2.4 GHz patch), mom-004 (Wilkinson), mom-005 (branch-line hybrid), mom-006 (Swanson hairpin BPF)

## Assumption being challenged

Phase 1.1.1.0 shipped a GPOF-fitted N-image DCIM kernel that approximates the spectral-domain Green's function $\tilde{G}^A(k_\rho)$ of a grounded dielectric slab as a sum of complex exponentials in $k_{z0}$. Phase 1.1.1.1 then refined the strip-discretisation mesh (edge-clustered, nz = 16) so mesh-error was not the bottleneck. Together they should have closed mom-002 to inside the Hammerstad-Jensen [35, 75] Ω corridor at 1 GHz on FR-4 (ε_r = 4.4, h = 1.6 mm).

They did not. Track AAAAA's measurement is unambiguous: $|Z_\text{in}|$ floors at ~2 kΩ and $\text{Im}(Z_\text{in}) \approx -2.1$ kΩ across mesh refinements from nz = 8 to nz = 32. Refinement does nothing because **the residual is not mesh-bound**. It is a **surface-wave pole signature**: the spectral denominator of the grounded slab has zeros on or near the real $k_\rho$ axis corresponding to the slow-wave bound modes of the substrate (TM₀, TE₁, …). At FR-4 / 1.6 mm / 1 GHz, the TM₀ pole sits at $k_\rho / k_0 \approx 1.6$ — squarely between $k_0$ (air) and $\sqrt{\varepsilon_r}\, k_0$ (slab bulk wavenumber), exactly where the bound surface-wave mode lives.

GPOF cannot fit a near-singular spectral function with smooth complex exponentials. The image-train converges in coefficient norm but never resolves the residue at the pole; the resulting space-domain Green's function is missing a Hankel-function term that decays only as $\rho^{-1/2}$ along the substrate, while the image-sum decays as $\rho^{-1}$. At the geometry scales of the mom-002 microstrip (W ≈ h, L = 30 mm ≫ h), the missing surface-wave contribution dominates the input impedance.

The fix is **analytic pole subtraction**: locate the discrete surface-wave poles by Newton-Raphson root-find in the complex $k_\rho$ plane, extract their residues in closed form, evaluate the surface-wave contribution analytically (Hankel function of the second kind at each pole's $k_{\rho,p}$), and run GPOF on the **pole-subtracted** smooth remainder. The full Phase 1.1.1.3 Sommerfeld integration is still months away — this spec is the targeted accuracy fix that unblocks every multilayer benchmark currently gated on "loose tolerances".

## Approach

For the grounded dielectric slab, the spectral-domain reflection coefficient has independent TE / TM channels whose denominators are

$$D_\text{TE}(k_\rho) = k_{z0} + k_{zd}\cot(k_{zd} h) \quad\text{(TE pole condition)}$$
$$D_\text{TM}(k_\rho) = \varepsilon_r k_{z0} + k_{zd}\cot(k_{zd} h) \quad\text{(TM pole condition)}$$

(Pozar §3.7, eq. 3.196–3.199; equivalent transverse-resonance forms in Felsen & Marcuvitz §5.) Surface-wave poles are zeros of these denominators on the proper Riemann sheet. For a lossless slab the poles lie just below the real $k_\rho$ axis on the sheet of $\text{Im}(k_{z0}) < 0$; mathematically they are real-axis poles indented into the lower half plane by the standard radiation-condition $\text{Im}(\omega) \to 0^-$ limit.

**Step 1 — Pole search.** Newton-Raphson root-find in the complex $k_\rho$ plane:

$$k_\rho^{(n+1)} = k_\rho^{(n)} - D(k_\rho^{(n)}) / D'(k_\rho^{(n)})$$

Initial guess from the quasi-static effective-permittivity approximation:

$$k_{\rho,0} \approx k_0 \sqrt{(\varepsilon_r + 1)/2}$$

which for FR-4 at 1 GHz gives $k_{\rho,0}/k_0 \approx 1.64$. The Jacobian $D'(k_\rho)$ is available in closed form via implicit differentiation of $k_{z0}(k_\rho) = \sqrt{k_0^2 - k_\rho^2}$ and $k_{zd}(k_\rho) = \sqrt{\varepsilon_r k_0^2 - k_\rho^2}$. Convergence tolerance: $|D(k_\rho)| < 10^{-12}$. Max iterations: 50 (typical converges in 5–10).

**Step 2 — Residue extraction.** At each converged pole $k_{\rho,p}$, the residue of the spectral Green's function is

$$\text{Res}_p\big[\tilde{G}^A(k_\rho)\big] = \frac{N(k_{\rho,p})}{D'(k_{\rho,p})}$$

where $N$ is the numerator of the spectral Green's function (also closed-form per Pozar §3.7). For the TM₀ pole on a grounded slab this reduces to the well-known Felsen-Marcuvitz §5 form involving the modal field profile evaluated at the slab top.

**Step 3 — Pole subtraction + GPOF fit on the smooth remainder.**

$$\tilde{G}^A_\text{reg}(k_\rho) = \tilde{G}^A(k_\rho) - \sum_p \frac{\text{Res}_p}{k_\rho - k_{\rho,p}}$$

is smooth everywhere on the integration contour; the existing OOOO GPOF machinery in `crates/yee-mom/src/gpof.rs` fits it cleanly with N = 5 images. The TE and TM channels are pole-subtracted independently using the corresponding $D_\text{TE}$ / $D_\text{TM}$ zero sets.

**Step 4 — Space-domain reconstruction.** The full space-domain Green's function is the sum of (a) the image-sum approximation of the regularised part, computed with the existing `MultilayerGreens::image_sum`, plus (b) the analytic surface-wave contribution:

$$G^A_\text{sw}(\rho, z, z') = \sum_p \text{Res}_p \cdot \frac{-j}{4} H_0^{(2)}(k_{\rho,p}\,\rho) \cdot \psi_p(z)\,\psi_p(z')$$

where $\psi_p(z)$ is the modal $z$-profile of the $p$-th surface wave (sinusoid inside the slab, exponentially decaying above). The Hankel $H_0^{(2)}$ asymptotic form gives the canonical $\rho^{-1/2}$ surface-wave decay rate that the pure image-sum cannot reproduce.

## Public API

Extend `MultilayerGreens` to optionally carry a pole list and modal-profile data. Construction-time pole search is one Newton solve per frequency per channel — milliseconds, negligible vs. the per-frequency MoM matrix fill.

```rust
pub struct MultilayerGreens {
    // ... existing OOOO fields (k0, eta0, eps_r, h, n_images,
    //     vector_images, scalar_images) ...

    /// Surface-wave poles of the TE channel, with attached residues
    /// and modal z-profile coefficients. Empty if pole extraction is
    /// disabled (i.e. `n_surface_wave_poles = 0`) or if Newton failed
    /// to converge for every initial guess.
    pub te_surface_waves: Vec<SurfaceWavePole>,
    /// Same for the TM channel. The dominant TM₀ pole lives here.
    pub tm_surface_waves: Vec<SurfaceWavePole>,
}

/// A single surface-wave pole: its complex `k_ρ` location, its
/// residue, and enough modal data to evaluate the closed-form
/// Hankel contribution at arbitrary (ρ, z, z').
pub struct SurfaceWavePole {
    /// Complex pole location in the `k_ρ` plane.
    pub k_rho: Complex64,
    /// Residue of the spectral Green's function at the pole.
    pub residue: Complex64,
    /// Modal `k_zd` inside the slab (for the modal `z`-profile).
    pub k_zd: Complex64,
}

impl MultilayerGreens {
    /// Build with N-image DCIM + pole-subtracted Sommerfeld
    /// extraction. `n_surface_wave_poles = 0` reproduces the
    /// Phase 1.1.1.0 (OOOO) GPOF-only path bit-for-bit.
    pub fn new_microstrip_sommerfeld(
        eps_r: f64,
        h: f64,
        freq_hz: f64,
        n_images: usize,
        n_surface_wave_poles: usize,
    ) -> Self;
}
```

`GreensSpec` gains a sibling variant; the OOOO `MicrostripDcim` variant stays for back-compat and as the fast path when surface waves are known absent (free-space, or high-frequency where the desired bandwidth is above the lowest-pole cutoff).

```rust
pub enum GreensSpec {
    FreeSpace,
    Microstrip { eps_r: f64, h_m: f64 },                       // Phase 1.1.0
    MicrostripDcim { eps_r: f64, h_m: f64, n_images: usize },  // Phase 1.1.1.0
    MicrostripSommerfeld {
        eps_r: f64,
        h_m: f64,
        n_images: usize,
        n_surface_wave_poles: usize,  // default 2: TM₀ + TE₁ envelopes
    },
}
```

Default `n_surface_wave_poles = 2` covers the FR-4 sweep up to ~10 GHz where only TM₀ is bound; the second slot is a placeholder that gracefully no-ops (Newton fails to find a second pole) if no second mode exists at the current frequency.

## Definition of done

1. Newton pole-search converges deterministically for FR-4 (ε_r = 4.4, h = 1.6 mm) at 1 GHz, 2.4 GHz, and 5 GHz. A unit test surfaces a sanity table:
   - 1 GHz: $k_{\rho,\text{TM₀}}/k_0 \in [1.55, 1.70]$, $|D| < 10^{-12}$, $\le 15$ iterations.
   - 2.4 GHz: $k_{\rho,\text{TM₀}}/k_0 \in [1.60, 1.75]$.
   - 5 GHz: $k_{\rho,\text{TM₀}}/k_0 \in [1.70, 1.90]$.
2. mom-002 `|Z_in|` at 1 GHz on the AAAAA nz = 16 edge-clustered mesh, with `MicrostripSommerfeld { eps_r: 4.4, h_m: 1.6e-3, n_images: 5, n_surface_wave_poles: 2 }`, lands in **[35, 75] Ω** — the Hammerstad-Jensen target. Tightens from the Phase 1.1.1.0 loose-tolerance gate.
3. mom-001 (free-space dipole, NEC-4 87 + j41 Ω) is **unchanged**: the surface-wave path is gated on `GreensSpec::MicrostripSommerfeld` only; `FreeSpace` and `Microstrip*` variants execute the OOOO codepath verbatim. CI test `dipole_z_at_resonance` stays green.
4. `MultilayerGreens::new_microstrip_sommerfeld(.., n_surface_wave_poles = 0)` produces bit-for-bit identical `(vector_images, scalar_images)` as `new_microstrip_with_n_images` — the OOOO back-compat tripwire.
5. Unit test on synthetic data: a hand-constructed spectral function with a known pole at $k_\rho = 1.6\, k_0$ + residue 1.0 + a smooth analytic remainder is recovered by the pole-search + subtraction pipeline to within $10^{-9}$ relative error.
6. `cargo build / clippy / test --release / fmt --check / doc --no-deps` on `-p yee-mom -p yee-validation`: all exit 0.

## Lane (impl)

`crates/yee-mom/**` (multilayer.rs additions; new `crates/yee-mom/src/sommerfeld.rs` for Newton + residue + Hankel; `GreensSpec` variant in lib.rs) + `crates/yee-validation/**` (mom-002 wire-up to the new variant; tolerance tightening). No other crates. The `bessel` Hankel-function evaluation pulls in `num_complex` only (Hankel asymptotic for large argument + small-argument series — same convention as Pozar, no new dep needed).

## Verification

```bash
cargo build  -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo test   -p yee-validation --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-mom
```

mom-001 must remain green (no surface-wave path executed); mom-002 must pass the new tight gate.

## Escape hatch

Two failure modes are anticipated.

**A — Newton fails to converge.** If the initial guess $k_{\rho,0} = k_0\sqrt{(\varepsilon_r + 1)/2}$ is far enough from the true pole that Newton overshoots into a different basin (or onto the wrong Riemann sheet), fall back to **Müller's method** (three-point parabolic extrapolation, no derivative needed; converges quadratically near a simple pole and tolerates a worse initial guess). If Müller also fails, surface as Phase 1.1.1.3 (full Sommerfeld integration with contour deformation around the pole locus) and degrade mom-002's gate back to the OOOO loose-tolerance ±50% corridor. Document the failure-mode signature in `validation/README.md`: which (ε_r, h, f) combo, which channel, what the initial-guess-vs-converged-pole distance was, where Newton's iterates landed.

**B — Pole found but residue extraction numerically unstable.** Symptom: $|D'(k_{\rho,p})| < 10^{-10}$ (double-pole proximity or near-cutoff degeneracy). Fall back to finite-differenced $D'$ with adaptive $h$-stepping, or — if that too fails — discard the pole and proceed with the OOOO image-only fit, documenting the loss. The TM₀ pole on FR-4 / 1.6 mm well above its cutoff is far from any degeneracy; this branch is anticipated only for edge cases (very thin substrates near the cutoff frequency of the next mode).

Blocked > 25 min → surface and stop.

## References

- D. M. Pozar, *Microwave Engineering*, 4th ed., §3.7, eq. 3.196–3.199 (grounded dielectric slab dispersion, TE/TM denominators).
- J. R. Mosig, "Integral equation technique for planar geometries," in *Numerical Techniques for Microwave and Millimeter-Wave Passive Structures*, T. Itoh (ed.), Wiley, 1989, Ch. 3.
- M. I. Aksun, "A robust approach for the derivation of closed-form Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44, no. 5, pp. 651–658, May 1996. (Aksun-1996 explicitly calls for **surface-wave pole extraction before GPOF**; this spec implements the missing step.)
- K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation formulations," *IEEE Trans. Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar 1997 — the canonical reference for surface-wave-augmented spectral-domain Green's functions in planar IE methods.
- L. B. Felsen and N. Marcuvitz, *Radiation and Scattering of Waves*, IEEE Press, 1994, Ch. 5 — rigorous Hankel-function asymptotic derivation of the surface-wave residue contribution.
