# Phase 1.1.1.2 — Sommerfeld surface-wave pole extraction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to drive this plan task-by-task.

**Goal:** Land Newton-Raphson surface-wave pole extraction inside `MultilayerGreens`, subtract the pole contribution from the spectral Green's function before the existing GPOF fit, add an analytic Hankel-function surface-wave term to the space-domain reconstruction, and tighten mom-002 from "loose tolerance" to the Hammerstad-Jensen `|Z_in| ∈ [35, 75] Ω` corridor at 1 GHz on FR-4.

**Companion spec:** `docs/superpowers/specs/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`

**Architecture:** One new `pub(crate) mod sommerfeld` under `crates/yee-mom/src/sommerfeld.rs` holding Newton + Müller-fallback pole search, residue extraction, and Hankel `H_0^{(2)}` evaluation. `MultilayerGreens` in `multilayer.rs` gains two new fields (`te_surface_waves`, `tm_surface_waves`) and one new constructor (`new_microstrip_sommerfeld`). `GreensSpec` in `lib.rs` gains a `MicrostripSommerfeld` variant. `yee-validation`'s mom-002 driver swaps to the new variant.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `crates/yee-mom/src/sommerfeld.rs` | Create | Newton/Müller pole search, residue, Hankel H₀⁽²⁾ evaluation |
| `crates/yee-mom/src/multilayer.rs` | Modify | Add `SurfaceWavePole`, `te_surface_waves`/`tm_surface_waves` fields, `new_microstrip_sommerfeld` constructor, pole-subtracted `fit_slab_dcim`, Hankel-term contribution in `scalar_vector` / `scalar_scalar` |
| `crates/yee-mom/src/lib.rs` | Modify | `GreensSpec::MicrostripSommerfeld` variant + `GreensSpec::build` arm |
| `crates/yee-mom/tests/sommerfeld_pole_search.rs` | Create | Step 1 verification: Newton convergence sanity table for FR-4 at 1 / 2.4 / 5 GHz |
| `crates/yee-mom/tests/sommerfeld_synthetic.rs` | Create | Step 2/3 verification: synthetic pole + residue recovery to 1e-9 |
| `crates/yee-validation/src/lib.rs` | Modify | mom-002 driver uses `GreensSpec::MicrostripSommerfeld`; tighten tolerance |
| `crates/yee-validation/validation/README.md` | Modify (if present) | mom-002 row updates from "loose" to "Hammerstad-Jensen [35, 75] Ω" |

---

## Step 0 — Re-derive the spectral denominators symbolically

**Files:** none yet — pen-and-paper / scratch markdown.

Before writing Newton code, write out $D_\text{TE}(k_\rho)$ and $D_\text{TM}(k_\rho)$ from first principles and cross-check against Pozar §3.7 (eq. 3.196–3.199). The transverse-resonance form on a grounded slab (PEC short at $z = -h$, half-space radiation at $z \to +\infty$) is:

- $k_{z0}(k_\rho) = \sqrt{k_0^2 - k_\rho^2}$, principal branch with $\text{Im}\,k_{z0} \le 0$ (outgoing wave above slab).
- $k_{zd}(k_\rho) = \sqrt{\varepsilon_r k_0^2 - k_\rho^2}$, real & positive when $k_\rho < \sqrt{\varepsilon_r}\,k_0$ (bound mode regime).
- $D_\text{TE}(k_\rho) = k_{z0} + k_{zd}\cot(k_{zd} h) = 0$.
- $D_\text{TM}(k_\rho) = \varepsilon_r k_{z0} + k_{zd}\cot(k_{zd} h)$, **but note** that for grounded-slab the canonical TM dispersion uses $\tan$ not $\cot$ (Pozar 3.197). Verify the sign / cot-vs-tan convention against Pozar before coding.

Confirm against a published TM₀ cutoff sanity check: TM₀ has no cutoff (propagates down to DC); the quasi-static limit gives $k_{\rho,\text{TM₀}} \to k_0\sqrt{(\varepsilon_r+1)/2}$ as $k_0 h \to 0$. This is the initial-guess formula used in Step 1.

Estimated LOC: 0 (paper). Verification: agent writes the two denominator expressions into the top doc-comment of `sommerfeld.rs` with a citation to Pozar eq. 3.196 / 3.197.

---

## Step 1 — Newton-Raphson pole search in complex $k_\rho$ plane

**Files:** `crates/yee-mom/src/sommerfeld.rs` (create), `crates/yee-mom/src/lib.rs` (declare `mod sommerfeld`).

- [ ] Implement `fn d_te(k_rho: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64` and `fn d_tm(...)` mirroring the Pozar forms verified in Step 0. Use `Complex64::sqrt` (principal-branch) for $k_{z0}$ and $k_{zd}$.
- [ ] Implement `fn d_prime_te / d_prime_tm` — closed-form derivatives via implicit differentiation: $\partial k_{z0}/\partial k_\rho = -k_\rho/k_{z0}$, $\partial k_{zd}/\partial k_\rho = -k_\rho/k_{zd}$, then chain through the $\cot(k_{zd} h)$ term using $d/dx[\cot x] = -\csc^2 x = -1 - \cot^2 x$.
- [ ] Implement `fn newton_pole(d, d_prime, k_rho_0, tol = 1e-12, max_iter = 50) -> Result<(Complex64, usize), PoleSearchError>` — pure Newton, deterministic. Return `(pole, iters_used)` on convergence or `PoleSearchError::NoConvergence { last_value, last_residual }` on failure.
- [ ] Implement `fn quasi_static_guess(eps_r: f64, k0: f64) -> Complex64` returning $k_0\sqrt{(\varepsilon_r + 1)/2}$ as a real-axis seed.
- [ ] Unit test `pole_search_fr4_at_three_frequencies` (`tests/sommerfeld_pole_search.rs`):
  - For (4.4, 1.6e-3, 1e9), (4.4, 1.6e-3, 2.4e9), (4.4, 1.6e-3, 5e9):
    - Seed Newton at the quasi-static guess.
    - Assert convergence ($|D| < 10^{-12}$) in $\le 15$ iterations.
    - Assert $k_{\rho,\text{TM₀}} / k_0$ lies in the spec'd corridor (1 GHz: [1.55, 1.70]; 2.4 GHz: [1.60, 1.75]; 5 GHz: [1.70, 1.90]).
  - Surface a small printout (only when `RUST_LOG=debug`) of the iter count + final $|D|$ so the sanity table is reproducible.

Estimated LOC: ~250 (Newton + derivative + tests). Verification: `cargo test -p yee-mom --release sommerfeld_pole_search` exits 0.

---

## Step 2 — Residue extraction at each pole

**Files:** `crates/yee-mom/src/sommerfeld.rs`.

- [ ] Define `pub struct SurfaceWavePole { pub k_rho, pub residue, pub k_zd }` (3 `Complex64` fields). Lives in `multilayer.rs` per spec API, but the computation is in `sommerfeld.rs`.
- [ ] Implement `fn residue_te / residue_tm (pole: Complex64, eps_r: f64, h: f64, k0: f64) -> Complex64` using $\text{Res}_p = N(k_{\rho,p}) / D'(k_{\rho,p})$. The numerator $N$ comes from the spectral Green's function: for the TE channel of $\tilde G^A$, $N_\text{TE}(k_\rho) = 2 k_{zd} \csc(k_{zd} h) / (\text{denominator-derivative-free part})$ — derive the exact form against Michalski-Mosig 1997 eq. (16)-(19) and verify dimensions before coding.
- [ ] Unit test `residue_smoke_fr4_1ghz`: residue is finite, non-zero, and finite-imaginary (the lossless slab gives a purely real residue on the proper sheet; lossy substrates would give complex). Tolerance loose — this is a sanity check, not the validation gate.
- [ ] Define `pub enum PoleSearchError { NoConvergence { ... }, DegeneratePole { d_prime_norm: f64 } }`. The `DegeneratePole` arm triggers when $|D'(k_{\rho,p})| < 10^{-10}$ at the converged pole — the escape-hatch fallback to finite-differenced $D'$.

Estimated LOC: ~150. Verification: `cargo test -p yee-mom --release sommerfeld::tests::residue_smoke_fr4_1ghz` exits 0.

---

## Step 3 — Pole-subtracted GPOF fit + Hankel-term contribution

**Files:** `crates/yee-mom/src/sommerfeld.rs` (Hankel), `crates/yee-mom/src/multilayer.rs` (pole-subtracted sampling + space-domain reconstruction).

- [ ] Implement `fn hankel_h0_2(z: Complex64) -> Complex64` — the Hankel function $H_0^{(2)}(z) = J_0(z) - j Y_0(z)$. For mom-002 geometry $\rho \in [10^{-4}, 0.1]$ m and $k_{\rho,p} \approx 1.6 k_0 \approx 33$ rad/m at 1 GHz, the argument is $|z| \in [3.3 \times 10^{-3}, 3.3]$ — straddles the small-arg / large-arg threshold. Use:
  - Small-argument series ($|z| < 8$): Hankel asymptotic series with $\ln(z/2) + \gamma$ logarithmic-term handling.
  - Large-argument asymptotic ($|z| \ge 8$): $H_0^{(2)}(z) \approx \sqrt{2/(\pi z)} \exp(-j(z - \pi/4))$.
  - Cross-check at $z = 8$: both branches agree to $10^{-10}$ relative.
- [ ] Modify `fn fit_slab_dcim` in `multilayer.rs` to optionally accept a list of `(k_rho_p, residue_p)` pairs. If non-empty, subtract $\sum_p \text{Res}_p / (k_\rho(t) - k_{\rho,p})$ from the sampled reflection coefficient before passing to GPOF. The original signature stays as a `..._unsubtracted` private helper; the spec'd public path goes through the new pole-aware entry point.
- [ ] Implement `MultilayerGreens::surface_wave_contribution(&self, r1, r2) -> Complex64` summing the Hankel-function term across all poles in `te_surface_waves` / `tm_surface_waves` for the requested channel. The modal $z$-profile $\psi_p(z) \psi_p(z')$ uses $\cos(k_{zd} \cdot (z + h)) / \cos(k_{zd} h)$ inside the slab (PEC ground at $-h$, peaks at $z = 0$); above the slab it is $\exp(-\alpha_0(z+z'))$ with $\alpha_0 = \sqrt{k_{\rho,p}^2 - k_0^2}$. For mom-002 strips at $z = 0$ both points sit at the slab top, so $\psi_p(0) = 1$ and the term collapses.
- [ ] Unit test `synthetic_pole_recovery` (`tests/sommerfeld_synthetic.rs`): construct a hand-rolled spectral function $\tilde G(k_\rho) = R_0/(k_\rho - k_{\rho,p}^*) + \text{(smooth analytic)}$ with $k_{\rho,p}^* = 1.6 k_0$, $R_0 = 1.0$. Run Newton on this synthetic $D$; assert the recovered pole and residue match the planted values to $10^{-9}$ relative error.

Estimated LOC: ~300 (Hankel ~120, pole-subtraction wiring ~80, Hankel-term Greens ~50, tests ~50). Verification: `cargo test -p yee-mom --release sommerfeld_synthetic` exits 0.

---

## Step 4 — `GreensSpec::MicrostripSommerfeld` variant + constructor

**Files:** `crates/yee-mom/src/lib.rs`, `crates/yee-mom/src/multilayer.rs`.

- [ ] Add `GreensSpec::MicrostripSommerfeld { eps_r: f64, h_m: f64, n_images: usize, n_surface_wave_poles: usize }` variant. `#[derive(Debug, Clone, Copy)]` already covers it.
- [ ] Add `GreensSpec::microstrip_sommerfeld(eps_r, h_m, n_images, n_sw_poles)` convenience constructor mirroring the OOOO style.
- [ ] Extend `GreensSpec::build` with the new arm routing to `MultilayerGreens::new_microstrip_sommerfeld(eps_r, h_m, freq_hz, n_images, n_sw_poles)`.
- [ ] Implement `MultilayerGreens::new_microstrip_sommerfeld`:
  - For each channel (TE, TM), run Newton from the quasi-static guess.
  - On Newton success: extract residue, store in `te_surface_waves` / `tm_surface_waves`. Subtract the pole from the sampled spectral function before the existing GPOF path.
  - On Newton failure (or `n_surface_wave_poles == 0`): pole list empty, fall through to the OOOO `fit_slab_dcim` path unchanged. **Important:** logged-but-non-fatal — the kernel must still build.
  - When more than 1 pole is requested per channel, seed the second Newton run from a higher-mode quasi-static estimate ($k_0 \sqrt{\varepsilon_r}$ for the TE₁ cutoff) and accept the result only if it converges to a *distinct* pole (i.e. $|k_{\rho,2} - k_{\rho,1}| > 0.01 k_0$).
- [ ] Override `Greens::scalar_vector / scalar_scalar / *_smooth` impls to add the surface-wave Hankel contribution on top of the existing image-sum + free-space terms.
- [ ] Unit test `n_sw_poles_zero_matches_phase_1_1_1_0`: with `n_surface_wave_poles = 0`, the constructed `MultilayerGreens` produces bit-for-bit identical `vector_images` / `scalar_images` to `new_microstrip_with_n_images`. OOOO tripwire preserved.

Estimated LOC: ~250 (constructor + Greens-trait overrides + tripwire test). Verification: `cargo test -p yee-mom --release multilayer::tests` exits 0; `n_equals_one_matches_phase_1_1_0` (OOOO) stays green.

---

## Step 5 — mom-002 re-wire + tolerance tighten

**Files:** `crates/yee-validation/src/lib.rs`, `crates/yee-validation/validation/README.md` (if it exists; otherwise add a comment in the driver function's doc).

- [ ] Locate `run_mom_002` (or equivalent driver function) in `yee-validation`. Currently uses `GreensSpec::MicrostripDcim { eps_r: 4.4, h_m: 1.6e-3, n_images: 5 }` (Phase 1.1.1.0). Update to `GreensSpec::MicrostripSommerfeld { eps_r: 4.4, h_m: 1.6e-3, n_images: 5, n_surface_wave_poles: 2 }`.
- [ ] Tighten the `assert!`/`expect_within` gate on `|Z_in|` at 1 GHz: from the OOOO loose `[35, 75] Ω` *target with documented misses* to a real `[35, 75] Ω` *hard gate*. Hammerstad-Jensen `Z_0 ≈ 50 Ω` ± 50% covers manufacturing-tolerance microstrip realisations of W/h ≈ 1.
- [ ] Update or create the validation table row for mom-002 indicating Phase 1.1.1.2 and the hard tolerance.
- [ ] mom-001 driver is untouched (still `GreensSpec::FreeSpace`); CI re-runs `dipole_z_at_resonance` to confirm no regression.

Estimated LOC: ~40 (mostly text + a single GreensSpec swap). Verification: `cargo test -p yee-validation --release run_mom_002` exits 0 with the tightened gate; `cargo test -p yee-mom --release dipole_z_at_resonance` exits 0 (no regression).

---

## Final verification

```bash
cargo build  -p yee-mom -p yee-validation
cargo clippy -p yee-mom -p yee-validation --all-targets -- -D warnings
cargo test   -p yee-mom --release
cargo test   -p yee-validation --release
cargo fmt    --check --all
cargo doc    --no-deps -p yee-mom
```

All six must exit 0. The OOOO tripwire `n_equals_one_matches_phase_1_1_0` and the mom-001 NEC-4 gate `dipole_z_at_resonance` both stay green — surface-wave path is opt-in, free-space path is unchanged.

---

## Estimated total

- LOC: ~990 (Newton + residue + Hankel ~400; multilayer.rs wiring ~250; GreensSpec + validation ~90; tests ~250).
- Wall-time per agent: 2–3 days. The high-risk step is Step 3's Hankel-function implementation — getting the small-argument logarithmic-term right is fiddly; a `cargo test` against a hand-table of $H_0^{(2)}$ values from Abramowitz & Stegun Table 9.x is the fastest sanity check.
- Risk concentration: Step 1's Newton can land on the wrong Riemann sheet if $k_{z0}$'s branch cut is misplaced. The spec's escape hatch (Müller's method + degrade tolerances) is the standard fallback. The OOOO tripwire ensures any regression of the image-only fast path is caught loudly.
