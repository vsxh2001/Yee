# fdtd-206 series-LC resonant frequency gate — Design Spec

**Phase:** 2.fdtd.6.1  
**ADR:** ADR-0080  
**Date:** 2026-05-29  
**Status:** Accepted

---

## 1. Goal

Close the ADR-0017 Phase 2.fdtd.6.1 gap: ship a **quantitative validation
gate for the `series_rlc` ODE integration path** in `LumpedRlcPort`.
The Phase 2.fdtd.6 energy-dissipation gate only verifies the pure-resistor
path; the series-RLC "compiles and self-tests" but has no published-benchmark
accuracy gate against an analytic reference.

The gate: simulate a lumped series-LC port (R small, L=1 nH, C chosen so
f₀=1 GHz exactly) inside a tiny PEC box, excite with a broadband impulse,
extract the ring-down frequency via DFT, and assert that the measured f₀ is
within ±2 % of the analytic 1/(2π√LC).

---

## 2. Physics

For a series R-L-C circuit (R → 0, underdamped), the natural resonant
frequency is

```
f₀ = 1 / (2π √(LC))
```

(Pozar, "Microwave Engineering," 4th ed., §2.4; Hayt & Kemmerly,
"Engineering Circuit Analysis," §14.1).  With L=1 nH and

```
C = 1 / (4π² f₀² L)  ≈  25.330 pF   (for f₀ = 1 GHz)
```

The FDTD integrates the series-RLC ODE alongside the Yee update (Taflove &
Hagness §15.10).  The centred-difference time-stepping introduces a relative
frequency error of order (dt ω₀)² / 24 ≈ 6×10⁻⁴ % — far below the 2 % gate
tolerance.

### Q factor

R = 1.0 Ω → Q = √(L/C)/R = 6.28/1.0 ≈ 6.28 (underdamped).
Amplitude decays by factor e⁻π ≈ 0.043 over Q=6.28 cycles ≈ 6.28 ns.
With a 5 000-step window (≈9.6 ns), the signal decays to ~0.8 % —
enough for a clean DFT peak.

### DFT peak width and resolution

- Time window: 5 000 × 1.925 ps ≈ 9.625 ns → frequency resolution 103.9 MHz.
- Peak FWHM (Lorentzian): f₀/Q ≈ 159 MHz.
- 1 000-bin DFT from 0.5 GHz to 1.5 GHz → 2 MHz bin spacing.
- Peak centre identifiable to ± bin-width / 2 = 1 MHz < 0.1 % of f₀.

---

## 3. Geometry

```
Grid : NX=5, NY=5, NZ=40, DX=1 mm   (5 × 5 × 40 mm, PEC on all faces)
Port : cell (2, 2, 20), z-axis
dt   : DX / (c √3) ≈ 1.925 ps
```

Lowest cavity mode of the 5×5×40 mm PEC box:

```
f_101 = (c/2) √((1/0.005)² + (1/0.040)²) ≈ 30.3 GHz >> f₀
```

So no cavity mode overlaps the LC resonance; the only sub-10 GHz dynamics
are the lumped LC oscillation.

### Stability checks

1. **Courant:** dt = DX/(c√3) < DX/c → standard 3-D Courant is satisfied.
2. **LC coupling coefficient** α_LC = dt²/(L ε₀ DX):
   = (1.925e-12)² / (1e-9 × 8.854e-12 × 1e-3) = 0.419 < 2 → stable.
3. **Nyquist:** f₀ = 1 GHz < f_Nyquist = 1/(2 dt) ≈ 260 GHz → satisfied.

---

## 4. Driver (Rust, in `yee-validation/src/lib.rs`)

```rust
pub struct LcResonanceResult {
    pub f_measured_hz: f64,
    pub f_analytic_hz: f64,
    pub rel_err: f64,
    pub passed: bool,          // |rel_err| < TOL_REL
}

pub fn fdtd206_run() -> LcResonanceResult { ... }
```

Steps:
1. Build 5×5×40 `YeeGrid::vacuum`.
2. Construct `LumpedRlcPort::series_rlc(cell, R, L, C, SourceWaveform::None)`.
3. **Kick phase** (N_KICK = 30 steps): inject a narrow Gaussian pulse into
   `grid.ez[(2,2,20)]` (broadband, sigma_t = 4 dt, centred at step 10).
4. **Ring-down phase** (N_RING = 5 000 steps): run with PEC boundary + LC
   correction; record `port.inductor_current()` each step.
5. **DFT scan**: 1 000 bins from 0.5 GHz to 1.5 GHz (2 MHz spacing);
   find bin with maximum amplitude.
6. Return `LcResonanceResult`.

### Gate parameters

```rust
FDTD206_L_H      : f64 = 1.0e-9;
FDTD206_F0_HZ    : f64 = 1.0e9;
FDTD206_C_F      : f64 = 1.0 / (4.0 * PI * PI * FDTD206_F0_HZ.powi(2) * FDTD206_L_H);
FDTD206_R_OHM    : f64 = 1.0;
FDTD206_TOL_F0   : f64 = 0.02;   // ±2 %
```

---

## 5. Integration test (`crates/yee-fdtd/tests/lumped_lc_resonance.rs`)

Self-contained test (no `yee-validation` dependency) that reimplements the
same physics.  Pattern: `ohmic_skin_depth.rs` (self-contained + detailed
diagnostic `eprintln!` + gate asserts).

---

## 6. Python wrapper (Phase 2.fdtd.py.6)

```python
from yee import run_lc_resonance
result = run_lc_resonance()
assert result.passed
print(result.f_measured_hz, result.f_analytic_hz, result.rel_err)
```

Pattern: same as `run_skin_depth()` in `yee-py/src/fdtd.rs`.

---

## 7. Tutorial 16

`docs/src/tutorials/16-fdtd-lumped-lc-resonance-from-python.md`

---

## 8. Aggregator registration

`run_fdtd_206_lumped_lc_resonance()` registered in `Report::run_all()`.
Expected status: **Passed** (< 0.1 s, NOT `#[ignore]`-gated).

---

## 9. DoD

- [ ] `fdtd206_run()` returns `LcResonanceResult` with `rel_err < 0.02`.
- [ ] `run_fdtd_206_lumped_lc_resonance()` in `run_all()` with status `Passed`.
- [ ] `cargo test -p yee-fdtd --test lumped_lc_resonance -- --nocapture` passes.
- [ ] `cargo test -p yee-validation -- fdtd_206 --nocapture` passes.
- [ ] `from yee import run_lc_resonance; assert run_lc_resonance().passed` passes.
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes.
- [ ] `cargo fmt --check --all` passes.
- [ ] ADR-0080 + spec + plan in repo; SUMMARY.md updated.
- [ ] Tutorial `16-fdtd-lumped-lc-resonance-from-python.md` in SUMMARY.md.

---

## 10. Risks

- **Frequency extraction accuracy**: DFT peak may be off by ±1 bin (1 MHz).
  Mitigation: use quadratic interpolation around the peak, or increase bin
  count.  At 2 MHz resolution the maximum bin-edge offset is 1 MHz = 0.1 %,
  well within the 2 % gate.
- **LC coupling instability**: α_LC = 0.419 < 2, so the explicit coupled scheme
  is stable (verified analytically in §3).  If divergence is observed, reduce
  DX or increase L.
- **Q too low for DFT**: with Q≈6.28 the signal decays 50× over the measurement
  window; the DFT SNR is still >> 1 for the kick amplitude used.

---

## References

- Pozar, "Microwave Engineering," 4th ed., §2.4 (series resonator).
- Hayt & Kemmerly, "Engineering Circuit Analysis," §14.1 (RLC natural
  frequency).
- Taflove & Hagness, "Computational Electrodynamics," 3rd ed., §15.10
  (lumped element FDTD update).
- ADR-0017 (Phase 2.fdtd.6 scope and explicit 2.fdtd.6.1 deferral).
