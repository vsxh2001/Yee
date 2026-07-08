//! Gate `engine-scale-001` (FS.2b, ADR-0207): the NTFF **absolute scale**
//! against the analytic Hertzian dipole — the one reference where the
//! source's current moment is known *exactly* from the injection
//! semantics, independent of any port model.
//!
//! A soft `E_z += s(t)` add per step is the equivalent current density
//! `J_z = −ε₀·s(t)/dt` over one cell, i.e. a Hertzian moment
//! `I·dl(ω) = ε₀·dx³·S(ω)/dt` (S = continuous-transform DFT of s). The
//! analytic far-field pattern amplitude (same e^{−jkr}/r-dropped
//! convention as `NtffState`):
//!
//! ```text
//! |F(θ)| = ω·μ₀·|I·dl(ω)|·sin θ / 4π
//! ```
//!
//! Every earlier NTFF gate (E.5, A.2, FS.1a.2, FS.1b.2) certified pattern
//! SHAPE (ratios); this is the first absolute-magnitude pin. **GREEN,
//! measured across three (dx, f) configs**: ratio 1.048 at θ = 90° and
//! 1.029 at 45° (dx ∈ {1.5, 2.0, 2.5} mm × f ∈ {1.8, 2.45, 3.2} GHz,
//! reproducible to 3 decimal places) — the free-space NTFF absolute
//! scale is right to ~3–5 % (single-cell discrete moment + box
//! truncation). Two measured lessons on the way: (1) a BASEBAND Gaussian
//! source's near-DC content survives the CPML and leaks into the
//! single-bin DFT (±40 % direction-dependent scatter measured) — hence
//! the zero-DC `GaussianPulseEz`; (2) the first gain run's +13 dB excess
//! is therefore NOT this transform — the patch fixture's equivalence box
//! intersects the whole-domain substrate slab (see engine-gain-001).
//!
//! ```bash
//! cargo test -p yee-engine --release --test ntff_absolute_scale -- --ignored --nocapture
//! ```

use yee_engine::sparams::single_bin_dft;
use yee_engine::{BackendChoice, BoundarySpec, JobEvent, JobSpec, NtffSpec, SourceSpec};

const N: usize = 90;
const DX_M: f64 = 2.0e-3;
const F_PROBE_HZ: f64 = 2.45e9;
const N_STEPS: usize = 3000;
const T0_STEPS: usize = 300;
const BW_HZ: f64 = 1.5e9;
const EPS0: f64 = 8.854_187_812_8e-12;
const MU0: f64 = 4.0e-7 * std::f64::consts::PI;

#[test]
#[ignore = "slow: one ~1 min release FDTD run + NTFF; engine-scale-001 gate (FS.2b) — run with --release --ignored"]
fn ntff_magnitude_matches_the_analytic_hertzian_dipole() {
    // Sweep knobs for the constancy forensic (defaults = the pinned run).
    let dx_m = std::env::var("YEE_SCALE_DX_MM")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .map(|mm| mm * 1e-3)
        .unwrap_or(DX_M);
    let f_probe = std::env::var("YEE_SCALE_F_GHZ")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .map(|g| g * 1e9)
        .unwrap_or(F_PROBE_HZ);
    let c = (N / 2, N / 2, N / 2);
    let deg = |d: f64| d * std::f64::consts::PI / 180.0;
    let directions = vec![(deg(90.0), 0.0), (deg(90.0), deg(90.0)), (deg(45.0), 0.0)];

    let spec = JobSpec {
        nx: N,
        ny: N,
        nz: N,
        dx_m,
        n_steps: N_STEPS,
        boundary: BoundarySpec::Cpml {
            npml: 10,
            axes: [true, true, true],
            faces: None,
        },
        // Modulated (zero-DC) pulse: a baseband Gaussian's near-DC
        // content survives the CPML and leaks into the single-bin DFT —
        // measured as ±40 % direction-dependent scatter in this very
        // gate before the switch.
        sources: vec![SourceSpec::GaussianPulseEz {
            cell: c,
            v0: 1.0,
            f0_hz: f_probe,
            bw_hz: BW_HZ,
            t0_steps: T0_STEPS,
        }],
        ports: vec![],
        aperture_ports: vec![],
        probes: vec![],
        slice: None,
        ntff: Some(NtffSpec {
            f_hz: f_probe,
            margin_cells: 15,
            k_min: None,
            directions: directions.clone(),
        }),
        materials: None,
        dt_s: None,
        backend: BackendChoice::Cpu,
    };

    let handle = yee_engine::submit(spec);
    let result = handle
        .events()
        .find_map(|e| match e {
            JobEvent::Done { result } => Some(result),
            JobEvent::Error { message } => panic!("job failed: {message}"),
            _ => None,
        })
        .expect("no Done event");
    let dt = result.dt_s;
    let ff = result.far_field.expect("no far field");

    // The injected waveform, reconstructed exactly
    // (Waveform::GaussianPulse: v0·exp(−((t−t0)/τ)²)·sin(2πf0(t−t0)),
    // τ = √(2ln2)/(π·bw)).
    let t0 = T0_STEPS as f64 * dt;
    let tau = (2.0_f64 * std::f64::consts::LN_2).sqrt() / (std::f64::consts::PI * BW_HZ);
    let s_series: Vec<f64> = (0..N_STEPS)
        .map(|n| {
            let t = n as f64 * dt;
            let arg = (t - t0) / tau;
            (-arg * arg).exp() * (std::f64::consts::TAU * f_probe * (t - t0)).sin()
        })
        .collect();
    let (sr, si) = single_bin_dft(&s_series, dt, f_probe);
    let s_mag = (sr * sr + si * si).sqrt() * dt; // continuous-transform scale

    let omega = std::f64::consts::TAU * f_probe;
    let idl = EPS0 * dx_m.powi(3) * s_mag / dt;
    eprintln!("engine-scale-001: |S(f)| = {s_mag:.4e} V·s, |I·dl| = {idl:.4e} A·m, dt = {dt:.3e}");

    for ((theta, phi), &f_ntff) in directions.iter().zip(&ff) {
        let f_analytic = omega * MU0 * idl * theta.sin() / (4.0 * std::f64::consts::PI);
        let ratio = f_ntff / f_analytic;
        eprintln!(
            "  θ = {:>3.0}°, φ = {:>3.0}°: NTFF {f_ntff:.4e}, analytic {f_analytic:.4e}, \
             ratio = {ratio:.4}",
            theta.to_degrees(),
            phi.to_degrees(),
        );
        // Pinned from the 3-config sweep (1.028–1.049 measured).
        assert!(
            (0.9..=1.15).contains(&ratio),
            "engine-scale-001 FAILED: NTFF/analytic ratio {ratio:.4} at θ = {:.0}° \
             outside [0.8, 1.25] — the NTFF absolute scale is off",
            theta.to_degrees(),
        );
    }
}
