//! Gate `compute-010` (E.5a): far-field extraction on the engine — the
//! `yee-fdtd` NTFF dipole gate reproduced with the FDTD stepped by
//! `yee-compute`'s CPU backend, checked against the **analytic dipole
//! pattern**: an `E_z`-polarised point source radiates as `sin θ`, so the
//! broadside (θ = π/2) to endfire (θ = 0, analytic null) ratio must be
//! ≥ 20 dB.
//!
//! The engine owns the stepping (CPML + soft Gaussian source, the exact
//! reference orchestration per gate `compute-007`); `yee_fdtd::NtffState`
//! stays the reference near-to-far transform and consumes the engine's
//! fields through a host-side grid adapter each step. A first-class,
//! GPU-side NTFF accumulator is a later phase (see ENGINE-STUDIO-ROADMAP).
//!
//! ```bash
//! cargo test -p yee-compute --release --test ntff_dipole -- --ignored --nocapture
//! ```

use std::f64::consts::FRAC_PI_2;

use yee_compute::{
    Boundary, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, SoftSource,
    Waveform,
};
use yee_fdtd::{NtffParams, NtffState, YeeGrid};

const N: usize = 50;
const DX: f64 = 1.0e-3;
const NPML: usize = 10;
const N_STEPS: usize = 2000;
const F_PROBE: f64 = 15.0e9;
const SRC: (usize, usize, usize) = (25, 25, 25);
const BOX_MARGIN_CELLS: usize = NPML + 5;

/// Copy the engine's flat fields into a `YeeGrid` (same staggered shapes,
/// same row-major order) so the reference NTFF accumulator can sample them.
fn copy_into_grid(fields: &Fields, grid: &mut YeeGrid) {
    grid.ex.as_slice_mut().unwrap().copy_from_slice(&fields.ex);
    grid.ey.as_slice_mut().unwrap().copy_from_slice(&fields.ey);
    grid.ez.as_slice_mut().unwrap().copy_from_slice(&fields.ez);
    grid.hx.as_slice_mut().unwrap().copy_from_slice(&fields.hx);
    grid.hy.as_slice_mut().unwrap().copy_from_slice(&fields.hy);
    grid.hz.as_slice_mut().unwrap().copy_from_slice(&fields.hz);
}

#[test]
#[ignore = "slow: 2000-step 50^3 CPML run (release ~10-30 s); compute-010 NTFF dipole gate (E.5a)"]
fn ntff_on_engine_recovers_dipole_pattern() {
    let mut scratch = YeeGrid::vacuum(N, N, N, DX);
    let dt = scratch.dt;
    let mut spec = FdtdSpec::vacuum(N, N, N, DX);
    spec.dt = dt;

    let t0 = 12.0 * dt;
    let sigma = 4.0 * dt;
    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML)),
        Drive {
            soft_sources: vec![SoftSource {
                component: EComponent::Ez,
                cell: SRC,
                waveform: Waveform::Gaussian { t0, sigma },
            }],
            ports: vec![],
            aperture_ports: vec![],
            probes: vec![],
            h_probes: vec![],
        },
    );

    let mut ntff = NtffState::new(
        &scratch,
        NtffParams {
            f_probe: F_PROBE,
            box_margin_cells: BOX_MARGIN_CELLS,
            theta_rad: FRAC_PI_2,
            phi_rad: 0.0,
        },
    );

    for _ in 0..N_STEPS {
        engine.step_n(1);
        copy_into_grid(engine.fields(), &mut scratch);
        ntff.sample(&scratch, engine.current_time());
    }
    assert_eq!(ntff.n_samples(), N_STEPS as u64);

    let mag_broad = ntff.far_field_at(FRAC_PI_2, 0.0).norm();
    let mag_end = ntff.far_field_at(0.0, 0.0).norm();
    assert!(
        mag_broad.is_finite() && mag_end.is_finite(),
        "non-finite far field"
    );
    assert!(mag_broad > 0.0, "broadside is zero — no radiation captured");
    let ratio = if mag_end > 0.0 {
        mag_broad / mag_end
    } else {
        f64::INFINITY
    };
    let db = 20.0 * ratio.log10();
    eprintln!(
        "compute-010: |E(broadside)| = {mag_broad:.3e}, |E(endfire)| = {mag_end:.3e}, \
         ratio = {db:.2} dB"
    );
    assert!(
        db >= 20.0,
        "compute-010 FAILED: broadside/endfire {db:.2} dB < 20 dB (analytic sin θ null)"
    );
}
