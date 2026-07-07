//! Gate `compute-007`: the E.2 driven CPU step — soft sources + resistive
//! port + probes — is **bit-exact** against the reference orchestration
//! (`WalkingSkeletonSolver` sub-step helpers + `LumpedRlcPort::pure_resistor`
//! with a `GaussianPulse` EMF), including the probe series.

use yee_compute::{
    Boundary, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, Probe, ResistivePort,
    SoftSource, Waveform,
};
use yee_fdtd::{FdtdSolver, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid};

const NX: usize = 16;
const NY: usize = 12;
const NZ: usize = 14;
const DX: f64 = 1e-3;
const STEPS: usize = 25;

const SOFT_CELL: (usize, usize, usize) = (4, 6, 7); // E_y
const PORT_CELL: (usize, usize, usize) = (10, 6, 7); // E_z
const PROBE_CELL: (usize, usize, usize) = (12, 6, 7); // E_z

#[test]
fn driven_step_is_bit_exact_against_reference() {
    let grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid.dt;
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = dt;

    let t0 = 8.0 * dt;
    let sigma = 3.0 * dt;
    let (v0, f0) = (1.0, 5.0e9);
    let bw = 0.8 * f0;
    let t0_steps = 6usize;
    let resistance = 50.0;

    // ---- reference: manual step body (the cavity / line-eeff pattern) ----
    let mut solver = WalkingSkeletonSolver::new(grid);
    let mut port = LumpedRlcPort::pure_resistor(
        PORT_CELL,
        resistance,
        SourceWaveform::GaussianPulse {
            v0,
            f0,
            bw,
            t0_steps,
        },
    );
    let mut ref_probe = Vec::with_capacity(STEPS);
    for n in 0..STEPS {
        let t = solver.current_time();
        solver.update_h_only();
        solver.apply_cpml_h();
        {
            let arg = (t - t0) / sigma;
            solver.grid_mut().ey[SOFT_CELL] += (-arg * arg).exp();
        }
        solver.update_e_only();
        solver.apply_cpml_e();
        port.correct_e(solver.grid_mut(), n, dt);
        solver.advance_clock();
        ref_probe.push(solver.grid().ez[PROBE_CELL]);
    }

    // ---- candidate: CpuFdtd::with_drive ----
    let drive = Drive {
        soft_sources: vec![SoftSource {
            component: EComponent::Ey,
            cell: SOFT_CELL,
            waveform: Waveform::Gaussian { t0, sigma },
        }],
        ports: vec![ResistivePort {
            cell: PORT_CELL,
            resistance,
            waveform: Waveform::GaussianPulse {
                v0,
                f0,
                bw,
                t0_steps,
            },
        }],
        aperture_ports: vec![],
        probes: vec![Probe {
            component: EComponent::Ez,
            cell: PROBE_CELL,
        }],
    };
    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        Boundary::PecBox,
        drive,
    );
    engine.step_n(STEPS);

    // Probe series bit-exact.
    let cand_probe = &engine.probe_series()[0];
    assert_eq!(cand_probe.len(), STEPS);
    for (n, (r, c)) in ref_probe.iter().zip(cand_probe).enumerate() {
        assert_eq!(
            r, c,
            "compute-007: probe diverged at step {n} ({r:e} vs {c:e})"
        );
    }

    // Full field state bit-exact.
    let grid = solver.grid();
    let fields = engine.fields();
    let pairs: [(&[f64], &[f64], &str); 6] = [
        (grid.ex.as_slice().unwrap(), &fields.ex, "ex"),
        (grid.ey.as_slice().unwrap(), &fields.ey, "ey"),
        (grid.ez.as_slice().unwrap(), &fields.ez, "ez"),
        (grid.hx.as_slice().unwrap(), &fields.hx, "hx"),
        (grid.hy.as_slice().unwrap(), &fields.hy, "hy"),
        (grid.hz.as_slice().unwrap(), &fields.hz, "hz"),
    ];
    for (reference, candidate, name) in pairs {
        let diff = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert_eq!(
            diff, 0.0,
            "compute-007: {name} diverged from the reference (max |Δ| = {diff:e})"
        );
    }
    assert!(
        ref_probe.iter().any(|v| *v != 0.0),
        "compute-007: probe never saw the drive — scenario broken"
    );
}
