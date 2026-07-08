//! Gate `compute-014` (S.10, ADR-0187): the engine's multi-cell **aperture
//! port** — modal `V = ∫E_z·dz`, aggregate branch current, sheet-current
//! back-action over the physical area — is **bit-exact** against the
//! validated reference `yee_fdtd::LumpedRlcPort::aperture` /
//! `correct_e_aperture` (Phase 2.fdtd.6.9, ADR-0125), pure-R arm, driven
//! by a `GaussianPulse` EMF, including the full final field state.

use yee_compute::{
    AperturePort, Boundary, CpuFdtd, Drive, EComponent, FdtdSpec, Fields, Materials, Probe,
    Waveform,
};
use yee_fdtd::{ApertureSpec, LumpedRlcPort, SourceWaveform, WalkingSkeletonSolver, YeeGrid};

const NX: usize = 18;
const NY: usize = 14;
const NZ: usize = 12;
const DX: f64 = 1e-3;
const STEPS: usize = 30;

const PORT_I: usize = 6;
const J_LO: usize = 5; // width band [5, 8) → 3 columns
const J_HI: usize = 8;
const K_TOP: usize = 4; // substrate column k = 0..4
const PROBE_CELL: (usize, usize, usize) = (12, 6, 2);

#[test]
fn aperture_port_is_bit_exact_against_reference() {
    let grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let dt = grid.dt;
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = dt;

    let cells: Vec<(usize, usize, usize)> = (J_LO..J_HI)
        .flat_map(|j| (0..K_TOP).map(move |k| (PORT_I, j, k)))
        .collect();
    let n_columns = J_HI - J_LO;
    let height = K_TOP as f64 * DX;
    let area = n_columns as f64 * DX * height;
    let (v0, f0) = (1.0, 5.0e9);
    let bw = 0.8 * f0;
    let t0_steps = 6usize;
    let resistance = 50.0;

    // ---- reference: WalkingSkeletonSolver + LumpedRlcPort::aperture ----
    let mut solver = WalkingSkeletonSolver::new(grid);
    let mut port = LumpedRlcPort::aperture(
        ApertureSpec {
            cells: cells.clone(),
            n_columns,
            area,
            height,
        },
        resistance,
        0.0,           // L = 0 (pure R)
        f64::INFINITY, // C = ∞ (no DC block)
        SourceWaveform::GaussianPulse {
            v0,
            f0,
            bw,
            t0_steps,
        },
    );
    let mut ref_probe = Vec::with_capacity(STEPS);
    for n in 0..STEPS {
        solver.update_h_only();
        solver.apply_cpml_h();
        solver.update_e_only();
        solver.apply_cpml_e();
        port.correct_e_aperture(solver.grid_mut(), n, dt);
        solver.advance_clock();
        ref_probe.push(solver.grid().ez[PROBE_CELL]);
    }

    // ---- candidate: CpuFdtd with Drive::aperture_ports ----
    let drive = Drive {
        soft_sources: vec![],
        ports: vec![],
        aperture_ports: vec![AperturePort {
            cells,
            n_columns,
            area,
            height,
            resistance,
            waveform: Waveform::GaussianPulse {
                v0,
                f0,
                bw,
                t0_steps,
            },
            record: false,
        }],
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
            "compute-014: probe diverged at step {n} ({r:e} vs {c:e})"
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
            "compute-014: {name} diverged from the reference (max |Δ| = {diff:e})"
        );
    }
    assert!(
        ref_probe.iter().any(|v| *v != 0.0),
        "compute-014: probe never saw the drive — scenario broken"
    );
}
