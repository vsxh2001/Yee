//! Gate `compute-018` (FS.0b.0, ADR-0208): constant [`GradedSpacings`]
//! arrays equal to the scalar `dx/dy/dz` must reproduce the uniform kernel
//! **bit-exactly** — probe series and every field component compared with
//! exact `f64` equality (the compute-007 idiom, no tolerance).
//!
//! The claim this pins: the kernel divides curl differences by primal/dual
//! *spacing* values (never precomputed inverses), and the uniform fill makes
//! every divisor bit-equal to the scalar, so the graded and uniform paths
//! execute identical FP operations. Scenarios cover both the CPML boundary
//! (bulk + ψ corrections) and the PEC box, with a soft source, a resistive
//! port, an aperture port, and probes, on asymmetric dims to catch any
//! axis mix-up in the spacing indexing.

use yee_compute::{
    AperturePort, Boundary, CpmlConfig, CpuFdtd, Drive, EComponent, FdtdSpec, Fields,
    GradedSpacings, Materials, Probe, ResistivePort, SoftSource, Waveform,
};

const NX: usize = 20;
const NY: usize = 16;
const NZ: usize = 12;
const DX: f64 = 1.0e-3;
const N_STEPS: usize = 150;

fn drive() -> Drive {
    let mut drive = Drive::default();
    drive.soft_sources.push(SoftSource {
        component: EComponent::Ez,
        cell: (10, 8, 6),
        waveform: Waveform::Gaussian {
            t0: 20.0 * 1.0e-12,
            sigma: 6.0 * 1.0e-12,
        },
    });
    drive.ports.push(ResistivePort {
        cell: (6, 8, 6),
        resistance: 50.0,
        waveform: Waveform::GaussianPulse {
            v0: 1.0,
            f0: 10.0e9,
            bw: 8.0e9,
            t0_steps: 30,
        },
    });
    drive.aperture_ports.push(AperturePort {
        cells: (7..10)
            .flat_map(|j| (0..3).map(move |k| (14, j, k)))
            .collect(),
        n_columns: 3,
        area: 3.0 * DX * (3.0 * DX),
        height: 3.0 * DX,
        resistance: 50.0,
        waveform: Waveform::GaussianPulse {
            v0: 0.5,
            f0: 12.0e9,
            bw: 10.0e9,
            t0_steps: 25,
        },
    });
    for cell in [(8, 8, 6), (12, 10, 4)] {
        drive.probes.push(Probe {
            component: EComponent::Ez,
            cell,
        });
    }
    drive
}

fn run(boundary: Boundary, graded: bool) -> CpuFdtd {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let mut engine = CpuFdtd::with_drive(
        spec,
        Fields::zero(&spec),
        Materials::default(),
        boundary,
        drive(),
    );
    if graded {
        engine.set_spacings(&GradedSpacings {
            dx: vec![DX; NX],
            dy: vec![DX; NY],
            dz: vec![DX; NZ],
        });
    }
    engine.step_n(N_STEPS);
    engine
}

fn assert_bit_identical(label: &str, uniform: &CpuFdtd, graded: &CpuFdtd) {
    for (name, a, b) in [
        (
            "probe 0",
            &uniform.probe_series()[0],
            &graded.probe_series()[0],
        ),
        (
            "probe 1",
            &uniform.probe_series()[1],
            &graded.probe_series()[1],
        ),
    ] {
        assert_eq!(a.len(), N_STEPS, "{label}/{name}: wrong series length");
        assert!(
            a.iter().any(|v| *v != 0.0),
            "{label}/{name}: probe stayed silent (test logic error)"
        );
        for (n, (u, g)) in a.iter().zip(b.iter()).enumerate() {
            assert!(
                u == g,
                "{label}/{name}: step {n} diverged (uniform {u:e}, graded {g:e})"
            );
        }
    }
    let (fu, fg) = (uniform.fields(), graded.fields());
    for (name, a, b) in [
        ("ex", &fu.ex, &fg.ex),
        ("ey", &fu.ey, &fg.ey),
        ("ez", &fu.ez, &fg.ez),
        ("hx", &fu.hx, &fg.hx),
        ("hy", &fu.hy, &fg.hy),
        ("hz", &fu.hz, &fg.hz),
    ] {
        let max_delta = a
            .iter()
            .zip(b.iter())
            .map(|(u, g)| (u - g).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_delta == 0.0,
            "{label}/{name}: max |Δ| = {max_delta:e} (must be exactly 0)"
        );
    }
}

#[test]
fn constant_spacings_are_bit_exact_under_cpml() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let boundary = Boundary::Cpml(CpmlConfig::for_spec(&spec, 4));
    let uniform = run(boundary.clone(), false);
    let graded = run(boundary, true);
    assert_bit_identical("cpml", &uniform, &graded);
}

#[test]
fn constant_spacings_are_bit_exact_under_pec_box() {
    let uniform = run(Boundary::PecBox, false);
    let graded = run(Boundary::PecBox, true);
    assert_bit_identical("pec", &uniform, &graded);
}
