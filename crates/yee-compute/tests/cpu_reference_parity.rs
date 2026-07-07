//! Gate `compute-001`: the rayon CPU backend is **bit-exact** against
//! `yee-fdtd`'s scalar FP64 reference kernels (uniform lossless vacuum,
//! PEC box).
//!
//! Every cell's half-step update is independent, so slab-parallelization must
//! not change a single bit — the assertion is max |Δ| == 0.0, not a
//! tolerance. Deliberately non-cubic dims so any index-order swap in the
//! flat-buffer port shows up as a mismatch.

use yee_compute::{CpuFdtd, FdtdSpec, Fields};
use yee_fdtd::YeeGrid;
use yee_fdtd::update::{update_e, update_h};

const NX: usize = 24;
const NY: usize = 20;
const NZ: usize = 22;
const DX: f64 = 1e-3;
const STEPS: usize = 25;

fn max_abs_diff(reference: &[f64], candidate: &[f64], name: &str) -> f64 {
    assert_eq!(reference.len(), candidate.len(), "{name} length mismatch");
    reference
        .iter()
        .zip(candidate)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max)
}

#[test]
fn cpu_backend_is_bit_exact_against_reference() {
    let mut grid = YeeGrid::vacuum(NX, NY, NZ, DX);

    // Same problem description; dt copied from the grid so both sides step
    // with the identical value regardless of how each derives its Courant
    // factor.
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = grid.dt;

    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.5);

    // Inject the identical Gaussian ball into the reference grid.
    let ezd = spec.ez_dims();
    for i in 0..ezd.0 {
        for j in 0..ezd.1 {
            for k in 0..ezd.2 {
                grid.ez[(i, j, k)] = init.ez[(i * ezd.1 + j) * ezd.2 + k];
            }
        }
    }

    let mut engine = CpuFdtd::new(spec, init);

    for _ in 0..STEPS {
        update_h(&mut grid);
        update_e(&mut grid);
    }
    engine.step_n(STEPS);

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
        let diff = max_abs_diff(reference, candidate, name);
        assert_eq!(
            diff, 0.0,
            "compute-001: {name} diverged from the scalar reference (max |Δ| = {diff:e})"
        );
    }

    // Anti-triviality: the pulse must actually have propagated — an
    // all-zeros implementation must not be able to pass the diff check.
    assert!(
        fields.hx.iter().any(|v| *v != 0.0),
        "compute-001: H_x never became non-zero; no propagation happened"
    );
    let center = fields.ez[((NX / 2) * ezd.1 + NY / 2) * ezd.2 + NZ / 2];
    assert!(
        (center - 1.0).abs() > 1e-6,
        "compute-001: the E_z pulse centre never evolved away from its initial value"
    );
}
