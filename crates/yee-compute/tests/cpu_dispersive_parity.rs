//! Gate `compute-011` (E.5c): the ADE dispersive E update on the CPU
//! backend is **bit-exact** against `yee_fdtd::dispersive::DispersiveState`
//! on a grid carrying all four material arms at once — a Drude block, a
//! Lorentz block, a Debye block, and vacuum elsewhere — over 25 steps from
//! a Gaussian E_z initial condition.

use yee_compute::{CpuFdtd, DispersiveMap, DispersiveMaterial, FdtdSpec, Fields};
use yee_fdtd::update::update_h;
use yee_fdtd::{DispersiveState, Material, MaterialMap, YeeGrid};

const NX: usize = 18;
const NY: usize = 16;
const NZ: usize = 14;
const DX: f64 = 1e-3;
const STEPS: usize = 25;

const DRUDE: (f64, f64, f64) = (1.0, 2.0e10 * std::f64::consts::TAU, 1.0e9);
const LORENTZ: (f64, f64, f64, f64) = (2.0, 1.5, 1.0e10 * std::f64::consts::TAU, 2.0e8);
const DEBYE: (f64, f64, f64) = (1.5, 8.0, 5.0e-11);

#[test]
fn dispersive_update_is_bit_exact_against_reference() {
    // ---- reference: YeeGrid + MaterialMap + DispersiveState ----
    let mut grid = YeeGrid::vacuum(NX, NY, NZ, DX);
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = grid.dt;
    let init = Fields::with_gaussian_ez(&spec, (NX / 2, NY / 2, NZ / 2), 2.0);
    let ezd = spec.ez_dims();
    for i in 0..ezd.0 {
        for j in 0..ezd.1 {
            for k in 0..ezd.2 {
                grid.ez[(i, j, k)] = init.ez[(i * ezd.1 + j) * ezd.2 + k];
            }
        }
    }

    let mut ref_map = MaterialMap::vacuum(NX, NY, NZ);
    let (di, dp, dg) = DRUDE;
    ref_map.set_box(
        2,
        7,
        2,
        7,
        2,
        7,
        Material::Drude {
            eps_inf: di,
            omega_p: dp,
            gamma: dg,
        },
    );
    let (li, ld, lo, lde) = LORENTZ;
    ref_map.set_box(
        10,
        16,
        3,
        9,
        4,
        10,
        Material::Lorentz {
            eps_inf: li,
            delta_eps: ld,
            omega_0: lo,
            delta: lde,
        },
    );
    let (bi, bd, bt) = DEBYE;
    ref_map.set_box(
        4,
        12,
        10,
        15,
        6,
        12,
        Material::Debye {
            eps_inf: bi,
            delta_eps: bd,
            tau: bt,
        },
    );
    let mut ref_state = DispersiveState::new(&ref_map);
    for _ in 0..STEPS {
        update_h(&mut grid);
        ref_state.update_e_with_dispersion(&mut grid, &ref_map);
    }

    // ---- candidate: CpuFdtd + DispersiveMap (identical regions) ----
    let mut map = DispersiveMap::vacuum(&spec);
    map.set_box(
        &spec,
        2,
        7,
        2,
        7,
        2,
        7,
        DispersiveMaterial::Drude {
            eps_inf: di,
            omega_p: dp,
            gamma: dg,
        },
    );
    map.set_box(
        &spec,
        10,
        16,
        3,
        9,
        4,
        10,
        DispersiveMaterial::Lorentz {
            eps_inf: li,
            delta_eps: ld,
            omega_0: lo,
            delta: lde,
        },
    );
    map.set_box(
        &spec,
        4,
        12,
        10,
        15,
        6,
        12,
        DispersiveMaterial::Debye {
            eps_inf: bi,
            delta_eps: bd,
            tau: bt,
        },
    );
    let mut engine = CpuFdtd::new(spec, init);
    engine.set_dispersive(map);
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
        let diff = reference
            .iter()
            .zip(candidate)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert_eq!(
            diff, 0.0,
            "compute-011: {name} diverged from the dispersive reference (max |Δ| = {diff:e})"
        );
    }
    assert!(
        fields.hx.iter().any(|v| *v != 0.0),
        "compute-011: no propagation happened"
    );
}
