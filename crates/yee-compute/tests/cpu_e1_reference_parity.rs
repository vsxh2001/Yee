//! Gate `compute-003`: the E.1 CPU backend is **bit-exact** against
//! `yee_fdtd::WalkingSkeletonSolver` on a heterogeneous scenario exercising
//! every new arm at once:
//!
//! - per-cell ε_r (dielectric slab), σ (lossy block, CA/CB path), μ_r
//!   (magnetic region),
//! - an interior PEC sheet with a slot (per-component masks),
//! - a driven Gaussian soft source,
//! - both boundary modes: CPML (Roden–Gedney) and the legacy PEC box.
//!
//! As with `compute-001`, the assertion is max |Δ| == 0.0 on all six field
//! components — the ports are line-for-line and parallelized only across
//! independent cells.

use ndarray::Array3;
use yee_compute::{Boundary, CpmlConfig, CpuFdtd, FdtdSpec, Fields, Materials};
use yee_fdtd::{CpmlParams, WalkingSkeletonSolver, YeeGrid};

const NX: usize = 24;
const NY: usize = 20;
const NZ: usize = 22;
const DX: f64 = 1e-3;
const NPML: usize = 5;
const STEPS: usize = 30;
const SOURCE: (usize, usize, usize) = (12, 10, 11);

/// Heterogeneous material description, produced once as flat vectors and
/// converted to `Array3` for the reference grid so both sides see the exact
/// same values.
struct Scenario {
    eps_r: Vec<f64>,
    mu_r: Vec<f64>,
    sigma: Vec<f64>,
    mask_ex: Vec<bool>,
    mask_ey: Vec<bool>,
}

fn scenario() -> Scenario {
    let celld = (NX + 1, NY + 1, NZ + 1);
    let n_cells = celld.0 * celld.1 * celld.2;
    let cell = |i: usize, j: usize, k: usize| (i * celld.1 + j) * celld.2 + k;

    let mut eps_r = vec![1.0; n_cells];
    let mut mu_r = vec![1.0; n_cells];
    let mut sigma = vec![0.0; n_cells];
    for i in 0..celld.0 {
        for j in 0..celld.1 {
            for k in 0..celld.2 {
                // Dielectric slab (FR-4-ish) across the lower z range.
                if (3..8).contains(&k) {
                    eps_r[cell(i, j, k)] = 4.3;
                }
                // Lossy block in the middle of the domain.
                if (5..10).contains(&i) && (5..10).contains(&j) && (9..13).contains(&k) {
                    sigma[cell(i, j, k)] = 0.5;
                }
                // Magnetic region across a y band.
                if (12..16).contains(&j) {
                    mu_r[cell(i, j, k)] = 2.0;
                }
            }
        }
    }

    // Interior PEC sheet at k = 15 clamping in-plane E (E_x, E_y), with a
    // slot at j ∈ [8, 12) left open.
    let exd = (NX, NY + 1, NZ + 1);
    let eyd = (NX + 1, NY, NZ + 1);
    let mut mask_ex = vec![false; exd.0 * exd.1 * exd.2];
    let mut mask_ey = vec![false; eyd.0 * eyd.1 * eyd.2];
    for i in 0..exd.0 {
        for j in 0..exd.1 {
            if !(8..12).contains(&j) {
                mask_ex[(i * exd.1 + j) * exd.2 + 15] = true;
            }
        }
    }
    for i in 0..eyd.0 {
        for j in 0..eyd.1 {
            if !(8..12).contains(&j) {
                mask_ey[(i * eyd.1 + j) * eyd.2 + 15] = true;
            }
        }
    }

    Scenario {
        eps_r,
        mu_r,
        sigma,
        mask_ex,
        mask_ey,
    }
}

fn reference_grid(sc: &Scenario) -> YeeGrid {
    let celld = (NX + 1, NY + 1, NZ + 1);
    let exd = (NX, NY + 1, NZ + 1);
    let eyd = (NX + 1, NY, NZ + 1);
    YeeGrid::vacuum(NX, NY, NZ, DX)
        .with_eps_r_cells(Array3::from_shape_vec(celld, sc.eps_r.clone()).unwrap())
        .with_mu_r_cells(Array3::from_shape_vec(celld, sc.mu_r.clone()).unwrap())
        .with_sigma_cells(Array3::from_shape_vec(celld, sc.sigma.clone()).unwrap())
        .with_pec_mask_ex(Array3::from_shape_vec(exd, sc.mask_ex.clone()).unwrap())
        .with_pec_mask_ey(Array3::from_shape_vec(eyd, sc.mask_ey.clone()).unwrap())
}

fn candidate_materials(sc: &Scenario) -> Materials {
    Materials {
        eps_r_cells: Some(sc.eps_r.clone()),
        mu_r_cells: Some(sc.mu_r.clone()),
        sigma_cells: Some(sc.sigma.clone()),
        pec_mask_ex: Some(sc.mask_ex.clone()),
        pec_mask_ey: Some(sc.mask_ey.clone()),
        pec_mask_ez: None,
        sheet_r_ohm: None,
    }
}

/// Drive both sides for `STEPS` steps with the identical Gaussian source
/// and assert bit-exact agreement on all six components.
fn assert_bit_exact(mut solver: WalkingSkeletonSolver, mut candidate: CpuFdtd, label: &str) {
    let dt = solver.dt();
    let t0 = 10.0 * dt;
    let sigma = 4.0 * dt;
    for _ in 0..STEPS {
        solver.step_with_source(SOURCE.0, SOURCE.1, SOURCE.2, t0, sigma);
        candidate.step_with_gaussian_ez(SOURCE, t0, sigma);
    }

    let grid = solver.grid();
    let fields = candidate.fields();
    let pairs: [(&[f64], &[f64], &str); 6] = [
        (grid.ex.as_slice().unwrap(), &fields.ex, "ex"),
        (grid.ey.as_slice().unwrap(), &fields.ey, "ey"),
        (grid.ez.as_slice().unwrap(), &fields.ez, "ez"),
        (grid.hx.as_slice().unwrap(), &fields.hx, "hx"),
        (grid.hy.as_slice().unwrap(), &fields.hy, "hy"),
        (grid.hz.as_slice().unwrap(), &fields.hz, "hz"),
    ];
    for (reference, cand, name) in pairs {
        assert_eq!(
            reference.len(),
            cand.len(),
            "{label}/{name} length mismatch"
        );
        let diff = reference
            .iter()
            .zip(cand)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert_eq!(
            diff, 0.0,
            "compute-003 [{label}]: {name} diverged from the reference (max |Δ| = {diff:e})"
        );
    }

    // Anti-triviality: the source must actually have propagated.
    assert!(
        fields.hx.iter().any(|v| *v != 0.0),
        "compute-003 [{label}]: H_x never became non-zero"
    );
}

#[test]
fn cpml_heterogeneous_is_bit_exact_against_reference() {
    let sc = scenario();
    let grid = reference_grid(&sc);
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = grid.dt;

    let params = CpmlParams::for_grid(&grid, NPML);
    let solver = WalkingSkeletonSolver::with_cpml(grid, params);
    let candidate = CpuFdtd::with_config(
        spec,
        Fields::zero(&spec),
        candidate_materials(&sc),
        Boundary::Cpml(CpmlConfig::for_spec(&spec, NPML)),
    );
    assert_bit_exact(solver, candidate, "cpml");
}

#[test]
fn pec_box_heterogeneous_is_bit_exact_against_reference() {
    let sc = scenario();
    let grid = reference_grid(&sc);
    let mut spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    spec.dt = grid.dt;

    let solver = WalkingSkeletonSolver::new(grid);
    let candidate = CpuFdtd::with_config(
        spec,
        Fields::zero(&spec),
        candidate_materials(&sc),
        Boundary::PecBox,
    );
    assert_bit_exact(solver, candidate, "pec-box");
}
