//! CPU H-probe unit test (FS.4.2a, deliverable 1): H sampling on the E.0
//! vacuum Gaussian scenario.
//!
//! `with_gaussian_ez` seeds an isotropic Gaussian on `E_z`
//! (`exp(-r²/2σ²)`), so `E_z(ci, cj+m, ck) == E_z(ci, cj-m, ck)` bit-for-bit
//! (`r_sq` is built from `dj*dj`, and `(-m)² == m²` exactly). `H_x` starts at
//! zero and its FIRST update is `H_x += coeff·(dE_y/dz − dE_z/dy)` with
//! `E_y == 0`, so `H_x(ci, j, ck) = −coeff·(E_z(ci, j+1, ck) −
//! E_z(ci, j, ck))` after exactly one step. Substituting the mirror
//! identity above gives an EXACT (not approximate) bit-negation:
//! `H_x(ci, cj−1, ck) == −H_x(ci, cj, ck)` at step 1 — IEEE 754 subtraction
//! and scalar multiplication are both sign-symmetric, so no floating-point
//! slop enters. This is the "known analytic situation" the FS.4.2a task
//! calls for; it does not depend on any interpretation of the stripline
//! gate physics.
//!
//! **E-probe regression** (existing consumer stays green, byte-identical):
//! `cpu_drive_parity.rs::driven_step_is_bit_exact_against_reference` is
//! untouched functionally by this change (its `Drive` literal only grew an
//! empty `h_probes: vec![]`) and still asserts the E-probe series matches
//! the `yee_fdtd` reference bit-for-bit — that is the regression evidence.

use yee_compute::{Boundary, CpuFdtd, Drive, FdtdSpec, Fields, HComponent, HProbe, Materials};

const NX: usize = 16;
const NY: usize = 16;
const NZ: usize = 16;
const DX: f64 = 1e-3;
const CENTER: (usize, usize, usize) = (8, 8, 8);
const SIGMA_CELLS: f64 = 2.0;

#[test]
fn h_probe_exact_bit_antisymmetry_about_gaussian_source_first_step() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, CENTER, SIGMA_CELLS);

    let drive = Drive {
        soft_sources: vec![],
        ports: vec![],
        aperture_ports: vec![],
        thin_wires: vec![],
        probes: vec![],
        h_probes: vec![
            // (ci, cj-1, ck): one cell "below" the source plane.
            HProbe {
                component: HComponent::Hx,
                cell: (CENTER.0, CENTER.1 - 1, CENTER.2),
            },
            // (ci, cj, ck): the source plane itself.
            HProbe {
                component: HComponent::Hx,
                cell: CENTER,
            },
        ],
    };
    let mut engine = CpuFdtd::with_drive(spec, init, Materials::default(), Boundary::None, drive);
    engine.step_n(1);

    let series = engine.h_probe_series();
    assert_eq!(series.len(), 2);
    assert_eq!(series[0].len(), 1);
    assert_eq!(series[1].len(), 1);
    let (h_below, h_at) = (series[0][0], series[1][0]);

    assert_ne!(
        h_below, 0.0,
        "h-probe never saw the source — scenario broken"
    );
    assert_eq!(
        h_below, -h_at,
        "H_x(cj-1) should be the EXACT bit-negation of H_x(cj) after the \
         first update_h on a mirror-symmetric E_z Gaussian (h_below={h_below:e}, h_at={h_at:e})"
    );
}

/// Longer run: the recorded H-probe stream stays finite and non-trivial
/// across multiple steps (the "stream nonzero" half of the deliverable),
/// once the field has left the single-step regime the exact identity above
/// was derived for.
#[test]
fn h_probe_stream_nonzero_and_finite_over_several_steps() {
    let spec = FdtdSpec::vacuum(NX, NY, NZ, DX);
    let init = Fields::with_gaussian_ez(&spec, CENTER, SIGMA_CELLS);
    const STEPS: usize = 5;

    let drive = Drive {
        soft_sources: vec![],
        ports: vec![],
        aperture_ports: vec![],
        thin_wires: vec![],
        probes: vec![],
        h_probes: vec![HProbe {
            component: HComponent::Hx,
            cell: (CENTER.0, CENTER.1 - 1, CENTER.2),
        }],
    };
    let mut engine = CpuFdtd::with_drive(spec, init, Materials::default(), Boundary::None, drive);
    engine.step_n(STEPS);

    let series = &engine.h_probe_series()[0];
    assert_eq!(series.len(), STEPS);
    assert!(series.iter().all(|v| v.is_finite()));
    assert!(
        series.iter().any(|v| *v != 0.0),
        "h-probe stream never left zero across {STEPS} steps"
    );
}
